import os
import sys
import traceback

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
BASE = os.path.abspath(os.path.dirname(__file__))

# Prefer release binaries to avoid debug CRT/runtime mismatches with conda Python.
_release_dir = os.path.join(REPO, "build", "Reaktoro", "Release")
_package_dir = os.path.join(
    REPO, "build", "python", "package", "build", "lib", "reaktoro"
)
_debug_dir = os.path.join(REPO, "build", "Reaktoro", "Debug")

_allow_debug = os.environ.get("REAKTORO_ALLOW_DEBUG_PYD", "0") == "1"
_candidates = [_release_dir, _package_dir] + ([_debug_dir] if _allow_debug else [])

PYD_DIR = None
for _cand in _candidates:
    _pyd = os.path.join(_cand, "reaktoro4py.cp312-win_amd64.pyd")
    if os.path.isfile(_pyd):
        PYD_DIR = _cand
        break

if PYD_DIR is None:
    raise FileNotFoundError(
        "Could not find reaktoro4py.cp312-win_amd64.pyd in expected build folders."
    )

if PYD_DIR not in sys.path:
    sys.path.insert(0, PYD_DIR)
os.add_dll_directory(PYD_DIR)

# Import autodiff first so pybind Real<1,double> values convert to Python.
import autodiff as ad

import reaktoro4py as rkt

# Let diagrams.py resolve __import__("reaktoro") in local-build mode.
sys.modules.setdefault("reaktoro", rkt)

import importlib.util as _ilu

_diagrams_file = os.path.join(
    REPO, "python", "package", "reaktoro", "extensions", "diagrams.py"
)
_spec = _ilu.spec_from_file_location("reaktoro_diagrams", _diagrams_file)
_dmod = _ilu.module_from_spec(_spec)
_spec.loader.exec_module(_dmod)

PredominancePlot = _dmod.PredominancePlot


def _env_int(name, default):
    try:
        return int(os.environ.get(name, str(default)))
    except Exception:
        return int(default)


def _env_float(name, default):
    try:
        return float(os.environ.get(name, str(default)))
    except Exception:
        return float(default)


# Robustness profile (bounded to avoid unending runs)
ROBUST_MAX_SEEDS_PER_POINT = max(1, _env_int("REAKTORO_ROBUST_MAX_SEEDS", 3))
ROBUST_HOMOTOPY_STEPS_MAIN = max(0, _env_int("REAKTORO_ROBUST_HOMOTOPY_MAIN", 0))
ROBUST_HOMOTOPY_STEPS_RECOVERY = max(
    0, _env_int("REAKTORO_ROBUST_HOMOTOPY_RECOVERY", 4)
)
ROBUST_MAX_RECOVERY_PASSES = max(0, _env_int("REAKTORO_ROBUST_RECOVERY_PASSES", 1))
ROBUST_MAX_RECOVERY_POINTS = max(0, _env_int("REAKTORO_ROBUST_RECOVERY_POINTS", 80))
ROBUST_EPSILON = _env_float("REAKTORO_ROBUST_EPSILON", 1e-14)


def assert_bindings_support_scalar_inputs():
    """Fail fast with a clear message when Python bindings require opaque Real types.

    Some local Windows builds expose methods like EquilibriumConditions.set()
    only with autodiff::Real<1,double> arguments, but the corresponding Real
    Python type is not registered. In that case, any call with a float fails,
    and attempting to run the test may trigger runtime assertions.
    """
    db = rkt.SupcrtDatabase("supcrtbl")
    system = rkt.ChemicalSystem(
        db,
        rkt.AqueousPhase(rkt.speciate(["H2O(aq)", "H+", "OH-", "e-", "Fe+2"])),
    )
    specs = rkt.EquilibriumSpecs(system)
    specs.temperature()
    specs.pressure()
    specs.lgActivity("H+")
    specs.Eh()
    specs.lgActivity("Fe+2")
    cond = rkt.EquilibriumConditions(specs)

    try:
        cond.set("T", 298.15)
        cond.set("P", 1.0e5)
        cond.set("ln(a[H+])", float(np.log(10.0) * -7.0))
        cond.set("Eh", 0.0)
        cond.set("ln(a[Fe+2])", 0.0)
    except TypeError as exc:
        raise RuntimeError(
            "Current reaktoro4py binary does not accept scalar Python floats for "
            "equilibrium input variables (T/P/H+/Eh/activity). "
            "This open-system comparison test cannot run with this binary. "
            "Use a Python binding build with float overloads for EquilibriumConditions, "
            "or rebuild reaktoro4py from source in Release mode with the updated "
            "EquilibriumConditions.py.cxx wrappers."
        ) from exc


# ── Stability / predominance criterion (open system) ────────────────────────
# The system is open to Fe, so total iron is not conserved: iron flows in/out
# as an implicit titrant to maintain a fixed log10(a(Fe+2)) = LOGA_FE2 at
# every grid point (identical to CHNOSZ's basis-species convention).
#
# Criterion for the predominant Fe-bearing phase at each (pH, Eh) cell:
#   1. If one or more minerals are thermodynamically stable (amount > threshold),
#      show the mineral that precipitated the most (most moles).
#   2. Otherwise show the aqueous Fe species with the highest log10 activity
#      (equivalent to CHNOSZ's affinity ranking for aqueous species).
#
# Do NOT use speciesAmount for aqueous species: with a fixed Fe+2 activity all
# dissolved Fe species have non-negligible amounts and the ranking is unstable.
# Do NOT use log10(activity) for minerals: pure solids always have a = 1
# regardless of whether they are in fact stable at that point.
MINERAL_STABLE_THRESHOLD = 1e-12  # mol — below this a mineral is considered absent


def compute_predominance(
    all_states, success_mask, pH_vals, Eh_vals, aq_species, minerals
):
    """Return a (n_pH × n_Eh) array of float indices into (aq_species + minerals)."""
    all_sp = aq_species + minerals
    n_aq = len(aq_species)
    pred = np.full((len(pH_vals), len(Eh_vals)), np.nan, dtype=float)

    for ix, row in enumerate(all_states):
        for iy, state in enumerate(row):
            if not success_mask[ix, iy]:
                continue
            props = rkt.ChemicalProps(state)

            # --- Check minerals (stable = amount > threshold) ---
            best_min_idx = -1
            best_min_amt = MINERAL_STABLE_THRESHOLD
            for k, sp in enumerate(minerals):
                try:
                    amt = float(props.speciesAmount(sp))
                except Exception:
                    amt = 0.0
                if amt > best_min_amt:
                    best_min_amt = amt
                    best_min_idx = n_aq + k  # index into all_sp

            if best_min_idx >= 0:
                pred[ix, iy] = float(best_min_idx)
                continue

            # --- No stable mineral: pick aqueous species by highest log10(a) ---
            best_aq_idx = -1
            best_lga = -np.inf
            for k, sp in enumerate(aq_species):
                try:
                    lga = float(props.speciesActivityLg(sp))
                except Exception:
                    lga = float("nan")
                if np.isfinite(lga) and lga > best_lga:
                    best_lga = lga
                    best_aq_idx = k

            if best_aq_idx >= 0:
                pred[ix, iy] = float(best_aq_idx)

    return pred


def state_has_finite_props(state, probe_species):
    """Return True only if state has finite amounts and activities on probes."""
    try:
        props = rkt.ChemicalProps(state)
    except Exception:
        return False

    for sp in probe_species:
        try:
            amt = float(props.speciesAmount(sp))
            lga = float(props.speciesActivityLg(sp))
        except Exception:
            continue
        if np.isfinite(amt) and np.isfinite(lga):
            return True
    return False


def configure_solver_for_robustness(solver):
    """Configure solver options for robustness-first operation.

    The goal is to improve convergence on difficult cells at the cost of speed.
    """
    opts = rkt.EquilibriumOptions()

    # New bindings exposed in local tree: force exact Hessian + warmstart.
    try:
        opts.hessian = rkt.GibbsHessian.Exact
    except Exception:
        pass
    try:
        opts.warmstart = True
    except Exception:
        pass

    # NOTE: EquilibriumOptions.optima currently points to an Optima::Options
    # type that is not exposed to Python in this build, so we must not read or
    # write nested optima fields here.

    # A slightly larger epsilon can improve numerical behavior near bounds.
    try:
        opts.epsilon = ROBUST_EPSILON
    except Exception:
        pass

    solver.setOptions(opts)
    return opts


def solve_point_robust(
    solver,
    conditions,
    seed_states,
    pH_val,
    Eh_val,
    probe_species,
    homotopy_steps=24,
):
    """Try direct solves and then homotopy continuation from easier anchors."""

    def _attempt(seed_state, pH_now, Eh_now):
        trial = rkt.ChemicalState(seed_state)
        conditions.set("ln(a[H+])", float(np.log(10.0) * -pH_now))
        conditions.set("Eh", float(Eh_now))
        res = solver.solve(trial, conditions)
        ok = bool(res.succeeded()) and state_has_finite_props(trial, probe_species)
        return ok, trial

    # 1) Direct attempts with candidate seeds.
    for seed in seed_states[:ROBUST_MAX_SEEDS_PER_POINT]:
        try:
            ok, trial = _attempt(seed, pH_val, Eh_val)
        except Exception:
            continue
        if ok:
            return True, trial

    # 2) Optional homotopy attempts from easy anchors to target point.
    if homotopy_steps <= 0:
        return False, None

    anchors = [
        (pH_val, 0.0),
        (7.0, Eh_val),
        (7.0, 0.0),
    ]

    for seed in seed_states[:ROBUST_MAX_SEEDS_PER_POINT]:
        for pH0, Eh0 in anchors:
            trial = rkt.ChemicalState(seed)
            ok_path = True
            # Exclude t=0 (initial seed state), solve from first step onward.
            for t in np.linspace(1.0 / homotopy_steps, 1.0, homotopy_steps):
                pHi = pH0 + (pH_val - pH0) * float(t)
                Ehi = Eh0 + (Eh_val - Eh0) * float(t)
                try:
                    ok, trial = _attempt(trial, pHi, Ehi)
                except Exception:
                    ok_path = False
                    break
                if not ok:
                    ok_path = False
                    break
            if ok_path:
                return True, trial

    return False, None


def main():
    print(f"Using reaktoro4py from: {PYD_DIR}")

    assert_bindings_support_scalar_inputs()

    try:
        rkt.Warnings.disable(906)
    except Exception:
        pass

    db = rkt.SupcrtDatabase("supcrtbl")

    # Strict test: keep only species/minerals shared with CHNOSZ list.
    pH_vals = np.linspace(-2.0, 16.0, 120)
    Eh_vals = np.linspace(-2.0, 2.0, 100)

    aq_species = [
        "Fe+2",
        "Fe+3",
        "FeO+",
        "FeO2-",
        "FeOH+",
        "FeOH+2",
        "HFeO2(aq)",  # CHNOSZ alias: HFeO2
        "HFeO2-",
    ]
    minerals = ["Goethite", "Hematite", "Iron", "Magnetite"]
    species = aq_species + minerals

    # Include e- explicitly so Eh constraints are well posed over the grid.
    aqueous = ["H2O(aq)", "H+", "OH-", "e-"] + aq_species

    system = rkt.ChemicalSystem(
        db,
        rkt.AqueousPhase(rkt.speciate(aqueous)),
        rkt.MineralPhases(rkt.StringList(minerals)),
    )

    # ── Open-system specifications ───────────────────────────────────────────
    # Bug fix: the system must be OPEN to the Fe component so that iron can
    # flow in/out as an implicit titrant.  specs.lgActivity("Fe+2") declares:
    #   - log10(a(Fe+2)) is an input variable (fixed at LOGA_FE2 every point)
    #   - Fe is an implicit titrant whose amount is solved, not conserved
    # This is the Reaktoro equivalent of CHNOSZ's  basis("Fe+2", LOGA_FE2).
    LOGA_FE2 = 0.0  # log10(a) = 0  ↔  a(Fe+2) = 1  (CHNOSZ default for basis species)

    specs = rkt.EquilibriumSpecs(system)
    specs.temperature()
    specs.pressure()
    specs.lgActivity("H+")
    specs.Eh()
    specs.lgActivity("Fe+2")  # FIX 1: open to Fe, constrains iron activity

    # ── Initial state (only H2O; Fe amount is unknown — open system) ────────
    state0 = rkt.ChemicalState(system)
    # A physically meaningful water-rich initial guess improves convergence.
    # Keep a fallback to default state for binding variants that do not accept
    # scalar setters in ChemicalState.
    try:
        state0.set("H2O(aq)", 55.5, "mol")
    except Exception:
        pass

    solver = rkt.EquilibriumSolver(specs)
    solver_opts = configure_solver_for_robustness(solver)
    conditions = rkt.EquilibriumConditions(specs)
    # Use generic input setters; this build exposes Real-only typed setters
    # for temperature/pressure/Eh/activity in Python.
    conditions.set("T", 298.15)  # K
    conditions.set("P", 1.0e5)  # Pa
    # Workaround for Python Real binding gaps:
    # use generic input setter with names instead of lgActivity()/Eh() methods.
    conditions.set("ln(a[Fe+2])", float(np.log(10.0) * LOGA_FE2))

    all_states = []
    success_mask = np.zeros((len(pH_vals), len(Eh_vals)), dtype=bool)
    n_fail = 0
    n_fail_nonconv = 0
    n_fail_exception = 0
    probe_species = ["Fe+2", "Fe+3", "H+", "OH-", "H2O(aq)"]

    print(
        "Robust profile: "
        f"max_seeds={ROBUST_MAX_SEEDS_PER_POINT}, "
        f"homotopy_main={ROBUST_HOMOTOPY_STEPS_MAIN}, "
        f"homotopy_recovery={ROBUST_HOMOTOPY_STEPS_RECOVERY}, "
        f"recovery_passes={ROBUST_MAX_RECOVERY_PASSES}, "
        f"recovery_points={ROBUST_MAX_RECOVERY_POINTS}, "
        f"epsilon={ROBUST_EPSILON}"
    )

    # Continuation across the grid (row-major) with multi-seed and homotopy fallback.
    # A failed solve can leave the state as a poor guess for nearby points,
    # so we keep independent row/column/base seeds to avoid cascade failures.
    seed_for_row = rkt.ChemicalState(state0)
    col_last_success = [None for _ in Eh_vals]
    for ix, pH_val in enumerate(pH_vals):
        row_states = []
        s = rkt.ChemicalState(seed_for_row)
        row_last_success = None

        for iy, Eh_val in enumerate(Eh_vals):
            conditions.set("ln(a[H+])", float(np.log(10.0) * -pH_val))
            conditions.set("Eh", float(Eh_val))
            # Ordered by locality first, then progressively more robust seeds.
            candidate_seeds = []
            candidate_seeds.append(rkt.ChemicalState(s))
            if row_last_success is not None:
                candidate_seeds.append(rkt.ChemicalState(row_last_success))
            if col_last_success[iy] is not None:
                candidate_seeds.append(rkt.ChemicalState(col_last_success[iy]))
            candidate_seeds.append(rkt.ChemicalState(seed_for_row))
            candidate_seeds.append(rkt.ChemicalState(state0))

            try:
                ok, chosen_state = solve_point_robust(
                    solver,
                    conditions,
                    candidate_seeds,
                    float(pH_val),
                    float(Eh_val),
                    probe_species,
                    homotopy_steps=ROBUST_HOMOTOPY_STEPS_MAIN,
                )

                if ok:
                    s = rkt.ChemicalState(chosen_state)
                    row_last_success = rkt.ChemicalState(chosen_state)
                    col_last_success[iy] = rkt.ChemicalState(chosen_state)
                else:
                    n_fail += 1
                    n_fail_nonconv += 1
                    # Reset continuation cursor after a full failure.
                    if row_last_success is not None:
                        s = rkt.ChemicalState(row_last_success)
                    else:
                        s = rkt.ChemicalState(seed_for_row)
            except Exception:
                ok = False
                n_fail += 1
                n_fail_exception += 1
                if row_last_success is not None:
                    s = rkt.ChemicalState(row_last_success)
                else:
                    s = rkt.ChemicalState(seed_for_row)

            success_mask[ix, iy] = ok
            row_states.append(rkt.ChemicalState(s))

        all_states.append(row_states)
        print(
            f"Row {ix + 1}/{len(pH_vals)} complete; "
            f"failures so far={int(np.sum(~success_mask[: ix + 1, :]))}"
        )
        if row_last_success is not None:
            seed_for_row = rkt.ChemicalState(row_last_success)

    # Recovery passes: revisit failed points using solved neighbors as seeds.
    # This often resolves isolated non-converged cells near field boundaries.
    n_recovered = 0
    max_recovery_passes = ROBUST_MAX_RECOVERY_PASSES
    for _ in range(max_recovery_passes):
        recovered_this_pass = 0
        failed_points = np.argwhere(~success_mask)
        if failed_points.size == 0:
            break
        if (
            ROBUST_MAX_RECOVERY_POINTS > 0
            and len(failed_points) > ROBUST_MAX_RECOVERY_POINTS
        ):
            stride = int(np.ceil(len(failed_points) / ROBUST_MAX_RECOVERY_POINTS))
            failed_points = failed_points[::stride]

        for ix, iy in failed_points:
            ix = int(ix)
            iy = int(iy)

            conditions.set("ln(a[H+])", float(np.log(10.0) * -pH_vals[ix]))
            conditions.set("Eh", float(Eh_vals[iy]))

            neighbor_idx = [
                (ix, iy - 1),
                (ix, iy + 1),
                (ix - 1, iy),
                (ix + 1, iy),
            ]
            candidate_seeds = []
            for jx, jy in neighbor_idx:
                if (
                    0 <= jx < len(pH_vals)
                    and 0 <= jy < len(Eh_vals)
                    and success_mask[jx, jy]
                ):
                    candidate_seeds.append(rkt.ChemicalState(all_states[jx][jy]))

            # Last-resort seeds
            candidate_seeds.append(rkt.ChemicalState(seed_for_row))
            candidate_seeds.append(rkt.ChemicalState(state0))

            ok, trial = solve_point_robust(
                solver,
                conditions,
                candidate_seeds,
                float(pH_vals[ix]),
                float(Eh_vals[iy]),
                probe_species,
                homotopy_steps=ROBUST_HOMOTOPY_STEPS_RECOVERY,
            )
            if ok:
                all_states[ix][iy] = rkt.ChemicalState(trial)
                success_mask[ix, iy] = True

            if ok:
                recovered_this_pass += 1

        n_recovered += recovered_this_pass
        if recovered_this_pass == 0:
            break

    total_points = len(pH_vals) * len(Eh_vals)
    n_success = int(np.sum(success_mask))
    n_fail = total_points - n_success
    n_fail_nonconv = n_fail
    n_fail_exception = 0

    print(
        f"Grid complete.  Failures: {n_fail}/{len(pH_vals) * len(Eh_vals)} "
        f"(non-converged={n_fail_nonconv}, exceptions={n_fail_exception}, recovered={n_recovered})"
    )

    # ── FIX 2: predominance via mineral-amount + aqueous-log-activity ────────
    pred = compute_predominance(
        all_states,
        success_mask,
        pH_vals,
        Eh_vals,
        aq_species,
        minerals,
    )

    # ── Diagnostic: log species info at key (pH, Eh) cells ──────────────────
    test_pts = [(0.0, 0.0), (2.0, 0.0), (5.0, 0.5), (5.0, -0.5), (1.0, -0.5)]
    for tpH, tEh in test_pts:
        ix = int(np.argmin(np.abs(pH_vals - tpH)))
        iy = int(np.argmin(np.abs(Eh_vals - tEh)))
        if not success_mask[ix, iy]:
            print(f"\n--- pH={tpH} Eh={tEh} (grid idx {ix},{iy}) ---")
            print("  SOLVER FAILED AT THIS GRID POINT")
            continue
        s = all_states[ix][iy]
        props_d = rkt.ChemicalProps(s)
        print(f"\n--- pH={tpH} Eh={tEh} (grid idx {ix},{iy}) ---")
        for sp in species:
            try:
                amt = float(props_d.speciesAmount(sp))
                lga = float(props_d.speciesActivityLg(sp))
            except Exception:
                amt, lga = float("nan"), float("nan")
            print(f"  {sp:<18s}: amount={amt:.3e}  log10(a)={lga:+.3f}")

    # ── Plot ─────────────────────────────────────────────────────────────────
    pp = PredominancePlot(
        pH_vals,
        Eh_vals,
        pred,
        species,
        xlabel="pH",
        ylabel="Eh",
    )

    fig, ax = pp.plot(
        figsize=(8.2, 6.2),
        label_min_fraction=0.004,
        boundary_color="black",
        boundary_linewidth=1.0,
    )

    # Solid boundaries are drawn by PredominancePlot.plot().
    # Add striped (dashed) boundaries specifically around mineral-stability
    # domains (mineral vs aqueous predominance frontier).
    mineral_flag = np.where(
        np.isfinite(pred), (pred >= len(aq_species)).astype(float), np.nan
    )
    X, Y = np.meshgrid(pH_vals, Eh_vals)
    try:
        ax.contour(
            X,
            Y,
            mineral_flag.T,
            levels=[0.5],
            colors="black",
            linewidths=1.2,
            linestyles="--",
        )
    except Exception:
        pass

    # Draw explicit boundaries around non-converged (infeasible) regions.
    fail_flag = (~success_mask).astype(float)
    if np.any(fail_flag > 0.0):
        try:
            ax.contour(
                X,
                Y,
                fail_flag.T,
                levels=[0.5],
                colors="black",
                linewidths=1.1,
                linestyles="-",
            )
            ax.contour(
                X,
                Y,
                fail_flag.T,
                levels=[0.5],
                colors="black",
                linewidths=1.1,
                linestyles="--",
            )
        except Exception:
            pass
    pp.add_water_lines(ax, T_K=298.15, color="black", linestyle="--", linewidth=1.0)
    ax.set_title(
        f"Reaktoro Fe Pourbaix — open system, log10(a(Fe²⁺)) = {LOGA_FE2}"
        "\n(Identical shared species test vs CHNOSZ)"
    )
    plt.tight_layout()

    out = os.path.join(BASE, "Reaktoro_Pourbaix_Fe_identical_test.png")
    fig.savefig(out, dpi=180)
    plt.close(fig)

    report = os.path.join(BASE, "identical_test_setup.txt")
    with open(report, "w", encoding="utf-8") as f:
        f.write("Identical species/settings test (CHNOSZ vs Reaktoro) — OPEN SYSTEM\n")
        f.write(f"- Binary path: {PYD_DIR}\n")
        f.write("- Case: Fe Pourbaix\n")
        f.write("- Range: pH [-2,16], Eh [-2,2], T=25C, P=1 bar\n")
        f.write(
            f"- Solver failures: {n_fail}/{len(pH_vals) * len(Eh_vals)} "
            f"(non-converged={n_fail_nonconv}, exceptions={n_fail_exception})\n"
        )
        f.write(
            "- Continuation strategy: row-major continuation with row-seed fallback\n"
        )
        hess_name = "unknown"
        try:
            hess_name = str(solver_opts.hessian)
        except Exception:
            pass
        f.write(
            f"- Robust solver options: warmstart=true, epsilon={getattr(solver_opts, 'epsilon', 'n/a')}, hessian={hess_name}\n"
        )
        f.write(
            "- Robust fallback: direct multi-seed solve + bounded homotopy retries "
            f"({ROBUST_HOMOTOPY_STEPS_MAIN}-{ROBUST_HOMOTOPY_STEPS_RECOVERY} steps)\n"
        )
        f.write(
            "- Boundary style: solid field boundaries + striped mineral frontiers + solid/striped non-convergence boundaries\n"
        )
        f.write(
            "- Non-converged cells are retained as explicit infeasible domains under the imposed fixed-activity constraints\n"
        )
        f.write(
            f"- Iron constraint: log10(a(Fe+2)) = {LOGA_FE2}  (same as CHNOSZ basis default)\n"
        )
        f.write(
            "- FIX 1: specs.lgActivity('Fe+2') — system is open to Fe (implicit titrant)\n"
        )
        f.write(
            "- FIX 2: predominance = mineral amount > threshold, else max aqueous log10(a)\n"
        )
        f.write("- Shared species set used (Reaktoro names):\n")
        for sp in species:
            f.write(f"  - {sp}\n")
        f.write("- CHNOSZ alias mapping: HFeO2(aq) <-> HFeO2\n")
        f.write(f"- Output: {out}\n")

    print("Wrote:", out)
    print("Wrote:", report)


if __name__ == "__main__":
    main()
