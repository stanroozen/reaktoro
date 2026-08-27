import os
import sys
import time

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
BASE = os.path.abspath(os.path.dirname(__file__))

# Prefer Release build; fall back to package build dir.
_release_dir = os.path.join(REPO, "build", "Reaktoro", "Release")
_package_dir = os.path.join(
    REPO, "build", "python", "package", "build", "lib", "reaktoro"
)
PYD_DIR = None
for _cand in [_release_dir, _package_dir]:
    if os.path.isfile(os.path.join(_cand, "reaktoro4py.cp312-win_amd64.pyd")):
        PYD_DIR = _cand
        break
if PYD_DIR is None:
    raise FileNotFoundError(
        "reaktoro4py.cp312-win_amd64.pyd not found in build folders."
    )
if PYD_DIR not in sys.path:
    sys.path.insert(0, PYD_DIR)
os.add_dll_directory(PYD_DIR)

import reaktoro4py as rkt
import autodiff as ad

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
MosaicPlot = _dmod.MosaicPlot

# ── Shared constants ─────────────────────────────────────────────────────────
MINERAL_STABLE_THRESHOLD = 1e-12  # mol — below this a mineral is absent


def _env_int(name, default):
    try:
        return int(os.environ.get(name, str(default)))
    except Exception:
        return int(default)


def _save_checkpoint(checkpoint_file, pred, next_row, n_ok, n_fail, label):
    if not checkpoint_file:
        return
    np.savez_compressed(
        checkpoint_file,
        pred=pred,
        next_row=int(next_row),
        n_ok=int(n_ok),
        n_fail=int(n_fail),
        label=str(label),
    )


def _predominant(state, aq_species, minerals):
    """Return index into (aq_species + minerals) for the dominant Fe species."""
    all_sp = aq_species + minerals
    n_aq = len(aq_species)

    # Minerals first: largest amount wins if above threshold.
    best_min_idx = -1
    best_min_amt = MINERAL_STABLE_THRESHOLD
    for k, sp in enumerate(minerals):
        try:
            amt = float(state.speciesAmount(sp))
        except Exception:
            amt = 0.0
        if amt > best_min_amt:
            best_min_amt = amt
            best_min_idx = n_aq + k

    if best_min_idx >= 0:
        return best_min_idx

    # Aqueous: highest log10(activity).
    props = rkt.ChemicalProps(state)
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

    return best_aq_idx


def _fresh_state(system, seed_amounts):
    """Build a clean initial state with given (species, amount, unit) triples."""
    state = rkt.ChemicalState(system)
    for sp, amt, unit in seed_amounts:
        try:
            state.set(sp, ad.real(amt), unit)
        except Exception:
            pass
    return state


def _solve_grid_ph_eh(
    system,
    specs,
    solver,
    conds,
    seed_amounts,
    pH_vals,
    Eh_vals,
    aq_species,
    minerals,
    label,
    checkpoint_file=None,
    checkpoint_every_rows=None,
):
    """Per-point fresh-state pH/Eh grid solve. Returns predominance array."""
    nx, ny = len(pH_vals), len(Eh_vals)
    pred = np.full((nx, ny), -1, dtype=int)
    n_ok = n_fail = 0
    start_row = 0

    log_every_rows = max(1, _env_int("REAKTORO_LOG_EVERY_ROWS", 10))
    if checkpoint_every_rows is None:
        checkpoint_every_rows = max(0, _env_int("REAKTORO_CHECKPOINT_EVERY_ROWS", 10))

    if checkpoint_file and os.path.isfile(checkpoint_file):
        try:
            ck = np.load(checkpoint_file)
            pred_ck = ck["pred"]
            if pred_ck.shape == pred.shape:
                pred = pred_ck.astype(int, copy=False)
                start_row = int(ck["next_row"])
                n_ok = int(ck["n_ok"])
                n_fail = int(ck["n_fail"])
                print(
                    f"  {label}: resumed from checkpoint row {start_row}/{nx}  ok={n_ok}  fail={n_fail}",
                    flush=True,
                )
        except Exception:
            pass

    t0 = time.time()

    for i in range(start_row, nx):
        pH = pH_vals[i]
        for j, Eh in enumerate(Eh_vals):
            state = _fresh_state(system, seed_amounts)
            try:
                conds.set("ln(a[H+])", float(np.log(10.0) * -pH))
                conds.set("Eh", float(Eh))
                res = solver.solve(state, conds)
                if res.succeeded():
                    idx = _predominant(state, aq_species, minerals)
                    if idx >= 0:
                        pred[i, j] = idx
                    n_ok += 1
                else:
                    n_fail += 1
            except Exception:
                n_fail += 1

        if checkpoint_every_rows > 0 and (
            ((i + 1) % checkpoint_every_rows == 0) or i == nx - 1
        ):
            _save_checkpoint(checkpoint_file, pred, i + 1, n_ok, n_fail, label)

        if ((i + 1) % log_every_rows == 0) or i == nx - 1:
            elapsed = time.time() - t0
            print(
                f"  {label}: row {i + 1}/{nx}  ok={n_ok}  fail={n_fail}  elapsed={elapsed:.1f}s",
                flush=True,
            )

    _save_checkpoint(checkpoint_file, pred, nx, n_ok, n_fail, label)
    print(f"  {label}: done  ok={n_ok}  fail={n_fail}", flush=True)
    return pred


def save_reaktoro_pourbaix(db):
    n_pH = max(20, _env_int("REAKTORO_POURBAIX_NPH", 120))
    n_Eh = max(20, _env_int("REAKTORO_POURBAIX_NEH", 100))
    pH_vals = np.linspace(-2.0, 16.0, n_pH)
    Eh_vals = np.linspace(-2.0, 2.0, n_Eh)

    aq_fe_oh = [
        "Fe+2",
        "Fe+3",
        "FeO(aq)",
        "FeO+",
        "FeO2-",
        "FeOH+",
        "FeOH+2",
        "HFeO2(aq)",
        "HFeO2-",
    ]
    min_fe_oh = ["Ferropericlase", "Goethite", "Hematite", "Iron", "Magnetite"]
    aqueous = ["H2O(aq)", "H+", "OH-", "e-"] + aq_fe_oh

    system = rkt.ChemicalSystem(
        db,
        rkt.AqueousPhase(rkt.speciate(aqueous)),
        rkt.MineralPhases(rkt.StringList(min_fe_oh)),
    )

    specs = rkt.EquilibriumSpecs(system)
    specs.temperature()
    specs.pressure()
    specs.lgActivity("H+")  # fixes pH via ln(a[H+])
    specs.Eh()

    conds = rkt.EquilibriumConditions(specs)
    conds.set("T", 298.15)
    conds.set("P", 1.0e5)

    solver = rkt.EquilibriumSolver(specs)

    seed_amounts = [
        ("H2O(aq)", 1.0, "kg"),
        ("Fe+2", 1e-6, "mol"),
    ]

    species = aq_fe_oh + min_fe_oh
    print(
        f"Pourbaix grid: {len(pH_vals)} x {len(Eh_vals)} = {len(pH_vals) * len(Eh_vals)} points"
    )
    ckpt_file = os.path.join(BASE, "Reaktoro_Pourbaix_Fe_CHNOSZmatched.checkpoint.npz")
    pred = _solve_grid_ph_eh(
        system,
        specs,
        solver,
        conds,
        seed_amounts,
        pH_vals,
        Eh_vals,
        aq_fe_oh,
        min_fe_oh,
        "Pourbaix",
        checkpoint_file=ckpt_file,
    )

    pp = PredominancePlot(
        pH_vals,
        Eh_vals,
        pred,
        species,
        xlabel="pH",
        ylabel="Eh",
        palette="tab20",
    )
    fig, ax = pp.plot(
        figsize=(8.2, 6.2),
        label_min_fraction=0.004,
        boundary_color="black",
        boundary_linewidth=1.0,
    )
    # Dashed frontier between mineral-dominant and aqueous-dominant cells.
    mineral_flag = np.where(pred >= 0, (pred >= len(aq_fe_oh)).astype(float), np.nan)
    if np.any(np.isfinite(mineral_flag)):
        X, Y = np.meshgrid(pH_vals, Eh_vals)
        try:
            ax.contour(
                X,
                Y,
                mineral_flag.T,
                levels=[0.5],
                colors="0.35",
                linewidths=1.2,
                linestyles="--",
            )
        except Exception:
            pass

    # Water-stability lines and neutral/reference guides.
    pp.add_water_lines(ax, T_K=298.15, color="black", linestyle="-.", linewidth=1.0)
    ax.axhline(0.0, color="0.55", linestyle=":", linewidth=1.0)
    ax.axvline(7.0, color="0.55", linestyle=":", linewidth=1.0)

    legend_handles = [
        Line2D(
            [0],
            [0],
            color="black",
            lw=1.4,
            linestyle="-",
            label="Predominance boundary (solid)",
        ),
        Line2D(
            [0],
            [0],
            color="0.35",
            lw=1.2,
            linestyle="--",
            label="Mineral-stability frontier (dashed)",
        ),
        Line2D(
            [0],
            [0],
            color="black",
            lw=1.0,
            linestyle="-.",
            label="Water-stability lines (dash-dot)",
        ),
        Line2D(
            [0],
            [0],
            color="0.55",
            lw=1.0,
            linestyle=":",
            label="Reference guides: Eh = 0, pH = 7",
        ),
    ]
    ax.legend(handles=legend_handles, loc="lower right", fontsize=8, framealpha=0.9)

    ax.set_title("Reaktoro Fe-O-H Pourbaix (CHNOSZ-matched ranges)\nT=25°C, P=1 bar")
    plt.tight_layout()

    out = os.path.join(BASE, "Reaktoro_Pourbaix_Fe_CHNOSZmatched.png")
    fig.savefig(out, dpi=180)
    plt.close(fig)
    return out


def save_reaktoro_mosaic(db):
    n_pH = max(20, _env_int("REAKTORO_MOSAIC_NPH", 100))
    n_Eh = max(20, _env_int("REAKTORO_MOSAIC_NEH", 80))
    pH_vals = np.linspace(0.0, 14.0, n_pH)
    Eh_vals = np.linspace(-1.0, 1.0, n_Eh)

    fe_aq = ["Fe+2", "Fe+3", "HFeO2-"]
    min_fe = ["Pyrite", "Pyrrhotite,trot", "Siderite", "Hematite", "Magnetite"]

    aqueous_full = [
        "H2O(aq)",
        "H+",
        "OH-",
        "e-",
        "Fe+2",
        "Fe+3",
        "HFeO2-",
        "SO4-2",
        "HSO4-",
        "HS-",
        "H2S(aq)",
        "CO3-2",
        "HCO3-",
        "CO2(aq)",
    ]

    system = rkt.ChemicalSystem(
        db,
        rkt.AqueousPhase(rkt.speciate(aqueous_full)),
        rkt.MineralPhases(rkt.StringList(min_fe)),
    )

    specs = rkt.EquilibriumSpecs(system)
    specs.temperature()
    specs.pressure()
    specs.lgActivity("H+")
    specs.Eh()

    conds = rkt.EquilibriumConditions(specs)
    conds.set("T", 298.15)
    conds.set("P", 1.0e5)

    solver = rkt.EquilibriumSolver(specs)

    seed_amounts = [
        ("H2O(aq)", 1.0, "kg"),
        ("Fe+2", 1e-6, "mol"),
        ("SO4-2", 1e-6, "mol"),
        ("CO3-2", 1.0, "mol"),
    ]

    print(
        f"Mosaic grid: {len(pH_vals)} x {len(Eh_vals)} = {len(pH_vals) * len(Eh_vals)} points"
    )
    ckpt_min = os.path.join(BASE, "Reaktoro_Mosaic_Min_CHNOSZmatched.checkpoint.npz")
    ckpt_aq = os.path.join(BASE, "Reaktoro_Mosaic_Aq_CHNOSZmatched.checkpoint.npz")

    # Mineral predominance
    min_pred = _solve_grid_ph_eh(
        system,
        specs,
        solver,
        conds,
        seed_amounts,
        pH_vals,
        Eh_vals,
        [],
        min_fe,
        "Mosaic-min",
        checkpoint_file=ckpt_min,
    )
    # Aqueous predominance (re-run with same fresh states)
    aq_pred = _solve_grid_ph_eh(
        system,
        specs,
        solver,
        conds,
        seed_amounts,
        pH_vals,
        Eh_vals,
        fe_aq,
        [],
        "Mosaic-aq",
        checkpoint_file=ckpt_aq,
    )

    layers = [
        {
            "species": min_fe,
            "predominance": min_pred,
            "palette": "Pastel1",
            "alpha": 1.0,
        },
        {"species": fe_aq, "predominance": aq_pred, "palette": "tab10", "alpha": 0.62},
    ]

    mp = MosaicPlot(pH_vals, Eh_vals, layers, xlabel="pH", ylabel="Eh")
    fig, ax = mp.plot(figsize=(8.2, 6.2), label_min_fraction=0.005)

    # Dashed line marks the mineral-layer active frontier.
    mineral_active = np.where(min_pred >= 0, 1.0, 0.0)
    X, Y = np.meshgrid(pH_vals, Eh_vals)
    try:
        ax.contour(
            X,
            Y,
            mineral_active.T,
            levels=[0.5],
            colors="0.35",
            linewidths=1.2,
            linestyles="--",
        )
    except Exception:
        pass

    mp.add_water_lines(ax, T_K=298.15, color="black", linestyle="-.", linewidth=1.0)
    ax.axhline(0.0, color="0.55", linestyle=":", linewidth=1.0)
    ax.axvline(7.0, color="0.55", linestyle=":", linewidth=1.0)

    legend_handles = [
        Line2D(
            [0],
            [0],
            color="black",
            lw=1.4,
            linestyle="-",
            label="Predominance boundary (solid)",
        ),
        Line2D(
            [0],
            [0],
            color="0.35",
            lw=1.2,
            linestyle="--",
            label="Mineral-stability frontier (dashed)",
        ),
        Line2D(
            [0],
            [0],
            color="black",
            lw=1.0,
            linestyle="-.",
            label="Water-stability lines (dash-dot)",
        ),
        Line2D(
            [0],
            [0],
            color="0.55",
            lw=1.0,
            linestyle=":",
            label="Reference guides: Eh = 0, pH = 7",
        ),
    ]
    ax.legend(handles=legend_handles, loc="lower right", fontsize=8, framealpha=0.9)

    ax.set_title("Reaktoro Fe-S-C-O-H Mosaic (CHNOSZ-matched ranges)\nT=25°C, P=1 bar")
    plt.tight_layout()

    out = os.path.join(BASE, "Reaktoro_Mosaic_Fe_CHNOSZmatched.png")
    fig.savefig(out, dpi=180)
    plt.close(fig)
    return out


def main():
    print(f"Using reaktoro4py from: {PYD_DIR}")
    try:
        rkt.Warnings.disable(906)
    except Exception:
        pass

    db = rkt.SupcrtDatabase("supcrtbl")

    pourbaix_nph = max(20, _env_int("REAKTORO_POURBAIX_NPH", 120))
    pourbaix_neh = max(20, _env_int("REAKTORO_POURBAIX_NEH", 100))
    mosaic_nph = max(20, _env_int("REAKTORO_MOSAIC_NPH", 100))
    mosaic_neh = max(20, _env_int("REAKTORO_MOSAIC_NEH", 80))

    out1 = save_reaktoro_pourbaix(db)
    out2 = save_reaktoro_mosaic(db)

    report = os.path.join(BASE, "Reaktoro_CHNOSZmatched_setup.txt")
    with open(report, "w", encoding="utf-8") as f:
        f.write(
            "Generated CHNOSZ-matched Reaktoro diagrams (per-point fresh-state solver)\n"
        )
        f.write(f"- Binary path: {PYD_DIR}\n")
        f.write(
            f"- Pourbaix: pH [-2,16] x Eh [-2,2], {pourbaix_nph}x{pourbaix_neh} points, T=25C, P=1 bar\n"
        )
        f.write(
            f"- Mosaic:   pH [0,14]  x Eh [-1,1], {mosaic_nph}x{mosaic_neh} points, T=25C, P=1 bar\n"
        )
        f.write("- Method: per-point EquilibriumSolver, fresh state each point\n")
        f.write(
            "- Boundary semantics: solid=predominance, dashed=minerals-vs-aqueous frontier, dash-dot=water stability, dotted=Eh0/pH7 guides\n"
        )
        f.write(f"- Output: {out1}\n")
        f.write(f"- Output: {out2}\n")

    print("Wrote:", out1)
    print("Wrote:", out2)
    print("Wrote:", report)


if __name__ == "__main__":
    main()
