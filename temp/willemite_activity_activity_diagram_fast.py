import importlib.util
import os
import sys
import optima
import gc

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REL = os.path.join(REPO, "build", "Reaktoro", "Release")
if REL not in sys.path:
    sys.path.insert(0, REL)

TUTORIAL = os.path.join(
    REPO,
    "DEW_Experimental_Benchmark",
    "Tutorial",
    "willemite_solubility_tutorial_dew17hp622_zn.py",
)
spec = importlib.util.spec_from_file_location("w", TUTORIAL)
w = importlib.util.module_from_spec(spec)
spec.loader.exec_module(w)

try:
    w.Warnings.disable(906)
except Exception:
    pass

w.USE_COMPETING_ZN_MINERALS = True
w.INCLUDE_COMPATIBLE_GANGUE_MINERALS = False
# Keep Zincite in the competing Zn mineral set and buffer with breccia gangue.
w.COMPETING_ZN_MINERALS = list(
    dict.fromkeys(w.COMPETING_ZN_MINERALS + ["hem", "ank", "q", "cc", "mag"])
)

extra = ["CO2,aq", "Fe+2", "Fe+3", "Ca+2", "Mg+2", "CaCO3,aq", "MgCO3,aq"]
w.AQUEOUS_SPECIES = list(dict.fromkeys(list(w.AQUEOUS_SPECIES) + extra))
# Favor low-sulfide hydrothermal chemistry for Willemite stability mapping.
w.AQUEOUS_SPECIES = [
    sp
    for sp in w.AQUEOUS_SPECIES
    if ("HS" not in sp and "SO" not in sp and "S2" not in sp)
]

# Re-seed minerals/aqueous species using the provided dolomitic breccia bulk trend
# (carbonate-rich, Zn-bearing, low sulfide) rather than equal 2 mol phase loads.
for m in w.COMPETING_ZN_MINERALS:
    w.INITIAL_SPECIES_AMOUNTS_MOL[m] = 1.0e-12

for m in ["hem", "ank", "q", "cc", "mag"]:
    w.INITIAL_SPECIES_AMOUNTS_MOL[m] = float(
        os.environ.get("WLM_GANGUE_SEED_MOL", "1.0e-6")
    )

# Keep tiny seeds for key Zn minerals to help solver nucleate alternatives.
w.INITIAL_SPECIES_AMOUNTS_MOL["Wlm"] = 1.0e-9
w.INITIAL_SPECIES_AMOUNTS_MOL["Znc"] = 1.0e-9

# Approximate 100 g bulk-rock mole trends (relative, not strict mass closure).
w.INITIAL_SPECIES_AMOUNTS_MOL["Zn2+"] = float(os.environ.get("WLM_ZN2_MOL", "2.29e-2"))
w.INITIAL_SPECIES_AMOUNTS_MOL["SiO2,aq"] = float(
    os.environ.get("WLM_SIO2_MOL", "3.24e-2")
)
w.INITIAL_SPECIES_AMOUNTS_MOL["Ca+2"] = float(os.environ.get("WLM_CA2_MOL", "5.08e-1"))
w.INITIAL_SPECIES_AMOUNTS_MOL["Mg+2"] = float(os.environ.get("WLM_MG2_MOL", "4.89e-1"))
w.INITIAL_SPECIES_AMOUNTS_MOL["Fe+3"] = float(os.environ.get("WLM_FE3_MOL", "3.5e-2"))
w.INITIAL_SPECIES_AMOUNTS_MOL["HCO3-"] = float(os.environ.get("WLM_HCO3_MOL", "1.04"))
w.INITIAL_SPECIES_AMOUNTS_MOL["CO3-2"] = float(os.environ.get("WLM_CO3_MOL", "2.0e-2"))
w.validate_user_inputs()

dew = w.Database.fromFile(w.PERPLEX_DATABASE_FILE)
sup = w.Database.fromFile(
    os.path.join(REPO, "embedded", "databases", "reaktoro", "supcrt07.yaml")
)
dew.addSpecies(sup.species("H2O(g)"))
dew.addSpecies(sup.species("CO2(g)"))

params = w.ActivityModelParamsPerplexDEW()
params.errorOnConflictingStandardState = False
aq = w.AqueousPhase(" ".join(w.AQUEOUS_SPECIES))
aq.setActivityModel(w.ActivityModelPerplexDEW(params))
mins = w.make_mineral_phases()
gas = w.GaseousPhase("H2O(g) CO2(g) O2")
gas.setActivityModel(w.ActivityModelPerplexGFSM(w.ActivityModelParamsPerplexGFSM()))

system = w.ChemicalSystem(dew, aq, mins, gas)
specs = w.EquilibriumSpecs(system)
specs.temperature()
specs.pressure()
specs.pH()
specs.fugacity("O2")
opts = w.EquilibriumOptions()
if hasattr(opts, "hessian") and hasattr(w, "GibbsHessian"):
    opts.hessian = w.GibbsHessian.Exact
if hasattr(opts, "warmstart"):
    opts.warmstart = True

# Ensure the iteration limit is applied through the nested Optima options.
MAX_ITERS = int(os.environ.get("WLM_MAX_ITERS", "1000"))
try:
    opt = optima.Options()
    opt.maxiters = MAX_ITERS
    opts.optima = opt
except Exception as exc:
    raise RuntimeError(
        "Failed to configure solver max iterations via opts.optima.maxiters"
    ) from exc


solver = w.make_equilibrium_solver(system, specs)
solver.setOptions(opts)

# High-effort rescue solver for stubborn points that fail under normal settings.
opts_rescue = w.EquilibriumOptions()
if hasattr(opts_rescue, "hessian") and hasattr(w, "GibbsHessian"):
    opts_rescue.hessian = w.GibbsHessian.Exact
if hasattr(opts_rescue, "warmstart"):
    opts_rescue.warmstart = True
try:
    opt_rescue = optima.Options()
    opt_rescue.maxiters = int(os.environ.get("WLM_RESCUE_MAX_ITERS", "5000"))
    opts_rescue.optima = opt_rescue
except Exception:
    pass

solver_rescue = w.make_equilibrium_solver(system, specs)
solver_rescue.setOptions(opts_rescue)


def make_seeded_state():
    st = w.make_base_state(system)
    st.set("H2O(g)", 1.0 - XCO2, "mol")
    st.set("CO2(g)", XCO2, "mol")
    st.set("CO2,aq", 1e-6 * XCO2, "mol")
    st.set("O2", 1e-20, "mol")
    return st


def solve_at(state, pH, logfO2, eqsolver):
    # Reuse solver instances so warmstart can carry information across retries.
    conds = w.EquilibriumConditions(specs)
    conds.temperature(T_C, "celsius")
    conds.pressure(P_BAR, "bar")
    conds.pH(float(pH))
    conds.fugacity("O2", 10.0 ** float(logfO2), "bar")
    return eqsolver.solve(state, conds)


# Broaden the grid around the stable anchor so the diagram includes the
# Willemite field boundary, not just its interior.
T_C = float(os.environ.get("WLM_T_C", "300.0"))
P_BAR = float(os.environ.get("WLM_P_BAR", "2000.0"))
XCO2 = float(os.environ.get("WLM_XCO2", "0.3"))
TH = 1e-8

# Hydrothermal Willemite window: neutral-basic, oxidized, low sulfide.
pH_min = float(os.environ.get("WLM_PH_MIN", "6.5"))
pH_max = float(os.environ.get("WLM_PH_MAX", "11.5"))
pH_n = int(os.environ.get("WLM_PH_N", "7"))
fo2_min = float(os.environ.get("WLM_LOGFO2_MIN", "-20.0"))
fo2_max = float(os.environ.get("WLM_LOGFO2_MAX", "-6.0"))
fo2_n = int(os.environ.get("WLM_LOGFO2_N", "7"))
pH_values = np.linspace(pH_min, pH_max, pH_n)
logfO2_values = np.linspace(fo2_min, fo2_max, fo2_n)
logaH_values = -pH_values

MINERAL_NAMES = list(dict.fromkeys(w.selected_mineral_names()))

field = np.full((len(logaH_values), len(logfO2_values)), -1, dtype=int)
wlm_amount = np.full((len(logaH_values), len(logfO2_values)), np.nan)
solved_states = [
    [None for _ in range(len(logfO2_values))] for _ in range(len(pH_values))
]

assemblage_labels = []
assemblage_to_id = {}


def classify_mineral_assemblage(state):
    stable = []
    for mineral in MINERAL_NAMES:
        try:
            amount = float(state.speciesAmount(mineral))
        except Exception:
            amount = 0.0
        if amount > TH:
            stable.append(mineral)

    return "+".join(stable) if stable else "(no mineral)"


def assemblage_id(label):
    if label not in assemblage_to_id:
        assemblage_to_id[label] = len(assemblage_labels)
        assemblage_labels.append(label)
    return assemblage_to_id[label]


converged_count = 0
failed_count = 0
current_retry_round = 0
MAX_SWEEPS = int(os.environ.get("WLM_MAX_SWEEPS", "6"))
ENABLE_GLOBAL_CONTINUATION = os.environ.get(
    "WLM_ENABLE_GLOBAL_CONTINUATION", "1"
).strip().lower() not in ("0", "false", "no")


def step_indices(start, stop):
    if start == stop:
        return []
    step = 1 if stop > start else -1
    return range(start + step, stop + step, step)


def iter_seed_states(i, j):
    seeds = [make_seeded_state()]
    seen = set()

    def push(ii, jj):
        if not (0 <= ii < len(pH_values) and 0 <= jj < len(logfO2_values)):
            return
        if (ii, jj) in seen:
            return
        seed_state = solved_states[ii][jj]
        if seed_state is not None:
            seeds.append(seed_state.clone())
            seen.add((ii, jj))

    neighbor_coords = [
        (i, j - 1),
        (i - 1, j),
        (i, j + 1),
        (i + 1, j),
        (i - 1, j - 1),
        (i - 1, j + 1),
        (i + 1, j - 1),
        (i + 1, j + 1),
    ]
    for ii, jj in neighbor_coords:
        push(ii, jj)

    # Expand to a small radius-2 neighborhood when immediate neighbors are missing.
    for di in range(-2, 3):
        for dj in range(-2, 3):
            if abs(di) <= 1 and abs(dj) <= 1:
                continue
            push(i + di, j + dj)

    # Final fallback: pull a few globally nearest solved states.
    if len(seeds) <= 2:
        candidates = []
        for ii in range(len(pH_values)):
            for jj in range(len(logfO2_values)):
                if solved_states[ii][jj] is not None and (ii, jj) not in seen:
                    dist = abs(ii - i) + abs(jj - j)
                    candidates.append((dist, ii, jj))
        candidates.sort(key=lambda t: t[0])
        for _, ii, jj in candidates[:8]:
            push(ii, jj)

    return seeds


def nearest_solved_indices(i, j, limit=8):
    candidates = []
    for ii in range(len(pH_values)):
        for jj in range(len(logfO2_values)):
            if solved_states[ii][jj] is not None:
                dist = abs(ii - i) + abs(jj - j)
                candidates.append((dist, ii, jj))
    candidates.sort(key=lambda t: t[0])
    return candidates[:limit]


def try_continuation_from_anchor(i, j, ai, aj, eqsolver):
    anchor = solved_states[ai][aj]
    if anchor is None:
        return False

    # Try two path orders to reduce path-dependence: fO2-then-pH and pH-then-fO2.
    for f_first in (True, False):
        st = anchor.clone()
        ci, cj = ai, aj
        ok = True

        if f_first:
            for jj in step_indices(cj, j):
                r = solve_at(
                    st, float(pH_values[ci]), float(logfO2_values[jj]), eqsolver
                )
                if not r.succeeded():
                    ok = False
                    break
                cj = jj
            if ok:
                for ii in step_indices(ci, i):
                    r = solve_at(
                        st, float(pH_values[ii]), float(logfO2_values[cj]), eqsolver
                    )
                    if not r.succeeded():
                        ok = False
                        break
                    ci = ii
        else:
            for ii in step_indices(ci, i):
                r = solve_at(
                    st, float(pH_values[ii]), float(logfO2_values[cj]), eqsolver
                )
                if not r.succeeded():
                    ok = False
                    break
                ci = ii
            if ok:
                for jj in step_indices(cj, j):
                    r = solve_at(
                        st, float(pH_values[ci]), float(logfO2_values[jj]), eqsolver
                    )
                    if not r.succeeded():
                        ok = False
                        break
                    cj = jj

        if ok and ci == i and cj == j:
            wlm = float(st.speciesAmount("Wlm"))
            wlm_amount[i, j] = wlm
            field[i, j] = assemblage_id(classify_mineral_assemblage(st))
            solved_states[i][j] = st.clone()
            return True

    return False


def solve_point(i, j):
    pH = float(pH_values[i])
    lf = float(logfO2_values[j])
    for seed in iter_seed_states(i, j):
        for eqsolver in (solver, solver_rescue):
            try:
                st = seed.clone()
                r = solve_at(st, pH, lf, eqsolver)
            except Exception as exc:
                print(
                    f"  Point ({pH:.2f}, {lf:.1f}): error - {type(exc).__name__}",
                    file=sys.stderr,
                )
                gc.collect()
                continue

            if r.succeeded():
                wlm = float(st.speciesAmount("Wlm"))
                wlm_amount[i, j] = wlm
                field[i, j] = assemblage_id(classify_mineral_assemblage(st))
                solved_states[i][j] = st.clone()
                return True

            gc.collect()

    # The global continuation fallback is expensive; use it only in later sweeps
    # when local warm-start retries have already filled most easy cells.
    if ENABLE_GLOBAL_CONTINUATION and current_retry_round >= 2:
        for _, ai, aj in nearest_solved_indices(i, j, limit=6):
            for eqsolver in (solver, solver_rescue):
                if try_continuation_from_anchor(i, j, ai, aj, eqsolver):
                    return True

    return False


for retry_round in range(MAX_SWEEPS):
    current_retry_round = retry_round
    round_failed = 0
    solved_before = int((field >= 0).sum())
    print(f"Retry sweep {retry_round + 1}/{MAX_SWEEPS}", file=sys.stderr)
    for i, pH in enumerate(pH_values):
        x_indices = (
            range(len(logfO2_values))
            if (i % 2 == 0)
            else range(len(logfO2_values) - 1, -1, -1)
        )
        for j in x_indices:
            if field[i, j] >= 0:
                continue
            if solve_point(i, j):
                pass
            else:
                round_failed += 1

    converged_count = int((field >= 0).sum())
    failed_count = int((field < 0).sum())
    solved_after = int((field >= 0).sum())
    print(
        f"  after sweep {retry_round + 1}: converged={converged_count}, failed={failed_count}",
        file=sys.stderr,
    )
    if round_failed == 0 or solved_after == solved_before:
        break

# Write CSV AFTER loop completes (outside loop structure)
csv_path = os.path.join(REPO, "temp", "willemite_activity_activity_field_slow.csv")
with open(csv_path, "w", encoding="utf-8") as f:
    f.write("pH,log_a_Hplus,logfO2,converged,wlm_mol,assemblage_id,assemblage\n")
    for i, pH in enumerate(pH_values):
        for j, lf in enumerate(logfO2_values):
            conv = int(field[i, j] >= 0)
            asm_id = int(field[i, j]) if conv else -1
            asm_lbl = assemblage_labels[asm_id] if conv else ""
            wlm = wlm_amount[i, j]
            wtxt = "nan" if np.isnan(wlm) else f"{wlm:.16e}"
            f.write(f"{pH:.4f},{-pH:.4f},{lf:.4f},{conv},{wtxt},{asm_id},{asm_lbl}\n")

print(
    f"Grid compute complete: converged={converged_count}, failed={failed_count}",
    file=sys.stderr,
)

fig, ax = plt.subplots(figsize=(8.2, 6.0))
extent = [
    logfO2_values.min(),
    logfO2_values.max(),
    logaH_values.min(),
    logaH_values.max(),
]

categorical = np.ma.masked_where(field < 0, field)
num_assemblages = max(len(assemblage_labels), 1)
cmap = plt.get_cmap("tab20", num_assemblages)
ax.imshow(
    categorical,
    origin="lower",
    extent=extent,
    aspect="auto",
    interpolation="nearest",
    cmap=cmap,
    vmin=0,
    vmax=max(num_assemblages - 1, 0),
    alpha=0.95,
)

nc_i, nc_j = np.where(field == -1)
if len(nc_i) > 0:
    ax.scatter(
        logfO2_values[nc_j],
        logaH_values[nc_i],
        marker="x",
        c="#d62728",
        s=14,
        linewidths=0.7,
    )

ax.set_xlabel("log10 fO2 (bar)")
ax.set_ylabel("log10 a(H+)")
ax.set_title("Mineral Assemblage Stability Fields (Zincite Included)")
ax.grid(True, alpha=0.2)

from matplotlib.patches import Patch

legend_handles = []
for k, label in enumerate(assemblage_labels):
    legend_handles.append(Patch(facecolor=cmap(k), alpha=0.9, label=label))
if np.any(field < 0):
    legend_handles.append(
        Patch(facecolor="#ffffff", alpha=0.0, label="x: not converged")
    )

ax.legend(
    handles=legend_handles,
    loc="best",
)

png_path = os.path.join(REPO, "temp", "willemite_activity_activity_diagram_slow.png")
fig.tight_layout()
fig.savefig(png_path, dpi=260)
plt.close(fig)
gc.collect()

print("MODEL=PerplexDEW+PerplexGFSM (Znc included, ank+cc+mag+hem+q, low sulfide)")
print(f"ANCHOR=T={T_C}C P={P_BAR}bar XCO2={XCO2}")
print(
    f"COUNTS total={field.size} converged={(field >= 0).sum()} assemblages={len(assemblage_labels)}"
)
for idx, label in enumerate(assemblage_labels):
    count = int((field == idx).sum())
    print(f"ASSEMBLAGE id={idx} count={count} label={label}")
print(f"CSV={csv_path}")
print(f"PNG={png_path}")
