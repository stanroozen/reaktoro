import importlib.util
import os
import sys

import matplotlib.pyplot as plt
import numpy as np

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

# Build system with zincite suppressed to expose possible willemite field.
w.USE_COMPETING_ZN_MINERALS = True
w.INCLUDE_COMPATIBLE_GANGUE_MINERALS = False
w.COMPETING_ZN_MINERALS = [m for m in w.COMPETING_ZN_MINERALS if m != "Znc"]
w.COMPETING_ZN_MINERALS = list(
    dict.fromkeys(w.COMPETING_ZN_MINERALS + ["hem", "cc", "q", "mag"])
)

extra = [
    "CO2,aq",
    "Fe+2",
    "Fe+3",
    "Ca+2",
    "Mg+2",
    "CaCO3,aq",
    "MgCO3,aq",
]
w.AQUEOUS_SPECIES = list(dict.fromkeys(list(w.AQUEOUS_SPECIES) + extra))

for m in ["hem", "cc", "q", "mag"]:
    w.INITIAL_SPECIES_AMOUNTS_MOL[m] = 2.0
w.INITIAL_SPECIES_AMOUNTS_MOL["Znc"] = 0.0

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

solver = w.make_equilibrium_solver(system, specs)
opts = w.EquilibriumOptions()
if hasattr(opts, "hessian") and hasattr(w, "GibbsHessian"):
    opts.hessian = w.GibbsHessian.Exact
if hasattr(opts, "warmstart"):
    opts.warmstart = True
solver.setOptions(opts)
conditions = w.EquilibriumConditions(specs)

# Fixed anchor conditions for activity-activity map.
T_C = 300.0
P_BAR = 2000.0
XCO2 = 0.3
TH = 1.0e-8

# Activity-activity axes: log a(H+) versus log fO2.
pH_values = np.linspace(2.5, 6.5, 25)
logfO2_values = np.linspace(-30.0, -18.0, 29)
logaH_values = -pH_values

# Grid values:
# -1: no convergence
#  0: converged, willemite absent
#  1: converged, willemite present
field = np.full((len(logaH_values), len(logfO2_values)), -1, dtype=int)
wlm_amount = np.full((len(logaH_values), len(logfO2_values)), np.nan, dtype=float)

for i, pH in enumerate(pH_values):
    state = w.make_base_state(system)
    state.set("H2O(g)", 1.0 - XCO2, "mol")
    state.set("CO2(g)", XCO2, "mol")
    state.set("CO2,aq", 1.0e-6 * XCO2, "mol")
    state.set("O2", 1.0e-20, "mol")

    for j, logfO2 in enumerate(logfO2_values):
        conditions.temperature(T_C, "celsius")
        conditions.pressure(P_BAR, "bar")
        conditions.pH(float(pH))
        conditions.fugacity("O2", 10.0 ** float(logfO2), "bar")

        res = solver.solve(state, conditions)
        if not res.succeeded():
            continue

        wlm = float(state.speciesAmount("Wlm"))
        wlm_amount[i, j] = wlm
        field[i, j] = 1 if wlm > TH else 0

# Save raw data table.
out_csv = os.path.join(REPO, "temp", "willemite_activity_activity_field.csv")
with open(out_csv, "w", encoding="utf-8") as f:
    f.write("pH,log_a_Hplus,logfO2,converged,wlm_mol,wlm_present\n")
    for i, pH in enumerate(pH_values):
        for j, logfO2 in enumerate(logfO2_values):
            converged = int(field[i, j] >= 0)
            present = int(field[i, j] == 1)
            wlm = wlm_amount[i, j]
            wlm_txt = "nan" if np.isnan(wlm) else f"{wlm:.16e}"
            f.write(
                f"{pH:.4f},{-pH:.4f},{logfO2:.4f},{converged},{wlm_txt},{present}\n"
            )

# Plot diagram.
fig, ax = plt.subplots(figsize=(8.8, 6.2))

# Plot converged stability field first.
extent = [
    logfO2_values.min(),
    logfO2_values.max(),
    logaH_values.min(),
    logaH_values.max(),
]

stable_mask = np.where(field == 1, 1.0, np.nan)
notstable_mask = np.where(field == 0, 1.0, np.nan)

ax.imshow(
    notstable_mask,
    origin="lower",
    extent=extent,
    aspect="auto",
    interpolation="nearest",
    cmap="Greys",
    alpha=0.35,
)
ax.imshow(
    stable_mask,
    origin="lower",
    extent=extent,
    aspect="auto",
    interpolation="nearest",
    cmap="Greens",
    alpha=0.85,
)

# Mark non-convergence points.
nc_i, nc_j = np.where(field == -1)
if len(nc_i) > 0:
    ax.scatter(
        logfO2_values[nc_j],
        logaH_values[nc_i],
        marker="x",
        s=12,
        c="#d62728",
        linewidths=0.7,
        label="No convergence",
    )

ax.set_xlabel("log10 fO2 (bar)")
ax.set_ylabel("log10 a(H+)")
ax.set_title("Willemite Activity-Activity Diagram (Zincite Suppressed)")
ax.grid(True, alpha=0.2)

from matplotlib.patches import Patch

legend_handles = [
    Patch(facecolor="#7f7f7f", alpha=0.35, label="Converged: Wlm absent"),
    Patch(facecolor="#2ca02c", alpha=0.85, label="Converged: Wlm stable"),
]
if len(nc_i) > 0:
    ax.legend(handles=legend_handles, loc="best")
else:
    ax.legend(handles=legend_handles, loc="best")

out_png = os.path.join(REPO, "temp", "willemite_activity_activity_diagram.png")
fig.tight_layout()
fig.savefig(out_png, dpi=260)

n_total = field.size
n_conv = int(np.sum(field >= 0))
n_stable = int(np.sum(field == 1))

print("MODEL=PerplexDEW+PerplexGFSM (Znc suppressed)")
print(f"ANCHOR=T={T_C}C P={P_BAR}bar XCO2={XCO2}")
print(
    f"GRID=pH[{pH_values.min():.2f},{pH_values.max():.2f}] x logfO2[{logfO2_values.min():.2f},{logfO2_values.max():.2f}]"
)
print(f"COUNTS total={n_total} converged={n_conv} wlm_stable={n_stable}")
print(f"CSV={out_csv}")
print(f"PNG={out_png}")
