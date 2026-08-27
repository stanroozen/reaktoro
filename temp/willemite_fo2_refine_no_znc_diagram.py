import importlib.util
import math
import os
import sys

import matplotlib.pyplot as plt

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

# Build with zincite suppressed.
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


def make_seeded_state(xco2):
    state = w.make_base_state(system)
    state.set("H2O(g)", 1.0 - xco2, "mol")
    state.set("CO2(g)", xco2, "mol")
    state.set("CO2,aq", max(1.0e-12, 1.0e-6 * xco2), "mol")
    state.set("O2", 1.0e-20, "mol")
    return state


def solve_at(state, anchor, logfO2):
    conditions.temperature(anchor["T"], "celsius")
    conditions.pressure(anchor["P"], "bar")
    conditions.pH(anchor["pH"])
    conditions.fugacity("O2", 10.0**logfO2, "bar")
    return solver.solve(state, conditions)


anchors = [
    {
        "name": "AnchorA_T240_P1000_pH3_XCO2_0.3",
        "T": 240.0,
        "P": 1000.0,
        "pH": 3.0,
        "XCO2": 0.3,
    },
    {
        "name": "AnchorB_T300_P2000_pH4.5_XCO2_0.3",
        "T": 300.0,
        "P": 2000.0,
        "pH": 4.5,
        "XCO2": 0.3,
    },
]

logfO2_values = [x / 4.0 for x in range(-96, -63)]  # -24.0 .. -15.75
TH = 1.0e-8

rows = []
summary = []

for anchor in anchors:
    xco2 = anchor["XCO2"]
    state = make_seeded_state(xco2)

    first_threshold = None

    for lf in logfO2_values:
        res = solve_at(state, anchor, lf)
        if not res.succeeded():
            state = make_seeded_state(xco2)
            res = solve_at(state, anchor, lf)
        if not res.succeeded():
            rows.append((anchor["name"], lf, 0, float("nan"), 0))
            continue

        wlm = float(state.speciesAmount("Wlm"))
        present = int(wlm > TH)
        rows.append((anchor["name"], lf, 1, wlm, present))

        if present and first_threshold is None:
            first_threshold = lf

    solved_count = sum(1 for r in rows if r[0] == anchor["name"] and r[2] == 1)
    hit_count = sum(1 for r in rows if r[0] == anchor["name"] and r[4] == 1)
    summary.append((anchor["name"], solved_count, hit_count, first_threshold))

out_csv = os.path.join(REPO, "temp", "willemite_fo2_refine_no_znc_data.csv")
with open(out_csv, "w", encoding="utf-8") as f:
    f.write("anchor,logfO2,converged,wlm_mol,wlm_present\n")
    for a, lf, conv, wlm, present in rows:
        wlm_txt = "nan" if math.isnan(wlm) else f"{wlm:.16e}"
        f.write(f"{a},{lf:.2f},{conv},{wlm_txt},{present}\n")

fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), sharey=True)
for i, anchor in enumerate(anchors):
    ax = axes[i]
    anchor_rows = [r for r in rows if r[0] == anchor["name"]]

    xs = [r[1] for r in anchor_rows if r[2] == 1]
    ys = [r[3] for r in anchor_rows if r[2] == 1]
    x_nc = [r[1] for r in anchor_rows if r[2] == 0]

    if xs:
        ax.plot(xs, ys, marker="o", linewidth=1.2, markersize=3, color="#1f77b4")
    if x_nc:
        ax.scatter(
            x_nc,
            [TH] * len(x_nc),
            marker="x",
            color="#d62728",
            s=20,
            label="No convergence",
        )

    ax.axhline(TH, color="#2ca02c", linestyle="--", linewidth=1.0)
    ax.set_yscale("log")
    ax.set_xlabel("log10 fO2 (bar)")
    ax.set_title(anchor["name"].replace("_", "\n", 1))
    ax.grid(True, alpha=0.25)

axes[0].set_ylabel("Willemite amount (mol)")
fig.suptitle("Willemite vs logfO2 (Zincite Suppressed)")
fig.tight_layout()

out_png = os.path.join(REPO, "temp", "willemite_fo2_refine_no_znc_diagram.png")
fig.savefig(out_png, dpi=220)

print("MODEL=PerplexDEW+PerplexGFSM (Znc suppressed)")
for name, solved_count, hit_count, first_threshold in summary:
    print(
        f"{name}: solved={solved_count} hits={hit_count} first_threshold={first_threshold}"
    )
print(f"CSV={out_csv}")
print(f"PNG={out_png}")
