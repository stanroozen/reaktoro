"""Debug non-convergence at specific pH/logfO2 points."""

import importlib.util
import os
import sys

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
w.COMPETING_ZN_MINERALS = [m for m in w.COMPETING_ZN_MINERALS if m != "Znc"]
w.COMPETING_ZN_MINERALS = list(
    dict.fromkeys(w.COMPETING_ZN_MINERALS + ["hem", "cc", "q", "mag"])
)

extra = ["CO2,aq", "Fe+2", "Fe+3", "Ca+2", "Mg+2", "CaCO3,aq", "MgCO3,aq"]
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

# Test a few failed points
failed_points = [
    (2.0, -10.0),  # pH=2, logfO2=-10
    (5.5, -24.0),  # pH=5.5, logfO2=-24
    (8.0, -14.0),  # pH=8, logfO2=-14
]

T_C = 300.0
P_BAR = 2000.0
XCO2 = 0.3

opts = w.EquilibriumOptions()
if hasattr(opts, "hessian") and hasattr(w, "GibbsHessian"):
    opts.hessian = w.GibbsHessian.Exact

print(f"Testing {len(failed_points)} previously failed points...\n")

for pH, logfO2 in failed_points:
    print(f"pH={pH:.1f}, logfO2={logfO2:.1f}:")

    # Make fresh solver and state
    solver = w.make_equilibrium_solver(system, specs)
    solver.setOptions(opts)

    st = w.make_base_state(system)
    st.set("H2O(g)", 1.0 - XCO2, "mol")
    st.set("CO2(g)", XCO2, "mol")
    st.set("CO2,aq", 1e-6 * XCO2, "mol")
    st.set("O2", 1e-20, "mol")

    conds = w.EquilibriumConditions(specs)
    conds.temperature(T_C, "celsius")
    conds.pressure(P_BAR, "bar")
    conds.pH(float(pH))
    conds.fugacity("O2", 10.0 ** float(logfO2), "bar")

    r = solver.solve(st, conds)

    print(f"  Succeeded: {r.succeeded()}")
    print(f"  Iterations: {r.iterations()}")

    if r.succeeded():
        wlm = float(st.speciesAmount("Wlm"))
        print(f"  Willemite: {wlm:.4e} mol")
        # Check some aqueous species
        print(f"  H+ conc: {float(st.speciesConcentration('H+')):e}")
        print(f"  Zn2+ conc: {float(st.speciesConcentration('Zn+2')):e}")
        print(f"  ZnO2-2 conc: {float(st.speciesConcentration('ZnO2-2')):e}")
    print()
