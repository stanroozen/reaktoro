"""Test if increasing iteration limit fixes non-convergence."""

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

# Test a problem point with different iteration limits
pH, logfO2 = 2.0, -10.0

T_C = 300.0
P_BAR = 2000.0
XCO2 = 0.3

print(f"Testing pH={pH}, logfO2={logfO2} with different iteration limits:\n")

for max_iters in [201, 500, 1000, 2000]:
    print(f"Max iterations: {max_iters}")

    opts = w.EquilibriumOptions()
    if hasattr(opts, "hessian") and hasattr(w, "GibbsHessian"):
        opts.hessian = w.GibbsHessian.Exact
    if hasattr(opts, "max_iterations"):
        opts.max_iterations = max_iters

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
    print()
