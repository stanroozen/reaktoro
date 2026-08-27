"""Test if warm-starting from neighbor helps convergence."""

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

T_C = 300.0
P_BAR = 2000.0
XCO2 = 0.3

opts = w.EquilibriumOptions()
if hasattr(opts, "hessian") and hasattr(w, "GibbsHessian"):
    opts.hessian = w.GibbsHessian.Exact


def make_fresh_state():
    st = w.make_base_state(system)
    st.set("H2O(g)", 1.0 - XCO2, "mol")
    st.set("CO2(g)", XCO2, "mol")
    st.set("CO2,aq", 1e-6 * XCO2, "mol")
    st.set("O2", 1e-20, "mol")
    return st


# Test 1: Try failed point with fresh state (should fail)
print("TEST 1: Failed point (pH=2.0, logfO2=-10) with FRESH state")
solver1 = w.make_equilibrium_solver(system, specs)
solver1.setOptions(opts)

st1 = make_fresh_state()
conds1 = w.EquilibriumConditions(specs)
conds1.temperature(T_C, "celsius")
conds1.pressure(P_BAR, "bar")
conds1.pH(2.0)
conds1.fugacity("O2", 10.0 ** (-10.0), "bar")

r1 = solver1.solve(st1, conds1)
print(f"  Succeeded: {r1.succeeded()}")
print(f"  Iterations: {r1.iterations()}\n")

# Test 2: Try successful neighbor (should converge)
print("TEST 2: Successful neighbor (pH=3.0, logfO2=-20) with FRESH state")
solver2 = w.make_equilibrium_solver(system, specs)
solver2.setOptions(opts)

st2 = make_fresh_state()
conds2 = w.EquilibriumConditions(specs)
conds2.temperature(T_C, "celsius")
conds2.pressure(P_BAR, "bar")
conds2.pH(3.0)
conds2.fugacity("O2", 10.0 ** (-20.0), "bar")

r2 = solver2.solve(st2, conds2)
print(f"  Succeeded: {r2.succeeded()}")
print(f"  Iterations: {r2.iterations()}\n")

# Test 3: Try failed point with WARM-START from neighbor
if r2.succeeded():
    print("TEST 3: Failed point (pH=2.0, logfO2=-10) with WARM-START from neighbor")
    solver3 = w.make_equilibrium_solver(system, specs)
    solver3.setOptions(opts)

    # Use st2 (converged state) as starting point
    conds3 = w.EquilibriumConditions(specs)
    conds3.temperature(T_C, "celsius")
    conds3.pressure(P_BAR, "bar")
    conds3.pH(2.0)
    conds3.fugacity("O2", 10.0 ** (-10.0), "bar")

    r3 = solver3.solve(st2, conds3)  # Note: st2, not a fresh state
    print(f"  Succeeded: {r3.succeeded()}")
    print(f"  Iterations: {r3.iterations()}\n")

    if r3.succeeded():
        print("CONCLUSION: Warm-starting fixes convergence!")
    else:
        print("CONCLUSION: Even warm-starting fails - deep thermodynamic issue")
