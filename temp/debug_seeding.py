"""Investigate what seeded state leads to in failed points."""

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

# Test at the problem point
pH, logfO2 = 2.0, -10.0
T_C = 300.0
P_BAR = 2000.0
XCO2 = 0.3

print(f"Analyzing initial state at pH={pH}, logfO2={logfO2}\n")

# Create seeded state
st = w.make_base_state(system)
st.set("H2O(g)", 1.0 - XCO2, "mol")
st.set("CO2(g)", XCO2, "mol")
st.set("CO2,aq", 1e-6 * XCO2, "mol")
st.set("O2", 1e-20, "mol")

print("Initial seed state (before solving):")
print(f"  CO2,aq amount: {float(st.speciesAmount('CO2,aq')):e} mol")
print(f"  Wlm amount: {float(st.speciesAmount('Wlm')):e} mol")
print(f"  H2O,aq amount: {float(st.speciesAmount('H2O,aq')):e} mol")

# Now test: try solving at a nearby CONVERGENT point
print("\n\nTesting at NEARBY CONVERGENT point (pH=2.5, logfO2=-20):")
pH2, logfO2_2 = 2.5, -20.0

opts = w.EquilibriumOptions()
if hasattr(opts, "hessian") and hasattr(w, "GibbsHessian"):
    opts.hessian = w.GibbsHessian.Exact

solver = w.make_equilibrium_solver(system, specs)
solver.setOptions(opts)

st2 = w.make_base_state(system)
st2.set("H2O(g)", 1.0 - XCO2, "mol")
st2.set("CO2(g)", XCO2, "mol")
st2.set("CO2,aq", 1e-6 * XCO2, "mol")
st2.set("O2", 1e-20, "mol")

conds2 = w.EquilibriumConditions(specs)
conds2.temperature(T_C, "celsius")
conds2.pressure(P_BAR, "bar")
conds2.pH(float(pH2))
conds2.fugacity("O2", 10.0 ** float(logfO2_2), "bar")

r2 = solver.solve(st2, conds2)
print(f"  Succeeded: {r2.succeeded()}")
print(f"  Iterations: {r2.iterations()}")
if r2.succeeded():
    print(f"  Willemite: {float(st2.speciesAmount('Wlm')):e} mol")

# Now try the failed point as CONTINUATION from the successful point
print(
    "\n\nRetrying FAILED point (pH=2.0, logfO2=-10) starting from successful neighbor:"
)

st_retry = st2  # Use converged state as starting point for the problem point

opts_retry = w.EquilibriumOptions()
if hasattr(opts_retry, "hessian") and hasattr(w, "GibbsHessian"):
    opts_retry.hessian = w.GibbsHessian.Exact

solver_retry = w.make_equilibrium_solver(system, specs)
solver_retry.setOptions(opts_retry)

conds_retry = w.EquilibriumConditions(specs)
conds_retry.temperature(T_C, "celsius")
conds_retry.pressure(P_BAR, "bar")
conds_retry.pH(float(pH))
conds_retry.fugacity("O2", 10.0 ** float(logfO2), "bar")

r_retry = solver_retry.solve(st_retry, conds_retry)
print(f"  Succeeded: {r_retry.succeeded()}")
print(f"  Iterations: {r_retry.iterations()}")
if r_retry.succeeded():
    print(f"  Willemite: {float(st_retry.speciesAmount('Wlm')):e} mol")
