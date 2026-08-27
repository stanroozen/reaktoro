import importlib.util
import math
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


# Model configuration: mixed fluid with H2O-CO2 plus O2 fugacity control.
w.USE_COMPETING_ZN_MINERALS = True
w.INCLUDE_COMPATIBLE_GANGUE_MINERALS = False
w.COMPETING_ZN_MINERALS = list(dict.fromkeys(w.COMPETING_ZN_MINERALS + ["hem", "cc", "q", "mag"]))

extra_aq = [
    "CO2,aq",
    "Ca+2", "CaCO3,aq", "Ca(HCO3)",
    "Mg+2", "MgCO3,aq", "Mg(HCO3)",
    "Fe+2", "Fe+3", "Fe(OH)+", "FeO,aq",
]
w.AQUEOUS_SPECIES = list(dict.fromkeys(list(w.AQUEOUS_SPECIES) + extra_aq))

for m in ["hem", "cc", "q", "mag"]:
    w.INITIAL_SPECIES_AMOUNTS_MOL[m] = 2.0

w.validate_user_inputs()


dew_db = w.Database.fromFile(w.PERPLEX_DATABASE_FILE)
supcrt = w.Database.fromFile(os.path.join(REPO, "embedded", "databases", "reaktoro", "supcrt07.yaml"))
dew_db.addSpecies(supcrt.species("H2O(g)"))
dew_db.addSpecies(supcrt.species("CO2(g)"))

# PerplexDEW in warn-only mode for known GFSM/HKF co-solvent conflicts.
params = w.ActivityModelParamsPerplexDEW()
params.errorOnConflictingStandardState = False
aq = w.AqueousPhase(" ".join(w.AQUEOUS_SPECIES))
aq.setActivityModel(w.ActivityModelPerplexDEW(params))
mins = w.make_mineral_phases()

# Keep O2 (DEW name) for fugacity control plus explicit H2O(g)/CO2(g) mixed fluid.
gas = w.GaseousPhase("H2O(g) CO2(g) O2")
gas.setActivityModel(w.ActivityModelPerplexGFSM(w.ActivityModelParamsPerplexGFSM()))

system = w.ChemicalSystem(dew_db, aq, mins, gas)

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

TH = 1.0e-8
N_GAS = 1.0

BASE = {
    "T": 240.0,
    "P": 1000.0,
    "pH": 3.0,
    "logfO2": -40.0,
    "XCO2": 0.3,
    "Fe": 0.0,
    "Ca": 0.0,
    "Mg": 0.0,
}


def solve_case(case):
    st = w.make_base_state(system)

    xco2 = max(1.0e-8, min(0.999999, float(case["XCO2"])))
    st.set("H2O(g)", (1.0 - xco2) * N_GAS, "mol")
    st.set("CO2(g)", xco2 * N_GAS, "mol")
    st.set("O2", 1.0e-20, "mol")
    st.set("CO2,aq", max(1.0e-12, 1.0e-6 * xco2), "mol")

    if case["Fe"] > 0.0:
        st.set("Fe+2", float(case["Fe"]), "mol")
    if case["Ca"] > 0.0:
        st.set("Ca+2", float(case["Ca"]), "mol")
    if case["Mg"] > 0.0:
        st.set("Mg+2", float(case["Mg"]), "mol")

    conditions.temperature(float(case["T"]), "celsius")
    conditions.pressure(float(case["P"]), "bar")
    conditions.pH(float(case["pH"]))
    conditions.fugacity("O2", 10.0 ** float(case["logfO2"]), "bar")

    res = solver.solve(st, conditions)
    if not res.succeeded():
        return False, None

    wlm = float(st.speciesAmount("Wlm"))
    return True, wlm


def first_threshold(var_name, values):
    found = None
    trace = []

    for v in values:
        case = dict(BASE)
        case[var_name] = float(v)
        ok, wlm = solve_case(case)
        trace.append((float(v), ok, None if wlm is None else float(wlm)))
        if ok and wlm is not None and wlm > TH and found is None:
            found = float(v)

    return found, trace


def print_trace(title, trace):
    print(f"\n[{title}]")
    for v, ok, wlm in trace:
        if not ok:
            print(f"  {v: .6g} -> NO_CONVERGENCE")
        else:
            print(f"  {v: .6g} -> Wlm={wlm:.3e} present={wlm > TH}")


print("MODEL=PerplexDEW(aq)+PerplexGFSM(H2O(g) CO2(g) O2)")
print(f"BASE= T={BASE['T']}C P={BASE['P']}bar pH={BASE['pH']} logfO2={BASE['logfO2']} XCO2={BASE['XCO2']} Fe={BASE['Fe']} Ca={BASE['Ca']} Mg={BASE['Mg']}")
print(f"THRESHOLD=Wlm>{TH:.1e} mol")

# 1) pH threshold
ph_values = [2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0]
ph_found, ph_trace = first_threshold("pH", ph_values)
print_trace("pH scan", ph_trace)
print(f"pH_first_threshold={ph_found}")

# 2) oxygen fugacity threshold (log10 bar)
fo2_values = [-60.0, -55.0, -50.0, -45.0, -40.0, -35.0, -30.0, -25.0, -20.0, -15.0, -10.0]
fo2_found, fo2_trace = first_threshold("logfO2", fo2_values)
print_trace("logfO2 scan", fo2_trace)
print(f"logfO2_first_threshold={fo2_found}")

# 3) XCO2 threshold
xco2_values = [0.02, 0.05, 0.1, 0.2, 0.35, 0.5, 0.7, 0.9]
xco2_found, xco2_trace = first_threshold("XCO2", xco2_values)
print_trace("XCO2 scan", xco2_trace)
print(f"XCO2_first_threshold={xco2_found}")

# 4) bulk composition thresholds (one-at-a-time)
comp_values = [0.0, 1.0e-8, 1.0e-7, 1.0e-6, 1.0e-5, 1.0e-4, 1.0e-3, 1.0e-2]
for elem in ["Fe", "Ca", "Mg"]:
    found, trace = first_threshold(elem, comp_values)
    print_trace(f"{elem} bulk scan (mol seed)", trace)
    print(f"{elem}_first_threshold={found}")

print("\nDONE")
