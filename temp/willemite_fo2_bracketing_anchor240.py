import importlib.util
import os
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REL = os.path.join(REPO, "build", "Reaktoro", "Release")
if REL not in sys.path:
    sys.path.insert(0, REL)

TUTORIAL = os.path.join(REPO, "DEW_Experimental_Benchmark", "Tutorial", "willemite_solubility_tutorial_dew17hp622_zn.py")
spec = importlib.util.spec_from_file_location("w", TUTORIAL)
w = importlib.util.module_from_spec(spec)
spec.loader.exec_module(w)

try:
    w.Warnings.disable(906)
except Exception:
    pass

w.USE_COMPETING_ZN_MINERALS = True
w.INCLUDE_COMPATIBLE_GANGUE_MINERALS = False
w.COMPETING_ZN_MINERALS = list(dict.fromkeys(w.COMPETING_ZN_MINERALS + ["hem", "cc", "q", "mag"]))

extra = ["CO2,aq", "Fe+2", "Ca+2", "Mg+2", "Fe+3", "CaCO3,aq", "MgCO3,aq"]
w.AQUEOUS_SPECIES = list(dict.fromkeys(list(w.AQUEOUS_SPECIES) + extra))
for m in ["hem", "cc", "q", "mag"]:
    w.INITIAL_SPECIES_AMOUNTS_MOL[m] = 2.0
w.validate_user_inputs()

dew = w.Database.fromFile(w.PERPLEX_DATABASE_FILE)
sup = w.Database.fromFile(os.path.join(REPO, "embedded", "databases", "reaktoro", "supcrt07.yaml"))
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
specs.temperature(); specs.pressure(); specs.pH(); specs.fugacity("O2")
solver = w.make_equilibrium_solver(system, specs)
opts = w.EquilibriumOptions()
if hasattr(opts, "hessian") and hasattr(w, "GibbsHessian"):
    opts.hessian = w.GibbsHessian.Exact
solver.setOptions(opts)
conds = w.EquilibriumConditions(specs)

# Requested anchor and reducing sweep.
T = 240.0
P = 1000.0
pH = 3.0
XCO2 = 0.3
TH = 1e-8
logfO2_values = [-70, -65, -60, -55, -50, -45, -40, -35, -30, -25, -20, -15, -10]

print("MODEL=PerplexDEW+PerplexGFSM")
print(f"ANCHOR=T={T}C P={P}bar pH={pH} XCO2={XCO2}")
print("SCAN=logfO2")

first = None
for lf in logfO2_values:
    st = w.make_base_state(system)
    st.set("H2O(g)", 1.0 - XCO2, "mol")
    st.set("CO2(g)", XCO2, "mol")
    st.set("CO2,aq", 1e-6 * XCO2, "mol")
    st.set("O2", 1e-20, "mol")

    conds.temperature(T, "celsius")
    conds.pressure(P, "bar")
    conds.pH(pH)
    conds.fugacity("O2", 10.0 ** float(lf), "bar")

    r = solver.solve(st, conds)
    if not r.succeeded():
        print(f"logfO2={lf: .1f} -> NO_CONVERGENCE")
        continue

    wlm = float(st.speciesAmount("Wlm"))
    present = wlm > TH
    print(f"logfO2={lf: .1f} -> Wlm={wlm:.3e} present={present}")
    if present and first is None:
        first = float(lf)

print(f"logfO2_first_threshold={first}")
