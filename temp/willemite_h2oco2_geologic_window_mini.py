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

w.USE_COMPETING_ZN_MINERALS = True
w.INCLUDE_COMPATIBLE_GANGUE_MINERALS = False
w.COMPETING_ZN_MINERALS = list(
    dict.fromkeys(w.COMPETING_ZN_MINERALS + ["hem", "cc", "q", "mag"])
)

extra = ["CO2,aq", "Ca+2", "Mg+2", "Fe+2", "Fe+3", "CaCO3,aq", "MgCO3,aq", "Fe(OH)+"]
w.AQUEOUS_SPECIES = list(dict.fromkeys(list(w.AQUEOUS_SPECIES) + extra))
for m in ["hem", "cc", "q", "mag"]:
    w.INITIAL_SPECIES_AMOUNTS_MOL[m] = 2.0
w.validate_user_inputs()

dew_db = w.Database.fromFile(w.PERPLEX_DATABASE_FILE)
supcrt = w.Database.fromFile(
    os.path.join(REPO, "embedded", "databases", "reaktoro", "supcrt07.yaml")
)
dew_db.addSpecies(supcrt.species("H2O(g)"))
dew_db.addSpecies(supcrt.species("CO2(g)"))

params = w.ActivityModelParamsPerplexDEW()
params.errorOnConflictingStandardState = False
aq = w.AqueousPhase(" ".join(w.AQUEOUS_SPECIES))
aq.setActivityModel(w.ActivityModelPerplexDEW(params))
mins = w.make_mineral_phases()
gas = w.GaseousPhase("H2O(g) CO2(g)")
gas.setActivityModel(w.ActivityModelPerplexGFSM(w.ActivityModelParamsPerplexGFSM()))
system = w.ChemicalSystem(dew_db, aq, mins, gas)

specs = w.EquilibriumSpecs(system)
specs.temperature()
specs.pressure()
specs.pH()
solver = w.make_equilibrium_solver(system, specs)
opts = w.EquilibriumOptions()
if hasattr(opts, "hessian") and hasattr(w, "GibbsHessian"):
    opts.hessian = w.GibbsHessian.Exact
solver.setOptions(opts)
conds = w.EquilibriumConditions(specs)

cases = [
    (240.0, 3.0, 1000.0, 0.1, 0.0, 0.0, 0.0),
    (320.0, 6.0, 3000.0, 0.5, 0.0, 0.0, 0.0),
    (240.0, 4.5, 3000.0, 0.5, 1e-3, 1e-3, 1e-3),
    (320.0, 3.0, 1000.0, 0.1, 1e-3, 0.0, 0.0),
]

TH = 1e-8
print("MODEL=PerplexDEW+PerplexGFSM mini geologic window")
for i, (T, pH, P, X, Fe, Ca, Mg) in enumerate(cases, start=1):
    st = w.make_base_state(system)
    st.set("H2O(g)", (1.0 - X), "mol")
    st.set("CO2(g)", X, "mol")
    st.set("CO2,aq", max(1e-12, 1e-6 * X), "mol")
    if Fe > 0:
        st.set("Fe+2", Fe, "mol")
    if Ca > 0:
        st.set("Ca+2", Ca, "mol")
    if Mg > 0:
        st.set("Mg+2", Mg, "mol")
    conds.temperature(T, "celsius")
    conds.pressure(P, "bar")
    conds.pH(pH)
    r = solver.solve(st, conds)
    if not r.succeeded():
        print(
            f"CASE {i}: NO_CONVERGENCE T={T} pH={pH} P={P} XCO2={X} Fe={Fe} Ca={Ca} Mg={Mg}"
        )
        continue
    wlm = float(st.speciesAmount("Wlm"))
    hem = float(st.speciesAmount("hem"))
    cc = float(st.speciesAmount("cc"))
    print(
        f"CASE {i}: OK T={T} pH={pH} P={P} XCO2={X} Fe={Fe} Ca={Ca} Mg={Mg} Wlm={wlm:.3e} hem={hem:.3e} cc={cc:.3e} WLM_PRESENT={wlm > TH}"
    )
