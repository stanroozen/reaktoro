import importlib.util
import os
import sys
import optima

repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
rel = os.path.join(repo, "build", "Reaktoro", "Release")
if rel not in sys.path:
    sys.path.insert(0, rel)

tutorial = os.path.join(
    repo,
    "DEW_Experimental_Benchmark",
    "Tutorial",
    "willemite_solubility_tutorial_dew17hp622_zn.py",
)
spec = importlib.util.spec_from_file_location("w", tutorial)
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
    os.path.join(repo, "embedded", "databases", "reaktoro", "supcrt07.yaml")
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
inner = optima.Options()
inner.maxiters = 1000
opts.optima = inner


def make_seeded_state():
    st = w.make_base_state(system)
    st.set("H2O(g)", 0.7, "mol")
    st.set("CO2(g)", 0.3, "mol")
    st.set("CO2,aq", 3e-7, "mol")
    st.set("O2", 1e-20, "mol")
    return st


pH_values = [2.0 + 0.5 * i for i in range(13)]
logfO2_values = [-36.0 + 2.0 * i for i in range(14)]
converged = 0
failed = 0

for pH in pH_values:
    for lf in logfO2_values:
        solver = w.make_equilibrium_solver(system, specs)
        solver.setOptions(opts)
        st = make_seeded_state()
        conds = w.EquilibriumConditions(specs)
        conds.temperature(300.0, "celsius")
        conds.pressure(2000.0, "bar")
        conds.pH(float(pH))
        conds.fugacity("O2", 10.0 ** float(lf), "bar")
        result = solver.solve(st, conds)
        if result.succeeded():
            converged += 1
        else:
            failed += 1
        print(
            f"POINT pH={pH:.1f} logfO2={lf:.1f} ok={result.succeeded()} iter={result.iterations()}",
            flush=True,
        )

print(f"FINAL converged={converged} failed={failed}", flush=True)
