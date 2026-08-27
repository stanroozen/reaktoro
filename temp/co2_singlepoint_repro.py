import importlib.util
import traceback

p = "DEW_Experimental_Benchmark/Tutorial/willemite_solubility_tutorial_dew17hp622_zn.py"
s = importlib.util.spec_from_file_location("m", p)
m = importlib.util.module_from_spec(s)
s.loader.exec_module(m)

m.USE_COMPETING_ZN_MINERALS = True
m.INCLUDE_COMPATIBLE_GANGUE_MINERALS = False
m.COMPETING_ZN_MINERALS = list(
    dict.fromkeys(m.COMPETING_ZN_MINERALS + ["hem", "cc", "q"])
)
m.validate_user_inputs()

db = m.Database.fromFile(m.PERPLEX_DATABASE_FILE)
aq = m.AqueousPhase(" ".join(m.AQUEOUS_SPECIES))
aq.setActivityModel(m.AQUEOUS_ACTIVITY_MODEL())
mins = m.make_mineral_phases()
gas = m.GaseousPhase("CO2")
gas.setActivityModel(m.ActivityModelIdealGas())
sysm = m.ChemicalSystem(db, aq, mins, gas)

specs = m.EquilibriumSpecs(sysm)
specs.temperature()
specs.pressure()
specs.pH()
specs.fugacity("CO2")
solver = m.EquilibriumSolver(specs)
conds = m.EquilibriumConditions(specs)

state = m.make_base_state(sysm)
state.set("CO2", 1.0e-20, "mol")
conds.temperature(300.0, "celsius")
conds.pressure(2000.0, "bar")
conds.pH(7.0)
conds.fugacity("CO2", 1.0e-3, "bar")

print("START_SOLVE")
try:
    result = solver.solve(state, conds)
    print("SOLVE_OK", bool(result.succeeded()))
except Exception as exc:
    print("PY_EXCEPTION", exc)
    traceback.print_exc()
print("END")
