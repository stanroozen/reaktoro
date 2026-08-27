import importlib.util
import json
import sys
from pathlib import Path

p = Path(
    r"DEW_Experimental_Benchmark/Tutorial/willemite_solubility_tutorial_dew17hp622_zn.py"
)
spec = importlib.util.spec_from_file_location("m", p)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

db_path = Path(m.PERPLEX_DATABASE_FILE)
with db_path.open("r", encoding="utf-8") as file:
    database_data = json.load(file)

available_gases = {
    name
    for name, entry in database_data.get("Species", {}).items()
    if entry.get("AggregateState") == "Gas"
}

combo = sys.argv[1] if len(sys.argv) > 1 else m.GFSM_GAS_PHASE_SPECIES
requested_gases = combo.split()
missing_gases = [name for name in requested_gases if name not in available_gases]
if missing_gases:
    raise RuntimeError(
        "Requested gas species not present in database: "
        + ", ".join(missing_gases)
        + ". Available gas species: "
        + ", ".join(sorted(available_gases))
    )

print("combo:", combo)
db = m.Database.fromFile(m.PERPLEX_DATABASE_FILE)
aq = m.AqueousPhase(" ".join(m.AQUEOUS_SPECIES))
aq.setActivityModel(m.AQUEOUS_ACTIVITY_MODEL())
mineral = m.MineralPhase(m.MINERAL_NAME)
gas = m.GaseousPhase(combo)
gas.setActivityModel(m.ActivityModelPerplexGFSM())
system = m.ChemicalSystem(db, aq, mineral, gas)

specs = m.EquilibriumSpecs(system)
specs.temperature()
specs.pressure()
specs.fugacity("O2")
solver = m.EquilibriumSolver(specs)
conds = m.EquilibriumConditions(specs)
conds.temperature(300.0, "celsius")
conds.pressure(2000.0, "bar")
conds.fugacity("O2", 1.0e-20, "bar")

st = m.make_base_state(system)
for s in combo.split():
    try:
        st.set(s, 1.0e-20, "mol")
    except Exception:
        pass

res = solver.solve(st, conds)
print("solve_succeeded:", bool(res.succeeded()))
