import importlib.util
from pathlib import Path

p = Path(
    r"DEW_Experimental_Benchmark/Tutorial/willemite_solubility_tutorial_dew17hp622_zn.py"
)
spec = importlib.util.spec_from_file_location("m", p)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

combos = [
    "CO2 H2 O2 SO2 H2S N2 NH3 CH4 CO",
    "H2O CO2 CO CH4 H2 H2S O2 SO2 N2 NH3",
    "O2 SO2",
]

for combo in combos:
    print("\n=== combo:", combo)
    try:
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
        conds.fugacity("O2", 1e-20, "bar")
        st = m.make_base_state(system)
        for s in combo.split():
            try:
                st.set(s, 1e-20, "mol")
            except Exception:
                pass
        res = solver.solve(st, conds)
        print("solve succeeded:", bool(res.succeeded()))
    except Exception as e:
        print("python exception:", e)
