import traceback
import importlib.util
from pathlib import Path
import numpy as np

try:
    p = Path(
        r"DEW_Experimental_Benchmark/Tutorial/willemite_solubility_tutorial_dew17hp622_zn.py"
    )
    spec = importlib.util.spec_from_file_location("m", p)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)

    db = m.Database.fromFile(m.PERPLEX_DATABASE_FILE)
    aq = m.AqueousPhase(" ".join(m.AQUEOUS_SPECIES))
    aq.setActivityModel(m.AQUEOUS_ACTIVITY_MODEL())
    mineral = m.MineralPhase(m.MINERAL_NAME)
    system = m.ChemicalSystem(db, aq, mineral)

    for gas, val in [("O2", 1e-20), ("SO2", 1e-8)]:
        specs = m.EquilibriumSpecs(system)
        specs.temperature()
        specs.pressure()
        specs.fugacity(gas)
        solver = m.EquilibriumSolver(specs)
        conds = m.EquilibriumConditions(specs)
        conds.temperature(300.0, "celsius")
        conds.pressure(2000.0, "bar")
        conds.fugacity(gas, val, "bar")
        st = m.make_base_state(system)
        res = solver.solve(st, conds)
        ok = bool(res.succeeded())
        mol = m.dissolved_element_molality(st, m.SOLVENT_SPECIES_NAME)
        print(gas, "solve", ok, "Zn_molality", float(mol) if np.isfinite(mol) else mol)
except Exception as e:
    print("ERROR", e)
    traceback.print_exc()
