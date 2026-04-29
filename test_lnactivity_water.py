"""Quick test: does lnActivity('WATER,AQ') constraint propagate to Si solubility?"""

import sys, os, importlib

if os.name == "nt":
    ep = sys.prefix
    env_paths = [
        ep,
        os.path.join(ep, "Library", "mingw-w64", "bin"),
        os.path.join(ep, "Library", "usr", "bin"),
        os.path.join(ep, "Library", "bin"),
        os.path.join(ep, "Scripts"),
        os.path.join(ep, "bin"),
    ]
    sr = os.environ.get("SystemRoot", r"C:\Windows")
    os.environ["PATH"] = ";".join(
        [p for p in env_paths + [os.path.join(sr, "System32"), sr] if os.path.isdir(p)]
    )

import numpy as np
import autodiff  # noqa

try:
    from reaktoro import *  # noqa
except ModuleNotFoundError:
    _pyd = r"C:\Users\stanroozen\Documents\Projects\reaktoro-dev\reaktoro\build-dew\Reaktoro\Release"
    sys.path.insert(0, _pyd)
    r = importlib.import_module("reaktoro4py")
    globals().update({k: getattr(r, k) for k in dir(r) if not k.startswith("_")})

print("Imports OK")

dew_db = DEWDatabase("dew2024-aqueous")
supcrt_db = SupcrtDatabase("supcrtbl")
mineral_species = supcrt_db.species("Quartz")
combined_db = Database(dew_db.species())
combined_db.addSpecies(mineral_species)
aqueous = AqueousPhase("WATER,AQ H+ OH- SiO2_aq HSiO3- Si2O4_aq Si3O6_aq")
aqueous.setActivityModel(ActivityModelPerplexDEW(ActivityDHModel.ExtendedDH))
mineral = MineralPhase("Quartz")
system = ChemicalSystem(combined_db, aqueous, mineral)
Warnings.disable(906)
solver = EquilibriumSolver(system)

# Use EquilibriumSpecs to register lnActivity("WATER,AQ") as a constraint
specs = EquilibriumSpecs(system)
specs.temperature()
specs.pressure()
specs.lnActivity("WATER,AQ")
solver2 = EquilibriumSolver(specs)

# Baseline at a_H2O = 1 (no constraint, use plain conditions)
conditions_base = EquilibriumConditions(system)
conditions_base.temperature(800.0, "celsius")
conditions_base.pressure(10000.0, "bar")
state0 = ChemicalState(system)
state0.set("WATER,AQ", 1.0, "kg")
state0.set("SiO2_aq", 1e-6, "mol")
state0.set("Quartz", 10.0, "mol")
result0 = solver.solve(state0, conditions_base)
aq0 = AqueousProps(state0)
m0 = (
    float(aq0.speciesMolality("SiO2_aq"))
    + float(aq0.speciesMolality("HSiO3-"))
    + 2 * float(aq0.speciesMolality("Si2O4_aq"))
    + 3 * float(aq0.speciesMolality("Si3O6_aq"))
)
print(f"a_H2O=1.00 (no constraint): total_Si={m0:.4f} mol/kg")

for a_h2o in [0.8, 0.6, 0.4, 0.2]:
    conditions = EquilibriumConditions(specs)
    conditions.temperature(800.0, "celsius")
    conditions.pressure(10000.0, "bar")
    conditions.lnActivity("WATER,AQ", np.log(a_h2o))
    state = ChemicalState(system)
    state.set("WATER,AQ", 1.0, "kg")
    state.set("SiO2_aq", 1e-6, "mol")
    state.set("Quartz", 10.0, "mol")
    result = solver2.solve(state, conditions)
    if result.succeeded():
        aq = AqueousProps(state)
        m_sio2 = float(aq.speciesMolality("SiO2_aq"))
        m_hsi = float(aq.speciesMolality("HSiO3-"))
        m_si2 = float(aq.speciesMolality("Si2O4_aq"))
        m_si3 = float(aq.speciesMolality("Si3O6_aq"))
        m_si_total = m_sio2 + m_hsi + 2 * m_si2 + 3 * m_si3
        print(
            f"a_H2O={a_h2o:.2f}: SiO2_aq={m_sio2:.4f} HSiO3-={m_hsi:.5f} Si2={m_si2:.4f} total_Si={m_si_total:.4f}"
        )
    else:
        print(f"a_H2O={a_h2o:.2f}: FAILED (convergence)")
