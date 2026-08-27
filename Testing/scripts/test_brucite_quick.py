"""
Quick diagnostic: test state.set and one equilibrium solve for brucite.
Run with: conda run -n reaktoro python temp/test_brucite_quick.py
"""

import sys, os

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(os.path.dirname(THIS_DIR))
BUILD_DIR = os.path.join(REPO_ROOT, "build", "Reaktoro", "Release")
sys.path.insert(0, BUILD_DIR)
import autodiff
from reaktoro4py import *

print("Imports OK", flush=True)

hp_db = Database.fromFile(
    os.path.join(
        REPO_ROOT, "embedded", "databases", "hollandpowell", "tc-ds62-reaktoro.json"
    )
)
dew_db = DEWDatabase("dew2024-aqueous")
combined_db = Database(dew_db.species())
combined_db.addSpecies(hp_db.species("br"))
print("Databases OK", flush=True)

aqueous = AqueousPhase("H2O(aq) H+(aq) OH-(aq) Mg+2(aq) MgOH+(aq)")
aqueous.setActivityModel(ActivityModelDEW())
mineral = MineralPhase("br")
system = ChemicalSystem(combined_db, aqueous, mineral)
print("System OK", flush=True)

state = ChemicalState(system)
print("Testing state.set with float...", flush=True)
try:
    state.set("H2O(aq)", 55.5, "mol")
    print("  float: OK", flush=True)
except TypeError as e:
    print(f"  float TypeError -> trying autodiff.real", flush=True)
    state.set("H2O(aq)", autodiff.real(55.5), "mol")
    print("  autodiff.real: OK", flush=True)

# Set rest of initial state
for name, val, unit in [
    ("H+(aq)", 1e-8, "mol"),
    ("OH-(aq)", 1e-8, "mol"),
    ("Mg+2(aq)", 1e-8, "mol"),
    ("br", 10.0, "mol"),
]:
    try:
        state.set(name, val, unit)
    except TypeError:
        state.set(name, autodiff.real(val), unit)
print("Initial state set OK", flush=True)

# Run single solve at 200C, 1kbar
specs = EquilibriumSpecs(system)
specs.temperature()
specs.pressure()
solver = EquilibriumSolver(specs)
cond = EquilibriumConditions(specs)
cond.temperature(200, "celsius")
cond.pressure(1000, "bar")
result = solver.solve(state, cond)
print(f"Solver succeeded: {result.succeeded()}", flush=True)
if result.succeeded():
    aq = AqueousProps(state)
    mg = float(aq.elementMolality("Mg"))
    print(f"Mg molality at 200C, 1kbar: {mg:.4e} mol/kg", flush=True)
else:
    print("FAILED to converge", flush=True)
