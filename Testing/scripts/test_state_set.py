import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "build" / "Reaktoro" / "Release"))
import autodiff
from reaktoro4py import *

db = DEWDatabase("dew2024-aqueous")
hp_db = Database.fromFile(
    str(
        REPO_ROOT / "embedded" / "databases" / "hollandpowell" / "tc-ds62-reaktoro.json"
    )
)
combined_db = Database(db.species())
combined_db.addSpecies(hp_db.species("br"))

aqueous = AqueousPhase("H2O(aq) H+(aq) OH-(aq) Mg+2(aq) MgOH+(aq)")
aqueous.setActivityModel(ActivityModelDEW())
mineral = MineralPhase("br")
system = ChemicalSystem(combined_db, aqueous, mineral)

state = ChemicalState(system)

# Test 1: plain float
print("Test 1: state.set with plain float...")
try:
    state.set("H2O(aq)", 55.5, "mol")
    print("  PASS: float accepted")
except TypeError as e:
    print(f"  FAIL TypeError: {e}")

# Test 2: autodiff.real
print("Test 2: state.set with autodiff.real...")
try:
    state.set("H2O(aq)", autodiff.real(55.5), "mol")
    print("  PASS: autodiff.real accepted")
except Exception as e:
    print(f"  FAIL {type(e).__name__}: {e}")

# Test 3: check autodiff types available
print("Test 3: autodiff types...")
print(f"  autodiff.real: {autodiff.real}")
try:
    v = autodiff.real(3.14)
    print(f"  autodiff.real(3.14) = {v}, type = {type(v)}")
except Exception as e:
    print(f"  Error: {e}")

# Test 4: full equilibrium solve at one T/P
print("Test 4: full equilibrium solve at 200C, 1kbar...")
state2 = ChemicalState(system)


def set_s(s, name, val, unit):
    try:
        s.set(name, float(val), unit)
    except TypeError:
        s.set(name, autodiff.real(float(val)), unit)


set_s(state2, "H2O(aq)", 55.5, "mol")
set_s(state2, "H+(aq)", 1e-8, "mol")
set_s(state2, "OH-(aq)", 1e-8, "mol")
set_s(state2, "Mg+2(aq)", 1e-8, "mol")
set_s(state2, "br", 10.0, "mol")

specs = EquilibriumSpecs(system)
specs.temperature()
specs.pressure()
solver = EquilibriumSolver(specs)
cond = EquilibriumConditions(specs)
cond.temperature(200, "celsius")
cond.pressure(1000, "bar")
result = solver.solve(state2, cond)
print(f"  Solver succeeded: {result.succeeded()}")
if result.succeeded():
    aq = AqueousProps(state2)
    mg = float(aq.elementMolality("Mg"))
    print(f"  Mg molality: {mg:.4e} mol/kg")
