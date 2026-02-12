"""
Check if H2O(aq) from supcrtbl is actually being used or if it's just a linker
"""

import sys
import os

# Setup paths
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(SCRIPT_DIR)
PYD_DIR = os.path.join(ROOT_DIR, "build-msvc", "Reaktoro", "Release")
if os.path.isdir(PYD_DIR) and PYD_DIR not in sys.path:
    sys.path.insert(0, PYD_DIR)

try:
    from reaktoro4py import *

    print("✓ Loaded reaktoro4py extension")
except Exception as e:
    print(f"✗ Failed to load reaktoro4py: {e}")
    sys.exit(1)

# Load databases
print("\n--- Loading databases ---")
dew_db = DEWDatabase("dew2019-aqueous")
supcrt_db = SupcrtDatabase("supcrtbl")

print(f"DEW species count: {len(dew_db.species())}")
print(f"SUPCRT species count: {len(supcrt_db.species())}")

# Check if H2O is in DEW
dew_species_names = [s.name() for s in dew_db.species()]
print(f"\nH2O in DEW: {'H2O(aq)' in dew_species_names}")

# Get H2O from both sources
print("\n--- H2O(aq) Sources ---")
try:
    h2o_supcrt = supcrt_db.species("H2O(aq)")
    print(f"✓ SUPCRT H2O(aq) found: {h2o_supcrt.name()}")
    print(f"  Formula: {h2o_supcrt.formula()}")
    print(f"  Aggregate state: {h2o_supcrt.aggregateState()}")
except Exception as e:
    print(f"✗ SUPCRT H2O(aq) error: {e}")

# Now test if changing to only DEW works
print("\n--- Test 1: Build system with SUPCRT H2O(aq) ---")
try:
    quartz_species = supcrt_db.species("Quartz")
    water_species = supcrt_db.species("H2O(aq)")

    combined_db = Database(dew_db.species())
    combined_db.addSpecies(quartz_species)
    combined_db.addSpecies(water_species)

    aqueous = AqueousPhase("H2O(aq) H+ OH- SiO2_aq")
    try:
        aqueous.setActivityModel(ActivityModelDEW())
    except:
        aqueous.setActivityModel(ActivityModelHKF())

    mineral = MineralPhase("Quartz")
    system1 = ChemicalSystem(combined_db, aqueous, mineral)
    print(f"✓ System 1 created with SUPCRT H2O(aq)")
    print(f"  Species in system: {len(system1.species())}")

    # Try to equilibrate and check if H2O thermodynamics matter
    solver = EquilibriumSolver(system1)
    state = ChemicalState(system1)
    state.set("H2O(aq)", 1.0, "kg")
    state.set("SiO2(aq)", 0.001, "mol")
    state.set("Quartz", 1.0, "kg")
    state.temperature(200, "celsius")
    state.pressure(1000, "bar")

    result = solver.solve(state)
    if result.solved():
        print(f"  ✓ Equilibrium solved")
        print(
            f"    H2O(aq) g0 = {state.speciesStandardGibbsEnergy('H2O(aq)')[0]:.2f} J/mol"
        )
    else:
        print(f"  ✗ Equilibrium failed")

except Exception as e:
    print(f"✗ System 1 error: {e}")
    import traceback

    traceback.print_exc()

print("\n--- Done ---")
