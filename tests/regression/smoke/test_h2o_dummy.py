"""
Test if H2O(aq) thermodynamics from SUPCRT are actually used by Reaktoro
or if it's just a dummy linker to the IAPWS-95 EOS.
"""

import sys
import os
import numpy as np

# Setup paths
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(SCRIPT_DIR)
PYD_DIR = os.path.join(ROOT_DIR, "build", "Reaktoro", "Release")
if os.path.isdir(PYD_DIR) and PYD_DIR not in sys.path:
    sys.path.insert(0, PYD_DIR)

from reaktoro4py import *

# Suppress warnings
try:
    Warnings.disable(906)
except:
    pass

# Load databases
dew_db = DEWDatabase("dew2019-aqueous")
supcrt_db = SupcrtDatabase("supcrtbl")

print("=" * 70)
print("Testing: Is SUPCRT H2O(aq) thermodynamics actually used?")
print("=" * 70)

# Build System 1: With SUPCRT H2O(aq)
print("\n[System 1] With SUPCRT H2O(aq) (current approach)")
quartz = supcrt_db.species("Quartz")
h2o_supcrt = supcrt_db.species("H2O(aq)")

db1 = Database(dew_db.species())
db1.addSpecies(quartz)
db1.addSpecies(h2o_supcrt)

aqueous1 = AqueousPhase("H2O(aq) H+ OH- SiO2_aq")
try:
    aqueous1.setActivityModel(ActivityModelDEW())
except:
    aqueous1.setActivityModel(ActivityModelHKF())

mineral1 = MineralPhase("Quartz")
system1 = ChemicalSystem(db1, aqueous1, mineral1)

# Build System 2: Without SUPCRT H2O(aq) - just DEW species + Quartz
# (testing if we can use a dummy or if we need SUPCRT entry)
print("\n[System 2] Attempting to use only DEW species + Quartz (NO H2O entry)")

try:
    db2 = Database(dew_db.species())
    db2.addSpecies(quartz)
    # NO h2o_supcrt added

    aqueous2 = AqueousPhase("H+ OH- SiO2_aq")  # removed H2O(aq)
    try:
        aqueous2.setActivityModel(ActivityModelDEW())
    except:
        aqueous2.setActivityModel(ActivityModelHKF())

    mineral2 = MineralPhase("Quartz")
    system2 = ChemicalSystem(db2, aqueous2, mineral2)
    print("âœ“ System 2 created without explicit H2O(aq) entry")
except Exception as e:
    print(f"âœ— System 2 error: {e}")
    system2 = None

# Test equilibration
print("\n" + "=" * 70)
print("Testing equilibration at 300Â°C, 10 MPa")
print("=" * 70)

T_C = 300
P_bar = 100

# System 1 test
print(f"\n[System 1] Equilibrating with SUPCRT H2O(aq)...")
try:
    solver1 = EquilibriumSolver(system1)
    state1 = ChemicalState(system1)
    state1.temperature(float(T_C), "celsius")
    state1.pressure(float(P_bar), "bar")

    # Set composition for open system (with solvent)
    state1.set("SiO2(aq)", 0.0001, "mol")
    state1.set("Quartz", 1.0, "kg")

    result1 = solver1.solve(state1)

    if result1.solved():
        print("âœ“ Equilibrium converged")
        g0_h2o = state1.speciesStandardGibbsEnergy("H2O(aq)")[0]
        T_K = state1.temperature()
        P_Pa = state1.pressure()
        print(f"  Temperature: {T_K:.2f} K ({T_K - 273.15:.1f}Â°C)")
        print(f"  Pressure: {P_Pa / 1e5:.2f} bar ({P_Pa / 1e8:.0f} MPa)")
        print(f"  H2O(aq) G0 = {g0_h2o:.2f} J/mol")

        species_amounts = []
        for i, species in enumerate(system1.species()):
            n = state1.speciesAmount(i)[0]
            if n > 1e-16:
                species_amounts.append((species.name(), float(n)))
        print(f"  Species amounts: {species_amounts[:5]}")
    else:
        print("âœ— Equilibrium failed to converge")
        print(f"  Message: {result1.message()}")

except Exception as e:
    print(f"âœ— System 1 equilibration error: {e}")
    import traceback

    traceback.print_exc()

# System 2 test (if created)
if system2:
    print(f"\n[System 2] Equilibrating without explicit H2O(aq) entry...")
    try:
        solver2 = EquilibriumSolver(system2)
        state2 = ChemicalState(system2)
        state2.temperature(float(T_C), "celsius")
        state2.pressure(float(P_bar), "bar")

        state2.set("SiO2(aq)", 0.0001, "mol")
        state2.set("Quartz", 1.0, "kg")

        result2 = solver2.solve(state2)

        if result2.solved():
            print("âœ“ Equilibrium converged")
            T_K = state2.temperature()
            P_Pa = state2.pressure()
            print(f"  Temperature: {T_K:.2f} K ({T_K - 273.15:.1f}Â°C)")
            print(f"  Pressure: {P_Pa / 1e5:.2f} bar ({P_Pa / 1e8:.0f} MPa)")

            species_amounts = []
            for i, species in enumerate(system2.species()):
                n = state2.speciesAmount(i)[0]
                if n > 1e-16:
                    species_amounts.append((species.name(), float(n)))
            print(f"  Species amounts: {species_amounts[:5]}")
        else:
            print("âœ— Equilibrium failed to converge")
            print(f"  Message: {result2.message()}")

    except Exception as e:
        print(f"âœ— System 2 equilibration error: {e}")

print("\n" + "=" * 70)
print("CONCLUSION")
print("=" * 70)
print("""
If System 1 and System 2 give similar results â†’ H2O(aq) is just a dummy linker
If System 1 works but System 2 fails â†’ H2O(aq) thermodynamics are being used

The key question: Do we NEED H2O(aq) in the database, or is it just a
placeholder that Reaktoro replaces with IAPWS-95 automatically?
""")

