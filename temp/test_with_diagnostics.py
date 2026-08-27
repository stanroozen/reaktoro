#!/usr/bin/env python3
import os
import sys

# Set diagnostic flags BEFORE importing reaktoro
os.environ["REAKTORO_FUGACITY_DIAGNOSTICS"] = "1"
os.environ["REAKTORO_FUGACITY_SKIP_STANDARD_THERMO"] = "1"

# Now import and run reproducer
sys.path.insert(
    0,
    r"C:\Users\stanroozen\Documents\Projects\reaktoro-dev\reaktoro\build\python\package\build\lib",
)
import reaktoro as m

print("Testing O2/SO2 gas mixture with diagnostic flags...")
try:
    gas = m.GaseousPhase("O2 SO2")
    print("GaseousPhase created successfully")

    db = m.Database("supcrt07.yaml")
    aq = m.AqueousPhase(m.speciate("H+ Cl-"))
    mineral = m.MineralPhase("Halite")

    system = m.ChemicalSystem(db, aq, mineral, gas)
    print("ChemicalSystem created successfully")

    specs = m.EquilibriumSpecs(system)
    specs.temperature()
    specs.pressure()
    print("EquilibriumSpecs created successfully")

except Exception as e:
    print(f"Error: {e}")
    import traceback

    traceback.print_exc()
