#!/usr/bin/env python
"""Test what constraint types are available in Reaktoro EquilibriumSpecs."""

import sys
import os

sys.path.insert(
    0, os.path.join(os.path.dirname(__file__), "build", "Reaktoro", "Release")
)

from reaktoro4py import *

# Create a simple test system
db = Database.fromFile(
    os.path.join(
        os.path.dirname(__file__),
        "embedded",
        "databases",
        "perplex",
        "DEW17HP622_Zn_2025-reaktoro.json",
    )
)

aq = AqueousPhase(["H2O", "H+", "OH-", "SiO2,aq", "HS-"])
aq.setActivityModel(ActivityModelPerplexDEW())

mineral = MineralPhase("Wlm")
gas = GaseousPhase(["O2"])
gas.setActivityModel(ActivityModelIdealGas())

sys = ChemicalSystem(db, aq, mineral, gas)

specs = EquilibriumSpecs(sys)

# List all constraint types by checking what methods exist
print("Available EquilibriumSpecs constraint methods:")
print("=" * 60)

constraint_methods = [
    m for m in dir(specs) if not m.startswith("_") and callable(getattr(specs, m))
]

# Filter for likely constraint methods
likely_constraints = [
    m
    for m in constraint_methods
    if any(
        x in m.lower()
        for x in [
            "temp",
            "press",
            "activity",
            "potential",
            "fugacity",
            "saturation",
            "affinity",
            "gibbs",
            "enthalpy",
            "entropy",
            "volume",
            "pH",
            "pe",
            "eh",
        ]
    )
]

for method in sorted(likely_constraints):
    print(f"  - {method}")

print("\n" + "=" * 60)
print("Testing specific constraints:")

# Test basic constraints
try:
    specs.temperature()
    print("✓ specs.temperature() works")
except Exception as e:
    print(f"✗ specs.temperature() failed: {e}")

try:
    specs.pressure()
    print("✓ specs.pressure() works")
except Exception as e:
    print(f"✗ specs.pressure() failed: {e}")

try:
    specs.pH()
    print("✓ specs.pH() works")
except Exception as e:
    print(f"✗ specs.pH() failed: {e}")

try:
    specs.lgActivity("H+")
    print("✓ specs.lgActivity() works")
except Exception as e:
    print(f"✗ specs.lgActivity() failed: {e}")

try:
    specs.fugacity("O2")
    print("✓ specs.fugacity() works")
except Exception as e:
    print(f"✗ specs.fugacity() failed: {e}")

try:
    specs.chemicalPotential("H+")
    print("✓ specs.chemicalPotential() works")
except Exception as e:
    print(f"✗ specs.chemicalPotential() failed: {e}")

try:
    specs.saturationIndex("Wlm")
    print("✓ specs.saturationIndex() works")
except Exception as e:
    print(f"✗ specs.saturationIndex() failed: {e}")

try:
    specs.affinity("Wlm")
    print("✓ specs.affinity() works")
except Exception as e:
    print(f"✗ specs.affinity() failed: {e}")

try:
    specs.activity("SiO2,aq")
    print("✓ specs.activity() works (without log)")
except Exception as e:
    print(f"✗ specs.activity() failed: {e}")
