#!/usr/bin/env python
"""Test individual constraint types in Reaktoro EquilibriumSpecs."""

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

constraints_to_test = [
    ("temperature", lambda s: s.temperature()),
    ("pressure", lambda s: s.pressure()),
    ("pH", lambda s: s.pH()),
    ("lgActivity (SiO2,aq)", lambda s: s.lgActivity("SiO2,aq")),
    ("lnActivity (SiO2,aq)", lambda s: s.lnActivity("SiO2,aq")),
    ("activity (SiO2,aq)", lambda s: s.activity("SiO2,aq")),
    ("fugacity (O2)", lambda s: s.fugacity("O2")),
    ("chemicalPotential (SiO2,aq)", lambda s: s.chemicalPotential("SiO2,aq")),
    ("Eh (redox potential)", lambda s: s.Eh()),
    ("pE (electrons)", lambda s: s.pE()),
    ("enthalpy", lambda s: s.enthalpy()),
    ("entropy", lambda s: s.entropy()),
    ("gibbsEnergy", lambda s: s.gibbsEnergy()),
    ("volume", lambda s: s.volume()),
    ("phaseVolume (Wlm)", lambda s: s.phaseVolume("Wlm")),
]

print("Individual Constraint Type Testing")
print("=" * 70)

for name, func in constraints_to_test:
    specs = EquilibriumSpecs(sys)
    specs.temperature()
    specs.pressure()

    try:
        func(specs)
        print(f"✓ {name:<40} WORKS")
    except Exception as e:
        error_msg = str(e)[:50] if len(str(e)) > 50 else str(e)
        print(f"✗ {name:<40} FAILED: {error_msg}")
