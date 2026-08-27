#!/usr/bin/env python
"""Debug script to check species names in DEW17HP622_Zn database."""

import sys
import os

sys.path.insert(
    0, os.path.join(os.path.dirname(__file__), "build", "Reaktoro", "Release")
)

from reaktoro4py import *

db = Database.fromFile(
    os.path.join(
        os.path.dirname(__file__),
        "embedded",
        "databases",
        "perplex",
        "DEW17HP622_Zn_2025-reaktoro.json",
    )
)

aq = AqueousPhase(["SiO2,aq", "HS-", "H+", "H2O"])
aq.setActivityModel(ActivityModelPerplexDEW())

mineral = MineralPhase("Wlm")

sys = ChemicalSystem(db, aq, mineral)

print("Species in system:")
for i in range(sys.species().size()):
    print(f'  {i}: "{sys.species(i).name()}"')
