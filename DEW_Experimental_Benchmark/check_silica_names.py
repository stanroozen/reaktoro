"""Find correct SiO2 species names in DEW database"""

import sys
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(SCRIPT_DIR)
PYD_DIR = os.path.join(ROOT_DIR, "build", "Reaktoro", "Release")
if os.path.isdir(PYD_DIR) and PYD_DIR not in sys.path:
    sys.path.insert(0, PYD_DIR)

from reaktoro4py import *

db = DEWDatabase("dew2019-aqueous")
species_names = [s.name() for s in db.species()]

# Find all silica species
sio2_species = [s for s in species_names if "SiO" in s or "Si" in s]
print(f"Silica-containing species in DEW ({len(sio2_species)}):")
for s in sorted(sio2_species):
    print(f"  - {s}")

