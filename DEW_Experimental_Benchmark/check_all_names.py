"""Find correct species names in DEW database"""

import sys
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(SCRIPT_DIR)
PYD_DIR = os.path.join(ROOT_DIR, "build-msvc", "Reaktoro", "Release")
if os.path.isdir(PYD_DIR) and PYD_DIR not in sys.path:
    sys.path.insert(0, PYD_DIR)

from reaktoro4py import *

db = DEWDatabase("dew2019-aqueous")
species_names = [s.name() for s in db.species()]

# Find specific species
searches = ["H2", "O2", "HSiO", "Si2O", "Si3O"]
for search in searches:
    matches = [s for s in species_names if search in s]
    print(f'Species with "{search}":')
    for s in sorted(matches):
        print(f"  - {s}")
    print()
