"""Check what the H2O species is actually called"""

import sys
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(SCRIPT_DIR)
PYD_DIR = os.path.join(ROOT_DIR, "build", "Reaktoro", "Release")
if os.path.isdir(PYD_DIR) and PYD_DIR not in sys.path:
    sys.path.insert(0, PYD_DIR)

from reaktoro4py import *

db = DEWDatabase("dew2019-aqueous")

# Try to get H2O_aq by its key
try:
    h2o = db.species("H2O_aq")
    print(f"âœ“ Found H2O_aq:")
    print(f"  - name(): {h2o.name()}")
    print(f"  - formula(): {h2o.formula()}")
except Exception as e:
    print(f"âœ— H2O_aq lookup failed: {e}")

# Also try WATER,AQ
try:
    h2o = db.species("WATER,AQ")
    print(f"âœ“ Found WATER,AQ:")
    print(f"  - name(): {h2o.name()}")
    print(f"  - formula(): {h2o.formula()}")
except Exception as e:
    print(f"âœ— WATER,AQ lookup failed: {e}")

# List first 5 species and their names
species_list = db.species()
print(f"\nFirst 5 species in database:")
for i, s in enumerate(species_list[:5]):
    print(f"  [{i}] name()={s.name()}, formula()={s.formula()}")

