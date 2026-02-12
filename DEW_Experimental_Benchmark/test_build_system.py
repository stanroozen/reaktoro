"""Quick test of build_system function"""

import sys
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(SCRIPT_DIR)
PYD_DIR = os.path.join(ROOT_DIR, "build-msvc", "Reaktoro", "Release")
if os.path.isdir(PYD_DIR) and PYD_DIR not in sys.path:
    sys.path.insert(0, PYD_DIR)

from reaktoro4py import *

print("Loading databases...")
dew_db = DEWDatabase("dew2019-aqueous")
supcrt_db = SupcrtDatabase("supcrtbl")

print("Building system...")
from quartz_solubility_analysis import build_system

try:
    system = build_system(dew_db, supcrt_db)
    print(f"[OK] System created with {len(system.species())} species")
    print(f"[OK] Test passed!")
except Exception as e:
    print(f"[FAIL] Error: {e}")
    import traceback

    traceback.print_exc()
