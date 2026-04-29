import sys
import os

# Add local build path
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BENCHMARK_DIR = os.path.dirname(SCRIPT_DIR)
ROOT_DIR = os.path.dirname(BENCHMARK_DIR)
PYD_DIR = os.path.join(ROOT_DIR, "build", "Reaktoro", "Release")
if os.path.isdir(PYD_DIR) and PYD_DIR not in sys.path:
    sys.path.insert(0, PYD_DIR)

from reaktoro4py import *

db_dew = DEWDatabase("dew2024-aqueous")

print("Aluminum species in DEW database:")
for species in db_dew.species():
    if "Al" in species.name():
        print(f"  {species.name()}")

