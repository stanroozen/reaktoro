import sys
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(SCRIPT_DIR)
PYD_DIR = os.path.join(ROOT_DIR, "build", "Reaktoro", "Release")
if os.path.isdir(PYD_DIR):
    sys.path.insert(0, PYD_DIR)

from reaktoro4py import *

print("=" * 80)
print("Available Silicate Species in DEW2019")
print("=" * 80)

db = DEWDatabase("dew2019-aqueous")
species_list = db.species()
silicate_species = sorted([sp.name() for sp in species_list if "Si" in sp.name()])

print(f"\nFound {len(silicate_species)} silicate species:\n")
for i, sp in enumerate(silicate_species, 1):
    print(f"  {i:2d}. {sp}")

print("\n" + "=" * 80)
print("Currently used in quartz_solubility_analysis:")
print("=" * 80)
current = [
    "WATER,AQ",
    "H+",
    "OH-",
    "SiO2_aq",
    "H2_aq",
    "O2_aq",
    "HO2-",
    "HSiO3-",
    "Si2O4_aq",
    "Si3O6_aq",
]
print("\nIncluded:")
for sp in current:
    if "Si" in sp:
        print(f"  âœ“ {sp}")

print("\nNot included:")
for sp in silicate_species:
    if sp not in current:
        print(f"  âœ— {sp}")

