import os
import sys

# Try to import Reaktoro; fall back to local extension module if not installed
try:
    from reaktoro import *
except ModuleNotFoundError:
    # Add local build path for reaktoro4py if available
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    PYD_DIR = os.path.join(SCRIPT_DIR, "build-msvc", "Reaktoro", "Release")
    if os.path.isdir(PYD_DIR) and PYD_DIR not in sys.path:
        sys.path.insert(0, PYD_DIR)
    try:
        from reaktoro4py import *

        print("Using local reaktoro4py extension from build-msvc.")
    except ModuleNotFoundError as e:
        raise ModuleNotFoundError(
            "Could not import Reaktoro. Install the 'reaktoro' package or ensure reaktoro4py is on PYTHONPATH."
        ) from e

# Load DEW2019 database
db = DEWDatabase("dew2019-aqueous")
species = db.species()

print(f"Total species in DEW2019: {len(species)}")
print("\nSilica-containing species:")

silica_species = []
for s in species:
    formula_str = str(s.formula())
    if "Si" in formula_str:
        silica_species.append(s)
        print(f"  Name: {s.name():<20} Formula: {formula_str}")

print(f"\nTotal silica species: {len(silica_species)}")

# Check for specific names
print("\nChecking for specific names:")
for name in ["SiO2,aq", "SiO2_aq", "H4SiO4(aq)", "HSiO3-", "Si2O4,aq", "Si3O6,aq"]:
    try:
        sp = species.get(name)
        print(f"  {name:<20} -> FOUND (formula: {sp.formula()})")
    except:
        print(f"  {name:<20} -> NOT FOUND")
