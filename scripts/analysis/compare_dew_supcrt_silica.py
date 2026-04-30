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

# Load SUPCRT database
print("=" * 70)
print("SUPCRT Database")
print("=" * 70)
supcrt_db = SupcrtDatabase("supcrtbl")
supcrt_species = supcrt_db.species()

print(f"Total species in SUPCRT: {len(supcrt_species)}")
print("\nSilica-containing aqueous species:")

silica_supcrt = []
for s in supcrt_species:
    formula_str = str(s.formula())
    if "Si" in formula_str and s.aggregateState() == AggregateState.Aqueous:
        silica_supcrt.append(s)
        print(f"  Name: {s.name():<25} Formula: {formula_str}")

print(f"\nTotal silica aqueous species in SUPCRT: {len(silica_supcrt)}")

# Check for specific names
print("\nChecking for specific names in SUPCRT:")
for name in ["SiO2(aq)", "H4SiO4(aq)", "HSiO3-", "H3SiO4-", "H2SiO4-2"]:
    try:
        sp = supcrt_species.get(name)
        print(f"  {name:<25} -> FOUND (formula: {str(sp.formula())})")
    except:
        print(f"  {name:<25} -> NOT FOUND")

print("\n" + "=" * 70)
print("DEW2019 Database")
print("=" * 70)

# Load DEW2019 database
dew_db = DEWDatabase("dew2019-aqueous")
dew_species = dew_db.species()

print(f"Total species in DEW2019: {len(dew_species)}")
print("\nSilica-containing species:")

silica_dew = []
for s in dew_species:
    formula_str = str(s.formula())
    if "Si" in formula_str:
        silica_dew.append(s)
        print(f"  Name: {s.name():<25} Formula: {formula_str}")

print(f"\nTotal silica species in DEW: {len(silica_dew)}")

# Check for specific names
print("\nChecking for specific names in DEW:")
for name in ["SiO2_aq", "H4SiO4(aq)", "HSiO3-", "Si2O4_aq", "Si3O6_aq"]:
    try:
        sp = dew_species.get(name)
        print(f"  {name:<25} -> FOUND (formula: {str(sp.formula())})")
    except:
        print(f"  {name:<25} -> NOT FOUND")

print("\n" + "=" * 70)
print("KEY FINDINGS")
print("=" * 70)
print("\nDEW uses different naming conventions:")
print("  - Commas are replaced with underscores: 'SiO2,aq' -> 'SiO2_aq'")
print("  - DEW has: SiO2_aq, HSiO3-, Si2O4_aq, Si3O6_aq")
print("  - DEW does NOT have: H4SiO4(aq)")
print("\nSUPCRT uses standard naming:")
print("  - SUPCRT has: H4SiO4(aq), SiO2(aq), HSiO3-, H3SiO4-, H2SiO4-2")
print("\nFor quartz solubility modeling with DEW:")
print("  - Use 'SiO2_aq' instead of 'H4SiO4(aq)'")
print("  - SiO2_aq is the neutral dissolved silica species")
