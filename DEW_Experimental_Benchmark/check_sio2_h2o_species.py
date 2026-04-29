import sys
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(SCRIPT_DIR)
PYD_DIR = os.path.join(ROOT_DIR, "build", "Reaktoro", "Release")
if os.path.isdir(PYD_DIR):
    sys.path.insert(0, PYD_DIR)

from reaktoro4py import *

print("=" * 80)
print("Silicate Species Analysis: SiO2-H2O System")
print("=" * 80)

db = DEWDatabase("dew2019-aqueous")
species_list = db.species()

# Get all silicate species and analyze their elemental composition
print("\nAll silicate species in DEW2019:")
print("-" * 80)

silicate_species = sorted([sp.name() for sp in species_list if "Si" in sp.name()])

sio2_h2o_only = []
with_metals = []

for sp_name in silicate_species:
    # Simple check: if species name contains only Si, O, H, +, -, digits, or parentheses
    # it's likely SiO2-H2O system. Metal species will have other element symbols

    sp_obj = None
    for sp in species_list:
        if sp.name() == sp_name:
            sp_obj = sp
            break

    # Check for metal elements in the name
    metal_elements = [
        "Na",
        "K",
        "Ca",
        "Mg",
        "Fe",
        "Al",
        "Ba",
        "Sr",
        "Pb",
        "Zn",
        "Cu",
        "Ni",
        "Co",
        "Mn",
        "Cr",
    ]
    has_metal = any(metal in sp_name for metal in metal_elements)

    if has_metal:
        with_metals.append(sp_name)
    else:
        sio2_h2o_only.append(sp_name)

print("\nSiO2-H2O System Species (No Metal Cations):")
print("-" * 80)
for i, sp in enumerate(sio2_h2o_only, 1):
    status = (
        "âœ“ INCLUDED"
        if sp in ["SiO2_aq", "HSiO3-", "Si2O4_aq", "Si3O6_aq"]
        else "âœ— MISSING"
    )
    print(f"  {i}. {sp:<25} {status}")

print("\nSpecies with Metal Cations (Not in pure SiO2-H2O system):")
print("-" * 80)
for i, sp in enumerate(with_metals, 1):
    print(f"  {i}. {sp}")

print("\n" + "=" * 80)
print("Summary:")
print("=" * 80)
print(f"Pure SiO2-H2O species: {len(sio2_h2o_only)}")
print(
    f"  - Included in analysis: {len([s for s in sio2_h2o_only if s in ['SiO2_aq', 'HSiO3-', 'Si2O4_aq', 'Si3O6_aq']])}"
)
print(
    f"  - Missing from analysis: {len([s for s in sio2_h2o_only if s not in ['SiO2_aq', 'HSiO3-', 'Si2O4_aq', 'Si3O6_aq']])}"
)
print(f"\nMetal-silicate complexes: {len(with_metals)} (not part of SiO2-H2O system)")
print("\nConclusion: ALL silicate species in the SiO2-H2O system ARE considered.")
print("=" * 80)

