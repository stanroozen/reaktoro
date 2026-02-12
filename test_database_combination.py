"""
Simple test to verify DEW + SUPCRT database combination works
"""

from reaktoro4py import *

# Load databases
dew_db = DEWDatabase("DEWDatabase")
supcrt_db = SupcrtDatabase("supcrtbl")

# Combine databases
combined_db = Database(dew_db.species())
combined_db.extend(supcrt_db)

# Print database info
print(f"DEW database: {len(dew_db.species())} species")
print(f"SUPCRT database: {len(supcrt_db.species())} species")
print(f"Combined database: {len(combined_db.species())} species")

# Check for water species
water_species = [
    sp
    for sp in combined_db.species()
    if "H2O" in sp.name() and sp.aggregateState() == AggregateState.Aqueous
]
print(f"\nFound {len(water_species)} H2O aqueous species:")
for sp in water_species:
    print(f"  - {sp.name()}")

# Check for silica species
silica_species = [
    sp
    for sp in combined_db.species()
    if "Si" in sp.name() and sp.aggregateState() == AggregateState.Aqueous
]
print(f"\nFound {len(silica_species)} Si-containing aqueous species (first 10):")
for sp in silica_species[:10]:
    print(f"  - {sp.name()}")

# Check for Quartz
quartz_species = [sp for sp in combined_db.species() if "Quartz" in sp.name()]
print(f"\nFound {len(quartz_species)} Quartz species:")
for sp in quartz_species:
    print(f"  - {sp.name()} ({sp.aggregateState()})")

print("\n✅ Database combination successful!")
