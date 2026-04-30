#!/usr/bin/env python
from sys import path

path.insert(0, "build-msvc/python")
from reaktoro4py import *

print("Loading DEW database...")
dew_db = DEWDatabase("DEWDatabase")

print(f"\nTotal species in DEW: {len(dew_db.species())}")
print("\nAqueous species containing 'H' or 'O' (first 30):")

count = 0
for sp in dew_db.species():
    if sp.aggregateState() == AggregateState.Aqueous and count < 30:
        print(f"  {sp.name()}")
        count += 1

print("\nSpecies containing 'water':")
for sp in dew_db.species():
    if "water" in sp.name().lower():
        print(f"  {sp.name()}")

print("\nSpecies containing 'H2O':")
for sp in dew_db.species():
    if "H2O" in sp.name():
        print(f"  {sp.name()}")
