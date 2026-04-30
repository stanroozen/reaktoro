from reaktoro4py import *

db = DEWDatabase("DEWDatabase")
water_species = [
    s.name() for s in db.species() if "H2O" in s.formula() or "H2O" in s.name()
]
print("Water species in DEW database:")
for name in water_species:
    print(f"  {name}")
