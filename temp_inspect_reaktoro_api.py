import os
import sys

sys.path.insert(
    0, os.path.join(os.path.dirname(__file__), "build", "Reaktoro", "Release")
)

import reaktoro4py as r

print("EquilibriumSpecs methods containing key words:")
specs_names = [
    n
    for n in dir(r.EquilibriumSpecs)
    if any(
        k in n.lower()
        for k in [
            "potential",
            "fugacity",
            "activity",
            "charge",
            "amount",
            "element",
            "phase",
            "constraint",
            "eh",
            "pe",
        ]
    )
]
for name in sorted(specs_names):
    print(name)

print("\nChemicalProps methods containing key words:")
props_names = [
    n
    for n in dir(r.ChemicalProps)
    if any(
        k in n.lower()
        for k in [
            "potential",
            "fugacity",
            "activity",
            "amount",
            "element",
            "phase",
            "species",
            "charge",
            "mu",
        ]
    )
]
for name in sorted(props_names):
    print(name)

print("\nAqueousProps methods containing key words:")
aq_names = [
    n
    for n in dir(r.AqueousProps)
    if any(
        k in n.lower()
        for k in [
            "eh",
            "pe",
            "activity",
            "fugacity",
            "potential",
            "species",
            "element",
            "charge",
        ]
    )
]
for name in sorted(aq_names):
    print(name)
