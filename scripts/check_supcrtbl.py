import sys, os

# The pyd lives in the same folder as Reaktoro.dll
PYD_DIR = r"C:\Users\stanroozen\Documents\Projects\reaktoro-dev\reaktoro\build\python\package\build\lib\reaktoro"
PY_PKG = r"C:\Users\stanroozen\Documents\Projects\reaktoro-dev\reaktoro\python\package"
sys.path.insert(
    0, os.path.dirname(PYD_DIR)
)  # .../build/lib  -- so 'reaktoro' is a package
sys.path.insert(0, PYD_DIR)  # .../reaktoro   -- direct import of reaktoro4py.pyd
sys.path.insert(0, PY_PKG)
from reaktoro4py import *

db = SupcrtDatabase("supcrtbl")


def aq(element_sym):
    return sorted(
        [
            s.name()
            for s in db.species()
            if element_sym in [e.symbol() for e in s.elements()]
            and s.aggregateState() == AggregateState.Aqueous
        ]
    )


def cr(element_sym):
    return sorted(
        [
            s.name()
            for s in db.species()
            if element_sym in [e.symbol() for e in s.elements()]
            and s.aggregateState() == AggregateState.Solid
        ]
    )


print("=== Fe aqueous ===")
print(aq("Fe")[:20])
print("=== Fe minerals ===")
print(cr("Fe")[:15])
print("=== Ca minerals ===")
print(cr("Ca")[:12])
print("=== Ca aqueous ===")
print(aq("Ca")[:12])
print("=== C aqueous (no Ca) ===")
c_aq = sorted(
    [
        s.name()
        for s in db.species()
        if "C" in [e.symbol() for e in s.elements()]
        and "Ca" not in [e.symbol() for e in s.elements()]
        and s.aggregateState() == AggregateState.Aqueous
    ]
)
print(c_aq[:20])
print("=== O2 gas ===")
gas = [
    s.name()
    for s in db.species()
    if s.aggregateState() == AggregateState.Gas
    and "O" in [e.symbol() for e in s.elements()]
]
print(sorted(gas)[:10])
