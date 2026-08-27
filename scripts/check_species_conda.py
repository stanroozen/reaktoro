from reaktoro import SupcrtDatabase, AggregateState

db = SupcrtDatabase("supcrtbl")


def elem_symbols(s):
    return [pair[0].symbol() for pair in s.elements()]


def aq(elem):
    return sorted(
        [
            s.name()
            for s in db.species()
            if elem in elem_symbols(s) and s.aggregateState() == AggregateState.Aqueous
        ]
    )


def cr(elem):
    return sorted(
        [
            s.name()
            for s in db.species()
            if elem in elem_symbols(s) and s.aggregateState() == AggregateState.Solid
        ]
    )


def gas():
    return sorted(
        [s.name() for s in db.species() if s.aggregateState() == AggregateState.Gas]
    )


print("=== Fe aqueous ===")
print(aq("Fe")[:20])
print("=== Fe minerals ===")
print(cr("Fe")[:15])
print("=== Ca minerals ===")
print(cr("Ca")[:12])
print("=== Ca aqueous ===")
print(aq("Ca")[:12])
c_aq = sorted(
    [
        s.name()
        for s in db.species()
        if "C" in elem_symbols(s)
        and "Ca" not in elem_symbols(s)
        and s.aggregateState() == AggregateState.Aqueous
    ]
)
print("=== C aqueous (no Ca) ===")
print(c_aq[:20])
print("=== Gases ===")
print(gas()[:15])
