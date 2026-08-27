"""Verify DEW + SUPCRT database combination works with current embedded names."""

import pytest

try:
    from reaktoro import *  # noqa: F401,F403
except Exception as exc:
    try:
        import reaktoro4py as _rkt

        globals().update(
            {
                name: getattr(_rkt, name)
                for name in dir(_rkt)
                if not name.startswith("_")
            }
        )
    except Exception:
        pytest.skip(f"reaktoro import failed: {exc}", allow_module_level=True)


def test_database_combination() -> None:
    dew_db = DEWDatabase("dew2024-aqueous")
    supcrt_db = SupcrtDatabase("supcrtbl")

    combined_db = Database(dew_db.species())
    combined_db.extend(supcrt_db)

    assert len(dew_db.species()) > 0
    assert len(supcrt_db.species()) > 0
    assert len(combined_db.species()) >= len(dew_db.species())

    water_species = [
        sp
        for sp in combined_db.species()
        if "H2O" in sp.name() and sp.aggregateState() == AggregateState.Aqueous
    ]
    assert len(water_species) > 0

    quartz_species = [sp for sp in combined_db.species() if "Quartz" in sp.name()]
    assert len(quartz_species) > 0
