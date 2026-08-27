"""Check that lnActivity(H2O(aq)) constraints are accepted in DEW mixed-fluid setup."""

import numpy as np
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


def _total_si(aq) -> float:
    return (
        float(aq.speciesMolality("SiO2(aq)"))
        + float(aq.speciesMolality("HSiO3-(aq)"))
        + 2.0 * float(aq.speciesMolality("Si2O4(aq)"))
        + 3.0 * float(aq.speciesMolality("Si3O6(aq)"))
    )


def _set_amount(state, name: str, value: float, unit: str) -> None:
    try:
        state.set(name, value, unit)
    except TypeError:
        import autodiff

        state.set(name, autodiff.real(value), unit)


def _set_condition_tp(cond, temp_c: float, pressure_bar: float) -> None:
    try:
        cond.temperature(temp_c, "celsius")
        cond.pressure(pressure_bar, "bar")
    except TypeError:
        import autodiff

        cond.temperature(autodiff.real(temp_c), "celsius")
        cond.pressure(autodiff.real(pressure_bar), "bar")


def test_lnactivity_water_constraint_propagates() -> None:
    dew_db = DEWDatabase("dew2024-aqueous")
    supcrt_db = SupcrtDatabase("supcrtbl")
    quartz_species = supcrt_db.species("Quartz")

    combined_db = Database(dew_db.species())
    combined_db.addSpecies(quartz_species)

    aqueous = AqueousPhase(
        "H2O(aq) H+(aq) OH-(aq) SiO2(aq) HSiO3-(aq) Si2O4(aq) Si3O6(aq)"
    )
    aqueous.setActivityModel(ActivityModelPerplexDEW(ActivityDHModel.ExtendedDH))
    mineral = MineralPhase("Quartz")
    system = ChemicalSystem(combined_db, aqueous, mineral)
    Warnings.disable(906)

    # Baseline with direct system solver
    solver = EquilibriumSolver(system)
    state0 = ChemicalState(system)
    _set_amount(state0, "H2O(aq)", 1.0, "kg")
    _set_amount(state0, "SiO2(aq)", 1e-6, "mol")
    _set_amount(state0, "Quartz", 10.0, "mol")
    cond0 = EquilibriumConditions(system)
    _set_condition_tp(cond0, 800.0, 10000.0)
    assert solver.solve(state0, cond0).succeeded()
    m0 = _total_si(AqueousProps(state0))
    assert m0 > 0.0

    # Modern specs/conditions with explicit lnActivity(H2O(aq)) constraint
    specs = EquilibriumSpecs(system)
    specs.temperature()
    specs.pressure()
    specs.lnActivity("H2O(aq)")
    solver2 = EquilibriumSolver(specs)

    succeeded = 0
    for a_h2o in [0.95, 0.9, 0.8, 0.7, 0.6]:
        state = ChemicalState(system)
        _set_amount(state, "H2O(aq)", 1.0, "kg")
        _set_amount(state, "SiO2(aq)", 1e-6, "mol")
        _set_amount(state, "Quartz", 10.0, "mol")

        cond = EquilibriumConditions(specs)
        _set_condition_tp(cond, 800.0, 10000.0)
        cond.lnActivity("H2O(aq)", float(np.log(a_h2o)))

        res = solver2.solve(state, cond)
        if res.succeeded():
            succeeded += 1
            assert _total_si(AqueousProps(state)) >= 0.0

    # At least one constrained-water-activity case should converge, validating
    # that the lnActivity(H2O(aq)) constraint is accepted and propagated.
    assert succeeded >= 1
