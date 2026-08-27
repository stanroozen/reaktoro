"""Verify legacy and modern EquilibriumSolver patterns yield comparable results."""

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


def test_modern_pattern_matches_legacy() -> None:
    dew_db = DEWDatabase("dew2024-aqueous")
    supcrt_db = SupcrtDatabase("supcrtbl")

    quartz_species = supcrt_db.species("Quartz")
    combined_db = Database(dew_db.species())
    combined_db.addSpecies(quartz_species)

    aqueous = AqueousPhase("H2O(aq) H+(aq) OH-(aq) SiO2(aq)")
    aqueous.setActivityModel(ActivityModelDEW())
    mineral = MineralPhase("Quartz")
    system = ChemicalSystem(combined_db, aqueous, mineral)

    solver_old = EquilibriumSolver(system)
    state_old = ChemicalState(system)
    _set_amount(state_old, "H2O(aq)", 1.0, "kg")
    _set_amount(state_old, "Quartz", 10.0, "mol")
    cond_old = EquilibriumConditions(system)
    _set_condition_tp(cond_old, 300.0, 1000.0)
    res_old = solver_old.solve(state_old, cond_old)
    assert res_old.succeeded()
    sol_old = float(AqueousProps(state_old).speciesMolality("SiO2(aq)"))

    specs = EquilibriumSpecs(system)
    specs.temperature()
    specs.pressure()
    solver_new = EquilibriumSolver(specs)
    cond_new = EquilibriumConditions(specs)
    _set_condition_tp(cond_new, 300.0, 1000.0)
    state_new = ChemicalState(system)
    _set_amount(state_new, "H2O(aq)", 1.0, "kg")
    _set_amount(state_new, "Quartz", 10.0, "mol")
    res_new = solver_new.solve(state_new, cond_new)
    assert res_new.succeeded()
    sol_new = float(AqueousProps(state_new).speciesMolality("SiO2(aq)"))

    assert sol_old > 0.0
    assert sol_new > 0.0
    rel_diff = abs(sol_old - sol_new) / max(sol_old, 1e-30)
    assert rel_diff < 1e-3
