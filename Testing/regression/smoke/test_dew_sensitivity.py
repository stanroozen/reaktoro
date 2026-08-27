"""
Smoke test: EquilibriumSensitivity produces finite, non-zero derivatives
(dndw with respect to temperature) for both ActivityModelDEW and
ActivityModelPerplexDEW.

Validates task 8 of the DEW/PerplexDEW workflow integration roadmap.
"""

import sys
import math
import pytest
from pathlib import Path


# Auto-discover reaktoro4py.pyd (conftest.py does this, but be explicit for standalone execution)
def _setup_path():
    testing_root = Path(__file__).parent.parent.parent
    repo_root = testing_root.parent if testing_root.name == "Testing" else testing_root

    search_dirs = [
        repo_root / "temp_build" / "build-dew" / "Reaktoro" / "Release",
        repo_root / "build-msvc" / "Reaktoro" / "Release",
        repo_root / "build" / "Reaktoro" / "Release",
    ]
    for d in search_dirs:
        if d.exists() and str(d) not in sys.path:
            sys.path.insert(0, str(d))
            break


_setup_path()


def _import():
    try:
        import reaktoro4py as rkt

        return rkt
    except ImportError as e:
        pytest.skip(f"reaktoro4py not available: {e}")


def _make_system(rkt, activity_model):
    """Minimal Mg-OH system with Brucite mineral."""
    dew = rkt.DEWDatabase("dew2024-aqueous")
    supcrt = rkt.SupcrtDatabase("supcrtbl")
    db = rkt.Database(dew.species())
    db.addSpecies(supcrt.species("Brucite"))

    # Species names must match the active DEW database naming in this build.
    aq = rkt.AqueousPhase("H2O(aq) H+(aq) OH-(aq) Mg+2(aq) MgOH+(aq)")
    aq.setActivityModel(activity_model)
    return rkt.ChemicalSystem(db, aq, rkt.MineralPhase("Brucite"))


def _make_state(rkt, system):
    state = rkt.ChemicalState(system)
    for name, val in [
        ("H2O(aq)", 55.5),
        ("H+(aq)", 1e-8),
        ("OH-(aq)", 1e-8),
        ("Mg+2(aq)", 1e-5),
        ("Brucite", 5.0),
    ]:
        try:
            state.set(name, val, "mol")
        except TypeError:
            import autodiff

            state.set(name, autodiff.real(val), "mol")
    return state


@pytest.mark.parametrize(
    "model_name,model_factory",
    [
        ("DEW", lambda rkt: rkt.ActivityModelDEW()),
        ("PerplexDEW", lambda rkt: rkt.ActivityModelPerplexDEW()),
    ],
    ids=["DEW", "PerplexDEW"],
)
def test_sensitivity_dndw_temperature_is_finite(model_name, model_factory):
    """
    EquilibriumSensitivity.dndw("T") must be a finite matrix (no NaN/Inf)
    after one solver call for both DEW and PerplexDEW activity models.
    """
    rkt = _import()
    system = _make_system(rkt, model_factory(rkt))

    specs = rkt.EquilibriumSpecs(system)
    specs.temperature()
    specs.pressure()
    solver = rkt.EquilibriumSolver(specs)

    conds = rkt.EquilibriumConditions(specs)
    conds.temperature(250.0, "celsius")
    conds.pressure(1000.0, "bar")

    state = _make_state(rkt, system)
    # Pre-equilibrate without sensitivity so the starting point is well-conditioned.
    pre = solver.solve(state, conds)
    assert pre.succeeded(), f"Pre-equilibration failed for {model_name}"

    sensitivity = rkt.EquilibriumSensitivity()
    result = solver.solve(state, sensitivity, conds)
    assert result.succeeded(), f"Sensitivity solve failed for {model_name}"

    dndw_T = sensitivity.dndw("T")
    # dndw_T is a numpy-compatible matrix: rows = species, col = input "T"
    n_rows = len(dndw_T)
    assert n_rows == system.species().size(), (
        f"dndw('T') row count {n_rows} != species count"
    )
    for i, val in enumerate(dndw_T):
        assert math.isfinite(float(val)), (
            f"{model_name}: dndw('T')[{i}] = {float(val)} is not finite"
        )


@pytest.mark.parametrize(
    "model_name,model_factory",
    [
        ("DEW", lambda rkt: rkt.ActivityModelDEW()),
        ("PerplexDEW", lambda rkt: rkt.ActivityModelPerplexDEW()),
    ],
    ids=["DEW", "PerplexDEW"],
)
def test_sensitivity_dndw_pressure_is_finite(model_name, model_factory):
    """
    EquilibriumSensitivity.dndw("P") must be finite for both DEW and PerplexDEW.
    """
    rkt = _import()
    system = _make_system(rkt, model_factory(rkt))

    specs = rkt.EquilibriumSpecs(system)
    specs.temperature()
    specs.pressure()
    solver = rkt.EquilibriumSolver(specs)

    conds = rkt.EquilibriumConditions(specs)
    conds.temperature(250.0, "celsius")
    conds.pressure(1000.0, "bar")

    state = _make_state(rkt, system)
    solver.solve(state, conds)  # pre-condition

    sensitivity = rkt.EquilibriumSensitivity()
    result = solver.solve(state, sensitivity, conds)
    assert result.succeeded(), f"Sensitivity solve failed for {model_name}"

    dndw_P = sensitivity.dndw("P")
    for i, val in enumerate(dndw_P):
        assert math.isfinite(float(val)), (
            f"{model_name}: dndw('P')[{i}] = {float(val)} is not finite"
        )
