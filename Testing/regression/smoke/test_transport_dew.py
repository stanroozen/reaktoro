"""
Smoke test: minimal operator-splitting reactive transport step completes
without NaN values and produces physically reasonable concentration gradients
for both ActivityModelDEW and ActivityModelPerplexDEW.

Since TransportSolver has no Python binding, this test implements a 3-cell
1D reactive-transport step by operator splitting:
  1. Transport (advection): shift fluid amounts one cell downstream.
  2. Reaction (equilibrium): re-equilibrate each cell at fixed T, P.

This exercises the ChemicalProps pipeline through DEW/PerplexDEW models in
exactly the same sequence a real transport solver would use.

Validates task 5 of the DEW/PerplexDEW workflow integration roadmap.
"""

import sys
import math
import copy
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
    """Simple Mg-OH aqueous system (no mineral) for clean transport test."""
    dew = rkt.DEWDatabase("dew2024-aqueous")
    # Species names must match the active DEW database naming in this build.
    aq = rkt.AqueousPhase("H2O(aq) H+(aq) OH-(aq) Mg+2(aq) MgOH+(aq)")
    aq.setActivityModel(activity_model)
    db = rkt.Database(dew.species())
    return rkt.ChemicalSystem(db, aq)


def _make_state(rkt, system, mg_mol: float):
    """Cell state with the given initial Mg amount (mol) at 200 C, 1 kbar."""
    state = rkt.ChemicalState(system)
    for name, val in [
        ("H2O(aq)", 55.5),
        ("H+(aq)", 1e-8),
        ("OH-(aq)", 1e-8),
        ("Mg+2(aq)", mg_mol),
        ("MgOH+(aq)", 1e-10),
    ]:
        try:
            state.set(name, val, "mol")
        except TypeError:
            import autodiff

            state.set(name, autodiff.real(val), "mol")
    return state


def _read_amounts(rkt, state, system) -> list:
    """Return a plain-float copy of the current species amount vector."""
    amounts = []
    for i in range(system.species().size()):
        amounts.append(float(state.speciesAmount(i)))
    return amounts


def _set_amounts(rkt, state, system, amounts: list):
    """Overwrite the species amount vector of a state."""
    for i, val in enumerate(amounts):
        state.setSpeciesAmount(i, max(val, 0.0))  # clamp negatives from transport


@pytest.mark.parametrize(
    "model_name,model_factory",
    [
        ("DEW", lambda rkt: rkt.ActivityModelDEW()),
        ("PerplexDEW", lambda rkt: rkt.ActivityModelPerplexDEW()),
    ],
    ids=["DEW", "PerplexDEW"],
)
def test_transport_no_nan(model_name, model_factory):
    """
    After 2 operator-splitting steps over 3 cells the species amounts must
    all be finite (no NaN or Inf from activity model evaluation).
    """
    rkt = _import()
    system = _make_system(rkt, model_factory(rkt))
    n_species = system.species().size()

    # Three cells with different initial Mg concentrations to create a gradient.
    mg_init = [1e-3, 1e-5, 1e-7]
    cells = [_make_state(rkt, system, mg) for mg in mg_init]

    specs = rkt.EquilibriumSpecs(system)
    specs.temperature()
    specs.pressure()
    solver = rkt.EquilibriumSolver(specs)
    conds = rkt.EquilibriumConditions(specs)
    conds.temperature(200.0, "celsius")
    conds.pressure(1000.0, "bar")

    for transport_step in range(2):
        # --- Reaction step: equilibrate each cell ---
        for cell_idx, cell in enumerate(cells):
            result = solver.solve(cell, conds)
            assert result.succeeded(), (
                f"{model_name}: step {transport_step + 1} cell {cell_idx} eq failed"
            )

        # --- Transport step: advect (shift amounts one cell downstream) ---
        amounts = [_read_amounts(rkt, c, system) for c in cells]
        alpha = 0.3  # Courant number (fraction of cell transferred per step)
        new_amounts = [list(a) for a in amounts]
        for cell_idx in range(1, len(cells)):
            for sp in range(n_species):
                flux = alpha * amounts[cell_idx - 1][sp]
                new_amounts[cell_idx][sp] += flux
                new_amounts[cell_idx - 1][sp] -= flux
        for cell_idx, cell in enumerate(cells):
            _set_amounts(rkt, cell, system, new_amounts[cell_idx])

    # Verify no NaN/Inf in any cell
    for cell_idx, cell in enumerate(cells):
        for sp in range(n_species):
            val = float(cell.speciesAmount(sp))
            assert math.isfinite(val), (
                f"{model_name}: cell {cell_idx} species {sp} = {val} not finite "
                f"after transport"
            )


@pytest.mark.parametrize(
    "model_name,model_factory",
    [
        ("DEW", lambda rkt: rkt.ActivityModelDEW()),
        ("PerplexDEW", lambda rkt: rkt.ActivityModelPerplexDEW()),
    ],
    ids=["DEW", "PerplexDEW"],
)
def test_transport_concentration_gradient_decreases(model_name, model_factory):
    """
    After operator-splitting transport the Mg concentration gradient between
    the inlet and outlet cell must be smaller than the initial gradient —
    confirming that advection is mixing the cells and the activity model
    does not produce pathological chemistry.
    """
    rkt = _import()
    system = _make_system(rkt, model_factory(rkt))

    mg_init = [1e-3, 1e-5, 1e-8]
    cells = [_make_state(rkt, system, mg) for mg in mg_init]

    specs = rkt.EquilibriumSpecs(system)
    specs.temperature()
    specs.pressure()
    solver = rkt.EquilibriumSolver(specs)
    conds = rkt.EquilibriumConditions(specs)
    conds.temperature(200.0, "celsius")
    conds.pressure(1000.0, "bar")

    def mg_molality(cell):
        props = rkt.AqueousProps(cell)
        return float(props.elementMolality("Mg"))

    # Initial gradient (before any transport)
    grad_before = mg_molality(cells[0]) - mg_molality(cells[-1])

    # 4 operator-splitting steps
    n_species = system.species().size()
    alpha = 0.3
    for _ in range(4):
        for cell in cells:
            solver.solve(cell, conds)
        amounts = [_read_amounts(rkt, c, system) for c in cells]
        new_amounts = [list(a) for a in amounts]
        for i in range(1, len(cells)):
            for sp in range(n_species):
                flux = alpha * amounts[i - 1][sp]
                new_amounts[i][sp] += flux
                new_amounts[i - 1][sp] -= flux
        for i, cell in enumerate(cells):
            _set_amounts(rkt, cell, system, new_amounts[i])

    grad_after = mg_molality(cells[0]) - mg_molality(cells[-1])

    assert grad_after < grad_before, (
        f"{model_name}: Mg gradient did not decrease after 4 transport steps "
        f"(before={grad_before:.2e}, after={grad_after:.2e})"
    )
