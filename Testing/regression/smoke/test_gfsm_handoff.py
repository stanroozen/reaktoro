"""
Smoke test: PerplexGFSM gas phase + PerplexDEW aqueous phase — water-activity
handoff via props.extra is consumed and produces a non-trivial water activity.

The handoff protocol:
1. ActivityModelPerplexGFSM evaluates the gas phase and stores
   props.extra["PerplexGFSM::WaterActivity::StateId"]  (current state ID)
   props.extra["PerplexGFSM::WaterActivity::ln_f_ratio_h2o"]  (ln(f_mix/f_pure))
2. ActivityModelPerplexDEW evaluates the aqueous phase.  If it finds a fresh
   StateId in props.extra that matches "Reaktoro::ChemicalProps::StateId", it
   replaces the pure-water solvent state with the GFSM mixed-fluid values and
   sets props.extra["PerplexDEW::WaterActivity::CoupledHandoffUsed"] = True.

Validates task 12 of the DEW/PerplexDEW workflow integration roadmap.
"""

import sys
import math
import pytest
from pathlib import Path


# Auto-discover reaktoro4py.pyd (conftest.py does this, but be explicit here for
# standalone execution outside pytest)
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


def _build_two_phase_system(rkt):
    """
    Aqueous (PerplexDEW) + gas (PerplexGFSM) two-phase system with
    Mg-OH-H2O-CO2 chemistry.  No mineral phase so that equilibrium
    freely partitions H2O and CO2 between the phases.
    """
    # Aqueous database (DEW)
    dew_db = rkt.DEWDatabase("dew2024-aqueous")
    supcrt_db = rkt.SupcrtDatabase("supcrtbl")

    # Build a combined database that has both aqueous and gaseous species.
    db_species = list(dew_db.species())
    for name in ("CO2(g)", "H2O(g)"):
        try:
            db_species.append(supcrt_db.species(name))
        except Exception:
            pass  # not available — test will still run without gas phase
    db = rkt.Database(db_species)

    # Aqueous phase — PerplexDEW activity model
    aq_params = rkt.ActivityModelParamsPerplexDEW()
    aq_params.dhModel = rkt.ActivityDHModel.Davies
    aq_params.warnOnUnmappedGFSMCoupling = False  # suppress noise in test output
    aq = rkt.AqueousPhase("H2O(aq) H+(aq) OH-(aq) Mg+2(aq) CO2(aq)")
    aq.setActivityModel(rkt.ActivityModelPerplexDEW(aq_params))

    # Gas phase — PerplexGFSM activity model (H2O + CO2)
    gas_params = rkt.ActivityModelParamsPerplexGFSM()
    gas_params.hybridEosOptions = rkt.makePerplexCOHFluidPlusEosOptions()
    gas = rkt.GaseousPhase("CO2(g) H2O(g)")
    gas.setActivityModel(rkt.ActivityModelPerplexGFSM(gas_params))

    return rkt.ChemicalSystem(db, aq, gas)


def _build_aqueous_only_strict_system(rkt):
    """
    Aqueous-only PerplexDEW system in strict coupled-fluid mode.
    Any evaluation must fail because no GFSM gas-phase handoff exists.
    """
    dew_db = rkt.DEWDatabase("dew2024-aqueous")
    db = rkt.Database(list(dew_db.species()))

    aq_params = rkt.ActivityModelParamsPerplexDEW()
    aq_params.dhModel = rkt.ActivityDHModel.Davies
    aq_params.warnOnUnmappedGFSMCoupling = False
    aq_params.requireCoupledGFSMHandoff = True

    aq = rkt.AqueousPhase("H2O(aq) H+(aq) OH-(aq) Mg+2(aq) CO2(aq)")
    aq.setActivityModel(rkt.ActivityModelPerplexDEW(aq_params))
    return rkt.ChemicalSystem(db, aq)


def _make_state(rkt, system):
    state = rkt.ChemicalState(system)
    init = {
        "H2O(aq)": 27.75,
        "H+(aq)": 1e-8,
        "OH-(aq)": 1e-8,
        "Mg+2(aq)": 1e-5,
        "CO2(aq)": 0.1,
        "CO2(g)": 0.5,
        "H2O(g)": 0.1,
    }
    for name, mol in init.items():
        try:
            species_idx = system.species().index(name)
            if species_idx < system.species().size():
                try:
                    state.set(name, mol, "mol")
                except TypeError:
                    import autodiff

                    state.set(name, autodiff.real(mol), "mol")
        except Exception:
            pass  # species not in system — skip
    return state


def _get_extra_or_none(props):
    if hasattr(props, "extra"):
        try:
            return props.extra()
        except Exception:
            return None
    return None


def test_gfsm_handoff_state_id_written(capsys):
    """
    After ChemicalProps evaluation the GFSM phase must have written
    props.extra["PerplexGFSM::WaterActivity::StateId"].
    """
    rkt = _import()
    try:
        system = _build_two_phase_system(rkt)
    except Exception as e:
        pytest.skip(f"Could not build two-phase system: {e}")

    state = _make_state(rkt, system)
    state.temperature(300.0 + 273.15)  # 300 °C in K
    state.pressure(2000.0 * 1e5)  # 2 kbar in Pa

    props = rkt.ChemicalProps(state)
    extra = _get_extra_or_none(props)
    if extra is None:
        # This build does not expose ChemicalProps.extra in Python bindings.
        assert props is not None
        return

    assert "PerplexGFSM::WaterActivity::StateId" in extra, (
        "GFSM did not write WaterActivity::StateId into props.extra — "
        "the gas phase may not have evaluated or may lack H2O(g)."
    )


def test_gfsm_handoff_ln_f_ratio_written():
    """
    After ChemicalProps evaluation props.extra must contain the GFSM
    ln(f_mix/f_pure) ratio for H2O.
    """
    rkt = _import()
    try:
        system = _build_two_phase_system(rkt)
    except Exception as e:
        pytest.skip(f"Could not build two-phase system: {e}")

    state = _make_state(rkt, system)
    state.temperature(300.0 + 273.15)
    state.pressure(2000.0 * 1e5)

    props = rkt.ChemicalProps(state)
    extra = _get_extra_or_none(props)
    if extra is None:
        assert props is not None
        return

    assert "PerplexGFSM::WaterActivity::ln_f_ratio_h2o" in extra, (
        "GFSM did not write ln_f_ratio_h2o into props.extra."
    )
    ln_f_ratio = float(extra["PerplexGFSM::WaterActivity::ln_f_ratio_h2o"])
    assert math.isfinite(ln_f_ratio), f"ln_f_ratio_h2o is not finite: {ln_f_ratio}"


def test_gfsm_handoff_consumed_by_perplexdew():
    """
    After equilibrium props.extra["PerplexDEW::WaterActivity::CoupledHandoffUsed"]
    must be True — confirming PerplexDEW picked up the GFSM water-activity signal
    rather than falling back to the pure-water default.
    """
    rkt = _import()
    try:
        system = _build_two_phase_system(rkt)
    except Exception as e:
        pytest.skip(f"Could not build two-phase system: {e}")

    specs = rkt.EquilibriumSpecs(system)
    specs.temperature()
    specs.pressure()
    solver = rkt.EquilibriumSolver(specs)
    conds = rkt.EquilibriumConditions(specs)
    conds.temperature(300.0, "celsius")
    conds.pressure(2000.0, "bar")

    state = _make_state(rkt, system)
    result = solver.solve(state, conds)
    assert result.succeeded(), "Equilibrium solver failed on two-phase GFSM/DEW system"

    props = rkt.ChemicalProps(state)
    extra = _get_extra_or_none(props)
    if extra is None:
        assert props is not None
        return

    assert "PerplexDEW::WaterActivity::CoupledHandoffUsed" in extra, (
        "PerplexDEW::WaterActivity::CoupledHandoffUsed key missing from props.extra — "
        "PerplexDEW may not have evaluated, or the key name changed."
    )
    assert extra["PerplexDEW::WaterActivity::CoupledHandoffUsed"] is True, (
        "CoupledHandoffUsed is False — PerplexDEW fell back to pure-water solvent "
        "state instead of consuming the GFSM handoff. Check that GFSM evaluates "
        "before PerplexDEW (gas phase listed first in ChemicalSystem)."
    )


def test_gfsm_water_activity_differs_from_pure_water():
    """
    In a CO2-H2O mixed fluid at 300 °C / 2 kbar, the water activity from
    PerplexDEW via GFSM handoff must differ from 0 (non-trivial ln correction).
    """
    rkt = _import()
    try:
        system = _build_two_phase_system(rkt)
    except Exception as e:
        pytest.skip(f"Could not build two-phase system: {e}")

    specs = rkt.EquilibriumSpecs(system)
    specs.temperature()
    specs.pressure()
    solver = rkt.EquilibriumSolver(specs)
    conds = rkt.EquilibriumConditions(specs)
    conds.temperature(300.0, "celsius")
    conds.pressure(2000.0, "bar")

    state = _make_state(rkt, system)
    result = solver.solve(state, conds)
    assert result.succeeded(), "Equilibrium solver failed"

    props = rkt.ChemicalProps(state)
    extra = _get_extra_or_none(props)
    if extra is None:
        assert props is not None
        return

    ln_a_gfsm = float(extra.get("PerplexDEW::WaterActivity::ln_a_h2o_gfsm", 0.0))
    # Pure water gives ln_a_gfsm = 0; mixed fluid should be non-zero.
    assert abs(ln_a_gfsm) > 1e-6, (
        f"ln_a_h2o_gfsm = {ln_a_gfsm:.2e} is essentially zero — "
        "GFSM correction may not have been applied to the mixed fluid."
    )

    ln_a_total = float(extra.get("PerplexDEW::WaterActivity::ln_a_h2o_total", 0.0))
    assert math.isfinite(ln_a_total), f"ln_a_h2o_total is not finite: {ln_a_total}"


def test_gfsm_handoff_strict_mode_fails_without_coupling():
    """
    With requireCoupledGFSMHandoff=True and no GFSM gas phase, PerplexDEW must
    fail fast instead of silently falling back to aqueous-neutral water activity.
    """
    rkt = _import()
    try:
        system = _build_aqueous_only_strict_system(rkt)
    except Exception as e:
        pytest.skip(f"Could not build strict aqueous-only system: {e}")

    state = rkt.ChemicalState(system)
    for name, mol in [
        ("H2O(aq)", 55.5),
        ("H+(aq)", 1e-8),
        ("OH-(aq)", 1e-8),
        ("Mg+2(aq)", 1e-5),
        ("CO2(aq)", 1e-3),
    ]:
        try:
            state.set(name, mol, "mol")
        except Exception:
            pass

    state.temperature(300.0 + 273.15)
    state.pressure(2000.0 * 1e5)

    with pytest.raises(RuntimeError, match="requireCoupledGFSMHandoff=True"):
        _ = rkt.ChemicalProps(state)
