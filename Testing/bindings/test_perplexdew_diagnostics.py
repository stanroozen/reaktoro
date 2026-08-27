"""
Regression tests for ActivityModelPerplexDEW diagnostic flags.

Tests:
  1. Default mode (errorOnConflictingStandardState=False): a species with both
     GFSM fluid index and PerplexDEW HKF model triggers a warning only —
     ChemicalSystem construction must NOT raise.
  2. Strict mode (errorOnConflictingStandardState=True): same species raises
     RuntimeError during ChemicalSystem construction.
  3. Unmapped GFSM coupling warning (warnOnUnmappedGFSMCoupling=True): a
     species whose formula has an aqueous suffix that, when stripped, matches a
     known GFSM key triggers a warning but does NOT raise.
  4. Unmapped GFSM coupling suppressed (warnOnUnmappedGFSMCoupling=False): same
     scenario emits no diagnostic (construction still succeeds).

Build requirements: reaktoro4py built in build-msvc/Reaktoro/Release.
"""

import pytest

try:
    import reaktoro4py as rkt
except Exception as e:
    pytest.skip(f"reaktoro4py not available: {e}", allow_module_level=True)

try:
    import autodiff
except Exception:
    autodiff = None


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _real(value: float):
    if autodiff is None:
        return value
    return autodiff.real(value)


def _make_perplexdew_hkf_model():
    """Return a StandardThermoModelPerplexDEW with minimal non-zero params."""
    params = rkt.StandardThermoModelParamsPerplexDEW()
    params.Gf = _real(-394359.0)  # J/mol  (CO2 approx)
    params.Sr = _real(213.79)  # J/(mol·K)
    params.wref = _real(1.23456e5)  # J/mol  Born coeff — triggers born_vec.active=True
    # HKF Helgeson power-series terms (zeros are fine for construction-time tests)
    params.a1 = _real(0.0)
    params.a2 = _real(0.0)
    params.a3 = _real(0.0)
    params.a4 = _real(0.0)
    params.c1 = _real(0.0)
    params.c2 = _real(0.0)
    return rkt.StandardThermoModelPerplexDEW(params)


def _make_conflict_species(name: str = "CO2hkf"):
    """
    Return a Species with a valid identifier and formula "CO2".
    The formula maps to a GFSM fluid key while the PerplexDEW HKF model is active,
    producing the intended conflict scenario.
    """
    return (
        rkt.Species("CO2")
        .withName(name)
        .withFormula("CO2")
        .withStandardThermoModel(_make_perplexdew_hkf_model())
    )


def _make_unmapped_species(name: str = "CO2aqsp"):
    """
    Return a neutral non-HKF species with formula CO2(aq).
    Uses withStandardGibbsEnergy only (no PerplexDEW HKF params), ensuring this is
    a neutral coupling-path diagnostic case.
    """
    return (
        rkt.Species("CO2")
        .withName(name)
        .withFormula("CO2(aq)")
        .withStandardGibbsEnergy(_real(-394359.0))
    )


def _make_test_db_conflict():
    """Database with H2O (from DEW embedded DB) + the conflict CO2_hkf species."""
    dew = rkt.DEWDatabase("dew2024-aqueous")
    water = dew.species("H2O(aq)")
    conflict_sp = _make_conflict_species()
    return rkt.Database([water, conflict_sp])


def _make_test_db_unmapped():
    """Database with H2O (from DEW embedded DB) + the unmapped-formula CO2_suffix."""
    dew = rkt.DEWDatabase("dew2024-aqueous")
    water = dew.species("H2O(aq)")
    unmapped_sp = _make_unmapped_species()
    return rkt.Database([water, unmapped_sp])


# ---------------------------------------------------------------------------
# Test 1 — conflict, default mode: warning only, no exception
# ---------------------------------------------------------------------------


def test_conflict_default_mode_no_raise():
    """
    When errorOnConflictingStandardState=False (default), a GFSM+HKF conflict
    emits a warning to stderr but ChemicalSystem construction succeeds.
    """
    db = _make_test_db_conflict()
    phase = rkt.AqueousPhase("H2O(aq) CO2hkf")
    act_params = rkt.ActivityModelParamsPerplexDEW()
    act_params.errorOnConflictingStandardState = False
    phase.setActivityModel(rkt.ActivityModelPerplexDEW(act_params))
    # Must NOT raise
    system = rkt.ChemicalSystem(db, phase)
    assert system is not None


# ---------------------------------------------------------------------------
# Test 2 — conflict, strict mode: RuntimeError raised
# ---------------------------------------------------------------------------


def test_conflict_strict_mode_raises():
    """
    When errorOnConflictingStandardState=True, a GFSM+HKF conflict raises
    RuntimeError during ChemicalSystem construction.
    """
    db = _make_test_db_conflict()
    phase = rkt.AqueousPhase("H2O(aq) CO2hkf")
    act_params = rkt.ActivityModelParamsPerplexDEW()
    act_params.errorOnConflictingStandardState = True
    phase.setActivityModel(rkt.ActivityModelPerplexDEW(act_params))
    with pytest.raises(RuntimeError):
        rkt.ChemicalSystem(db, phase)


# ---------------------------------------------------------------------------
# Test 3 — unmapped GFSM coupling, warn mode: no exception
# ---------------------------------------------------------------------------


def test_unmapped_gfsm_warn_no_raise():
    """
    When warnOnUnmappedGFSMCoupling=True (default), a species with formula
    'CO2_aq' (suffix mismatch) triggers a warning but construction succeeds.
    """
    db = _make_test_db_unmapped()
    phase = rkt.AqueousPhase("H2O(aq) CO2aqsp")
    act_params = rkt.ActivityModelParamsPerplexDEW()
    act_params.warnOnUnmappedGFSMCoupling = True
    phase.setActivityModel(rkt.ActivityModelPerplexDEW(act_params))
    system = rkt.ChemicalSystem(db, phase)
    assert system is not None


# ---------------------------------------------------------------------------
# Test 4 — unmapped GFSM coupling, suppressed: no exception, no warn
# ---------------------------------------------------------------------------


def test_unmapped_gfsm_suppressed_no_raise():
    """
    When warnOnUnmappedGFSMCoupling=False, the suffix-mismatch check is
    skipped entirely; construction still succeeds.
    """
    db = _make_test_db_unmapped()
    phase = rkt.AqueousPhase("H2O(aq) CO2aqsp")
    act_params = rkt.ActivityModelParamsPerplexDEW()
    act_params.warnOnUnmappedGFSMCoupling = False
    phase.setActivityModel(rkt.ActivityModelPerplexDEW(act_params))
    system = rkt.ChemicalSystem(db, phase)
    assert system is not None


# ---------------------------------------------------------------------------
# Test 5 — correct formula passes cleanly (no conflict, no mismatch)
# ---------------------------------------------------------------------------


def test_clean_species_no_diagnostic():
    """
    A neutral species with formula "N2" (GFSM pidx=10) but NO PerplexDEW HKF
    model produces no conflict and no unmapped warning — clean construction.
    """
    dew = rkt.DEWDatabase("dew2024-aqueous")
    water = dew.species("H2O(aq)")
    n2_sp = (
        rkt.Species("N2")
        .withName("N2test")
        .withFormula("N2")
        .withStandardGibbsEnergy(_real(0.0))
    )
    db = rkt.Database([water, n2_sp])
    phase = rkt.AqueousPhase("H2O(aq) N2test")
    act_params = rkt.ActivityModelParamsPerplexDEW()
    act_params.errorOnConflictingStandardState = True  # strict — must still pass
    act_params.warnOnUnmappedGFSMCoupling = True
    phase.setActivityModel(rkt.ActivityModelPerplexDEW(act_params))
    system = rkt.ChemicalSystem(db, phase)
    assert system is not None


# ---------------------------------------------------------------------------
# Test 6 — params field round-trip
# ---------------------------------------------------------------------------


def test_params_field_roundtrip():
    """ActivityModelParamsPerplexDEW fields are readable/writable from Python."""
    p = rkt.ActivityModelParamsPerplexDEW()
    assert p.errorOnConflictingStandardState is False
    assert p.warnOnUnmappedGFSMCoupling is True
    assert p.requireCoupledGFSMHandoff is False

    p.errorOnConflictingStandardState = True
    p.warnOnUnmappedGFSMCoupling = False
    p.requireCoupledGFSMHandoff = True
    assert p.errorOnConflictingStandardState is True
    assert p.warnOnUnmappedGFSMCoupling is False
    assert p.requireCoupledGFSMHandoff is True
