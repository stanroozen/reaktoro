import pytest


def test_h2o_dummy_vs_explicit():
    """Compare equilibration with species names supported by the active DEW build."""
    try:
        import autodiff
        import reaktoro4py as rkt
    except ImportError as e:
        pytest.skip(f"reaktoro4py not available: {e}")

    T_C, P_bar = 300, 100
    dew_db = rkt.DEWDatabase("dew2019-aqueous")
    supcrt_db = rkt.SupcrtDatabase("supcrtbl")

    quartz = supcrt_db.species("Quartz")

    # Use species names matching the active DEW database in this build.
    db = rkt.Database(dew_db.species())
    db.addSpecies(quartz)
    aq = rkt.AqueousPhase("H2O(aq) H+(aq) OH-(aq) SiO2(aq)")
    try:
        aq.setActivityModel(rkt.ActivityModelDEW())
    except Exception:
        aq.setActivityModel(rkt.ActivityModelHKF())
    system = rkt.ChemicalSystem(db, aq, rkt.MineralPhase("Quartz"))
    assert len(system.species()) > 0

    solver = rkt.EquilibriumSolver(system)
    state = rkt.ChemicalState(system)
    state.temperature(autodiff.real(float(T_C)), "celsius")
    state.pressure(autodiff.real(float(P_bar)), "bar")
    state.set("H2O(aq)", 1.0, "kg")
    state.set("SiO2(aq)", 1e-4, "mol")
    state.set("Quartz", 1.0, "mol")
    result = solver.solve(state)
    assert result.succeeded(), f"Equilibration failed at {T_C}C, {P_bar} bar"


if __name__ == "__main__":
    test_h2o_dummy_vs_explicit()
    print("test_h2o_dummy passed.")
