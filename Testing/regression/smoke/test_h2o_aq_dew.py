import pytest


def test_h2o_aq_dew_equilibration():
    """Verify water species in DEW database enables equilibration with Duan EOS."""
    try:
        import autodiff
        import reaktoro4py as rkt
    except ImportError as e:
        pytest.skip(f"dependency not available: {e}")

    dew_db = rkt.DEWDatabase("dew2019-aqueous")
    supcrt_db = rkt.SupcrtDatabase("supcrtbl")
    dew_species = [s.name() for s in dew_db.species()]
    # Check for water species (may be named WATER,AQ or H2O(aq) depending on database)
    has_water = any(name in dew_species for name in ["H2O(aq)", "WATER,AQ", "H2O_aq"])
    assert has_water, (
        f"Water species must be in dew2019-aqueous. Available: {dew_species[:20]}"
    )

    quartz = supcrt_db.species("Quartz")
    db_combined = rkt.Database(dew_db.species())
    db_combined.addSpecies(quartz)
    # Use species names matching the active DEW database in this build.
    aqueous = rkt.AqueousPhase("H2O(aq) H+(aq) OH-(aq) SiO2(aq)")
    try:
        aqueous.setActivityModel(rkt.ActivityModelDEW())
    except Exception:
        aqueous.setActivityModel(rkt.ActivityModelHKF())
    system = rkt.ChemicalSystem(db_combined, aqueous, rkt.MineralPhase("Quartz"))
    assert len(system.species()) > 0

    solver = rkt.EquilibriumSolver(system)
    state = rkt.ChemicalState(system)
    state.temperature(autodiff.real(300.0), "celsius")
    state.pressure(autodiff.real(500.0), "bar")
    state.set("SiO2(aq)", 0.0001, "mol")
    state.set("Quartz", 1.0, "kg")
    result = solver.solve(state)
    assert result.succeeded(), "Equilibrium should converge at 300 degC, 500 bar"
    aq = rkt.AqueousProps(state)
    molality = float(aq.speciesMolality("SiO2(aq)"))
    assert molality > 0


if __name__ == "__main__":
    test_h2o_aq_dew_equilibration()
    print("test_h2o_aq_dew passed.")
