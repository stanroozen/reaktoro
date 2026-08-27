import pytest


def test_perplex_conditions_nosilence():
    """PerplexDEW conditions sweep without silencing warnings."""
    try:
        import autodiff
        import reaktoro4py as rkt
    except ImportError as e:
        pytest.skip(f"dependency not available: {e}")

    dew_db = rkt.DEWDatabase("dew2024-aqueous")
    supcrt_db = rkt.SupcrtDatabase("supcrtbl")
    combined_db = rkt.Database(dew_db.species())
    combined_db.addSpecies(supcrt_db.species("Quartz"))

    aq = rkt.AqueousPhase(
        "H2O(aq) H+(aq) OH-(aq) SiO2(aq) H2(aq) O2(aq) HO2-(aq) HSiO3-(aq) Si2O4(aq) Si3O6(aq)"
    )
    aq.setActivityModel(rkt.ActivityModelPerplexDEW())
    sys_ = rkt.ChemicalSystem(combined_db, aq, rkt.MineralPhase("Quartz"))

    specs = rkt.EquilibriumSpecs(sys_)
    specs.temperature()
    specs.pressure()
    solver = rkt.EquilibriumSolver(specs)
    conds = rkt.EquilibriumConditions(specs)

    conditions = [(200, 500), (300, 1000), (400, 2000)]
    failures = []
    for T_C, P_bar in conditions:
        state = rkt.ChemicalState(sys_)
        state.set("H2O(aq)", autodiff.real(1.0), "kg")
        state.set("H+(aq)", autodiff.real(1e-8), "mol")
        state.set("OH-(aq)", autodiff.real(1e-8), "mol")
        state.set("SiO2(aq)", autodiff.real(1e-6), "mol")
        state.set("Quartz", autodiff.real(10.0), "mol")
        conds.temperature(float(T_C), "celsius")
        conds.pressure(float(P_bar), "bar")
        result = solver.solve(state, conds)
        if not result.succeeded():
            failures.append((T_C, P_bar))
    assert not failures, f"PerplexDEW solve failed at conditions: {failures}"


if __name__ == "__main__":
    test_perplex_conditions_nosilence()
    print("test_perplex_conditions_nosilence passed.")
