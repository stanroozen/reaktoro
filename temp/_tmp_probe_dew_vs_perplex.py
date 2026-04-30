import importlib.util

script = r"c:\Users\stanroozen\Documents\Projects\reaktoro-dev\reaktoro\DEW_Experimental_Benchmark\quartz_solubility_analysis_v2_dew24.py"
spec = importlib.util.spec_from_file_location("q", script)
q = importlib.util.module_from_spec(spec)
spec.loader.exec_module(q)

dew_db = q.DEWDatabase("dew2024-aqueous")
supcrt_db = q.SupcrtDatabase("supcrtbl")

def check_backend(backend):
    system = q.build_system(dew_db, supcrt_db, q.MINERAL_CONFIG, water_config=q.DEW_CONFIG, model_backend=backend)
    specs = q.EquilibriumSpecs(system)
    specs.temperature()
    specs.pressure()
    solver = q.EquilibriumSolver(specs)
    conds = q.EquilibriumConditions(specs)

    print(f"\nBackend={backend}")
    for P_kbar, T_C in [(1.0, 350), (0.25, 360), (0.000293842, 100), (10.0, 700)]:
        state = q.ChemicalState(system)
        state.set("WATER,AQ", q.to_real(1.0), "kg")
        state.set("H+", q.to_real(1e-8), "mol")
        state.set("OH-", q.to_real(1e-8), "mol")
        state.set(q.MINERAL_CONFIG["solute_species"], q.to_real(1e-6), "mol")
        state.set(q.MINERAL_CONFIG["mineral_name"], q.to_real(10.0), "mol")
        conds.temperature(float(T_C), "celsius")
        conds.pressure(float(P_kbar * 1000.0), "bar")
        r = solver.solve(state, conds)
        print(f"  P={P_kbar} T={T_C}: ok={bool(r.succeeded())}, iter={int(r.iterations())}")

check_backend("DEW")
check_backend("PerplexDEW")
