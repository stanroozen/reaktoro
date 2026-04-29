import importlib.util

script = r"c:\Users\stanroozen\Documents\Projects\reaktoro-dev\reaktoro\DEW_Experimental_Benchmark\quartz_solubility_analysis_v2_dew24.py"
spec = importlib.util.spec_from_file_location("q", script)
q = importlib.util.module_from_spec(spec)
spec.loader.exec_module(q)

print("module loaded, backend required symbols present:", all(name in q.__dict__ for name in q.PERPLEXDEW_REQUIRED_SYMBOLS))

dew_db = q.DEWDatabase("dew2024-aqueous")
supcrt_db = q.SupcrtDatabase("supcrtbl")
system = q.build_system(dew_db, supcrt_db, q.MINERAL_CONFIG, water_config=q.DEW_CONFIG, model_backend="PerplexDEW")

specs = q.EquilibriumSpecs(system)
specs.temperature()
specs.pressure()
solver = q.EquilibriumSolver(specs)
conds = q.EquilibriumConditions(specs)


def run(P_kbar, T_C):
    state = q.ChemicalState(system)
    state.set("WATER,AQ", q.to_real(1.0), "kg")
    state.set("H+", q.to_real(1e-8), "mol")
    state.set("OH-", q.to_real(1e-8), "mol")
    state.set(q.MINERAL_CONFIG["solute_species"], q.to_real(1e-6), "mol")
    state.set(q.MINERAL_CONFIG["mineral_name"], q.to_real(10.0), "mol")
    conds.temperature(float(T_C), "celsius")
    conds.pressure(float(P_kbar * 1000.0), "bar")
    res = solver.solve(state, conds)
    ok = bool(res.succeeded())
    print(f"P={P_kbar} kbar T={T_C}C succeeded={ok}")
    if ok:
        try:
            aq = q.AqueousProps(state)
            m = q.total_element_molality(aq, q.MINERAL_CONFIG, q.get_solute_species_list(q.MINERAL_CONFIG))
            print("  molality=", m)
        except Exception as e:
            print("  aqprops/total molality failed:", repr(e))
    else:
        for attr in ["message", "error", "status", "numIterations", "iterations"]:
            if hasattr(res, attr):
                try:
                    print(f"  {attr}=", getattr(res, attr)())
                except Exception:
                    try:
                        print(f"  {attr}=", getattr(res, attr))
                    except Exception:
                        pass
        print("  result type:", type(res))
        print("  result dir sample:", [n for n in dir(res) if any(k in n.lower() for k in ["error", "status", "iter", "message", "succeed"])])


for P, T in [
    (1.0, 350),
    (0.25, 360),
    (0.35, 390),
    (0.5, 500),
    (0.000293842, 100),
    (0.000293842, 20),
]:
    run(P, T)
