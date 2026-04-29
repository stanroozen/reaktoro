import importlib.util
import numpy as np

script = r"c:\Users\stanroozen\Documents\Projects\reaktoro-dev\reaktoro\DEW_Experimental_Benchmark\quartz_solubility_analysis_v2_dew24.py"
spec = importlib.util.spec_from_file_location("q", script)
q = importlib.util.module_from_spec(spec)
spec.loader.exec_module(q)

dew_db = q.DEWDatabase("dew2024-aqueous")
supcrt_db = q.SupcrtDatabase("supcrtbl")
system = q.build_system(dew_db, supcrt_db, q.MINERAL_CONFIG, water_config=q.DEW_CONFIG, model_backend="PerplexDEW")

specs = q.EquilibriumSpecs(system)
specs.temperature()
specs.pressure()
solver = q.EquilibriumSolver(specs)
conds = q.EquilibriumConditions(specs)


def init_state():
    s = q.ChemicalState(system)
    s.set("WATER,AQ", q.to_real(1.0), "kg")
    s.set("H+", q.to_real(1e-8), "mol")
    s.set("OH-", q.to_real(1e-8), "mol")
    s.set(q.MINERAL_CONFIG["solute_species"], q.to_real(1e-6), "mol")
    s.set(q.MINERAL_CONFIG["mineral_name"], q.to_real(10.0), "mol")
    return s


P_kbar = 1.0
T = np.linspace(45, 605, 20)

state = init_state()
ok_cont = 0
iters_cont = []
for t in T:
    conds.temperature(float(t), "celsius")
    conds.pressure(float(P_kbar * 1000), "bar")
    r = solver.solve(state, conds)
    it = int(r.iterations()) if hasattr(r, "iterations") else -1
    iters_cont.append(it)
    ok_cont += int(bool(r.succeeded()))
print("continuation: ok", ok_cont, "/", len(T), "iters", iters_cont)

ok_fresh = 0
iters_fresh = []
for t in T:
    state = init_state()
    conds.temperature(float(t), "celsius")
    conds.pressure(float(P_kbar * 1000), "bar")
    r = solver.solve(state, conds)
    it = int(r.iterations()) if hasattr(r, "iterations") else -1
    iters_fresh.append(it)
    ok_fresh += int(bool(r.succeeded()))
print("fresh-state: ok", ok_fresh, "/", len(T), "iters", iters_fresh)
