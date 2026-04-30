import importlib.util
import numpy as np
import pandas as pd

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

curves = pd.read_csv(r"c:\Users\stanroozen\Documents\Projects\reaktoro-dev\reaktoro\DEW_Experimental_Benchmark\quartz_curves_dew24_PerplexDEW.csv")

for P_kbar in [0.25, 0.3, 0.5]:
    sub = curves[(curves["curve_type"] == "isobar") & (np.isclose(curves["P_kbar"], P_kbar, rtol=0, atol=1e-9))]
    Tvals = sub["T_C"].to_numpy(dtype=float)

    succ = 0
    fail = 0
    aq_err = 0
    iters = []

    for T_C in Tvals:
        state = q.ChemicalState(system)
        state.set("WATER,AQ", q.to_real(1.0), "kg")
        state.set("H+", q.to_real(1e-8), "mol")
        state.set("OH-", q.to_real(1e-8), "mol")
        state.set(q.MINERAL_CONFIG["solute_species"], q.to_real(1e-6), "mol")
        state.set(q.MINERAL_CONFIG["mineral_name"], q.to_real(10.0), "mol")
        conds.temperature(float(T_C), "celsius")
        conds.pressure(float(P_kbar * 1000.0), "bar")
        r = solver.solve(state, conds)
        iters.append(int(r.iterations()))
        if r.succeeded():
            succ += 1
            try:
                aq = q.AqueousProps(state)
                _ = q.total_element_molality(aq, q.MINERAL_CONFIG, q.get_solute_species_list(q.MINERAL_CONFIG))
            except Exception:
                aq_err += 1
        else:
            fail += 1

    print(f"P={P_kbar} kbar: success={succ}, fail={fail}, aqprops_err={aq_err}, iters_unique={sorted(set(iters))[:6]}{'...' if len(set(iters))>6 else ''}")
