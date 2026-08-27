import importlib.util
import sys
import os

p = r"c:\Users\stanroozen\Documents\Projects\reaktoro-dev\reaktoro\DEW_Experimental_Benchmark\Tutorial\willemite_stability_diagram_suite.py"
spec = importlib.util.spec_from_file_location("ws", p)
ws = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ws)
ws.OUTPUT_T_PH = os.path.join(
    r"c:\Users\stanroozen\Documents\Projects\reaktoro-dev\reaktoro\temp",
    "willemite_stability_temperature_ph_fixed_sio2_bench.png",
)

m = ws.load_tutorial_module(ws.TUTORIAL_PATH)
try:
    m.Warnings.disable(906)
except Exception:
    pass
m.validate_user_inputs()
m.USE_COMPETING_ZN_MINERALS = True
cats = list(m.selected_mineral_names()) + ["AqueousOnly", "NoConvergence"]


def run_case(name, eps=None, barrier=None, ideal=None):
    old_make = m.make_equilibrium_solver

    def patched_make(system, specs):
        s = old_make(system, specs)
        try:
            o = m.EquilibriumOptions()
            if eps is not None and hasattr(o, "epsilon"):
                o.epsilon = float(eps)
            if barrier is not None and hasattr(o, "logarithm_barrier_factor"):
                o.logarithm_barrier_factor = float(barrier)
            if ideal is not None and hasattr(o, "use_ideal_activity_models"):
                o.use_ideal_activity_models = bool(ideal)
            s.setOptions(o)
        except Exception as exc:
            print(name, "option-set-warning", exc)
        return s

    m.make_equilibrium_solver = patched_make
    t = ws.save_temperature_ph_map(m, cats)
    m.make_equilibrium_solver = old_make
    total = t["converged"] + t["failed"]
    print(name, t["converged"], "/", total)


cases = [
    ("default", None, None, None),
    ("eps1e-13", 1e-13, None, None),
    ("eps1e-13_bar1e-2", 1e-13, 1e-2, None),
    ("eps1e-13_ideal", 1e-13, None, True),
]

for case in cases:
    run_case(*case)
