import importlib.util
import numpy as np

SCRIPT = r"c:\Users\stanroozen\Documents\Projects\reaktoro-dev\reaktoro\DEW_Experimental_Benchmark\Tutorial\willemite_activity_activity_phase_diagram.py"

spec = importlib.util.spec_from_file_location("diag", SCRIPT)
d = importlib.util.module_from_spec(spec)
spec.loader.exec_module(d)

m = d.load_tutorial_module(d.TUTORIAL_PATH)
try:
    m.Warnings.disable(906)
except Exception:
    pass
m.validate_user_inputs()
m.USE_COMPETING_ZN_MINERALS = True

system = m.build_tutorial_system()
specs = m.EquilibriumSpecs(system)
specs.temperature()
specs.pressure()
specs.lgActivity(d.ACTIVITY_SPECIES_X)
specs.lgActivity(d.ACTIVITY_SPECIES_Y)
solver = m.make_equilibrium_solver(system, specs)
cond = m.EquilibriumConditions(specs)
cond.temperature(float(d.TEMPERATURE_C), "celsius")
cond.pressure(float(d.PRESSURE_KBAR * 1000.0), "bar")

xs = np.linspace(d.LOG_A_X_MIN, d.LOG_A_X_MAX, 13)
ys = np.linspace(d.LOG_A_Y_MIN, d.LOG_A_Y_MAX, 13)

none_pts = []
noconv_pts = []

for y in ys:
    for x in xs:
        st = m.make_base_state(system)
        cond.lgActivity(d.ACTIVITY_SPECIES_X, float(x))
        cond.lgActivity(d.ACTIVITY_SPECIES_Y, float(y))
        res = solver.solve(st, cond)
        if not res.succeeded():
            noconv_pts.append((x, y))
            continue

        maxa = 0.0
        for nm in m.selected_mineral_names():
            try:
                a = float(st.speciesAmount(nm))
            except Exception:
                a = 0.0
            if a > maxa:
                maxa = a
        if maxa < 1.0e-14:
            none_pts.append((x, y, maxa))

print(f"sample grid: {len(xs)}x{len(ys)}")
print(f"none points: {len(none_pts)}")
print(f"noconv points: {len(noconv_pts)}")

if none_pts:
    arr = np.array(none_pts, dtype=float)
    print("none bounds:")
    print(f"  x: [{arr[:, 0].min():.3f}, {arr[:, 0].max():.3f}]")
    print(f"  y: [{arr[:, 1].min():.3f}, {arr[:, 1].max():.3f}]")
    print(f"  max solid amount range: [{arr[:, 2].min():.3e}, {arr[:, 2].max():.3e}]")

if noconv_pts:
    arr = np.array(noconv_pts, dtype=float)
    print("noconv bounds:")
    print(f"  x: [{arr[:, 0].min():.3f}, {arr[:, 0].max():.3f}]")
    print(f"  y: [{arr[:, 1].min():.3f}, {arr[:, 1].max():.3f}]")
