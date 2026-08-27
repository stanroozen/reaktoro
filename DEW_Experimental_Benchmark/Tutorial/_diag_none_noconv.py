import importlib.util
import numpy as np

DIAG_PATH = r"c:\Users\stanroozen\Documents\Projects\reaktoro-dev\reaktoro\DEW_Experimental_Benchmark\Tutorial\willemite_activity_activity_phase_diagram.py"


def run():
    spec = importlib.util.spec_from_file_location("diag", DIAG_PATH)
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

    xs = np.linspace(d.LOG_A_X_MIN, d.LOG_A_X_MAX, 21)
    ys = np.linspace(d.LOG_A_Y_MIN, d.LOG_A_Y_MAX, 21)
    minerals = list(m.selected_mineral_names())

    none_pts = []
    max_solid_succeeded = []
    nconv_pts = []
    total = 0
    succ = 0

    for y in ys:
        for x in xs:
            total += 1
            st = m.make_base_state(system)
            cond.lgActivity(d.ACTIVITY_SPECIES_X, float(x))
            cond.lgActivity(d.ACTIVITY_SPECIES_Y, float(y))
            res = solver.solve(st, cond)
            if not res.succeeded():
                nconv_pts.append((x, y))
                continue
            succ += 1

            amounts = []
            for nm in minerals:
                try:
                    a = float(st.speciesAmount(nm))
                except Exception:
                    a = 0.0
                amounts.append(a)

            maxa = max(amounts) if amounts else 0.0
            max_solid_succeeded.append(maxa)
            if maxa < 1.0e-14:
                none_pts.append((x, y, maxa, sum(amounts)))

    print(
        f"grid {len(xs)}x{len(ys)} total={total} succeeded={succ} noconv={len(nconv_pts)} none={len(none_pts)}"
    )

    if nconv_pts:
        arr = np.array(nconv_pts, dtype=float)
        print(
            "noconv bounds: "
            f"x[{arr[:, 0].min():.3f}, {arr[:, 0].max():.3f}], "
            f"y[{arr[:, 1].min():.3f}, {arr[:, 1].max():.3f}]"
        )

    if none_pts:
        arr = np.array([(p[0], p[1], p[2], p[3]) for p in none_pts], dtype=float)
        print(
            "none bounds: "
            f"x[{arr[:, 0].min():.3f}, {arr[:, 0].max():.3f}], "
            f"y[{arr[:, 1].min():.3f}, {arr[:, 1].max():.3f}]"
        )
        print(
            f"none max-solid amount: min={arr[:, 2].min():.3e} max={arr[:, 2].max():.3e}"
        )
        print(
            f"none total-solid amount: min={arr[:, 3].min():.3e} max={arr[:, 3].max():.3e}"
        )

    max_solid_succeeded = np.array(max_solid_succeeded, dtype=float)
    for t in [1e-14, 1e-16, 1e-20, 1e-30]:
        n = int(np.sum(max_solid_succeeded < t))
        print(f"none_count threshold {t:.0e} = {n}")


if __name__ == "__main__":
    run()
