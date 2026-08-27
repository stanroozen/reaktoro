import argparse
import importlib.util
import os
import sys

import numpy as np


def load_tutorial_module(repo_root):
    rel = os.path.join(repo_root, "build", "Reaktoro", "Release")
    if rel not in sys.path:
        sys.path.insert(0, rel)

    tutorial = os.path.join(
        repo_root,
        "DEW_Experimental_Benchmark",
        "Tutorial",
        "willemite_solubility_tutorial_dew17hp622_zn.py",
    )
    spec = importlib.util.spec_from_file_location("w", tutorial)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def present(state, name, threshold):
    try:
        return float(state.speciesAmount(name)) > threshold
    except Exception:
        return False


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gas", choices=["O2", "CO2"], required=True)
    args = parser.parse_args()

    repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    w = load_tutorial_module(repo)

    w.USE_COMPETING_ZN_MINERALS = True
    w.INCLUDE_COMPATIBLE_GANGUE_MINERALS = False
    w.COMPETING_ZN_MINERALS = list(
        dict.fromkeys(w.COMPETING_ZN_MINERALS + ["hem", "cc", "q"])
    )

    extra_aq = [
        "Ca+2",
        "Ca(OH)+",
        "CaCO3,aq",
        "Ca(HCO3)",
        "CaCl+",
        "CaCl2,aq",
        "CaSO4,aq",
        "Fe+2",
        "Fe+3",
        "Fe(OH)+",
        "FeO,aq",
        "HFeO2-",
        "FeCl+",
        "FeCl2,aq",
        "FeCl2+",
    ]
    w.AQUEOUS_SPECIES = list(dict.fromkeys(list(w.AQUEOUS_SPECIES) + extra_aq))

    for m in ["hem", "cc", "q"]:
        w.INITIAL_SPECIES_AMOUNTS_MOL[m] = 2.0

    w.validate_user_inputs()

    database = w.Database.fromFile(w.PERPLEX_DATABASE_FILE)
    aq = w.AqueousPhase(" ".join(w.AQUEOUS_SPECIES))
    aq.setActivityModel(w.AQUEOUS_ACTIVITY_MODEL())
    mins = w.make_mineral_phases()
    gas = w.GaseousPhase(args.gas)
    gas.setActivityModel(w.ActivityModelIdealGas())
    system = w.ChemicalSystem(database, aq, mins, gas)

    specs = w.EquilibriumSpecs(system)
    specs.temperature()
    specs.pressure()
    specs.pH()
    specs.fugacity(args.gas)

    solver = w.make_equilibrium_solver(system, specs)
    opts = w.EquilibriumOptions()
    if hasattr(opts, "hessian") and hasattr(w, "GibbsHessian"):
        opts.hessian = w.GibbsHessian.Exact
    if hasattr(opts, "warmstart"):
        opts.warmstart = True
    solver.setOptions(opts)
    conditions = w.EquilibriumConditions(specs)

    t_vals = [200.0, 250.0, 300.0, 350.0, 400.0]
    ph_vals = np.linspace(4.0, 10.0, 7)
    logf_vals = (
        [-30.0, -25.0, -20.0, -15.0, -10.0]
        if args.gas == "O2"
        else [-5.0, -4.0, -3.0, -2.0, -1.0]
    )
    p_bar = 2000.0
    th = 1.0e-8

    targets = {
        "Wlm+hem+cc": lambda st: (
            present(st, "Wlm", th) and present(st, "hem", th) and present(st, "cc", th)
        ),
        "Wlm+hem+cc+q": lambda st: (
            present(st, "Wlm", th)
            and present(st, "hem", th)
            and present(st, "cc", th)
            and present(st, "q", th)
        ),
        "Wlm+hem+cc+q+Znc": lambda st: (
            present(st, "Wlm", th)
            and present(st, "hem", th)
            and present(st, "cc", th)
            and present(st, "q", th)
            and present(st, "Znc", th)
        ),
    }

    hits = {k: [] for k in targets}
    solved = 0
    total = len(t_vals) * len(ph_vals) * len(logf_vals)

    for t in t_vals:
        state = w.make_base_state(system)
        try:
            state.set(args.gas, 1.0e-20, "mol")
        except Exception:
            pass

        for ph in ph_vals:
            for lf in logf_vals:
                conditions.temperature(float(t), "celsius")
                conditions.pressure(float(p_bar), "bar")
                conditions.pH(float(ph))
                conditions.fugacity(args.gas, float(10.0**lf), "bar")

                result = solver.solve(state, conditions)
                if not result.succeeded():
                    continue

                solved += 1
                point = (float(t), float(ph), float(lf))
                for name, check in targets.items():
                    if check(state):
                        hits[name].append(point)

    print(f"GAS={args.gas}")
    print(f"SOLVED={solved}/{total}")
    for name, pts in hits.items():
        print(f"\n{name}: {len(pts)} points")
        if not pts:
            continue

        ts = [p[0] for p in pts]
        phs = [p[1] for p in pts]
        lfs = [p[2] for p in pts]
        print(f"  T range (C): {min(ts):.1f} to {max(ts):.1f}")
        print(f"  pH range: {min(phs):.2f} to {max(phs):.2f}")
        print(f"  log10(f{args.gas}/bar) range: {min(lfs):.1f} to {max(lfs):.1f}")
        for p in pts[:12]:
            print(f"    T={p[0]:.1f} C, pH={p[1]:.2f}, logf{args.gas}={p[2]:.1f}")


if __name__ == "__main__":
    main()
