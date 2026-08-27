import importlib.util
import os
import sys

import numpy as np


REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REL = os.path.join(REPO, "build", "Reaktoro", "Release")
if REL not in sys.path:
    sys.path.insert(0, REL)

TUTORIAL = os.path.join(
    REPO,
    "DEW_Experimental_Benchmark",
    "Tutorial",
    "willemite_solubility_tutorial_dew17hp622_zn.py",
)
spec = importlib.util.spec_from_file_location("w", TUTORIAL)
w = importlib.util.module_from_spec(spec)
spec.loader.exec_module(w)


w.USE_COMPETING_ZN_MINERALS = True
w.INCLUDE_COMPATIBLE_GANGUE_MINERALS = False
w.COMPETING_ZN_MINERALS = list(
    dict.fromkeys(w.COMPETING_ZN_MINERALS + ["hem", "cc", "mag", "q"])
)

extra_aq = [
    "Ca+2",
    "Ca(OH)+",
    "CaCO3,aq",
    "Ca(HCO3)",
    "CaCl+",
    "CaCl2,aq",
    "CaSO4,aq",
    "Mg+2",
    "MgOH+",
    "MgCO3,aq",
    "Mg(HCO3)",
    "MgCl+",
    "MgSO4,aq",
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

for m in ["hem", "cc", "mag", "q"]:
    w.INITIAL_SPECIES_AMOUNTS_MOL[m] = 2.0

w.validate_user_inputs()


Tvals = [200.0, 250.0, 300.0, 350.0, 400.0]
pHvals = np.linspace(4.0, 10.0, 7)
Pbar = 2000.0
TH = 1.0e-8


def present(st, name):
    try:
        return float(st.speciesAmount(name)) > TH
    except Exception:
        return False


def summarize_hits(tag, hits, vars_fmt):
    print(f"\n[{tag}]")
    for name, pts in hits.items():
        print(f"{name}: {len(pts)} points")
        if not pts:
            continue
        Ts = [p[0] for p in pts]
        pHs = [p[1] for p in pts]
        print(f"  T range (C): {min(Ts):.1f} to {max(Ts):.1f}")
        print(f"  pH range: {min(pHs):.2f} to {max(pHs):.2f}")
        print("  first points:")
        for p in pts[:15]:
            print(f"    T={p[0]:.1f} C, pH={p[1]:.2f}, {vars_fmt(p)}")


# 1) Redox-buffered study: vary fO2 with O2 gas phase.
def run_redox_scan():
    database = w.Database.fromFile(w.PERPLEX_DATABASE_FILE)
    aq = w.AqueousPhase(" ".join(w.AQUEOUS_SPECIES))
    aq.setActivityModel(w.AQUEOUS_ACTIVITY_MODEL())
    mins = w.make_mineral_phases()
    gas = w.GaseousPhase("O2")
    gas.setActivityModel(w.ActivityModelIdealGas())
    system = w.ChemicalSystem(database, aq, mins, gas)

    specs = w.EquilibriumSpecs(system)
    specs.temperature()
    specs.pressure()
    specs.pH()
    specs.fugacity("O2")

    solver = w.make_equilibrium_solver(system, specs)
    opts = w.EquilibriumOptions()
    if hasattr(opts, "hessian") and hasattr(w, "GibbsHessian"):
        opts.hessian = w.GibbsHessian.Exact
    if hasattr(opts, "warmstart"):
        opts.warmstart = True
    solver.setOptions(opts)
    conditions = w.EquilibriumConditions(specs)

    logfO2_vals = [-30.0, -25.0, -20.0, -15.0, -10.0]
    total = len(Tvals) * len(pHvals) * len(logfO2_vals)

    targets = {
        "Wlm+hem+cc": lambda st: (
            present(st, "Wlm") and present(st, "hem") and present(st, "cc")
        ),
        "Wlm+hem+cc+q": lambda st: (
            present(st, "Wlm")
            and present(st, "hem")
            and present(st, "cc")
            and present(st, "q")
        ),
        "RELAXED_Wlm+hem+(cc|mag)+q": lambda st: (
            present(st, "Wlm")
            and present(st, "hem")
            and (present(st, "cc") or present(st, "mag"))
            and present(st, "q")
        ),
    }
    hits = {k: [] for k in targets}
    solved = 0

    for T in Tvals:
        state = w.make_base_state(system)
        try:
            state.set("O2", 1.0e-20, "mol")
        except Exception:
            pass

        for pH in pHvals:
            for lfO2 in logfO2_vals:
                conditions.temperature(float(T), "celsius")
                conditions.pressure(float(Pbar), "bar")
                conditions.pH(float(pH))
                conditions.fugacity("O2", float(10.0**lfO2), "bar")

                result = solver.solve(state, conditions)
                if not result.succeeded():
                    continue

                solved += 1
                point = (float(T), float(pH), float(lfO2))
                for name, check in targets.items():
                    if check(state):
                        hits[name].append(point)

    print(f"REDOX_SOLVED={solved}/{total}")
    summarize_hits("RedoxBuffered", hits, lambda p: f"logfO2={p[2]:.1f}")


# 2) Carbonate-buffered study: vary fCO2 with CO2 gas phase.
def run_carb_scan():
    database = w.Database.fromFile(w.PERPLEX_DATABASE_FILE)
    aq = w.AqueousPhase(" ".join(w.AQUEOUS_SPECIES))
    aq.setActivityModel(w.AQUEOUS_ACTIVITY_MODEL())
    mins = w.make_mineral_phases()
    gas = w.GaseousPhase("CO2")
    gas.setActivityModel(w.ActivityModelIdealGas())
    system = w.ChemicalSystem(database, aq, mins, gas)

    specs = w.EquilibriumSpecs(system)
    specs.temperature()
    specs.pressure()
    specs.pH()
    specs.fugacity("CO2")

    solver = w.make_equilibrium_solver(system, specs)
    opts = w.EquilibriumOptions()
    if hasattr(opts, "hessian") and hasattr(w, "GibbsHessian"):
        opts.hessian = w.GibbsHessian.Exact
    if hasattr(opts, "warmstart"):
        opts.warmstart = True
    solver.setOptions(opts)
    conditions = w.EquilibriumConditions(specs)

    logfCO2_vals = [-5.0, -4.0, -3.0, -2.0, -1.0]
    total = len(Tvals) * len(pHvals) * len(logfCO2_vals)

    targets = {
        "Wlm+hem+cc": lambda st: (
            present(st, "Wlm") and present(st, "hem") and present(st, "cc")
        ),
        "Wlm+hem+cc+q": lambda st: (
            present(st, "Wlm")
            and present(st, "hem")
            and present(st, "cc")
            and present(st, "q")
        ),
        "RELAXED_Wlm+hem+(cc|mag)+q": lambda st: (
            present(st, "Wlm")
            and present(st, "hem")
            and (present(st, "cc") or present(st, "mag"))
            and present(st, "q")
        ),
    }
    hits = {k: [] for k in targets}
    solved = 0

    for T in Tvals:
        state = w.make_base_state(system)
        try:
            state.set("CO2", 1.0e-20, "mol")
        except Exception:
            pass

        for pH in pHvals:
            for lfCO2 in logfCO2_vals:
                conditions.temperature(float(T), "celsius")
                conditions.pressure(float(Pbar), "bar")
                conditions.pH(float(pH))
                conditions.fugacity("CO2", float(10.0**lfCO2), "bar")

                result = solver.solve(state, conditions)
                if not result.succeeded():
                    continue

                solved += 1
                point = (float(T), float(pH), float(lfCO2))
                for name, check in targets.items():
                    if check(state):
                        hits[name].append(point)

    print(f"CARB_SOLVED={solved}/{total}")
    summarize_hits("CarbonateBuffered", hits, lambda p: f"logfCO2={p[2]:.1f}")


if __name__ == "__main__":
    run_redox_scan()
    run_carb_scan()
