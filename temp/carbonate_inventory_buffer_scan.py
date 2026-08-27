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
system = w.build_tutorial_system()

specs = w.EquilibriumSpecs(system)
specs.temperature()
specs.pressure()
specs.pH()

solver = w.make_equilibrium_solver(system, specs)
opts = w.EquilibriumOptions()
if hasattr(opts, "hessian") and hasattr(w, "GibbsHessian"):
    opts.hessian = w.GibbsHessian.Exact
if hasattr(opts, "warmstart"):
    opts.warmstart = True
solver.setOptions(opts)
conditions = w.EquilibriumConditions(specs)

Tvals = [200.0, 250.0, 300.0, 350.0, 400.0]
pHvals = np.linspace(4.0, 10.0, 7)
log_dic_vals = [-8.0, -7.0, -6.0, -5.0, -4.0, -3.0, -2.0]
Pbar = 2000.0
TH = 1.0e-8


def present(st, name):
    try:
        return float(st.speciesAmount(name)) > TH
    except Exception:
        return False


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
    "Wlm+hem+cc+q+Znc": lambda st: (
        present(st, "Wlm")
        and present(st, "hem")
        and present(st, "cc")
        and present(st, "q")
        and present(st, "Znc")
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
total = len(Tvals) * len(pHvals) * len(log_dic_vals)

for T in Tvals:
    for pH in pHvals:
        for ldic in log_dic_vals:
            state = w.make_base_state(system)

            # Carbonate inventory proxy for buffering in absence of stable CO2 fugacity solve path.
            dic = 10.0 ** float(ldic)
            state.set("HCO3-", dic, "mol")
            state.set("CO3-2", 0.1 * dic, "mol")

            conditions.temperature(float(T), "celsius")
            conditions.pressure(float(Pbar), "bar")
            conditions.pH(float(pH))

            result = solver.solve(state, conditions)
            if not result.succeeded():
                continue

            solved += 1
            point = (float(T), float(pH), float(ldic))
            for name, check in targets.items():
                if check(state):
                    hits[name].append(point)

print(f"SOLVED={solved}/{total}")
for name, pts in hits.items():
    print(f"\n{name}: {len(pts)} points")
    if not pts:
        continue

    Ts = [p[0] for p in pts]
    pHs = [p[1] for p in pts]
    lds = [p[2] for p in pts]
    print(f"  T range (C): {min(Ts):.1f} to {max(Ts):.1f}")
    print(f"  pH range: {min(pHs):.2f} to {max(pHs):.2f}")
    print(f"  log10(DIC mol) range: {min(lds):.1f} to {max(lds):.1f}")
    print("  first points:")
    for p in pts[:20]:
        print(f"    T={p[0]:.1f} C, pH={p[1]:.2f}, logDIC={p[2]:.1f}")
