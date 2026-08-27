import importlib.util
import os
import sys
import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
rel = os.path.join(REPO, "build", "Reaktoro", "Release")
if rel not in sys.path:
    sys.path.insert(0, rel)

TUTORIAL = os.path.join(
    REPO,
    "DEW_Experimental_Benchmark",
    "Tutorial",
    "willemite_solubility_tutorial_dew17hp622_zn.py",
)
spec = importlib.util.spec_from_file_location("w", TUTORIAL)
w = importlib.util.module_from_spec(spec)
spec.loader.exec_module(w)

# Build a stable system for T-pH scan without gas fugacity constraints.
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

# Seed gangue minerals so they can persist/disappear during equilibration.
w.INITIAL_SPECIES_AMOUNTS_MOL["hem"] = 2.0
w.INITIAL_SPECIES_AMOUNTS_MOL["cc"] = 2.0
w.INITIAL_SPECIES_AMOUNTS_MOL["q"] = 2.0

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
solver.setOptions(opts)
conditions = w.EquilibriumConditions(specs)

Tvals = np.linspace(100.0, 400.0, 16)
pHvals = np.linspace(0.0, 14.0, 29)
Pbar = 2000.0
TH = 1e-8


def present(st, name):
    try:
        return float(st.speciesAmount(name)) > TH
    except Exception:
        return False


combos = {
    "Wlm+hem+cc": lambda st: (
        present(st, "Wlm") and present(st, "hem") and present(st, "cc")
    ),
    "Wlm+hem+q": lambda st: (
        present(st, "Wlm") and present(st, "hem") and present(st, "q")
    ),
    "Wlm+cc+q": lambda st: (
        present(st, "Wlm") and present(st, "cc") and present(st, "q")
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
}

hits = {k: [] for k in combos}
solved = 0

for T in Tvals:
    for pH in pHvals:
        st = w.make_base_state(system)
        conditions.temperature(float(T), "celsius")
        conditions.pressure(float(Pbar), "bar")
        conditions.pH(float(pH))
        r = solver.solve(st, conditions)
        if not r.succeeded():
            continue
        solved += 1
        for name, check in combos.items():
            if check(st):
                hits[name].append((float(T), float(pH)))

print(f"SOLVED={solved}/{len(Tvals) * len(pHvals)}")
for name, pts in hits.items():
    print(f"\n{name}: {len(pts)} points")
    if pts:
        Ts = [p[0] for p in pts]
        pHs = [p[1] for p in pts]
        print(f"  T range: {min(Ts):.1f} to {max(Ts):.1f} C")
        print(f"  pH range: {min(pHs):.2f} to {max(pHs):.2f}")
        print("  first points:")
        for p in pts[:20]:
            print(f"    T={p[0]:.1f} C, pH={p[1]:.2f}")
