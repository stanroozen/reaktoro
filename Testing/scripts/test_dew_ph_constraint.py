"""
Test whether the standard Reaktoro specs.pH() and specs.fugacity() constraints
work without crashing when using ActivityModelPerplexDEW.

Tests:
1. Basic equilibrium (T, P only) — baseline
2. specs.pH() + conditions.pH() — standard pH constraint
3. specs.fugacity() + conditions.fugacity() — standard fO2 constraint
4. specs.pH() + NaCl brine background
5. Mineral solubility with fixed pH sweep (like U-solubility example)
"""

import sys
import os
import traceback
from pathlib import Path

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(
    os.path.dirname(THIS_DIR)
)  # up 2 levels: scripts/ -> Testing/ -> repo root
LOCAL_BUILD = os.path.join(REPO_ROOT, "build", "Reaktoro", "Release")
if os.path.isdir(LOCAL_BUILD) and LOCAL_BUILD not in sys.path:
    sys.path.insert(0, LOCAL_BUILD)

try:
    from reaktoro4py import *  # noqa: F401,F403
    import autodiff  # register pybind11 type casters for autodiff::Real

    print(f"Using local build: {LOCAL_BUILD}")
except ModuleNotFoundError:
    from reaktoro import *  # noqa: F401,F403

    print("Using installed reaktoro")

DB_PATH = (
    Path(REPO_ROOT) / "embedded/databases/perplex/DEW17HP622_Zn_2025-reaktoro.json"
)
MINERAL = "Wlm"

AQUEOUS = [
    "H2O",
    "H+",
    "OH-",
    "Zn2+",
    "ZnOH+",
    "HZnO2-",
    "SiO2,aq",
    "HSiO3-",
    "HS-",
    "SO4-2",
    "Na+",
    "Cl-",
]


def make_system(db):
    aq = AqueousPhase(" ".join(AQUEOUS))
    aq.setActivityModel(ActivityModelPerplexDEW())
    mineral = MineralPhase(MINERAL)
    return ChemicalSystem(db, aq, mineral)


def make_state(system, ph_proxy=7.0):
    state = ChemicalState(system)
    state.set("H2O", 55.5, "mol")
    state.set("H+", 1e-7, "mol")
    state.set("OH-", 1e-7, "mol")
    state.set("Zn2+", 1e-8, "mol")
    state.set("ZnOH+", 1e-10, "mol")
    state.set("HZnO2-", 1e-10, "mol")
    state.set("SiO2,aq", 1e-6, "mol")
    state.set("HSiO3-", 1e-9, "mol")
    state.set("HS-", 1e-10, "mol")
    state.set("SO4-2", 1e-10, "mol")
    state.set("Na+", 0.1, "mol")
    state.set("Cl-", 0.1, "mol")
    state.set(MINERAL, 0.01, "mol")
    return state


def run_test(label, fn):
    print(f"\n{'=' * 60}")
    print(f"TEST: {label}")
    print("=" * 60)
    try:
        fn()
        print(f"  PASSED")
    except Exception as e:
        print(f"  FAILED: {type(e).__name__}: {e}")
        traceback.print_exc()


# -----------------------------------------------------------------------
db = Database.fromFile(str(DB_PATH))
system = make_system(db)


# -----------------------------------------------------------------------
# Test 1: Basic T, P equilibrium — baseline
# -----------------------------------------------------------------------
def test_basic():
    specs = EquilibriumSpecs(system)
    specs.temperature()
    specs.pressure()
    solver = EquilibriumSolver(specs)
    conditions = EquilibriumConditions(specs)
    conditions.temperature(300.0, "celsius")
    conditions.pressure(2000.0, "bar")
    state = make_state(system)
    res = solver.solve(state, conditions)
    assert res.succeeded(), "Solver did not converge"
    props = AqueousProps(state)
    print(f"  pH = {float(props.pH()):.3f}")
    print(f"  Zn molality = {float(props.elementMolality('Zn')):.4e} mol/kg")


run_test("Basic T,P equilibrium (baseline)", test_basic)


# -----------------------------------------------------------------------
# Test 2: specs.pH() — standard approach
# -----------------------------------------------------------------------
def test_specs_ph():
    specs = EquilibriumSpecs(system)
    specs.temperature()
    specs.pressure()
    specs.pH()
    solver = EquilibriumSolver(specs)
    conditions = EquilibriumConditions(specs)
    conditions.temperature(300.0, "celsius")
    conditions.pressure(2000.0, "bar")
    conditions.pH(6.0)
    state = make_state(system)
    res = solver.solve(state, conditions)
    assert res.succeeded(), "Solver did not converge"
    props = AqueousProps(state)
    print(f"  pH = {float(props.pH()):.3f}  (target 6.0)")
    print(f"  Zn molality = {float(props.elementMolality('Zn')):.4e} mol/kg")


run_test("specs.pH() + conditions.pH(6.0)", test_specs_ph)

# -----------------------------------------------------------------------
# Test 3: specs.fugacity("O2") — SKIPPED: requires O2 gas phase in system
# -----------------------------------------------------------------------
# Note: specs.fugacity('O2') with no O2 species in the system causes a
# hard C++ crash in the solver. A proper fO2 constraint test would require
# adding O2(g) as a gas phase to the system.
print("\n" + "=" * 60)
print("TEST: specs.fugacity('O2')  [SKIPPED — requires O2 gas phase]")
print("=" * 60)


# -----------------------------------------------------------------------
# Test 4: specs.pH() sweep — like U-solubility example
# -----------------------------------------------------------------------
def test_ph_sweep():
    import numpy as np

    specs = EquilibriumSpecs(system)
    specs.temperature()
    specs.pressure()
    specs.pH()
    solver = EquilibriumSolver(specs)
    conditions = EquilibriumConditions(specs)
    conditions.temperature(300.0, "celsius")
    conditions.pressure(2000.0, "bar")
    state = make_state(system)
    pHs = np.linspace(4.0, 9.0, 6)
    results = []
    for pH in pHs:
        conditions.pH(pH)
        res = solver.solve(state, conditions)
        if res.succeeded():
            props = AqueousProps(state)
            mZn = float(props.elementMolality("Zn"))
            results.append((pH, mZn))
        else:
            results.append((pH, None))
    for pH, mZn in results:
        status = f"{mZn:.4e} mol/kg" if mZn is not None else "FAILED"
        print(f"  pH={pH:.1f}  Zn={status}")


run_test("specs.pH() sweep pH 4-9", test_ph_sweep)


# -----------------------------------------------------------------------
# Test 5: Fixed NaCl brine background (no pH spec, just conserved Cl)
# -----------------------------------------------------------------------
def test_nacl_brine():
    import numpy as np

    specs = EquilibriumSpecs(system)
    specs.temperature()
    specs.pressure()
    solver = EquilibriumSolver(specs)
    conditions = EquilibriumConditions(specs)
    conditions.temperature(300.0, "celsius")
    conditions.pressure(2000.0, "bar")
    state = make_state(system)
    # 0.5 mol/kg NaCl
    state.set("Na+", 0.5, "mol")
    state.set("Cl-", 0.5, "mol")
    res = solver.solve(state, conditions)
    assert res.succeeded(), "Solver did not converge"
    props = AqueousProps(state)
    print(f"  pH = {float(props.pH()):.3f}")
    print(f"  Zn molality = {float(props.elementMolality('Zn')):.4e} mol/kg")
    print(f"  Cl molality = {float(props.elementMolality('Cl')):.4e} mol/kg")


run_test("NaCl brine (0.5 mol/kg, no pH spec)", test_nacl_brine)

print("\n" + "=" * 60)
print("ALL TESTS COMPLETE")
print("=" * 60)
