"""
Test if silica proxy constraint is incompatible with high pH.
Try solving without fixing silica proxy.
"""

import importlib.util
import os
import sys

import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(os.path.dirname(SCRIPT_DIR))
TUTORIAL_PATH = os.path.join(
    SCRIPT_DIR,
    "willemite_solubility_tutorial_dew17hp622_zn.py",
)
LOCAL_REAKTORO_RELEASE_DIR = os.path.join(REPO_ROOT, "build", "Reaktoro", "Release")
LOCAL_REAKTORO_PYD_DIR = os.path.join(
    REPO_ROOT,
    "build",
    "python",
    "package",
    "build",
    "lib",
    "reaktoro",
)

TEMPERATURE_C = 300.0
PRESSURE_KBAR = 2.0


def load_tutorial_module(path):
    if sys.platform.startswith("win") and os.path.isdir(LOCAL_REAKTORO_RELEASE_DIR):
        if LOCAL_REAKTORO_RELEASE_DIR not in sys.path:
            sys.path.insert(0, LOCAL_REAKTORO_RELEASE_DIR)
    if os.path.isdir(LOCAL_REAKTORO_PYD_DIR) and LOCAL_REAKTORO_PYD_DIR not in sys.path:
        sys.path.append(LOCAL_REAKTORO_PYD_DIR)

    spec = importlib.util.spec_from_file_location("willemite_tutorial", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_no_silica_constraint(module, pH_range):
    """Test convergence WITHOUT fixing silica proxy."""
    print("\nWithout silica proxy constraint:")

    system = module.build_tutorial_system()
    specs = module.EquilibriumSpecs(system)
    specs.temperature()
    specs.pressure()
    specs.lgActivity("H+")
    # NOTE: NO SiO2 activity constraint

    solver = module.make_equilibrium_solver(system, specs)
    conditions = module.EquilibriumConditions(specs)
    pressure_bar = PRESSURE_KBAR * 1000.0

    converged = 0

    for pH_value in pH_range:
        state = module.make_base_state(system)
        conditions.temperature(float(TEMPERATURE_C), "celsius")
        conditions.pressure(float(pressure_bar), "bar")
        conditions.lgActivity("H+", float(-pH_value))
        # NOT setting silica inventory

        result = solver.solve(state, conditions)
        if result.succeeded():
            converged += 1
            print(f"  pH={pH_value:.1f}: ✓")
        else:
            print(f"  pH={pH_value:.1f}: ✗")

    pct = 100.0 * converged / len(pH_range)
    print(f"\nConvergence (no silica fix): {converged}/{len(pH_range)} ({pct:.1f}%)")
    return converged


def test_with_silica_proxy(module, log_sio2_proxy_mol, pH_range):
    """Test convergence WITH fixed silica proxy (current approach)."""
    print(f"\nWith silica proxy = 10^{log_sio2_proxy_mol} mol:")

    system = module.build_tutorial_system()
    specs = module.EquilibriumSpecs(system)
    specs.temperature()
    specs.pressure()
    specs.lgActivity("H+")

    solver = module.make_equilibrium_solver(system, specs)
    conditions = module.EquilibriumConditions(specs)
    pressure_bar = PRESSURE_KBAR * 1000.0

    converged = 0

    for pH_value in pH_range:
        state = module.make_base_state(system)
        silica_moles = float(10.0**log_sio2_proxy_mol)
        try:
            state.set(module.SIO2_PROXY_SPECIES, silica_moles, "mol")
        except Exception:
            pass

        conditions.temperature(float(TEMPERATURE_C), "celsius")
        conditions.pressure(float(pressure_bar), "bar")
        conditions.lgActivity("H+", float(-pH_value))

        result = solver.solve(state, conditions)
        if result.succeeded():
            converged += 1
            print(f"  pH={pH_value:.1f}: ✓")
        else:
            print(f"  pH={pH_value:.1f}: ✗")

    pct = 100.0 * converged / len(pH_range)
    print(f"\nConvergence (with silica fix): {converged}/{len(pH_range)} ({pct:.1f}%)")
    return converged


def main():
    module = load_tutorial_module(TUTORIAL_PATH)

    try:
        module.Warnings.disable(906)
    except Exception:
        pass

    module.validate_user_inputs()
    module.USE_COMPETING_ZN_MINERALS = True

    pH_range = np.linspace(3.0, 10.0, 8)

    print("=" * 70)
    print("Testing if Silica Proxy Constraint Causes High-pH Infeasibility")
    print("=" * 70)

    print("\nTest 1: High pH WITHOUT any silica constraint")
    conv_no_silica = test_no_silica_constraint(module, pH_range)

    print("\n" + "=" * 70)
    print("Test 2: High pH WITH current silica proxy")
    conv_with_silica = test_with_silica_proxy(module, -3.0, pH_range)

    print("\n" + "=" * 70)
    print("RESULT")
    print("=" * 70)
    if conv_no_silica > conv_with_silica:
        print(f"✓ Silica proxy IS the problem: {conv_no_silica} vs {conv_with_silica}")
        print("   Solution: Use softer silica control (free speciation vs fixed proxy)")
    elif conv_no_silica == conv_with_silica:
        print(f"✗ Silica proxy is NOT the issue")
        print("   The problem is the fundamental pH constraint or DEW model at high pH")
    else:
        print(
            f"? Unexpected: silica proxy improved convergence ({conv_with_silica} vs {conv_no_silica})"
        )


if __name__ == "__main__":
    main()
