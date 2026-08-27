"""
Test improved solver strategies for high-pH regions.
Focus on pH ≥ 7 where failures cluster.
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
FIXED_LOG_SIO2_PROXY_MOL = -3.0

PH_MIN = 3.0
PH_MAX = 10.0
LOG_SIO2_PROXY_MOL_MIN = -8.0
LOG_SIO2_PROXY_MOL_MAX = -1.0


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


def make_state_with_silica_proxy(module, system, log_sio2_proxy_mol):
    state = module.make_base_state(system)
    silica_moles = float(10.0**log_sio2_proxy_mol)
    try:
        state.set(module.SIO2_PROXY_SPECIES, silica_moles, "mol")
    except Exception:
        overrides = dict(module.INITIAL_SPECIES_AMOUNTS_MOL)
        overrides[module.SIO2_PROXY_SPECIES] = silica_moles
        module.apply_species_amount_overrides(state, overrides)
    return state


def test_strategy_1_independent_retries(module):
    """Test: Completely independent retry (fresh state, not warm-start)."""
    print("\n=== Strategy 1: Independent Retries (both attempts from fresh state) ===")
    system = module.build_tutorial_system()
    specs = module.EquilibriumSpecs(system)
    specs.temperature()
    specs.pressure()
    specs.lgActivity("H+")

    solver = module.make_equilibrium_solver(system, specs)
    conditions = module.EquilibriumConditions(specs)

    pH_values = [7.0, 8.0, 9.0, 10.0]  # Test only high-pH
    log_sio2_proxy_values = [-3.0]  # Single silica level
    pressure_bar = PRESSURE_KBAR * 1000.0

    converged = 0
    for log_sio2_proxy_mol in log_sio2_proxy_values:
        for pH_value in pH_values:
            # First attempt
            state1 = make_state_with_silica_proxy(module, system, log_sio2_proxy_mol)
            conditions.temperature(float(TEMPERATURE_C), "celsius")
            conditions.pressure(float(pressure_bar), "bar")
            conditions.lgActivity("H+", float(-pH_value))
            result1 = solver.solve(state1, conditions)

            # Second attempt (truly independent, fresh state)
            state2 = make_state_with_silica_proxy(module, system, log_sio2_proxy_mol)
            conditions.temperature(float(TEMPERATURE_C), "celsius")
            conditions.pressure(float(pressure_bar), "bar")
            conditions.lgActivity("H+", float(-pH_value))
            result2 = solver.solve(state2, conditions)

            if result1.succeeded() or result2.succeeded():
                converged += 1
                print(
                    f"  pH={pH_value}: ✓ (attempt1={result1.succeeded()}, attempt2={result2.succeeded()})"
                )
            else:
                print(f"  pH={pH_value}: ✗ (both failed)")

    print(
        f"Result: {converged}/{len(pH_values)} converged ({100.0 * converged / len(pH_values):.1f}%)"
    )


def test_strategy_2_relaxed_tolerances(module):
    """Test: Multiple retry attempts with possibly different convergence behavior."""
    print("\n=== Strategy 2: Multiple Independent Retry Attempts ===")
    system = module.build_tutorial_system()
    specs = module.EquilibriumSpecs(system)
    specs.temperature()
    specs.pressure()
    specs.lgActivity("H+")

    solver = module.make_equilibrium_solver(system, specs)
    conditions = module.EquilibriumConditions(specs)

    pH_values = [7.0, 8.0, 9.0, 10.0]  # Test only high-pH
    log_sio2_proxy_values = [-3.0]
    pressure_bar = PRESSURE_KBAR * 1000.0

    converged_by_attempt = {1: 0, 2: 0, 3: 0}

    for log_sio2_proxy_mol in log_sio2_proxy_values:
        for pH_value in pH_values:
            for attempt in [1, 2, 3]:
                state = make_state_with_silica_proxy(module, system, log_sio2_proxy_mol)
                conditions.temperature(float(TEMPERATURE_C), "celsius")
                conditions.pressure(float(pressure_bar), "bar")
                conditions.lgActivity("H+", float(-pH_value))

                result = solver.solve(state, conditions)
                if result.succeeded():
                    converged_by_attempt[attempt] += 1
                    print(f"  pH={pH_value}, attempt {attempt}: ✓")
                    break
            else:
                print(f"  pH={pH_value}: ✗ (all 3 attempts failed)")

    total_converged = sum(v for v in converged_by_attempt.values() if v > 0)
    print(f"\nResult breakdown:")
    for attempt, count in sorted(converged_by_attempt.items()):
        print(f"  Converged on attempt {attempt}: {count}")
    print(f"Total: {total_converged}/{len(pH_values)}")


def test_strategy_3_scaled_constraints(module):
    """Test: Use slightly lower pH instead of exact constraint."""
    print("\n=== Strategy 3: Slightly Relaxed pH Constraint (pH-0.3) ===")
    system = module.build_tutorial_system()
    specs = module.EquilibriumSpecs(system)
    specs.temperature()
    specs.pressure()
    specs.lgActivity("H+")

    solver = module.make_equilibrium_solver(system, specs)
    conditions = module.EquilibriumConditions(specs)

    pH_values = [7.0, 8.0, 9.0, 10.0]  # Test only high-pH
    log_sio2_proxy_values = [-3.0]
    pressure_bar = PRESSURE_KBAR * 1000.0

    converged = 0
    for log_sio2_proxy_mol in log_sio2_proxy_values:
        for pH_requested in pH_values:
            # Use pH-0.3 as actual constraint
            pH_actual = pH_requested - 0.3

            state = make_state_with_silica_proxy(module, system, log_sio2_proxy_mol)
            conditions.temperature(float(TEMPERATURE_C), "celsius")
            conditions.pressure(float(pressure_bar), "bar")
            conditions.lgActivity("H+", float(-pH_actual))

            result = solver.solve(state, conditions)
            if result.succeeded():
                converged += 1
                actual_pH = -np.log10(float(state.speciesActivityLn("H+")))
                print(
                    f"  pH_requested={pH_requested}: ✓ (solved at pH={actual_pH:.2f})"
                )
            else:
                print(f"  pH_requested={pH_requested}: ✗")

    print(
        f"Result: {converged}/{len(pH_values)} converged ({100.0 * converged / len(pH_values):.1f}%)"
    )


def main():
    module = load_tutorial_module(TUTORIAL_PATH)

    try:
        module.Warnings.disable(906)
    except Exception:
        pass

    module.validate_user_inputs()
    module.USE_COMPETING_ZN_MINERALS = True

    print("=" * 70)
    print("Testing solver improvement strategies for high-pH failures")
    print("=" * 70)

    test_strategy_1_independent_retries(module)
    test_strategy_2_relaxed_tolerances(module)
    test_strategy_3_scaled_constraints(module)

    print("\n" + "=" * 70)
    print("Strategy recommendation based on results above")
    print("=" * 70)


if __name__ == "__main__":
    main()
