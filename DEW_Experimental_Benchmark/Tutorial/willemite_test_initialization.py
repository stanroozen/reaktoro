"""
Test initialization strategies: better starting guesses for high-pH.
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


def make_state_default(module, system, log_sio2_proxy_mol):
    """Default initialization (aqueous-heavy)."""
    state = module.make_base_state(system)
    silica_moles = float(10.0**log_sio2_proxy_mol)
    try:
        state.set(module.SIO2_PROXY_SPECIES, silica_moles, "mol")
    except Exception:
        overrides = dict(module.INITIAL_SPECIES_AMOUNTS_MOL)
        overrides[module.SIO2_PROXY_SPECIES] = silica_moles
        module.apply_species_amount_overrides(state, overrides)
    return state


def make_state_precipitate_initialized(module, system, log_sio2_proxy_mol):
    """High-pH initialization: pre-precipitate most Zn."""
    state = make_state_default(module, system, log_sio2_proxy_mol)

    # Try to pre-set Zn as precipitate
    try:
        # Set most Zn as Znc (zincite)
        total_zn_moles = 20.0  # From tutorial
        state.set("Znc", total_zn_moles * 0.95, "mol")  # 95% as precipitate
        state.set("Zn2+", total_zn_moles * 0.01, "mol")  # 1% aqueous
        state.set("ZnOH+", total_zn_moles * 0.04, "mol")  # 4% as hydroxide complex
    except Exception:
        pass

    return state


def make_state_charge_balanced(module, system, log_sio2_proxy_mol, pH_value):
    """Initialize with charge-balanced aqueous speciation for target pH."""
    state = make_state_default(module, system, log_sio2_proxy_mol)

    # At high pH, ZnOH+ dominates, so pre-seed the state with that expectation
    try:
        if pH_value >= 7:
            # High pH: Zn exists mostly as ZnOH+ or even Zn(OH)3-
            state.set("Zn2+", 0.1, "mol")
            state.set("ZnOH+", 1.0, "mol")
            state.set("ZnO", 1.0, "mol")
            state.set("HZnO2-", 1.0, "mol")
    except Exception:
        pass

    return state


def test_initialization_strategy(module, strategy_name, init_func, pH_range):
    """Test convergence with a specific initialization strategy."""
    print(f"\n{strategy_name}:")

    system = module.build_tutorial_system()
    specs = module.EquilibriumSpecs(system)
    specs.temperature()
    specs.pressure()
    specs.lgActivity("H+")

    solver = module.make_equilibrium_solver(system, specs)
    conditions = module.EquilibriumConditions(specs)
    pressure_bar = PRESSURE_KBAR * 1000.0

    converged = 0
    failed_list = []

    for pH_value in pH_range:
        # Use custom initialization
        if init_func == "charge_balanced":
            state = make_state_charge_balanced(
                module, system, FIXED_LOG_SIO2_PROXY_MOL, pH_value
            )
        elif init_func == "precipitate":
            state = make_state_precipitate_initialized(
                module, system, FIXED_LOG_SIO2_PROXY_MOL
            )
        else:
            state = make_state_default(module, system, FIXED_LOG_SIO2_PROXY_MOL)

        conditions.temperature(float(TEMPERATURE_C), "celsius")
        conditions.pressure(float(pressure_bar), "bar")
        conditions.lgActivity("H+", float(-pH_value))

        result = solver.solve(state, conditions)
        if result.succeeded():
            converged += 1
        else:
            failed_list.append(pH_value)

    pct = 100.0 * converged / len(pH_range)
    print(f"  Convergence: {converged}/{len(pH_range)} ({pct:.1f}%)")
    if failed_list:
        print(f"  Failed at pH: {[f'{p:.1f}' for p in failed_list[:4]]}")

    return converged, len(pH_range)


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
    print("Testing Initialization Strategies for High-pH Convergence")
    print("=" * 70)

    strategies = {
        "DEFAULT (aqueous-heavy)": "default",
        "PRECIPITATE-INIT (Zn pre-precipitated)": "precipitate",
        "CHARGE-BALANCED (pH-aware speciation)": "charge_balanced",
    }

    results = {}
    for strategy_name, init_func in strategies.items():
        converged, total = test_initialization_strategy(
            module, strategy_name, init_func, pH_range
        )
        results[strategy_name] = (converged, total)

    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    for strategy_name, (converged, total) in results.items():
        pct = 100.0 * converged / total if total > 0 else 0
        print(f"{strategy_name:50s}: {converged:2d}/{total} ({pct:5.1f}%)")

    print("\n" + "=" * 70)
    print("INTERPRETATION")
    print("=" * 70)
    print("If all strategies fail equally, the issue is thermodynamic infeasibility.")
    print(
        "If CHARGE-BALANCED works best, initialization matters—we can use pH-aware preconditioning."
    )
    print(
        "If PRECIPITATE-INIT works best, we need pre-equilibration at near-solution pH before solving."
    )


if __name__ == "__main__":
    main()
