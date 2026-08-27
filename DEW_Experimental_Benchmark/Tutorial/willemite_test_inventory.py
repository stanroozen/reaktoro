"""
Test if reducing Zn inventory helps high-pH convergence.
At high pH, maybe 20 mol Zn is thermodynamically infeasible.
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


def make_state_with_inventory(module, system, log_sio2_proxy_mol, zn_inventory_factor):
    """Create state with scaled Zn inventory."""
    state = module.make_base_state(system)
    silica_moles = float(10.0**log_sio2_proxy_mol)
    try:
        state.set(module.SIO2_PROXY_SPECIES, silica_moles, "mol")
    except Exception:
        overrides = dict(module.INITIAL_SPECIES_AMOUNTS_MOL)
        overrides[module.SIO2_PROXY_SPECIES] = silica_moles
        module.apply_species_amount_overrides(state, overrides)

    # Scale Zn minerals
    for mineral_name in module.selected_mineral_names():
        try:
            current = float(state.speciesAmount(mineral_name, "mol"))
            if current > 0:
                state.set(mineral_name, current * zn_inventory_factor, "mol")
        except Exception:
            pass

    return state


def test_inventory(module, zn_factor, zn_factor_label, pH_range):
    """Test convergence with scaled Zn inventory."""
    print(f"\n{zn_factor_label}:")

    system = module.build_tutorial_system()
    specs = module.EquilibriumSpecs(system)
    specs.temperature()
    specs.pressure()
    specs.lgActivity("H+")

    solver = module.make_equilibrium_solver(system, specs)
    conditions = module.EquilibriumConditions(specs)
    pressure_bar = PRESSURE_KBAR * 1000.0

    converged = 0
    converged_phs = []
    failed_phs = []

    for pH_value in pH_range:
        state = make_state_with_inventory(
            module, system, FIXED_LOG_SIO2_PROXY_MOL, zn_factor
        )
        conditions.temperature(float(TEMPERATURE_C), "celsius")
        conditions.pressure(float(pressure_bar), "bar")
        conditions.lgActivity("H+", float(-pH_value))

        result = solver.solve(state, conditions)
        if result.succeeded():
            converged += 1
            converged_phs.append(pH_value)
        else:
            failed_phs.append(pH_value)

    pct = 100.0 * converged / len(pH_range)
    print(f"  Convergence: {converged}/{len(pH_range)} ({pct:.1f}%)")
    if converged_phs:
        print(f"  Converged at pH: {[f'{p:.1f}' for p in converged_phs]}")
    if failed_phs:
        print(f"  Failed at pH: {[f'{p:.1f}' for p in failed_phs]}")

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
    print("Testing Zn Inventory Scaling Effect on High-pH Convergence")
    print("=" * 70)
    print(f"Default Zn inventory (Wlm): 10 mol")

    inventories = {
        1.0: "1.0x (DEFAULT - 10 mol Zn)",
        0.5: "0.5x (5 mol Zn)",
        0.1: "0.1x (1 mol Zn)",
        0.01: "0.01x (0.1 mol Zn)",
        0.001: "0.001x (0.01 mol Zn)",
    }

    results = {}
    for factor, label in inventories.items():
        converged, total = test_inventory(module, factor, label, pH_range)
        results[label] = (converged, total)

    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    for label, (converged, total) in results.items():
        pct = 100.0 * converged / total if total > 0 else 0
        print(f"{label:35s}: {converged:2d}/{total} ({pct:5.1f}%)")

    print("\n" + "=" * 70)
    print("ANALYSIS")
    print("=" * 70)
    print("If convergence improves dramatically at lower Zn inventory, the issue is")
    print("saturation/precipitation creating a degenerate or ill-posed system.")
    print("If convergence stays constant, the issue is fundamental (model/database).")


if __name__ == "__main__":
    main()
