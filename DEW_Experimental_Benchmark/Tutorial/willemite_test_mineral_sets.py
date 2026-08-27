"""
Test reduced mineral sets to improve high-pH convergence.
Compare full 13-mineral set vs curated subsets.
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


def test_mineral_set(module, mineral_set_name, mineral_list, pH_range):
    """Test convergence with a specific mineral set."""
    print(f"\n{mineral_set_name}:")
    print(f"  Minerals: {mineral_list}")

    # Modify global to select this mineral set
    original_minerals = module.COMPETING_ZN_MINERALS
    module.COMPETING_ZN_MINERALS = mineral_list

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
        state = make_state_with_silica_proxy(module, system, FIXED_LOG_SIO2_PROXY_MOL)
        conditions.temperature(float(TEMPERATURE_C), "celsius")
        conditions.pressure(float(pressure_bar), "bar")
        conditions.lgActivity("H+", float(-pH_value))

        result = solver.solve(state, conditions)
        if result.succeeded():
            converged += 1
        else:
            failed_list.append(pH_value)

    # Restore original
    module.COMPETING_ZN_MINERALS = original_minerals

    pct = 100.0 * converged / len(pH_range)
    print(f"  Convergence: {converged}/{len(pH_range)} ({pct:.1f}%)")
    if failed_list:
        print(f"  Failed at pH: {[f'{p:.1f}' for p in failed_list[:5]]}")
        if len(failed_list) > 5:
            print(f"            ... and {len(failed_list) - 5} more")

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
    print("Testing Reduced Mineral Sets for High-pH Convergence")
    print("=" * 70)

    # Test different mineral subsets
    mineral_sets = {
        "FULL (13 minerals)": [
            "Wlm",
            "Sph",
            "Smth",
            "Znc",
            "Znks",
            "Ghn",
            "Frk",
            "Hrds",
            "ZnSp",
            "HZnc",
            "Zn",
            "Wrt",
            "Zn-St",
        ],
        "CORE (5 minerals)": [
            "Wlm",  # Willemite - primary focus
            "Sph",  # Sphalerite - common sulfide
            "Smth",  # Smithsonite - common carbonate
            "Znc",  # Zincite - oxide
            "HZnc",  # Hydrozincite - hydroxide
        ],
        "MINIMAL (3 minerals)": [
            "Wlm",  # Willemite
            "Znc",  # Zincite - most stable at high pH
            "HZnc",  # Hydrozincite - forms at high pH + H2O
        ],
        "OXIDE-ONLY (2 minerals)": [
            "Wlm",  # Willemite
            "Znc",  # Zincite
        ],
    }

    results_summary = {}
    for set_name, mineral_list in mineral_sets.items():
        converged, total = test_mineral_set(module, set_name, mineral_list, pH_range)
        results_summary[set_name] = (converged, total)

    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    for set_name, (converged, total) in results_summary.items():
        pct = 100.0 * converged / total
        print(f"{set_name:30s}: {converged:2d}/{total} ({pct:5.1f}%)")

    print("\n" + "=" * 70)
    print("RECOMMENDATION")
    print("=" * 70)
    best_set = max(results_summary.items(), key=lambda x: x[1][0])
    print(f"Best convergence: {best_set[0]} ({best_set[1][0]}/{best_set[1][1]})")
    print(
        "\nTo apply: Update COMPETING_ZN_MINERALS in tutorial module to use the best set."
    )


if __name__ == "__main__":
    main()
