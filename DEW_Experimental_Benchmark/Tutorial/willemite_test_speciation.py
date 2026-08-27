"""
Test if aqueous speciation complexity is the root cause.
Test with increasingly simplified aqueous species lists.
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


def make_state_with_silica_proxy(
    module, system, log_sio2_proxy_mol, aqueous_species_list
):
    # Create state with specified aqueous species
    state = module.create_equilibrium_state(system)

    # Set base species
    state.set("H2O", 55.5, "mol")

    # Set silica proxy
    silica_moles = float(10.0**log_sio2_proxy_mol)
    state.set(module.SIO2_PROXY_SPECIES, silica_moles, "mol")

    # Set minimal Zn seed
    for sp in aqueous_species_list:
        if "Zn" in sp or "Wlm" in sp:
            state.set(sp, 1.0e-20, "mol")

    return state


def test_speciation_set(module, spec_name, aqueous_species, pH_range):
    """Test convergence with a specific aqueous speciation set."""
    print(f"\n{spec_name}:")
    print(f"  Aqueous species count: {len(aqueous_species)}")
    print(f"  Zn species: {[s for s in aqueous_species if 'Zn' in s]}")

    # Backup originals
    original_aqueous = module.AQUEOUS_SPECIES
    original_zn_minerals = module.COMPETING_ZN_MINERALS

    # Set new speciation
    module.AQUEOUS_SPECIES = aqueous_species
    module.COMPETING_ZN_MINERALS = ["Wlm", "Znc"]  # Keep minimal mineral set

    try:
        system = module.build_tutorial_system()
    except Exception as e:
        print(f"  ✗ Failed to build system: {e}")
        module.AQUEOUS_SPECIES = original_aqueous
        module.COMPETING_ZN_MINERALS = original_zn_minerals
        return 0, len(pH_range)

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
        state = make_state_with_silica_proxy(
            module, system, FIXED_LOG_SIO2_PROXY_MOL, aqueous_species
        )
        conditions.temperature(float(TEMPERATURE_C), "celsius")
        conditions.pressure(float(pressure_bar), "bar")
        conditions.lgActivity("H+", float(-pH_value))

        result = solver.solve(state, conditions)
        if result.succeeded():
            converged += 1
        else:
            failed_list.append(pH_value)

    # Restore
    module.AQUEOUS_SPECIES = original_aqueous
    module.COMPETING_ZN_MINERALS = original_zn_minerals

    pct = 100.0 * converged / len(pH_range)
    print(f"  Convergence: {converged}/{len(pH_range)} ({pct:.1f}%)")
    if failed_list:
        print(f"  Failed at pH: {[f'{p:.1f}' for p in failed_list[:3]]}")

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
    print("Testing Aqueous Speciation Complexity Impact on Convergence")
    print("=" * 70)

    # Define speciation sets
    speciation_sets = {
        "FULL (current)": module.AQUEOUS_SPECIES,
        "MINIMAL (only essential Zn + Si/S)": [
            "H2O",
            "H+",
            "OH-",
            "Na+",
            "Zn2+",  # Only main Zn species
            "ZnOH+",  # And its main complex
            "SiO2,aq",  # Silicon
            "HSiO3-",
            "HS-",  # Sulfur
            "Cl-",  # Counter-ion
        ],
        "VERY_MINIMAL (only Zn2+ + essentials)": [
            "H2O",
            "H+",
            "OH-",
            "Na+",
            "Zn2+",  # Only divalent Zn
            "SiO2,aq",  # Silicon
            "HSiO3-",
            "HS-",  # Sulfur
            "Cl-",  # Counter-ion
        ],
    }

    results = {}
    for spec_name, aqueous_list in speciation_sets.items():
        converged, total = test_speciation_set(
            module, spec_name, aqueous_list, pH_range
        )
        results[spec_name] = (converged, total)

    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    for spec_name, (converged, total) in results.items():
        pct = 100.0 * converged / total if total > 0 else 0
        print(f"{spec_name:40s}: {converged:2d}/{total} ({pct:5.1f}%)")


if __name__ == "__main__":
    main()
