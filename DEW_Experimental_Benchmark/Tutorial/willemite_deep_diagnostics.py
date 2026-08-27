"""
Deep diagnostic: Compare converged (low pH) vs failed (high pH) system states.
Check mineral stability, charge balance, and convergence diagnostics.
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


def print_state_diagnostics(module, state, label):
    """Print detailed state information for debugging."""
    print(f"\n{label}")
    print("-" * 70)

    # Mineral amounts
    print("Minerals (amounts in mol):")
    for mineral_name in module.selected_mineral_names():
        try:
            amount = float(state.speciesAmount(mineral_name))
            if amount > 0:
                print(f"  {mineral_name}: {amount:.6e}")
        except Exception:
            pass

    # Total dissolved Zn
    try:
        solvent_name = module.infer_solvent_species_name(
            module.AQUEOUS_SPECIES,
            module.SOLVENT_SPECIES_NAME,
            module.INITIAL_SPECIES_AMOUNTS_MOL,
        )
        zn_molality = module.dissolved_element_molality(state, solvent_name)
        print(f"\nDissolved Zn molality: {zn_molality:.6e} mol/kg")
    except Exception as e:
        print(f"\nDissolved Zn: (error: {e})")

    # Key species amounts (aqueous)
    print("\nKey aqueous species (amounts in mol):")
    key_species = [
        "H+",
        "OH-",
        "H2O,aq",
        "SiO2,aq",
        "H2SiO3,aq",
        "HSiO3-",
        "ZnOH+",
        "Zn++",
        "ZnOH2,aq",
        "ZnO,aq",
        "Zn(OH)3-",
        "Zn(OH)4--",
        "NaCl,aq",
        "Na+",
        "Cl-",
    ]
    for sp_name in key_species:
        try:
            amount = float(state.speciesAmount(sp_name))
            if amount > 0:
                print(f"  {sp_name}: {amount:.6e}")
        except Exception:
            pass

    # Activity values
    print("\nKey species activities (log10):")
    key_species_act = ["H+", "OH-", "Zn++", "SiO2,aq"]
    for sp_name in key_species_act:
        try:
            lnact = float(state.speciesActivityLn(sp_name))
            log10act = lnact / np.log(10.0)
            print(f"  {sp_name}: {log10act:.4f}")
        except Exception:
            pass


def main():
    module = load_tutorial_module(TUTORIAL_PATH)

    try:
        module.Warnings.disable(906)
    except Exception:
        pass

    module.validate_user_inputs()
    module.USE_COMPETING_ZN_MINERALS = True

    system = module.build_tutorial_system()
    specs = module.EquilibriumSpecs(system)
    specs.temperature()
    specs.pressure()
    specs.lgActivity("H+")

    solver = module.make_equilibrium_solver(system, specs)
    conditions = module.EquilibriumConditions(specs)
    pressure_bar = PRESSURE_KBAR * 1000.0

    print("=" * 70)
    print("Deep System Diagnostics: Converged vs Failed Points")
    print("=" * 70)

    # Test converged point (low pH)
    print("\n>>> CONVERGED POINT (pH=3.0) <<<")
    state_converged = make_state_with_silica_proxy(
        module, system, FIXED_LOG_SIO2_PROXY_MOL
    )
    conditions.temperature(float(TEMPERATURE_C), "celsius")
    conditions.pressure(float(pressure_bar), "bar")
    conditions.lgActivity("H+", float(-3.0))  # pH=3
    result_low = solver.solve(state_converged, conditions)
    print(f"Solve result: {result_low.succeeded()}")
    print_state_diagnostics(module, state_converged, "State after solve (pH=3):")

    # Test failed point (high pH)
    print("\n\n>>> FAILED POINT (pH=10.0) <<<")
    state_failed = make_state_with_silica_proxy(
        module, system, FIXED_LOG_SIO2_PROXY_MOL
    )
    conditions.temperature(float(TEMPERATURE_C), "celsius")
    conditions.pressure(float(pressure_bar), "bar")
    conditions.lgActivity("H+", float(-10.0))  # pH=10
    result_high = solver.solve(state_failed, conditions)
    print(f"Solve result: {result_high.succeeded()}")
    print_state_diagnostics(
        module, state_failed, "State after attempted solve (pH=10):"
    )

    # Try intermediate pH
    print("\n\n>>> INTERMEDIATE POINT (pH=6.5) <<<")
    state_intermediate = make_state_with_silica_proxy(
        module, system, FIXED_LOG_SIO2_PROXY_MOL
    )
    conditions.temperature(float(TEMPERATURE_C), "celsius")
    conditions.pressure(float(pressure_bar), "bar")
    conditions.lgActivity("H+", float(-6.5))  # pH=6.5
    result_mid = solver.solve(state_intermediate, conditions)
    print(f"Solve result: {result_mid.succeeded()}")
    print_state_diagnostics(
        module, state_intermediate, "State after attempted solve (pH=6.5):"
    )

    # Check which minerals are available
    print("\n\n" + "=" * 70)
    print("Available Minerals in System:")
    print("=" * 70)
    try:
        for mineral in module.selected_mineral_names():
            print(f"  - {mineral}")
    except Exception as e:
        print(f"Error listing minerals: {e}")

    # Check USE_COMPETING_ZN_MINERALS flag
    print(f"\nUSE_COMPETING_ZN_MINERALS = {module.USE_COMPETING_ZN_MINERALS}")


if __name__ == "__main__":
    main()
