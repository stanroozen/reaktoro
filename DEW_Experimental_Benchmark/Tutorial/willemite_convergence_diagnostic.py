"""
Diagnostic script to analyze nonconvergence regions in Willemite diagrams.
Maps where solves fail and tests solver options to improve convergence.
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
TEMP_MIN_C = 50.0
TEMP_MAX_C = 500.0
LOG_FH2S_MIN = -20.0
LOG_FH2S_MAX = -6.0

GRID_NX = 8
GRID_NY = 8


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


def test_ph_sio2_convergence(module):
    """Test pH vs silica diagram and report failure locations."""
    print("\n=== pH vs Silica Convergence Diagnostic ===")
    system = module.build_tutorial_system()
    specs = module.EquilibriumSpecs(system)
    specs.temperature()
    specs.pressure()
    specs.lgActivity("H+")

    solver = module.make_equilibrium_solver(system, specs)
    conditions = module.EquilibriumConditions(specs)

    pH_values = np.linspace(PH_MIN, PH_MAX, GRID_NX)
    log_sio2_proxy_values = np.linspace(
        LOG_SIO2_PROXY_MOL_MIN, LOG_SIO2_PROXY_MOL_MAX, GRID_NY
    )
    pressure_bar = PRESSURE_KBAR * 1000.0

    converged = 0
    failed = 0
    failed_points = []
    convergence_grid = np.full((len(log_sio2_proxy_values), len(pH_values)), 0)

    for j, log_sio2_proxy_mol in enumerate(log_sio2_proxy_values):
        for i, pH_value in enumerate(pH_values):
            state = make_state_with_silica_proxy(module, system, log_sio2_proxy_mol)
            conditions.temperature(float(TEMPERATURE_C), "celsius")
            conditions.pressure(float(pressure_bar), "bar")
            conditions.lgActivity("H+", float(-pH_value))

            result = solver.solve(state, conditions)
            if result.succeeded():
                converged += 1
                convergence_grid[j, i] = 1
            else:
                failed += 1
                failed_points.append((pH_value, log_sio2_proxy_mol))

    print(
        f"Converged: {converged}/{converged + failed} ({100.0 * converged / (converged + failed):.1f}%)"
    )

    # Analyze failure pattern
    if failed_points:
        print(f"\nFailed points ({len(failed_points)} total):")

        # Group by pH
        ph_groups = {}
        for pH, log_sio2 in failed_points:
            pH_bin = f"{float(pH):.1f}"
            if pH_bin not in ph_groups:
                ph_groups[pH_bin] = []
            ph_groups[pH_bin].append(log_sio2)

        print("  Failures by pH:")
        for pH_bin in sorted(ph_groups.keys()):
            log_sio2s = ph_groups[pH_bin]
            print(
                f"    pH={pH_bin}: {len(log_sio2s)} failures at log10(SiO2)={sorted(log_sio2s)}"
            )

        # Group by silica
        sio2_groups = {}
        for pH, log_sio2 in failed_points:
            sio2_bin = f"{float(log_sio2):.1f}"
            if sio2_bin not in sio2_groups:
                sio2_groups[sio2_bin] = []
            sio2_groups[sio2_bin].append(pH)

        print("\n  Failures by silica:")
        for sio2_bin in sorted(sio2_groups.keys()):
            phs = sio2_groups[sio2_bin]
            print(f"    log10(SiO2)={sio2_bin}: pH range {min(phs):.1f}-{max(phs):.1f}")

    # Print convergence grid visualization
    print("\n  Convergence pattern (1=converged, 0=failed):")
    print("    pH:", " ".join(f"{p:.1f}" for p in pH_values))
    for j in range(len(log_sio2_proxy_values) - 1, -1, -1):
        row_str = f"    log10(SiO2)={log_sio2_proxy_values[j]:.1f}: " + "".join(
            str(convergence_grid[j, i]) for i in range(len(pH_values))
        )
        print(row_str)


def test_temperature_ph_convergence(module):
    """Test Temperature vs pH diagram and report failure locations."""
    print("\n=== Temperature vs pH Convergence Diagnostic ===")
    system = module.build_tutorial_system()
    specs = module.EquilibriumSpecs(system)
    specs.temperature()
    specs.pressure()
    specs.lgActivity("H+")

    solver = module.make_equilibrium_solver(system, specs)
    conditions = module.EquilibriumConditions(specs)

    pH_values = np.linspace(PH_MIN, PH_MAX, GRID_NX)
    temperature_values = np.linspace(TEMP_MIN_C, TEMP_MAX_C, GRID_NY)
    pressure_bar = PRESSURE_KBAR * 1000.0

    converged = 0
    failed = 0
    failed_points = []
    convergence_grid = np.full((len(temperature_values), len(pH_values)), 0)

    for j, temperature_c in enumerate(temperature_values):
        for i, pH_value in enumerate(pH_values):
            state = make_state_with_silica_proxy(
                module, system, FIXED_LOG_SIO2_PROXY_MOL
            )
            conditions.temperature(float(temperature_c), "celsius")
            conditions.pressure(float(pressure_bar), "bar")
            conditions.lgActivity("H+", float(-pH_value))

            result = solver.solve(state, conditions)
            if result.succeeded():
                converged += 1
                convergence_grid[j, i] = 1
            else:
                failed += 1
                failed_points.append((temperature_c, pH_value))

    print(
        f"Converged: {converged}/{converged + failed} ({100.0 * converged / (converged + failed):.1f}%)"
    )

    if failed_points:
        print(f"Failed points ({len(failed_points)} total):")

        # Group by pH
        ph_groups = {}
        for T, pH in failed_points:
            pH_bin = f"{float(pH):.1f}"
            if pH_bin not in ph_groups:
                ph_groups[pH_bin] = []
            ph_groups[pH_bin].append(T)

        print("  Failures by pH:")
        for pH_bin in sorted(ph_groups.keys()):
            temps = ph_groups[pH_bin]
            print(
                f"    pH={pH_bin}: {len(temps)} failures at T={min(temps):.0f}-{max(temps):.0f}°C"
            )

        # Group by temperature
        temp_groups = {}
        for T, pH in failed_points:
            T_bin = f"{int(T)}"
            if T_bin not in temp_groups:
                temp_groups[T_bin] = []
            temp_groups[T_bin].append(pH)

        print("\n  Failures by temperature:")
        for T_bin in sorted(temp_groups.keys(), key=int):
            phs = temp_groups[T_bin]
            print(f"    T={T_bin}°C: pH range {min(phs):.1f}-{max(phs):.1f}")


def test_ph_fh2s_convergence(module):
    """Test pH vs fH2S diagram and report failure locations."""
    print("\n=== pH vs fH2S Convergence Diagnostic ===")
    system = module.build_system_with_h2s_gas_phase()
    specs = module.EquilibriumSpecs(system)
    specs.temperature()
    specs.pressure()
    specs.lgActivity("H+")
    specs.fugacity(module.SULFUR_FUGACITY_SPECIES)

    solver = module.make_equilibrium_solver(system, specs)
    conditions = module.EquilibriumConditions(specs)

    pH_values = np.linspace(PH_MIN, PH_MAX, GRID_NX)
    logf_values = np.linspace(LOG_FH2S_MIN, LOG_FH2S_MAX, GRID_NY)
    pressure_bar = PRESSURE_KBAR * 1000.0

    converged = 0
    failed = 0
    failed_points = []
    convergence_grid = np.full((len(logf_values), len(pH_values)), 0)

    for j, logf in enumerate(logf_values):
        for i, pH_value in enumerate(pH_values):
            state = make_state_with_silica_proxy(
                module, system, FIXED_LOG_SIO2_PROXY_MOL
            )
            try:
                state.set(module.SULFUR_FUGACITY_SPECIES, 1.0e-20, "mol")
            except Exception:
                pass

            conditions.temperature(float(TEMPERATURE_C), "celsius")
            conditions.pressure(float(pressure_bar), "bar")
            conditions.lgActivity("H+", float(-pH_value))
            conditions.fugacity(
                module.SULFUR_FUGACITY_SPECIES, float(10.0**logf), "bar"
            )

            result = solver.solve(state, conditions)
            if result.succeeded():
                converged += 1
                convergence_grid[j, i] = 1
            else:
                failed += 1
                failed_points.append((pH_value, logf))

    print(
        f"Converged: {converged}/{converged + failed} ({100.0 * converged / (converged + failed):.1f}%)"
    )

    if failed_points:
        print(f"Failed points ({len(failed_points)} total):")

        # Group by fH2S
        fh2s_groups = {}
        for pH, logf in failed_points:
            logf_bin = f"{float(logf):.1f}"
            if logf_bin not in fh2s_groups:
                fh2s_groups[logf_bin] = []
            fh2s_groups[logf_bin].append(pH)

        print("  Failures by log10(fH2S):")
        for logf_bin in sorted(fh2s_groups.keys()):
            phs = fh2s_groups[logf_bin]
            print(f"    log10(fH2S)={logf_bin}: pH range {min(phs):.1f}-{max(phs):.1f}")

        # Group by pH
        ph_groups = {}
        for pH, logf in failed_points:
            pH_bin = f"{float(pH):.1f}"
            if pH_bin not in ph_groups:
                ph_groups[pH_bin] = []
            ph_groups[pH_bin].append(logf)

        print("\n  Failures by pH:")
        for pH_bin in sorted(ph_groups.keys()):
            logfs = ph_groups[pH_bin]
            print(
                f"    pH={pH_bin}: log10(fH2S) range {min(logfs):.1f}-{max(logfs):.1f}"
            )


def main():
    module = load_tutorial_module(TUTORIAL_PATH)

    try:
        module.Warnings.disable(906)
    except Exception:
        pass

    module.validate_user_inputs()
    module.USE_COMPETING_ZN_MINERALS = True

    test_ph_sio2_convergence(module)
    test_temperature_ph_convergence(module)
    test_ph_fh2s_convergence(module)

    print("\n=== Summary ===")
    print("Diagnostics complete. Check failure patterns above to identify:")
    print("1. Which parameter ranges fail most (pH, T, fH2S, silica)")
    print("2. Whether higher max_iterations helps")
    print("3. Whether failures cluster at phase boundaries or extreme conditions")


if __name__ == "__main__":
    main()
