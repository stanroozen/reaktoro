"""
Mineral Solubility Calculation Utilities
=========================================

This module provides reusable, mineral-agnostic functions for calculating
mineral solubility using the Reaktoro library. Functions are designed to work
with any mineral/aqueous species combination and follow software design best practices.

Can be used for:
- Quartz solubility
- Calcite solubility
- Any other mineral + aqueous species system

Author: Reaktoro Tutorial
License: MIT
"""

import os
import sys
import numpy as np
import pandas as pd
import warnings

# Import Reaktoro with fallback to local build
try:
    from reaktoro import *
except ModuleNotFoundError:
    # Fallback to local build
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    ROOT_DIR = os.path.dirname(SCRIPT_DIR)
    PYD_DIR = os.path.join(ROOT_DIR, "build", "Reaktoro", "Release")
    if os.path.isdir(PYD_DIR) and PYD_DIR not in sys.path:
        sys.path.insert(0, PYD_DIR)
    from reaktoro4py import *

# Suppress convergence warnings
try:
    Warnings.disable(906)
except Exception:
    pass


# =============================================================================
# SATURATION PRESSURE CALCULATION
# =============================================================================


def calculate_saturation_pressure(temperature_celsius):
    """
    Calculate water saturation pressure using the Antoine equation.

    Parameters
    ----------
    temperature_celsius : float
        Temperature in Celsius (valid range: 0-374Â°C)

    Returns
    -------
    float
        Saturation pressure in bar

    Notes
    -----
    Uses Antoine equation coefficients valid for 0-374Â°C range.
    Returns NaN for temperatures outside valid range.
    """
    if temperature_celsius < 0 or temperature_celsius > 374:
        return np.nan

    T_kelvin = temperature_celsius + 273.15
    A, B, C = 5.40221, 1838.675, -31.737
    log10_P = A - B / (T_kelvin + C)
    return 10**log10_P


def calculate_saturation_pressure_kbar(temperature_celsius):
    """
    Calculate water saturation pressure in kilobars.

    Parameters
    ----------
    temperature_celsius : float
        Temperature in Celsius

    Returns
    -------
    float
        Saturation pressure in kbar
    """
    P_bar = calculate_saturation_pressure(temperature_celsius)
    return P_bar / 1000.0 if not np.isnan(P_bar) else np.nan


# =============================================================================
# CHEMICAL SYSTEM SETUP
# =============================================================================


def build_chemical_system(
    aqueous_db,
    mineral_db,
    aqueous_species_string,
    mineral_species_list,
    activity_model="DEW",
    dew_config=None,
):
    """
    Build a generic ChemicalSystem for mineral solubility calculations.

    This function is mineral-agnostic and can be used for any mineral + aqueous
    species combination (quartz, calcite, etc.)

    Parameters
    ----------
    aqueous_db : Database
        Database for aqueous species (e.g., DEW2024())
    mineral_db : Database
        Database for minerals (e.g., SUPCRTBL())
    aqueous_species_string : str
        Space-separated list of aqueous species names
        Example: "WATER,AQ H+ OH- SiO2_aq HSiO3- Na+ Cl-"
    mineral_species_list : list of str
        List of mineral species names to include
        Example: ["Quartz"] or ["Calcite"]
    activity_model : str
        Activity model for aqueous phase: "DEW" or "HKF" (default: "DEW")
    dew_config : dict, optional
        Configuration for DEW water model

    Returns
    -------
    ChemicalSystem
        Ready-to-use chemical system with aqueous and mineral phases

    Examples
    --------
    >>> dew_db = DEW2024()
    >>> supcrt_db = SUPCRTBL()
    >>> system = build_chemical_system(
    ...     dew_db, supcrt_db,
    ...     "WATER,AQ H+ OH- SiO2_aq HSiO3-",
    ...     ["Quartz"],
    ...     activity_model="DEW"
    ... )
    """
    if dew_config is None:
        dew_config = {}

    # Start with aqueous database
    combined_db = Database(aqueous_db.species())

    # Add mineral species from mineral database
    for mineral_name in mineral_species_list:
        try:
            mineral_species = mineral_db.species(mineral_name)
            combined_db.addSpecies(mineral_species)
        except Exception as e:
            print(f"Warning: Could not add {mineral_name}: {e}")

    # Create aqueous phase
    aqueous = AqueousPhase(aqueous_species_string)

    # Configure activity model
    try:
        if activity_model.upper() == "DEW":
            aqueous.setActivityModel(ActivityModelDEW())
            print(f"âœ“ Aqueous phase configured with DEW activity model")
        else:
            aqueous.setActivityModel(ActivityModelHKF())
            print(f"âœ“ Aqueous phase configured with HKF activity model")
    except Exception as e:
        print(f"Warning: Could not configure {activity_model} model: {e}")
        try:
            aqueous.setActivityModel(ActivityModelHKF())
        except Exception:
            pass

    # Create mineral phases
    mineral_phases = [MineralPhase(name) for name in mineral_species_list]

    # Create system with all phases
    return ChemicalSystem(combined_db, aqueous, *mineral_phases)


# Convenience function for backward compatibility
def build_chemical_system_dew(dew_database, supcrtbl_database, dew_config=None):
    """
    Convenience wrapper for quartz solubility with DEW databases.

    This is a simplified interface that maintains backward compatibility
    while using the generic build_chemical_system function.

    Parameters
    ----------
    dew_database : Database
        DEW2024 database
    supcrtbl_database : Database
        SUPCRTBL database
    dew_config : dict, optional
        Configuration for DEW water model

    Returns
    -------
    ChemicalSystem
        Chemical system configured for quartz solubility
    """
    return build_chemical_system(
        aqueous_db=dew_database,
        mineral_db=supcrtbl_database,
        aqueous_species_string="WATER,AQ H+ OH- SiO2_aq H2_aq O2_aq HO2- HSiO3- Si2O4_aq Si3O6_aq",
        mineral_species_list=["Quartz"],
        activity_model="DEW",
        dew_config=dew_config,
    )


# =============================================================================
# EQUILIBRIUM SOLVER
# =============================================================================


def solve_quartz_equilibrium(system, temperature_celsius, pressure_bar, verbose=False):
    """
    Solve equilibrium for quartz + water system at given T, P.

    Parameters
    ----------
    system : ChemicalSystem
        Chemical system (from build_chemical_system_dew)
    temperature_celsius : float
        Temperature in Celsius
    pressure_bar : float
        Pressure in bar
    verbose : bool
        Print diagnostic information if True

    Returns
    -------
    tuple
        (success: bool, state: ChemicalState or None)
        If successful, returns (True, state) with converged state.
        If failed, returns (False, None).

    Examples
    --------
    >>> success, state = solve_quartz_equilibrium(system, 250, 500)
    >>> if success:
    ...     print("Converged successfully")
    """
    try:
        solver = EquilibriumSolver(system)
        state = ChemicalState(system)

        # Initialize with 1 kg water, quartz, and trace silica species
        state.set("WATER,AQ", 1.0, "kg")
        state.set("H+", 1e-8, "mol")
        state.set("OH-", 1e-8, "mol")
        state.set("SiO2_aq", 1e-6, "mol")
        state.set("Quartz", 10.0, "mol")

        # Set pressure and temperature
        state.pressure(float(pressure_bar), "bar")
        state.temperature(float(temperature_celsius), "celsius")

        # Solve
        result = solver.solve(state)

        if result and result.succeeded():
            if verbose:
                print(
                    f"  âœ“ Converged at T={temperature_celsius}Â°C, P={pressure_bar:.0f} bar"
                )
            return True, state
        else:
            if verbose:
                print(
                    f"  âœ— Failed to converge at T={temperature_celsius}Â°C, P={pressure_bar:.0f} bar"
                )
            return False, None

    except Exception as e:
        if verbose:
            print(
                f"  âœ— Exception at T={temperature_celsius}Â°C, P={pressure_bar:.0f} bar: {e}"
            )
        return False, None


def get_silica_molality(state):
    """
    Extract SiO2_aq molality from converged equilibrium state.

    Parameters
    ----------
    state : ChemicalState
        Converged chemical state

    Returns
    -------
    float
        Molality (mol/kg solvent) of dissolved silica, or NaN if unavailable
    """
    try:
        aqprops = AqueousProps(state)
        molality = float(aqprops.speciesMolality("SiO2_aq"))
        return molality
    except Exception:
        return np.nan


# =============================================================================
# SOLUBILITY CURVE CALCULATION
# =============================================================================


def calculate_solubility_curve(system, temperature_range_celsius, pressure_bar):
    """
    Calculate quartz solubility curve across a temperature range at fixed pressure.

    Parameters
    ----------
    system : ChemicalSystem
        Chemical system
    temperature_range_celsius : array-like
        Temperature values in Celsius
    pressure_bar : float
        Pressure in bar (fixed)

    Returns
    -------
    dict
        Contains 'T_C' (temperature array) and 'molality' (solubility array).
        NaN values indicate convergence failures.

    Examples
    --------
    >>> T_range = np.linspace(100, 300, 50)
    >>> curve = calculate_solubility_curve(system, T_range, 500)
    >>> print(f"Converged {np.sum(~np.isnan(curve['molality']))} points")
    """
    molalities = []

    for T_C in temperature_range_celsius:
        success, state = solve_quartz_equilibrium(system, T_C, pressure_bar)

        if success:
            molality = get_silica_molality(state)
        else:
            molality = np.nan

        molalities.append(molality)

    return {
        "T_C": np.array(temperature_range_celsius),
        "molality": np.array(molalities),
    }


def calculate_saturation_curve(system, temperature_range_celsius, n_points=100):
    """
    Calculate quartz solubility along water saturation curve (Psat).

    Parameters
    ----------
    system : ChemicalSystem
        Chemical system
    temperature_range_celsius : array-like
        Temperature values in Celsius (0-374Â°C valid)
    n_points : int
        Number of points for saturation curve

    Returns
    -------
    dict
        Contains 'T_C', 'P_kbar', and 'molality' arrays

    Examples
    --------
    >>> T_range = np.linspace(100, 350, 100)
    >>> psat_curve = calculate_saturation_curve(system, T_range)
    """
    molalities = []
    valid_pressures = []
    valid_temps = []

    for T_C in temperature_range_celsius:
        P_kbar = calculate_saturation_pressure_kbar(T_C)

        if np.isnan(P_kbar):
            molalities.append(np.nan)
            continue

        P_bar = P_kbar * 1000.0
        success, state = solve_quartz_equilibrium(system, T_C, P_bar)

        if success:
            molality = get_silica_molality(state)
            valid_temps.append(T_C)
            valid_pressures.append(P_kbar)
            molalities.append(molality)
        else:
            molalities.append(np.nan)

    return {
        "T_C": np.array(temperature_range_celsius),
        "P_kbar": np.array(
            [calculate_saturation_pressure_kbar(T) for T in temperature_range_celsius]
        ),
        "molality": np.array(molalities),
    }


# =============================================================================
# DATA LOADING
# =============================================================================


def load_experimental_data(csv_file):
    """
    Load and organize experimental quartz solubility data from CSV.

    Parameters
    ----------
    csv_file : str
        Path to CSV file with columns: T_C, P_kbar, molality_m, reference, experiment_type

    Returns
    -------
    pd.DataFrame
        Organized with separate flags for Psat vs fixed-pressure experiments
    """
    df = pd.read_csv(csv_file)
    df = df[["T_C", "P_kbar", "molality_m", "reference", "experiment_type"]].copy()
    df["P_bar"] = df["P_kbar"] * 1000.0
    df["experiment_id"] = df["reference"] + " (" + df["experiment_type"] + ")"

    # Identify saturation pressure vs controlled pressure experiments
    kennedy_controlled = {0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.75}
    is_kennedy = df["reference"].str.contains("Kennedy", case=False, na=False)
    is_controlled = (
        is_kennedy & df["P_kbar"].notna() & df["P_kbar"].isin(kennedy_controlled)
    )
    df["is_psat"] = ~is_controlled
    df = df.sort_values(["is_psat", "P_kbar", "T_C"]).reset_index(drop=True)

    return df


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================


def count_valid_points(molality_array):
    """Count non-NaN points in molality array."""
    return np.sum(~np.isnan(molality_array))


def get_valid_temperature_range(molality_array, temperature_array):
    """
    Get min/max temperatures where solubility was successfully calculated.

    Returns
    -------
    tuple
        (T_min, T_max) or (None, None) if no valid points
    """
    valid_idx = np.where(~np.isnan(molality_array))[0]
    if len(valid_idx) == 0:
        return None, None
    return temperature_array[valid_idx[0]], temperature_array[valid_idx[-1]]


# =============================================================================
# HIGH-LEVEL CALCULATION FUNCTIONS (FOR NOTEBOOK USE)
# =============================================================================


def calculate_mineral_solubility_curves(
    system,
    specs,
    conditions,
    temperature_range,
    pressure_list_kbar,
    mineral_name,
    solubility_species,
    include_psat=True,
    n0_mineral=10.0,
    verbose=True,
):
    """
    High-level function to calculate solubility curves for multiple pressures.

    This is the main function called from notebooks - it handles all the details.

    Parameters
    ----------
    system : ChemicalSystem
        The chemical system
    specs : EquilibriumSpecs
        Equilibrium specifications
    conditions : EquilibriumConditions
        Equilibrium conditions
    temperature_range : array-like
        Temperature values in Celsius
    pressure_list_kbar : array-like
        Pressure values in kbar to calculate
    mineral_name : str
        Name of mineral species (e.g., "Quartz", "Calcite")
    solubility_species : str
        Aqueous species to measure (e.g., "SiO2_aq", "Ca+2")
    include_psat : bool
        Whether to include saturation pressure curve
    n0_mineral : float
        Initial amount of mineral (mol)
    verbose : bool
        Print progress messages

    Returns
    -------
    dict
        Dictionary with pressure as keys, each containing 'T_C' and 'molality' arrays
    """
    import matplotlib.pyplot as plt

    solver = EquilibriumSolver(specs)
    aprops = AqueousProps(system)

    solubility_curves = {}

    if verbose:
        print(f"Calculating {mineral_name} solubility curves...")
        print("=" * 60)

    # Calculate for each fixed pressure
    for P_kbar in pressure_list_kbar:
        P_bar = P_kbar * 1000.0
        if verbose:
            print(f"  P = {P_kbar:.3f} kbar ({P_bar:.0f} bar)...")

        molalities = []
        for T_C in temperature_range:
            # Create initial state
            state = ChemicalState(system)
            state.set("WATER,AQ", 1.0, "kg")
            state.set("H+", 1e-8, "mol")
            state.set("OH-", 1e-8, "mol")
            state.set(solubility_species, 1e-6, "mol")
            state.set(mineral_name, n0_mineral, "mol")

            # Set conditions
            conditions.temperature(T_C, "celsius")
            conditions.pressure(P_bar, "bar")

            # Solve
            try:
                result = solver.solve(state, conditions)
                if result.succeeded():
                    aprops.update(state)
                    molality = float(aprops.speciesMolality(solubility_species))
                else:
                    molality = np.nan
            except Exception:
                molality = np.nan

            molalities.append(molality)

        solubility_curves[P_kbar] = {
            "T_C": np.array(temperature_range),
            "molality": np.array(molalities),
        }

        if verbose:
            n_valid = np.sum(~np.isnan(molalities))
            print(f"    âœ“ Converged {n_valid}/{len(temperature_range)} points")

    # Calculate saturation curve if requested
    if include_psat:
        if verbose:
            print(f"  Calculating saturation curve...")

        T_psat_range = np.linspace(100, 374, len(temperature_range))
        psat_molalities = []

        for T_C in T_psat_range:
            P_kbar = calculate_saturation_pressure_kbar(T_C)

            if np.isnan(P_kbar):
                psat_molalities.append(np.nan)
                continue

            P_bar = P_kbar * 1000.0

            # Create initial state
            state = ChemicalState(system)
            state.set("WATER,AQ", 1.0, "kg")
            state.set("H+", 1e-8, "mol")
            state.set("OH-", 1e-8, "mol")
            state.set(solubility_species, 1e-6, "mol")
            state.set(mineral_name, n0_mineral, "mol")

            # Set conditions
            conditions.temperature(T_C, "celsius")
            conditions.pressure(P_bar, "bar")

            # Solve
            try:
                result = solver.solve(state, conditions)
                if result.succeeded():
                    aprops.update(state)
                    molality = float(aprops.speciesMolality(solubility_species))
                else:
                    molality = np.nan
            except Exception:
                molality = np.nan

            psat_molalities.append(molality)

        solubility_curves["Psat"] = {
            "T_C": T_psat_range,
            "P_kbar": np.array(
                [calculate_saturation_pressure_kbar(T) for T in T_psat_range]
            ),
            "molality": np.array(psat_molalities),
        }

        if verbose:
            n_valid = np.sum(~np.isnan(psat_molalities))
            print(f"    âœ“ Saturation: {n_valid}/{len(T_psat_range)} points")

    if verbose:
        print("=" * 60)
        print("âœ“ Calculations complete")

    return solubility_curves


def plot_mineral_solubility_with_experiments(
    solubility_curves,
    exp_data,
    mineral_name,
    mineral_formula,
    output_prefix=None,
    low_P_threshold=1.0,
    show_plots=True,
):
    """
    Create publication-quality plots comparing calculated curves with experiments.

    Generates two plots:
    - Low pressure (<threshold)
    - High pressure (>=threshold)

    Parameters
    ----------
    solubility_curves : dict
        Output from calculate_mineral_solubility_curves()
    exp_data : pd.DataFrame
        Experimental data with columns: T_C, P_kbar, molality_m, reference, is_psat
    mineral_name : str
        Name for labels (e.g., "Quartz")
    mineral_formula : str
        Chemical formula (e.g., "SiOâ‚‚")
    output_prefix : str, optional
        Prefix for output filenames. If None, uses mineral_name.lower()
    low_P_threshold : float
        Pressure (kbar) separating low/high plots (default: 1.0)
    show_plots : bool
        Whether to display plots interactively

    Returns
    -------
    tuple
        (low_P_filename, high_P_filename) - paths to saved plot files
    """
    import matplotlib.pyplot as plt

    if output_prefix is None:
        output_prefix = mineral_name.lower()

    # Separate experimental data
    if len(exp_data) > 0:
        non_psat_data = exp_data[~exp_data["is_psat"]]
        psat_data = exp_data[exp_data["is_psat"]]

        low_P_data = non_psat_data[non_psat_data["P_kbar"] < low_P_threshold]
        high_P_data = non_psat_data[non_psat_data["P_kbar"] >= low_P_threshold]

        low_P_pressures = sorted(low_P_data["P_kbar"].unique()) if len(low_P_data) > 0 else []
        high_P_pressures = sorted(high_P_data["P_kbar"].unique()) if len(high_P_data) > 0 else []
    else:
        low_P_pressures = []
        high_P_pressures = []
        low_P_data = pd.DataFrame()
        high_P_data = pd.DataFrame()
        psat_data = pd.DataFrame()

    # Author marker mapping
    author_markers = {
        "Kennedy_1950": "o",
        "Hemley_1980": "^",
        "Morey_Fournier_Rowe_1962": "s",
        "Walther_Orville_1983": "D",
        "Manning_1994": "v",
        "Newton_Manning_2000": "p",
    }

    # =========================================================================
    # LOW PRESSURE PLOT
    # =========================================================================
    low_P_file = f"{output_prefix}_solubility_low_pressure.png"

    if len(low_P_pressures) > 0 or 'Psat' in solubility_curves:
        fig, ax = plt.subplots(figsize=(14, 8))

        n_low = max(len(low_P_pressures), 1)
        colors_low = plt.cm.viridis(np.linspace(0, 0.9, n_low))
        P_to_color_low = {P: colors_low[i] for i, P in enumerate(low_P_pressures)}

        # Plot experimental data
        for P_kbar in low_P_pressures:
            subset = low_P_data[low_P_data["P_kbar"] == P_kbar]
            for author in subset["reference"].unique():
                author_subset = subset[subset["reference"] == author]
                marker = author_markers.get(author, "o")
                ax.scatter(
                    author_subset["T_C"], author_subset["molality_m"],
                    c=[P_to_color_low[P_kbar]], marker=marker, s=70,
                    alpha=0.7, edgecolors="black", linewidths=0.4,
                    label=f"Exp P={P_kbar:.2f} kbar ({author})", zorder=10,
                )

        # Plot Psat experiments
        if len(psat_data) > 0:
            low_psat_data = psat_data[(psat_data["P_kbar"] < low_P_threshold) | (psat_data["P_kbar"].isna())]
            for author in low_psat_data["reference"].unique():
                author_psat = low_psat_data[low_psat_data["reference"] == author]
                marker = author_markers.get(author, "s")
                ax.scatter(
                    author_psat["T_C"], author_psat["molality_m"],
                    c="purple", marker=marker, s=80, alpha=0.8,
                    edgecolors="darkviolet", linewidths=0.5,
                    label=f"Exp P=Psat ({author})", zorder=11,
                )

        # Plot calculated curves
        for P_kbar in low_P_pressures:
            if P_kbar in solubility_curves:
                curve = solubility_curves[P_kbar]
                valid = ~np.isnan(curve["molality"])
                ax.plot(
                    curve["T_C"][valid], curve["molality"][valid],
                    color=P_to_color_low[P_kbar], linewidth=2.0,
                    label=f"Calc P={P_kbar:.2f} kbar", zorder=5,
                )

        # Plot Psat curve
        if "Psat" in solubility_curves:
            curve = solubility_curves["Psat"]
            valid = ~np.isnan(curve["molality"])
            ax.plot(
                curve["T_C"][valid], curve["molality"][valid],
                color="purple", linewidth=3.0, linestyle="-",
                label="Calc P=Psat", zorder=6, alpha=0.9,
            )

        ax.set_yscale("log")
        ax.set_ylim(1e-4, 1e-1)
        ax.set_xlabel("Temperature (Â°C)", fontsize=14, fontweight="bold")
        ax.set_ylabel(f"{mineral_name} Solubility (mol/kg-Hâ‚‚O)", fontsize=14, fontweight="bold")
        ax.set_title(f"{mineral_name} Solubility: Low Pressure (<{low_P_threshold} kbar)",
                    fontsize=16, fontweight="bold", pad=20)
        ax.grid(True, which="both", alpha=0.3, linestyle="--")

        info_text = f"DEW2024 (aqueous) + SUPCRTBL ({mineral_name})"
        ax.text(0.02, 0.98, info_text, transform=ax.transAxes, fontsize=8,
               verticalalignment="top", bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))

        ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=9, framealpha=0.9, ncol=1)
        fig.tight_layout()
        plt.savefig(low_P_file, dpi=300, bbox_inches="tight")
        print(f"âœ“ Low-pressure plot saved: {low_P_file}")
        if show_plots:
            plt.show()
        else:
            plt.close(fig)

    # =========================================================================
    # HIGH PRESSURE PLOT
    # =========================================================================
    high_P_file = f"{output_prefix}_solubility_high_pressure.png"

    if len(high_P_pressures) > 0:
        # Combine all high-P experiments
        high_P_all_data = pd.concat([
            psat_data[(psat_data["P_kbar"] >= low_P_threshold) & (psat_data["P_kbar"].notna())] if len(psat_data) > 0 else pd.DataFrame(),
            non_psat_data[non_psat_data["P_kbar"] >= low_P_threshold] if len(non_psat_data) > 0 else pd.DataFrame(),
        ], ignore_index=True)

        if len(high_P_all_data) > 0:
            high_P_all_pressures = sorted(high_P_all_data["P_kbar"].unique())

            fig, ax = plt.subplots(figsize=(14, 8))

            n_high = max(len(high_P_all_pressures), 1)
            colors_high = plt.cm.plasma(np.linspace(0, 0.9, n_high))
            P_to_color_high = {P: colors_high[i] for i, P in enumerate(high_P_all_pressures)}

            # Plot experimental data
            for P_kbar in high_P_all_pressures:
                subset = high_P_all_data[high_P_all_data["P_kbar"] == P_kbar]
                for author in subset["reference"].unique():
                    author_subset = subset[subset["reference"] == author]
                    marker = author_markers.get(author, "o")
                    ax.scatter(
                        author_subset["T_C"], author_subset["molality_m"],
                        c=[P_to_color_high[P_kbar]], marker=marker, s=70,
                        alpha=0.7, edgecolors="black", linewidths=0.4,
                        label=f"Exp P={P_kbar:.2f} kbar ({author})", zorder=10,
                    )

            # Plot calculated curves
            for P_kbar in high_P_all_pressures:
                if P_kbar in solubility_curves:
                    curve = solubility_curves[P_kbar]
                    valid = ~np.isnan(curve["molality"])
                    ax.plot(
                        curve["T_C"][valid], curve["molality"][valid],
                        color=P_to_color_high[P_kbar], linewidth=2.0,
                        label=f"Calc P={P_kbar:.2f} kbar", zorder=5,
                    )

            ax.set_yscale("log")
            ax.set_xlabel("Temperature (Â°C)", fontsize=14, fontweight="bold")
            ax.set_ylabel(f"{mineral_name} Solubility (mol/kg-Hâ‚‚O)", fontsize=14, fontweight="bold")
            ax.set_title(f"{mineral_name} Solubility: High Pressure (â‰¥{low_P_threshold} kbar)",
                        fontsize=16, fontweight="bold", pad=20)
            ax.grid(True, which="both", alpha=0.3, linestyle="--")

            info_text = f"DEW2024 (aqueous) + SUPCRTBL ({mineral_name})"
            ax.text(0.02, 0.98, info_text, transform=ax.transAxes, fontsize=8,
                   verticalalignment="top", bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))

            ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=9, framealpha=0.9, ncol=1)
            fig.tight_layout()
            plt.savefig(high_P_file, dpi=300, bbox_inches="tight")
            print(f"âœ“ High-pressure plot saved: {high_P_file}")
            if show_plots:
                plt.show()
            else:
                plt.close(fig)

    return low_P_file, high_P_file
        Experimental data with columns: T_C, P_kbar, molality_m, reference, is_psat
    mineral_name : str
        Name for labels (e.g., "Quartz")
    mineral_formula : str
        Chemical formula (e.g., "SiOâ‚‚")
    output_prefix : str, optional
        Prefix for output filenames. If None, uses mineral_name.lower()
    low_P_threshold : float
        Pressure (kbar) separating low/high plots (default: 1.0)
    show_plots : bool
        Whether to display plots interactively

    Returns
    -------
    tuple
        (T_min, T_max) or (None, None) if no valid points
    """
    valid_idx = np.where(~np.isnan(molality_array))[0]
    if len(valid_idx) == 0:
        return None, None
    return temperature_array[valid_idx[0]], temperature_array[valid_idx[-1]]

