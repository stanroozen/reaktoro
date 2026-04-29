"""
Quartz Solubility Analysis using Reaktoro with DEW2024
Compares calculated solubilities with experimental data
"""

import autodiff
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import sys

# Try to import Reaktoro; fall back to local extension module if not installed
try:
    from reaktoro import *  # noqa: F401,F403
except ModuleNotFoundError:
    # Add local build path for reaktoro4py if available
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    ROOT_DIR = os.path.dirname(SCRIPT_DIR)
    PYD_DIR = os.path.join(ROOT_DIR, "build", "Reaktoro", "Release")
    if os.path.isdir(PYD_DIR) and PYD_DIR not in sys.path:
        sys.path.insert(0, PYD_DIR)
    try:
        from reaktoro4py import *  # noqa: F401,F403

        print("Using local reaktoro4py extension from build.")
    except ModuleNotFoundError as e:
        raise ModuleNotFoundError(
            "Could not import Reaktoro. Install the 'reaktoro' package or ensure reaktoro4py is on PYTHONPATH."
        ) from e

# Silence repeated non-convergence warnings while we handle failures manually
try:
    Warnings.disable(906)
except Exception:
    pass

# =============================================================================
# Configuration
# =============================================================================

# Get script directory for relative paths
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CSV_FILE = os.path.join(SCRIPT_DIR, "quartz_DEW_testset.csv")
OUTPUT_PLOT = os.path.join(SCRIPT_DIR, "quartz_solubility_comparison.png")

# Temperature range for solubility curves (Â°C)
T_MIN, T_MAX = 150, 550
N_POINTS = 100

# =============================================================================
# DEW Water Model Configuration
# =============================================================================
# Default: Duan & Zang 2005 EOS, Power Function dielectric, Volume integration for Gibbs
# Can be customized by modifying these values or by passing to build_system()

DEW_CONFIG = {
    "eos_model": "ZhangDuan2005",  # Options: WagnerPruss, HGK, ZhangDuan2005, ZhangDuan2009
    "dielectric_model": "PowerFunction",  # Options: PowerFunction, JohnsonNorton1991
    "gibbs_model": "DewIntegral",  # Options: DewIntegral, DelaneyHelgeson1978
    "born_model": "Shock92Dew",  # Options: Shock92Dew, Shock92 (for neutral species)
}

# =============================================================================
# Helper Functions
# =============================================================================


def build_system(dew_db, supcrt_db, water_config=None):
    """
    Build and return a ChemicalSystem combining DEW aqueous species (including H2O_aq) with Quartz from SUPCRTBL.

    Parameters:
    -----------
    dew_db : DEWDatabase
        DEW database containing aqueous species
    supcrt_db : SupcrtDatabase
        SUPCRT database containing minerals
    water_config : dict, optional
        Custom water model configuration. If None, uses DEW_CONFIG defaults.
        Supported keys:
        - eos_model: WagnerPruss, HGK, ZhangDuan2005 (default), ZhangDuan2009
        - dielectric_model: PowerFunction (default), JohnsonNorton1991
        - gibbs_model: DewIntegral (default), DelaneyHelgeson1978
        - born_model: Shock92Dew (default), Shock92
    """
    if water_config is None:
        water_config = DEW_CONFIG

    # Cherry-pick Quartz from SUPCRT
    quartz_species = supcrt_db.species("Quartz")

    # Create combined database: all DEW species (includes H2O_aq) + Quartz
    combined_db = Database(dew_db.species())
    combined_db.addSpecies(quartz_species)

    # Define the aqueous phase with WATER,AQ from DEW database
    # (WATER,AQ is the species name, derived from H2O_aq key in YAML)
    aqueous = AqueousPhase(
        "WATER,AQ H+ OH- SiO2_aq H2_aq O2_aq HO2- HSiO3- Si2O4_aq Si3O6_aq"
    )

    # Apply DEW activity model with configurable water thermodynamics
    try:
        # Create DEW model with custom water options
        params = StandardThermoModelParamsDEW()

        # Map string config to enum values
        eos_map = {
            "WagnerPruss": WaterEosModel.WagnerPruss,
            "HGK": WaterEosModel.HGK,
            "ZhangDuan2005": WaterEosModel.ZhangDuan2005,
            "ZhangDuan2009": WaterEosModel.ZhangDuan2009,
        }

        dielectric_map = {
            "PowerFunction": WaterDielectricModel.PowerFunction,
            "JohnsonNorton1991": WaterDielectricModel.JohnsonNorton1991,
        }

        gibbs_map = {
            "DewIntegral": WaterGibbsModel.DewIntegral,
            "DelaneyHelgeson1978": WaterGibbsModel.DelaneyHelgeson1978,
        }

        born_map = {
            "Shock92Dew": WaterBornModel.Shock92Dew,
        }

        # Set water model options from config
        params.waterOptions.eosModel = eos_map.get(
            water_config.get("eos_model", "ZhangDuan2005"), WaterEosModel.ZhangDuan2005
        )
        params.waterOptions.dielectricModel = dielectric_map.get(
            water_config.get("dielectric_model", "PowerFunction"),
            WaterDielectricModel.PowerFunction,
        )
        params.waterOptions.gibbsModel = gibbs_map.get(
            water_config.get("gibbs_model", "DewIntegral"), WaterGibbsModel.DewIntegral
        )
        params.waterOptions.bornModel = born_map.get(
            water_config.get("born_model", "Shock92Dew"), WaterBornModel.Shock92Dew
        )

        dew_model = StandardThermoModelDEW(params)
        aqueous.setActivityModel(dew_model)

        # Log water model configuration
        eos_name = water_config.get("eos_model", "ZhangDuan2005")
        diel_name = water_config.get("dielectric_model", "PowerFunction")
        gibbs_name = water_config.get("gibbs_model", "DewIntegral")
        print(
            f"âœ“ DEW configured: EOS={eos_name}, Dielectric={diel_name}, Gibbs={gibbs_name}"
        )

    except Exception as e:
        print(f"Warning: Could not configure DEW with custom water options: {e}")
        print("  Falling back to default ActivityModelDEW()")
        try:
            aqueous.setActivityModel(ActivityModelDEW())
        except NameError:
            print("  ActivityModelDEW not available, using HKF fallback")
            aqueous.setActivityModel(ActivityModelHKF())

    # Define mineral phase (Quartz from SUPCRT database)
    mineral = MineralPhase("Quartz")

    # Create system from combined database and phases
    system = ChemicalSystem(combined_db, aqueous, mineral)
    return system


def calculate_quartz_solubility(T_C, P_bar, dew_db, supcrt_db, water_config=None):
    """
    Calculate quartz solubility at given T and P using Reaktoro with DEW water model.
    Combines aqueous species from DEW with minerals from SUPCRTBL.

    Parameters:
    -----------
    T_C : float
        Temperature in Â°C
    P_bar : float
        Pressure in bar
    dew_db : DEWDatabase
        DEW database for aqueous species
    supcrt_db : SupcrtDatabase
        SUPCRTBL database for minerals
    water_config : dict, optional
        Custom water model configuration (see build_system for details)

    Returns:
    --------
    float : Quartz molality (mol/kg-H2O), or NaN if calculation fails
    """
    try:
        # Ensure Reaktoro-friendly scalar types
        T_C = float(T_C)
        P_bar = float(P_bar)
        T_real = autodiff.real(T_C)
        P_real = autodiff.real(P_bar)

        # Build system with configurable water model
        system = build_system(dew_db, supcrt_db, water_config=water_config)
        # Build system with configurable water model
        system = build_system(dew_db, supcrt_db, water_config=water_config)

        def build_state():
            """Create a fresh ChemicalState with consistent seeds."""
            s = ChemicalState(system)
            s.temperature(T_real, "celsius")
            s.pressure(P_real, "bar")
            s.set(
                "WATER,AQ", 1.0, "kg"
            )  # Solvent; DEW computes properties via Duan EOS
            s.set("H+", 1e-8, "mol")
            s.set("OH-", 1e-8, "mol")
            s.set("SiO2_aq", 1e-6, "mol")
            s.set("Quartz", 10.0, "mol")
            return s

        # Equilibrate using EquilibriumSolver with a fallback cold start on failure
        solver = EquilibriumSolver(system)
        state = build_state()
        result = solver.solve(state)

        if not result.succeeded():
            # Retry with a fresh solver/state as a cold start (no special options)
            solver = EquilibriumSolver(system)
            state = build_state()
            result = solver.solve(state)

        if not result.succeeded():
            return np.nan

        # Get dissolved silica molality from aqueous properties
        # DEW database uses SiO2_aq for dissolved silica
        try:
            aqprops = AqueousProps(state)
            molality = float(aqprops.speciesMolality("SiO2_aq"))
            return molality
        except Exception as e:
            print(
                f"Warning: Could not extract SiO2(aq) molality at T={T_C}Â°C, P={P_bar} bar: {e}"
            )
            return np.nan

    except Exception as e:
        print(f"Warning: Failed at T={T_C}Â°C, P={P_bar} bar: {e}")
        return np.nan


def load_experimental_data(csv_file):
    """Load and organize experimental data from CSV."""
    df = pd.read_csv(csv_file)

    # Extract relevant columns
    df = df[["T_C", "P_kbar", "molality_m", "reference", "experiment_type"]].copy()

    # Convert P from kbar to bar
    df["P_bar"] = df["P_kbar"] * 1000.0

    # Handle NaN pressures (LV curve data) - use saturation pressure approximation
    # For simplicity, we'll skip these or use a placeholder
    df = df.dropna(subset=["P_kbar"])

    # Create experiment identifier (reference + experiment type)
    df["experiment_id"] = df["reference"] + " (" + df["experiment_type"] + ")"

    # Add kbar category: nearest integer kbar (Â±0.5 tolerance by definition)
    df["kbar_cat"] = df["P_kbar"].round().astype(int)
    df = df.sort_values(["kbar_cat", "T_C"]).reset_index(drop=True)

    return df


# =============================================================================
# Main Script
# =============================================================================


def main():
    print("=" * 80)
    print("Quartz Solubility Analysis - Reaktoro DEW2024")
    print("=" * 80)

    # Load experimental data
    print("\n[1] Loading experimental data...")
    if not os.path.exists(CSV_FILE):
        print(f"    WARNING: Experimental data file not found: {CSV_FILE}")
        print("    Proceeding with calculated curves only (no experimental comparison)")
        exp_data = pd.DataFrame()
    else:
        exp_data = load_experimental_data(CSV_FILE)
        print(f"    Loaded {len(exp_data)} experimental data points")
        if len(exp_data) > 0:
            print(
                f"    Temperature range: {exp_data['T_C'].min():.0f} - {exp_data['T_C'].max():.0f} Â°C"
            )
            print(
                f"    Pressure range: {exp_data['P_kbar'].min():.3f} - {exp_data['P_kbar'].max():.3f} kbar"
            )

    # Get unique experiments and pressures
    if len(exp_data) > 0:
        experiments = exp_data["experiment_id"].unique()
        pressures_kbar = sorted(exp_data["P_kbar"].unique())
        kbar_categories = sorted(exp_data["kbar_cat"].unique())
        print(f"    Experiments: {len(experiments)}")
        print(f"    Kbar categories present: {kbar_categories}")
    else:
        experiments = []
        pressures_kbar = [0.5, 1.0, 2.0, 5.0]  # Default pressures
        kbar_categories = [1, 2, 5]
        print(f"    Using default pressure conditions: {pressures_kbar} kbar")

    # Initialize databases
    print("\n[2] Initializing Reaktoro databases...")
    try:
        dew_db = DEWDatabase("dew2019-aqueous")
        print("    Successfully loaded DEW2019 aqueous database")
    except Exception as e:
        print(f"    ERROR: Failed to load DEW database: {e}")
        raise

    try:
        supcrt_db = SupcrtDatabase("supcrtbl")
        print("    Successfully loaded SUPCRTBL mineral database")
    except Exception as e:
        print(f"    ERROR: Failed to load SUPCRTBL database: {e}")
        raise

    print(f"    Pressure conditions: {len(pressures_kbar)}")

    # Build system once (independent of pressure)
    system = build_system(dew_db, supcrt_db)

    # Calculate solubility curves for each pressure
    print("\n[3] Calculating quartz solubility curves...")
    T_range = np.linspace(T_MIN, T_MAX, N_POINTS)

    solubility_curves = {}
    for P_kbar in pressures_kbar:
        P_bar = P_kbar * 1000.0
        print(f"    P = {P_kbar:.2f} kbar ({P_bar:.0f} bar)...")

        # Initialize solver and state; reuse state across temperatures for robust convergence
        solver = EquilibriumSolver(system)
        state = ChemicalState(system)
        state.set("WATER,AQ", 1.0, "kg")
        state.set("H+", 1e-8, "mol")
        state.set("OH-", 1e-8, "mol")
        state.set("SiO2_aq", 1e-6, "mol")
        state.set("Quartz", 10.0, "mol")
        state.pressure(float(P_bar), "bar")

        molalities = []
        for T_C in T_range:
            state.temperature(float(T_C), "celsius")
            result = solver.solve(state)

            if result.succeeded():
                try:
                    aqprops = AqueousProps(state)
                    molality = float(aqprops.speciesMolality("SiO2_aq"))
                except Exception as e:
                    print(
                        f"Warning: Could not extract SiO2_aq molality at T={T_C:.1f}Â°C, P={P_bar:.1f} bar: {e}"
                    )
                    molality = np.nan
            else:
                molality = np.nan

            molalities.append(molality)

        solubility_curves[P_kbar] = {"T_C": T_range, "molality": np.array(molalities)}

        valid_points = np.sum(~np.isnan(molalities))
        if valid_points > 0:
            valid_idx = np.where(~np.isnan(molalities))[0]
            first_idx = valid_idx[0]
            last_idx = valid_idx[-1]
            first_T = T_range[first_idx]
            last_T = T_range[last_idx]
            first_m = molalities[first_idx]
            last_m = molalities[last_idx]
            print(
                f"       Calculated {valid_points}/{N_POINTS} points successfully "
                f"(first valid: T={first_T:.1f} Â°C, m={first_m:.3e}; "
                f"last valid: T={last_T:.1f} Â°C, m={last_m:.3e})"
            )
        else:
            print(f"       Calculated {valid_points}/{N_POINTS} points successfully")

    # Plotting
    print("\n[4] Creating plots...")

    fig, ax = plt.subplots(figsize=(12, 8))

    # Color schemes
    # Assign consistent colors per kbar category
    colors_cat = plt.cm.tab20(np.linspace(0, 0.9, max(len(kbar_categories), 1)))
    cat_to_color = {
        cat: colors_cat[i % len(colors_cat)] for i, cat in enumerate(kbar_categories)
    }

    # Marker styles
    markers = ["o", "s", "^", "v", "D", "<", ">", "p", "*", "h"]

    # Plot experimental data grouped by kbar category (consistent colors per category)
    if len(exp_data) > 0:
        for cat in kbar_categories:
            subset = exp_data[exp_data["kbar_cat"] == cat]
            if len(subset) == 0:
                continue
            ax.scatter(
                subset["T_C"],
                subset["molality_m"],
                c=[cat_to_color[cat]],
                marker="o",
                s=70,
                alpha=0.7,
                edgecolors="black",
                linewidths=0.4,
                label=f"Exp Pâ‰ˆ{cat} kbar",
                zorder=10,
            )

    # Plot calculated solubility curves per integer kbar category (matching experimental grouping)
    for cat in kbar_categories:
        P_kbar = float(cat)
        if P_kbar not in solubility_curves:
            # If exact integer kbar curve not computed, skip
            continue
        curve = solubility_curves[P_kbar]
        valid = ~np.isnan(curve["molality"])
        ax.plot(
            curve["T_C"][valid],
            curve["molality"][valid],
            color=cat_to_color[cat],
            linewidth=2.5,
            linestyle="-",
            label=f"Calc P={cat} kbar",
            zorder=5,
        )

    # Set logarithmic y-axis
    ax.set_yscale("log")

    # Labels and formatting
    ax.set_xlabel("Temperature (Â°C)", fontsize=14, fontweight="bold")
    ax.set_ylabel("Quartz Solubility (mol/kg-Hâ‚‚O)", fontsize=14, fontweight="bold")
    ax.set_title(
        "Quartz Solubility: DEW2024 Calculations vs Experimental Data",
        fontsize=16,
        fontweight="bold",
        pad=20,
    )

    # Grid
    ax.grid(True, which="both", alpha=0.3, linestyle="--")

    # Legend outside frame (right side)
    ax.legend(
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        fontsize=9,
        framealpha=0.9,
        ncol=1,
    )

    # Adjust layout to prevent legend cutoff
    plt.tight_layout()

    # Save figure
    plt.savefig(OUTPUT_PLOT, dpi=300, bbox_inches="tight")
    print(f"    Plot saved to: {OUTPUT_PLOT}")

    # Compute predicted values at experimental points and plot residuals
    if len(exp_data) > 0:

        def interp_pred(P_kbar, T_C):
            curve = solubility_curves.get(float(P_kbar))
            if curve is None:
                return np.nan
            T = curve["T_C"]
            m = curve["molality"]
            valid = ~np.isnan(m)
            if np.sum(valid) < 2:
                return np.nan
            return float(np.interp(T_C, T[valid], m[valid]))

        exp_data["predicted_m"] = exp_data.apply(
            lambda r: interp_pred(r["P_kbar"], r["T_C"]), axis=1
        )
        exp_data["residual_m"] = exp_data["molality_m"] - exp_data["predicted_m"]

        fig2, ax2 = plt.subplots(figsize=(12, 6))
        for cat in kbar_categories:
            subset = exp_data[exp_data["kbar_cat"] == cat]
            if len(subset) == 0:
                continue
            ax2.scatter(
                subset["T_C"],
                subset["residual_m"],
                c=[cat_to_color[cat]],
                marker="o",
                s=60,
                alpha=0.8,
                edgecolors="black",
                linewidths=0.4,
                label=f"Pâ‰ˆ{cat} kbar",
            )
        ax2.axhline(0.0, color="gray", linestyle="--", linewidth=1.0)
        ax2.set_xlabel("Temperature (Â°C)", fontsize=12, fontweight="bold")
        ax2.set_ylabel(
            "Residual (measured âˆ’ predicted) mol/kg-Hâ‚‚O", fontsize=12, fontweight="bold"
        )
        ax2.set_title(
            "Quartz Solubility Residuals by Pressure Category",
            fontsize=14,
            fontweight="bold",
        )
        ax2.grid(True, alpha=0.3, linestyle="--")
        ax2.legend(
            loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=9, framealpha=0.9
        )
        plt.tight_layout()
        residuals_path = os.path.join(SCRIPT_DIR, "quartz_solubility_residuals.png")
        plt.savefig(residuals_path, dpi=300, bbox_inches="tight")
        print(f"    Residuals plot saved to: {residuals_path}")

    # Show plots
    plt.show()

    print("\n" + "=" * 80)
    print("Analysis complete!")
    print("=" * 80)


if __name__ == "__main__":
    main()

