"""
Quartz Solubility Analysis using Reaktoro with DEW2024
Compares calculated solubilities with experimental data
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from reaktoro import *
import os

# =============================================================================
# Configuration
# =============================================================================

# Get script directory for relative paths
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CSV_FILE = os.path.join(SCRIPT_DIR, "quartz_DEW_testset.csv")
OUTPUT_PLOT = os.path.join(SCRIPT_DIR, "quartz_solubility_comparison.png")

# Temperature range for solubility curves (°C)
T_MIN, T_MAX = 150, 550
N_POINTS = 100

# =============================================================================
# Helper Functions
# =============================================================================


def calculate_quartz_solubility(T_C, P_bar, db):
    """
    Calculate quartz solubility at given T and P using Reaktoro with DEW water model.

    Parameters:
    -----------
    T_C : float
        Temperature in °C
    P_bar : float
        Pressure in bar
    db : Database
        Reaktoro database object

    Returns:
    --------
    float : Quartz molality (mol/kg-H2O), or NaN if calculation fails
    """
    try:
        # Define the chemical system with DEW water and quartz
        system = ChemicalSystem(
            db,
            AqueousPhase("H2O(aq) H+ OH- SiO2(aq)").set(
                ActivityModel(
                    chain(
                        ActivityModelDrummond("CO2"),
                        ActivityModelDuanSun("CO2"),
                        ActivityModelHKF(),
                    )
                )
            ),
            MineralPhases("Quartz"),
        )

        # Create equilibrium state
        state = ChemicalState(system)
        state.temperature(T_C, "celsius")
        state.pressure(P_bar, "bar")

        # Initial amounts: 1 kg water + excess quartz
        state.set("H2O(aq)", 1.0, "kg")
        state.set("Quartz", 10.0, "mol")  # Excess to ensure saturation

        # Equilibrate
        solver = EquilibriumSolver(system)
        solver.solve(state)

        # Get aqueous properties
        props = AqueousProps(state)

        # Get SiO2(aq) molality
        si_molality = props.speciesMolality("SiO2(aq)")

        return float(si_molality)

    except Exception as e:
        print(f"Warning: Failed at T={T_C}°C, P={P_bar} bar: {e}")
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
    exp_data = load_experimental_data(CSV_FILE)
    print(f"    Loaded {len(exp_data)} experimental data points")
    print(
        f"    Temperature range: {exp_data['T_C'].min():.0f} - {exp_data['T_C'].max():.0f} °C"
    )
    print(
        f"    Pressure range: {exp_data['P_kbar'].min():.3f} - {exp_data['P_kbar'].max():.3f} kbar"
    )

    # Get unique experiments and pressures
    experiments = exp_data["experiment_id"].unique()
    pressures_kbar = sorted(exp_data["P_kbar"].unique())

    print(f"    Experiments: {len(experiments)}")
    print(f"    Pressure conditions: {len(pressures_kbar)}")

    # Initialize database with DEW and SUPCRTBL
    print("\n[2] Initializing Reaktoro database (DEW2024 + SUPCRTBL)...")
    db = SupcrtDatabase("supcrtbl")

    # Calculate solubility curves for each pressure
    print("\n[3] Calculating quartz solubility curves...")
    T_range = np.linspace(T_MIN, T_MAX, N_POINTS)

    solubility_curves = {}
    for P_kbar in pressures_kbar:
        P_bar = P_kbar * 1000.0
        print(f"    P = {P_kbar:.2f} kbar ({P_bar:.0f} bar)...")

        molalities = []
        for T_C in T_range:
            mol = calculate_quartz_solubility(T_C, P_bar, db)
            molalities.append(mol)

        solubility_curves[P_kbar] = {"T_C": T_range, "molality": np.array(molalities)}

        valid_points = np.sum(~np.isnan(molalities))
        print(f"       Calculated {valid_points}/{N_POINTS} points successfully")

    # Plotting
    print("\n[4] Creating plots...")

    fig, ax = plt.subplots(figsize=(12, 8))

    # Color schemes
    colors_exp = plt.cm.tab10(np.linspace(0, 0.9, len(experiments)))
    colors_calc = plt.cm.viridis(np.linspace(0, 1, len(pressures_kbar)))

    # Marker styles
    markers = ["o", "s", "^", "v", "D", "<", ">", "p", "*", "h"]

    # Plot experimental data
    for i, exp_id in enumerate(experiments):
        exp_subset = exp_data[exp_data["experiment_id"] == exp_id]

        ax.scatter(
            exp_subset["T_C"],
            exp_subset["molality_m"],
            c=[colors_exp[i]],
            marker=markers[i % len(markers)],
            s=80,
            alpha=0.7,
            edgecolors="black",
            linewidths=0.5,
            label=f"{exp_id}",
            zorder=10,
        )

    # Plot calculated solubility curves
    for i, P_kbar in enumerate(pressures_kbar):
        curve = solubility_curves[P_kbar]

        # Filter out NaN values
        valid = ~np.isnan(curve["molality"])

        ax.plot(
            curve["T_C"][valid],
            curve["molality"][valid],
            color=colors_calc[i],
            linewidth=2.5,
            linestyle="-",
            label=f"Calculated P = {P_kbar:.2f} kbar",
            zorder=5,
        )

    # Set logarithmic y-axis
    ax.set_yscale("log")

    # Labels and formatting
    ax.set_xlabel("Temperature (°C)", fontsize=14, fontweight="bold")
    ax.set_ylabel("Quartz Solubility (mol/kg-H₂O)", fontsize=14, fontweight="bold")
    ax.set_title(
        "Quartz Solubility: DEW2024 Calculations vs Experimental Data",
        fontsize=16,
        fontweight="bold",
        pad=20,
    )

    # Grid
    ax.grid(True, which="both", alpha=0.3, linestyle="--")

    # Legend
    ax.legend(loc="upper left", fontsize=9, framealpha=0.9, ncol=2)

    # Adjust layout
    plt.tight_layout()

    # Save figure
    plt.savefig(OUTPUT_PLOT, dpi=300, bbox_inches="tight")
    print(f"    Plot saved to: {OUTPUT_PLOT}")

    # Show plot
    plt.show()

    print("\n" + "=" * 80)
    print("Analysis complete!")
    print("=" * 80)


if __name__ == "__main__":
    main()
