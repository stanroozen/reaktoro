"""
Quartz Solubility Analysis using Reaktoro with ThermoFun (MINES16 + AQ17)
Compares calculated solubilities with experimental data
Exact same calculation as dew24 script, but databases:
- Minerals + quartz + water EOS: MINES16 (thermofun)
- Aqueous silica species and water: AQ17 (thermofun)
"""

import pandas as pd
import numpy as np
import matplotlib

matplotlib.use("Agg")
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
    from reaktoro4py import *  # noqa: F401,F403

    print("Using local reaktoro4py extension from build.")

# Silence repeated non-convergence warnings
try:
    Warnings.disable(906)
except Exception:
    pass

# =============================================================================
# Configuration
# =============================================================================

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CSV_FILE = os.path.join(SCRIPT_DIR, "..", "quartz_DEW", "quartz_DEW_testset.csv")
OUTPUT_PLOT_LOW = os.path.join(
    SCRIPT_DIR, "quartz_solubility_comparison_low_P_mines16_aq17.png"
)
OUTPUT_PLOT_HIGH = os.path.join(
    SCRIPT_DIR, "quartz_solubility_comparison_high_P_mines16_aq17.png"
)
RESIDUALS_PLOT = os.path.join(
    SCRIPT_DIR, "quartz_solubility_residuals_mines16_aq17.png"
)

T_MIN, T_MAX = 150, 550
N_POINTS = 100
DEFAULT_PRESSURES = [0.5, 1.0, 2.0, 5.0, 10.0]

# =============================================================================
# Saturation Pressure Functions
# =============================================================================


def psat_bar(T_C):
    """Calculate water saturation pressure (bar) using Antoine equation."""
    if T_C < 0 or T_C > 374:
        return np.nan
    T_K = T_C + 273.15
    A, B, C = 5.40221, 1838.675, -31.737
    log10_P = A - B / (T_K + C)
    return 10**log10_P


def psat_kbar(T_C):
    """Calculate water saturation pressure (kbar)."""
    P = psat_bar(T_C)
    return P / 1000.0 if not np.isnan(P) else np.nan


# =============================================================================
# Helper Functions
# =============================================================================


def build_system_thermofun(aq_db, minerals_db):
    """Build ChemicalSystem with AQ17 aqueous species and Quartz from MINES16."""
    # Select explicit aqueous species by name from AQ17
    aq_species_names = [
        "H2O@",  # water solvent in AQ17 (aqueous)
        "H+",
        "OH-",
        "SiO2@",
        "HSiO3-",
        "NaHSiO3@",
        "MgSiO3@",
        "CaSiO3@",
    ]

    combined_db = Database()
    for name in aq_species_names:
        try:
            combined_db.addSpecies(aq_db.species(name))
        except Exception as e:
            print(f"    WARNING: Skipping missing AQ17 species '{name}': {e}")

    # Add Quartz from MINES16 (avoid AQ17 Quartz to use MINES16 thermodynamics)
    quartz_species = minerals_db.species("Quartz")
    combined_db.addSpecies(quartz_species)

    aqueous = AqueousPhase("H2O@ H+ OH- SiO2@ HSiO3- NaHSiO3@ MgSiO3@ CaSiO3@")

    # Use HKF activity model appropriate for ThermoFun datasets
    try:
        aqueous.setActivityModel(ActivityModelHKF())
        print("âœ“ ThermoFun configured: phase activity=ActivityModelHKF")
    except Exception as e:
        print(f"Warning: Could not configure HKF activity model: {e}")
        aqueous.setActivityModel(ActivityModelDebyeHuckel())

    mineral = MineralPhase("Quartz")
    system = ChemicalSystem(combined_db, aqueous, mineral)
    return system


def load_experimental_data(csv_file):
    """Load and organize experimental data from CSV."""
    df = pd.read_csv(csv_file)
    df = df[["T_C", "P_kbar", "molality_m", "reference", "experiment_type"]].copy()
    df["P_bar"] = df["P_kbar"] * 1000.0
    df["experiment_id"] = df["reference"] + " (" + df["experiment_type"] + ")"

    kennedy_controlled_pressures = {0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.75}
    is_kennedy = df["reference"].str.contains("Kennedy", case=False, na=False)
    is_kennedy_controlled = (
        is_kennedy
        & df["P_kbar"].notna()
        & df["P_kbar"].isin(kennedy_controlled_pressures)
    )
    df["is_psat"] = ~is_kennedy_controlled
    df = df.sort_values(["is_psat", "P_kbar", "T_C"]).reset_index(drop=True)
    return df


def compute_residuals(exp_data, curves):
    """Compute per-point residuals between experiments and calculated curves."""
    residuals = []
    non_psat_data = exp_data[~exp_data["is_psat"]]
    for _, row in non_psat_data.iterrows():
        P_kbar, T_C, m_exp = row["P_kbar"], row["T_C"], row["molality_m"]
        curve = curves.get(P_kbar)
        if curve is None:
            continue
        T_range = curve["T_C"]
        m_calc = curve["molality"]
        if len(T_range) == 0:
            continue
        i = np.abs(T_range - T_C).argmin()
        m_calc_at_T = m_calc[i]
        residuals.append((P_kbar, T_C, m_exp, m_calc_at_T, m_calc_at_T - m_exp))

    res_df = pd.DataFrame(
        residuals,
        columns=["P_kbar", "T_C", "molality_exp", "molality_calc", "residual"],
    )
    return res_df


def main():
    print("Quartz Solubility Analysis - ThermoFun MINES16 + AQ17")
    print("\n[1] Loading experimental dataset...")
    exp_data = load_experimental_data(CSV_FILE)
    non_psat_data = exp_data[~exp_data["is_psat"]]
    psat_data = exp_data[exp_data["is_psat"]]

    # Initialize databases
    print("\n[2] Initializing ThermoFun databases...")
    try:
        aq_db = ThermoFunDatabase("aq17")
        print("    Successfully loaded AQ17 (ThermoFun) aqueous database")
    except Exception as e:
        print(f"    ERROR: Failed to load AQ17 ThermoFun database: {e}")
        raise

    try:
        minerals_db = ThermoFunDatabase("mines16")
        print("    Successfully loaded MINES16 (ThermoFun) minerals database")
    except Exception as e:
        print(f"    ERROR: Failed to load MINES16 ThermoFun database: {e}")
        raise

    system = build_system_thermofun(aq_db, minerals_db)

    # Calculate solubility curves for each experimental pressure (drops NaN)
    print("\n[3] Calculating quartz solubility curves...")
    solubility_curves = {}

    pressures_for_curves = sorted(non_psat_data["P_kbar"].dropna().unique())

    for P_kbar in pressures_for_curves:
        P_bar = P_kbar * 1000.0
        print(f"    P = {P_kbar:.3f} kbar ({P_bar:.0f} bar)...")

        # Determine T range from experiments at this pressure (Â±5%)
        P_tol = 0.05 * P_kbar
        exp_at_P = non_psat_data[
            (non_psat_data["P_kbar"] >= P_kbar - P_tol)
            & (non_psat_data["P_kbar"] <= P_kbar + P_tol)
        ]
        if len(exp_at_P) > 0:
            T_min_cat = exp_at_P["T_C"].min()
            T_max_cat = exp_at_P["T_C"].max()
            T_span = T_max_cat - T_min_cat
            T_min = max(25, T_min_cat - 0.05 * T_span) if T_span > 0 else T_min_cat - 50
            T_max = (
                min(1000, T_max_cat + 0.05 * T_span) if T_span > 0 else T_max_cat + 50
            )
        else:
            T_min, T_max = T_MIN, T_MAX

        T_range = np.linspace(T_min, T_max, N_POINTS)
        solver = EquilibriumSolver(system)
        state = ChemicalState(system)
        state.set("H2O@", 1.0, "kg")
        state.set("H+", 1e-8, "mol")
        state.set("OH-", 1e-8, "mol")
        state.set("SiO2@", 1e-6, "mol")
        state.set("Quartz", 10.0, "mol")
        state.pressure(float(P_bar), "bar")

        molalities = []
        for T_C in T_range:
            state.temperature(float(T_C), "celsius")
            try:
                result = solver.solve(state)
            except Exception:
                result = None
            if result and result.succeeded():
                try:
                    aqprops = AqueousProps(state)
                    molality = float(aqprops.speciesMolality("SiO2@"))
                except Exception:
                    molality = np.nan
            else:
                molality = np.nan
            molalities.append(molality)

        solubility_curves[P_kbar] = {"T_C": T_range, "molality": np.array(molalities)}
        valid_points = np.sum(~np.isnan(molalities))
        if valid_points > 0:
            valid_idx = np.where(~np.isnan(molalities))[0]
            first_T = T_range[valid_idx[0]]
            last_T = T_range[valid_idx[-1]]
            print(
                f"       Calculated {valid_points}/{N_POINTS} points (T: {first_T:.0f}-{last_T:.0f}Â°C)"
            )

    # Calculate Psat solubility curve
    print("    P = Psat curve...")
    if len(psat_data) > 0:
        T_psat_min = max(25, psat_data["T_C"].min() - 25)
        T_psat_max = min(374, psat_data["T_C"].max() + 25)
    else:
        T_psat_min, T_psat_max = 100, 374
    T_psat_range = np.linspace(T_psat_min, T_psat_max, N_POINTS)
    P_psat_values = np.array([psat_kbar(T) for T in T_psat_range])
    valid_temps = ~np.isnan(P_psat_values)

    psat_molalities = []
    for i, T_C in enumerate(T_psat_range):
        if not valid_temps[i]:
            psat_molalities.append(np.nan)
            continue
        P_bar = P_psat_values[i] * 1000.0
        solver = EquilibriumSolver(system)
        state = ChemicalState(system)
        state.set("H2O@", 1.0, "kg")
        state.set("H+", 1e-8, "mol")
        state.set("OH-", 1e-8, "mol")
        state.set("SiO2@", 1e-6, "mol")
        state.set("Quartz", 10.0, "mol")
        state.pressure(float(P_bar), "bar")
        state.temperature(float(T_C), "celsius")
        try:
            result = solver.solve(state)
        except Exception:
            result = None
        if result and result.succeeded():
            try:
                aqprops = AqueousProps(state)
                molality = float(aqprops.speciesMolality("SiO2@"))
            except Exception:
                molality = np.nan
        else:
            molality = np.nan
        psat_molalities.append(molality)

    solubility_curves["Psat"] = {
        "T_C": T_psat_range,
        "P_kbar": P_psat_values,
        "molality": np.array(psat_molalities),
    }
    print(
        f"       Calculated {np.sum(~np.isnan(psat_molalities))}/{N_POINTS} points along Psat curve"
    )

    print("\n[4] Creating plots...")
    low_P_threshold = 1.0
    non_psat_data_local = exp_data[~exp_data["is_psat"]]
    low_P_data = non_psat_data_local[non_psat_data_local["P_kbar"] < low_P_threshold]
    high_P_data = non_psat_data_local[non_psat_data_local["P_kbar"] >= low_P_threshold]

    low_P_pressures = (
        sorted(low_P_data["P_kbar"].unique()) if len(low_P_data) > 0 else []
    )
    high_P_pressures = (
        sorted(high_P_data["P_kbar"].unique()) if len(high_P_data) > 0 else []
    )

    print("    Creating low-pressure plot (<1 kbar)...")
    fig1, ax1 = plt.subplots(figsize=(14, 8))
    n_low = max(len(low_P_pressures), 1)
    colors_low = plt.cm.viridis(np.linspace(0, 0.9, n_low))
    P_to_color_low = {P: colors_low[i] for i, P in enumerate(low_P_pressures)}

    for P_kbar in low_P_pressures:
        data_at_P = low_P_data[low_P_data["P_kbar"] == P_kbar]
        ax1.scatter(
            data_at_P["T_C"],
            data_at_P["molality_m"],
            label=f"Exp {P_kbar:.2f} kbar",
            color=P_to_color_low[P_kbar],
            s=45,
            edgecolors="black",
        )
        if P_kbar in solubility_curves:
            curve = solubility_curves[P_kbar]
            ax1.plot(
                curve["T_C"],
                curve["molality"],
                color=P_to_color_low[P_kbar],
                linewidth=2,
                linestyle="-",
            )

    if "Psat" in solubility_curves:
        curve = solubility_curves["Psat"]
        ax1.plot(
            curve["T_C"],
            curve["molality"],
            color="gray",
            linewidth=2.5,
            linestyle="--",
            label="Psat curve",
        )

    ax1.set_xlabel("Temperature (Â°C)", fontsize=14, fontweight="bold")
    ax1.set_ylabel("Quartz Solubility (mol/kg-Hâ‚‚O)", fontsize=14, fontweight="bold")
    ax1.set_title(
        "Quartz Solubility: Low Pressure (<1 kbar)", fontsize=16, fontweight="bold"
    )
    ax1.legend(fontsize=12, loc="upper left", frameon=True)
    ax1.grid(True, which="both", linestyle="--", alpha=0.4)
    plt.tight_layout()
    fig1.savefig(OUTPUT_PLOT_LOW, dpi=300)

    print("    Creating high-pressure plot (>=1 kbar)...")
    fig2, ax2 = plt.subplots(figsize=(14, 8))
    n_high = max(len(high_P_pressures), 1)
    colors_high = plt.cm.plasma(np.linspace(0.1, 0.95, n_high))
    P_to_color_high = {P: colors_high[i] for i, P in enumerate(high_P_pressures)}

    for P_kbar in high_P_pressures:
        data_at_P = high_P_data[high_P_data["P_kbar"] == P_kbar]
        ax2.scatter(
            data_at_P["T_C"],
            data_at_P["molality_m"],
            label=f"Exp {P_kbar:.2f} kbar",
            color=P_to_color_high[P_kbar],
            s=45,
            edgecolors="black",
        )
        if P_kbar in solubility_curves:
            curve = solubility_curves[P_kbar]
            ax2.plot(
                curve["T_C"],
                curve["molality"],
                color=P_to_color_high[P_kbar],
                linewidth=2,
                linestyle="-",
            )

    ax2.set_xlabel("Temperature (Â°C)", fontsize=14, fontweight="bold")
    ax2.set_ylabel("Quartz Solubility (mol/kg-Hâ‚‚O)", fontsize=14, fontweight="bold")
    ax2.set_title(
        "Quartz Solubility: High Pressure (>=1 kbar)", fontsize=16, fontweight="bold"
    )
    ax2.legend(fontsize=12, loc="upper left", frameon=True)
    ax2.grid(True, which="both", linestyle="--", alpha=0.4)
    plt.tight_layout()
    fig2.savefig(OUTPUT_PLOT_HIGH, dpi=300)

    print("    Creating residual plots...")
    residuals_df = compute_residuals(exp_data, solubility_curves)
    fig3, ax3 = plt.subplots(figsize=(14, 8))
    for P_kbar in sorted(residuals_df["P_kbar"].unique()):
        sub = residuals_df[residuals_df["P_kbar"] == P_kbar]
        ax3.scatter(sub["T_C"], sub["residual"], label=f"{P_kbar:.2f} kbar", s=45)
    ax3.axhline(0, color="black", linewidth=1)
    ax3.set_xlabel("Temperature (Â°C)", fontsize=14, fontweight="bold")
    ax3.set_ylabel("Residual (calc - exp) mol/kg-Hâ‚‚O", fontsize=14, fontweight="bold")
    ax3.set_title(
        "Quartz Solubility Residuals (MINES16 + AQ17)", fontsize=16, fontweight="bold"
    )
    ax3.legend(fontsize=12, loc="upper left", frameon=True)
    ax3.grid(True, which="both", linestyle="--", alpha=0.4)
    plt.tight_layout()
    fig3.savefig(RESIDUALS_PLOT, dpi=300)

    print(
        "\n================================================================================"
    )
    print("Analysis complete!")
    print(
        "================================================================================"
    )


if __name__ == "__main__":
    main()
