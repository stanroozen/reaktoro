"""
Quartz Solubility Analysis using Reaktoro with DEW2019 + SUPCRT07
Compares calculated solubilities with experimental data
Uses Quartz thermodynamics from SUPCRT07 database
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
    PYD_DIR = os.path.join(ROOT_DIR, "build-msvc", "Reaktoro", "Release")
    if os.path.isdir(PYD_DIR) and PYD_DIR not in sys.path:
        sys.path.insert(0, PYD_DIR)
    try:
        from reaktoro4py import *  # noqa: F401,F403

        print("Using local reaktoro4py extension from build-msvc.")
    except ModuleNotFoundError as e:
        raise ModuleNotFoundError(
            "Could not import Reaktoro. Install the 'reaktoro' package or ensure reaktoro4py is on PYTHONPATH."
        ) from e

# Silence repeated non-convergence warnings
try:
    Warnings.disable(906)
except Exception:
    pass

# =============================================================================
# Configuration
# =============================================================================

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CSV_FILE = os.path.join(SCRIPT_DIR, "quartz_DEW_testset.csv")
OUTPUT_PLOT_LOW = os.path.join(
    SCRIPT_DIR, "quartz_solubility_comparison_low_P_supcrt07.png"
)
OUTPUT_PLOT_HIGH = os.path.join(
    SCRIPT_DIR, "quartz_solubility_comparison_high_P_supcrt07.png"
)
RESIDUALS_PLOT = os.path.join(SCRIPT_DIR, "quartz_solubility_residuals_supcrt07.png")

T_MIN, T_MAX = 150, 550
N_POINTS = 100
DEFAULT_PRESSURES = [0.5, 1.0, 2.0, 5.0, 10.0]

DEW_CONFIG = {
    "eos_model": "ZhangDuan2005",
    "dielectric_model": "PowerFunction",
    "gibbs_model": "DewIntegral",
    "born_model": "Shock92Dew",
}

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


def build_system(dew_db, supcrt_db, water_config=None):
    """Build ChemicalSystem combining DEW aqueous species with Quartz from SUPCRT07."""
    if water_config is None:
        water_config = DEW_CONFIG

    quartz_species = supcrt_db.species("Quartz")
    combined_db = Database(dew_db.species())
    combined_db.addSpecies(quartz_species)

    aqueous = AqueousPhase(
        "WATER,AQ H+ OH- SiO2_aq H2_aq O2_aq HO2- HSiO3- Si2O4_aq Si3O6_aq"
    )

    try:
        params = StandardThermoModelParamsDEW()
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
        aqueous.setActivityModel(ActivityModelDEW())
        eos_name = water_config.get("eos_model", "ZhangDuan2005")
        print(f"✓ DEW configured: EOS={eos_name}; phase activity=ActivityModelDEW")

    except Exception as e:
        print(f"Warning: Could not configure DEW: {e}")
        print("  Falling back to default ActivityModelDEW()")
        try:
            aqueous.setActivityModel(ActivityModelDEW())
        except NameError:
            aqueous.setActivityModel(ActivityModelHKF())

    mineral = MineralPhase("Quartz")
    system = ChemicalSystem(combined_db, aqueous, mineral)
    return system


def load_experimental_data(csv_file):
    """Load and organize experimental data from CSV."""
    df = pd.read_csv(csv_file)
    df = df[["T_C", "P_kbar", "molality_m", "reference", "experiment_type"]].copy()
    df["P_bar"] = df["P_kbar"] * 1000.0
    # Keep NaN pressures for Hemley and other saturation curve experiments
    df["experiment_id"] = df["reference"] + " (" + df["experiment_type"] + ")"

    # Mark experiments ON the saturation pressure curve
    # Kennedy controlled-pressure experiments (0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.75 kbar) are non-Psat
    # All other experiments are Psat
    kennedy_controlled_pressures = {0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.75}
    is_kennedy = df["reference"].str.contains("Kennedy", case=False, na=False)
    is_kennedy_controlled = (
        is_kennedy
        & df["P_kbar"].notna()
        & df["P_kbar"].isin(kennedy_controlled_pressures)
    )

    # Everything that's NOT Kennedy controlled-pressure is Psat
    df["is_psat"] = ~is_kennedy_controlled

    df = df.sort_values(["is_psat", "P_kbar", "T_C"]).reset_index(drop=True)

    return df


# =============================================================================
# Main Script
# =============================================================================


def main():
    print("=" * 80)
    print("Quartz Solubility Analysis - DEW2019 + SUPCRT07")
    print("=" * 80)

    # Load experimental data
    print("\n[1] Loading experimental data...")
    if not os.path.exists(CSV_FILE):
        print(f"    WARNING: Experimental data file not found: {CSV_FILE}")
        exp_data = pd.DataFrame()
    else:
        exp_data = load_experimental_data(CSV_FILE)
        print(f"    Loaded {len(exp_data)} experimental data points")
        if len(exp_data) > 0:
            print(
                f"    Temperature range: {exp_data['T_C'].min():.0f} - {exp_data['T_C'].max():.0f} °C"
            )
            print(
                f"    Pressure range: {exp_data['P_kbar'].min():.3f} - {exp_data['P_kbar'].max():.3f} kbar"
            )

    if len(exp_data) > 0:
        experiments = exp_data["experiment_id"].unique()
        # Separate Psat and non-Psat experiments
        psat_data = exp_data[exp_data["is_psat"]]
        non_psat_data = exp_data[~exp_data["is_psat"]]

        # For non-Psat experiments, use only Kennedy controlled pressures
        kennedy_controlled_pressures = [0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.75]
        pressures_kbar = sorted(
            [
                p
                for p in kennedy_controlled_pressures
                if p in non_psat_data["P_kbar"].values
            ]
        )

        print(f"    Experiments: {len(experiments)}")
        print(f"    Psat experiments: {len(psat_data)} points")
        print(
            f"    Non-Psat experiments: {len(non_psat_data)} points at {len(pressures_kbar)} pressures"
        )
    else:
        experiments = []
        pressures_kbar = DEFAULT_PRESSURES
        psat_data = pd.DataFrame()
        non_psat_data = pd.DataFrame()

    # Initialize databases
    print("\n[2] Initializing Reaktoro databases...")
    try:
        dew_db = DEWDatabase("dew2019-aqueous")
        print("    Successfully loaded DEW2019 aqueous database")
    except Exception as e:
        print(f"    ERROR: Failed to load DEW database: {e}")
        raise

    try:
        supcrt_db = SupcrtDatabase("supcrt07")
        print("    Successfully loaded SUPCRT07 mineral database")
    except Exception as e:
        print(f"    ERROR: Failed to load SUPCRT07 database: {e}")
        raise

    system = build_system(dew_db, supcrt_db)

    # Calculate solubility curves for each experimental pressure (drops NaN)
    print("\n[3] Calculating quartz solubility curves...")
    solubility_curves = {}

    pressures_for_curves = sorted(exp_data["P_kbar"].dropna().unique())

    for P_kbar in pressures_for_curves:
        P_bar = P_kbar * 1000.0
        print(f"    P = {P_kbar:.3f} kbar ({P_bar:.0f} bar)...")

        # Determine T range from experiments at this pressure (±5%)
        P_tol = 0.05 * P_kbar
        exp_at_P = exp_data[
            (exp_data["P_kbar"] >= P_kbar - P_tol)
            & (exp_data["P_kbar"] <= P_kbar + P_tol)
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
                f"       Calculated {valid_points}/{N_POINTS} points (T: {first_T:.0f}-{last_T:.0f}°C)"
            )

    # Calculate Psat solubility curve (following saturation pressure)
    print("    P = Psat curve...")
    # Limit to valid range where Psat is defined (below critical point)
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
        state.set("WATER,AQ", 1.0, "kg")
        state.set("H+", 1e-8, "mol")
        state.set("OH-", 1e-8, "mol")
        state.set("SiO2_aq", 1e-6, "mol")
        state.set("Quartz", 10.0, "mol")
        state.pressure(float(P_bar), "bar")
        state.temperature(float(T_C), "celsius")

        result = solver.solve(state)
        if result.succeeded():
            try:
                aqprops = AqueousProps(state)
                molality = float(aqprops.speciesMolality("SiO2_aq"))
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

    valid_psat_points = np.sum(~np.isnan(psat_molalities))
    print(f"       Calculated {valid_psat_points}/{N_POINTS} points along Psat curve")

    # Plotting
    print("\n[4] Creating plots...")

    # Separate data into low (<1 kbar) and high (≥1 kbar) pressure ranges
    low_P_threshold = 1.0

    low_P_data = non_psat_data[non_psat_data["P_kbar"] < low_P_threshold]
    high_P_data = non_psat_data[non_psat_data["P_kbar"] >= low_P_threshold]

    low_P_pressures = (
        sorted(low_P_data["P_kbar"].unique()) if len(low_P_data) > 0 else []
    )
    high_P_pressures = (
        sorted(high_P_data["P_kbar"].unique()) if len(high_P_data) > 0 else []
    )

    # =========================================================================
    # PLOT 1: Low Pressure (<1 kbar) with Psat
    # =========================================================================
    print("    Creating low-pressure plot (<1 kbar)...")
    fig1, ax1 = plt.subplots(figsize=(14, 8))

    # Generate colors for low-pressure experiments
    n_low = max(len(low_P_pressures), 1)
    colors_low = plt.cm.viridis(np.linspace(0, 0.9, n_low))
    P_to_color_low = {
        P: colors_low[i % len(colors_low)] for i, P in enumerate(low_P_pressures)
    }

    # Map authors to different marker shapes
    author_markers = {
        "Kennedy_1950": "o",
        "Hemley_1980": "^",
        "Morey_Fournier_Rowe_1962": "s",
        "Walther_Orville_1983": "D",
        "Manning_1994": "v",
        "Newton_Manning_2000": "p",
    }

    # Plot low-pressure experimental data
    for P_kbar in low_P_pressures:
        P_tol = 0.05 * P_kbar if P_kbar > 0.1 else 0.01
        subset = low_P_data[
            (low_P_data["P_kbar"] >= P_kbar - P_tol)
            & (low_P_data["P_kbar"] <= P_kbar + P_tol)
        ]
        if len(subset) == 0:
            continue

        # Plot by author within this pressure group
        for author in subset["reference"].unique():
            author_subset = subset[subset["reference"] == author]
            marker = author_markers.get(author, "o")
            ax1.scatter(
                author_subset["T_C"],
                author_subset["molality_m"],
                c=[P_to_color_low[P_kbar]],
                marker=marker,
                s=70,
                alpha=0.7,
                edgecolors="black",
                linewidths=0.4,
                label=f"Exp P={P_kbar:.2f} kbar ({author})",
                zorder=10,
            )

    # Plot Psat experiments (low-pressure only on this plot)
    if len(psat_data) > 0:
        # Include psat with P < 1.0 kbar OR P is NaN (saturation curve experiments without explicit pressure)
        low_psat_data = psat_data[
            (psat_data["P_kbar"] < 1.0) | (psat_data["P_kbar"].isna())
        ]
        for author in low_psat_data["reference"].unique():
            author_psat = low_psat_data[low_psat_data["reference"] == author]
            marker = author_markers.get(author, "s")
            ax1.scatter(
                author_psat["T_C"],
                author_psat["molality_m"],
                c="purple",
                marker=marker,
                s=80,
                alpha=0.8,
                edgecolors="darkviolet",
                linewidths=0.5,
                label=f"Exp P=Psat ({author})",
                zorder=11,
            )

    # Plot calculated curves for low pressures
    for P_kbar in low_P_pressures:
        if P_kbar not in solubility_curves:
            continue
        curve = solubility_curves[P_kbar]
        valid = ~np.isnan(curve["molality"])
        ax1.plot(
            curve["T_C"][valid],
            curve["molality"][valid],
            color=P_to_color_low[P_kbar],
            linewidth=2.0,
            linestyle="-",
            label=f"Calc P={P_kbar:.2f} kbar",
            zorder=5,
        )

    # Plot Psat curve
    if "Psat" in solubility_curves:
        curve_psat = solubility_curves["Psat"]
        valid_psat_m = ~np.isnan(curve_psat["molality"])
        ax1.plot(
            curve_psat["T_C"][valid_psat_m],
            curve_psat["molality"][valid_psat_m],
            color="purple",
            linewidth=3.0,
            linestyle="-",
            label="Calc P=Psat",
            zorder=6,
            alpha=0.9,
        )

    ax1.set_yscale("log")
    ax1.set_ylim(1e-4, 1e-1)
    ax1.set_xlabel("Temperature (°C)", fontsize=14, fontweight="bold")
    ax1.set_ylabel("Quartz Solubility (mol/kg-H₂O)", fontsize=14, fontweight="bold")
    ax1.set_title(
        "Quartz Solubility: Low Pressure (<1 kbar)",
        fontsize=16,
        fontweight="bold",
        pad=20,
    )
    ax1.grid(True, which="both", alpha=0.3, linestyle="--")

    # Add database/model info annotation
    info_text = "DEW19 (species) + SUPCRT07 (quartz) + Zhang-Duan 2005 EOS (H₂O)"
    ax1.text(
        0.02,
        0.98,
        info_text,
        transform=ax1.transAxes,
        fontsize=8,
        verticalalignment="top",
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5),
    )

    ax1.legend(
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        fontsize=9,
        framealpha=0.9,
        ncol=1,
    )
    plt.tight_layout()
    plt.savefig(OUTPUT_PLOT_LOW, dpi=300, bbox_inches="tight")
    print(f"    Low-P plot saved to: {OUTPUT_PLOT_LOW}")
    plt.close(fig1)

    # =========================================================================
    # PLOT 2: High Pressure (>=1 kbar)
    # =========================================================================
    print("    Creating high-pressure plot (>=1 kbar)...")
    fig2, ax2 = plt.subplots(figsize=(14, 8))

    # Collect all experiments >= 1.0 kbar (from both psat and non-psat data)
    high_P_all_data = pd.concat(
        [
            psat_data[(psat_data["P_kbar"] >= 1.0) & (psat_data["P_kbar"].notna())],
            non_psat_data[non_psat_data["P_kbar"] >= 1.0],
        ],
        ignore_index=True,
    )

    # Get all unique pressures for these high-P experiments
    high_P_all_pressures = sorted(high_P_all_data["P_kbar"].unique())

    # Generate colors for high-pressure experiments
    n_high = max(len(high_P_all_pressures), 1)
    colors_high = plt.cm.plasma(np.linspace(0, 0.9, n_high))
    P_to_color_high = {
        P: colors_high[i % len(colors_high)] for i, P in enumerate(high_P_all_pressures)
    }

    # Plot high-pressure experimental data by actual pressure
    for P_kbar in high_P_all_pressures:
        P_tol = 0.05 * P_kbar
        subset = high_P_all_data[
            (high_P_all_data["P_kbar"] >= P_kbar - P_tol)
            & (high_P_all_data["P_kbar"] <= P_kbar + P_tol)
        ]
        if len(subset) == 0:
            continue

        # Plot by author within this pressure group
        for author in subset["reference"].unique():
            author_subset = subset[subset["reference"] == author]
            marker = author_markers.get(author, "o")
            ax2.scatter(
                author_subset["T_C"],
                author_subset["molality_m"],
                c=[P_to_color_high[P_kbar]],
                marker=marker,
                s=70,
                alpha=0.7,
                edgecolors="black",
                linewidths=0.4,
                label=f"Exp P={P_kbar:.2f} kbar ({author})",
                zorder=10,
            )

    # Plot calculated curves for high pressures
    for P_kbar in high_P_all_pressures:
        if P_kbar not in solubility_curves:
            continue
        curve = solubility_curves[P_kbar]
        valid = ~np.isnan(curve["molality"])
        ax2.plot(
            curve["T_C"][valid],
            curve["molality"][valid],
            color=P_to_color_high[P_kbar],
            linewidth=2.0,
            linestyle="-",
            label=f"Calc P={P_kbar:.2f} kbar",
            zorder=5,
        )

    ax2.set_yscale("log")
    ax2.set_xlabel("Temperature (°C)", fontsize=14, fontweight="bold")
    ax2.set_ylabel("Quartz Solubility (mol/kg-H₂O)", fontsize=14, fontweight="bold")
    ax2.set_title(
        "Quartz Solubility: High Pressure (>=1 kbar)",
        fontsize=16,
        fontweight="bold",
        pad=20,
    )
    ax2.grid(True, which="both", alpha=0.3, linestyle="--")

    # Add database/model info annotation
    info_text = "DEW19 (species) + SUPCRT07 (quartz) + Zhang-Duan 2005 EOS (H₂O)"
    ax2.text(
        0.02,
        0.98,
        info_text,
        transform=ax2.transAxes,
        fontsize=8,
        verticalalignment="top",
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5),
    )

    ax2.legend(
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        fontsize=9,
        framealpha=0.9,
        ncol=1,
    )
    plt.tight_layout()
    plt.savefig(OUTPUT_PLOT_HIGH, dpi=300, bbox_inches="tight")
    print(f"    High-P plot saved to: {OUTPUT_PLOT_HIGH}")
    plt.close(fig2)

    # =========================================================================
    # PLOT 3: Absolute and Relative Differences (measured - predicted)
    # =========================================================================
    print("    Creating residual plots...")

    # Predict molality at each experimental (T, P)
    predicted = []
    solver_resid = EquilibriumSolver(system)
    for _, row in exp_data.iterrows():
        T_C = float(row["T_C"])
        P_kbar_row = row["P_kbar"]

        # For NaN pressures (on Psat curve), compute Psat at that temperature
        if pd.isna(P_kbar_row):
            P_kbar_row = psat_kbar(T_C)

        if pd.isna(P_kbar_row):
            predicted.append(np.nan)
            continue

        P_bar = float(P_kbar_row) * 1000.0

        state = ChemicalState(system)
        state.set("WATER,AQ", 1.0, "kg")
        state.set("H+", 1e-8, "mol")
        state.set("OH-", 1e-8, "mol")
        state.set("SiO2_aq", 1e-6, "mol")
        state.set("Quartz", 10.0, "mol")
        state.pressure(P_bar, "bar")
        state.temperature(float(T_C), "celsius")

        result = solver_resid.solve(state)
        if result.succeeded():
            try:
                aqprops = AqueousProps(state)
                molality = float(aqprops.speciesMolality("SiO2_aq"))
            except Exception:
                molality = np.nan
        else:
            molality = np.nan

        predicted.append(molality)

    exp_resid = exp_data.copy()
    exp_resid["predicted_m"] = predicted
    exp_resid["abs_diff"] = exp_resid["molality_m"] - exp_resid["predicted_m"]
    exp_resid["rel_diff_pct"] = np.where(
        exp_resid["predicted_m"] > 0,
        100.0 * exp_resid["abs_diff"] / exp_resid["predicted_m"],
        np.nan,
    )

    # Color by pressure; Psat (NaN pressure) in purple
    finite_pressures = exp_resid["P_kbar"].dropna()
    if len(finite_pressures) > 0:
        norm = plt.Normalize(vmin=finite_pressures.min(), vmax=finite_pressures.max())
    else:
        norm = plt.Normalize(vmin=0.0, vmax=1.0)
    cmap = plt.cm.plasma

    fig_res, (ax_abs, ax_rel) = plt.subplots(2, 1, figsize=(12, 10), sharex=True)

    for _, row in exp_resid.iterrows():
        T_C = row["T_C"]
        P_kbar_row = row["P_kbar"]
        author = row["reference"]
        marker = author_markers.get(author, "o")

        if pd.isna(P_kbar_row):
            color = "purple"
        else:
            color = cmap(norm(P_kbar_row))

        ax_abs.scatter(
            T_C,
            row["abs_diff"],
            c=[color],
            marker=marker,
            s=50,
            alpha=0.7,
            edgecolors="black",
            linewidths=0.4,
        )

        ax_rel.scatter(
            T_C,
            row["rel_diff_pct"],
            c=[color],
            marker=marker,
            s=50,
            alpha=0.7,
            edgecolors="black",
            linewidths=0.4,
        )

    ax_abs.axhline(0, color="gray", linestyle="--", linewidth=1.0)
    ax_abs.set_ylabel("Measured - Predicted (mol/kg)", fontsize=12, fontweight="bold")
    ax_abs.grid(True, alpha=0.3, linestyle="--")

    ax_rel.axhline(0, color="gray", linestyle="--", linewidth=1.0)
    ax_rel.set_ylabel("Relative Diff (%)", fontsize=12, fontweight="bold")
    ax_rel.set_xlabel("Temperature (°C)", fontsize=12, fontweight="bold")
    ax_rel.grid(True, alpha=0.3, linestyle="--")

    plt.tight_layout()
    plt.savefig(RESIDUALS_PLOT, dpi=300, bbox_inches="tight")
    print(f"    Residual plots saved to: {RESIDUALS_PLOT}")

    print("\n" + "=" * 80)
    print("Analysis complete!")
    print("=" * 80)


if __name__ == "__main__":
    main()
