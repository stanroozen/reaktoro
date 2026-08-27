"""
Mineral Solubility Analysis using Reaktoro with DEW2024
Generic framework for comparing calculated solubilities with experimental data
Includes per-kbar-category temperature ranges, Psat curve, and uncertainty analysis

Easily adaptable for different minerals by changing MINERAL_CONFIG section
"""

try:
    import autodiff  # noqa: F401
except ModuleNotFoundError:
    autodiff = None
import pandas as pd
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os
import sys
import re

# Try to import Reaktoro; fall back to local extension module if not installed
try:
    from reaktoro import *  # noqa: F401,F403
except ModuleNotFoundError:
    # Add local build path for reaktoro4py if available
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    BENCHMARK_DIR = os.path.dirname(SCRIPT_DIR)
    ROOT_DIR = os.path.dirname(BENCHMARK_DIR)
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

# Silence repeated non-convergence warnings
try:
    Warnings.disable(906)
except Exception:
    pass

# =============================================================================
# Configuration
# =============================================================================

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_MINERAL_SOL_DIR = os.path.dirname(SCRIPT_DIR)
_DEW_BENCHMARK_DIR = os.path.dirname(_MINERAL_SOL_DIR)
_ROOT_DIR = os.path.dirname(_DEW_BENCHMARK_DIR)
SUPCRTBL_DATABASE_FILE = os.path.join(
    _ROOT_DIR, "embedded", "databases", "reaktoro", "supcrtbl.yaml"
)

# ============================================================
# MINERAL CONFIGURATION - Change these for different minerals
# ============================================================
MINERAL_CONFIG = {
    # Mineral identification
    "mineral_name": "Wollastonite",  # Name in database
    "mineral_formula": "CaSiO3",  # Chemical formula
    "target_element": "Si",  # Element to report total dissolved molality
    "solute_species": "SiO2(aq)",  # Primary aqueous species (DEW name)
    # Aqueous species to include (besides water, H+, OH-)
    "include_elements": ["Ca", "Si", "O", "H"],
    "aqueous_species": "",
    # File paths
    "csv_file": "wollastonite_DEW_testset.csv",
    "output_prefix": "wollastonite",  # Prefix for output files
    # Plot settings
    "plot_title": "Wollastonite Solubility",
    "y_label": "Total Si molality (mol/kg-H₂O)",
}
# ============================================================

# Generate file paths from mineral config
CSV_FILE = os.path.join(SCRIPT_DIR, MINERAL_CONFIG["csv_file"])
OUTPUT_PLOT_HIGH = os.path.join(
    SCRIPT_DIR,
    f"{MINERAL_CONFIG['output_prefix']}_solubility_comparison_high_P_dew24.png",
)
RESIDUALS_PLOT = os.path.join(
    SCRIPT_DIR, f"{MINERAL_CONFIG['output_prefix']}_solubility_residuals_dew24.png"
)

# Calculation settings
T_MIN, T_MAX = 150, 550
N_POINTS = 100
DEFAULT_PRESSURES = [0.5, 1.0, 2.0, 5.0, 10.0]

# Water model configuration
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


def build_author_markers(references):
    """Assign a unique marker shape to each author/reference."""
    marker_cycle = [
        "o",
        "s",
        "^",
        "D",
        "v",
        "p",
        "P",
        "X",
        "*",
        "h",
        "H",
        "<",
        ">",
        "8",
    ]
    unique_authors = sorted({str(ref) for ref in references if pd.notna(ref)})
    return {
        author: marker_cycle[i % len(marker_cycle)]
        for i, author in enumerate(unique_authors)
    }


def get_solute_species_list(mineral_config):
    """Return a unique, ordered list of aqueous solute species to sum."""
    species = []
    for key in ("solute_species", "aqueous_species"):
        value = mineral_config.get(key, "")
        if not value:
            continue
        for name in value.split():
            if name and name not in species:
                species.append(name)
    return species


def total_solute_molality(aqprops, solute_species_list):
    """Sum molalities for all aqueous solute species in the list."""
    total = 0.0
    for name in solute_species_list:
        total += float(aqprops.speciesMolality(name))
    return total


def total_element_molality(aqprops, mineral_config, solute_species_list):
    """Return total dissolved element molality (uses element stoichiometry)."""
    element = mineral_config.get("target_element")
    if element:
        element_species = mineral_config.get("_element_species_map", {}).get(
            element, []
        )
        if element_species:
            total = 0.0
            for name, coeff in element_species:
                total += float(aqprops.speciesMolality(name)) * float(coeff)
            return total
        return float(aqprops.elementMolality(element))
    return total_solute_molality(aqprops, solute_species_list)


def validate_aqueous_species(dew_db, aqueous_species_str):
    """Check that all aqueous species exist in the DEW database."""
    missing = []
    for name in aqueous_species_str.split():
        try:
            dew_db.species(name)
        except Exception:
            missing.append(name)
    if missing:
        print(
            "    WARNING: Missing aqueous species in DEW database: "
            + ", ".join(missing)
        )
    else:
        print("    All requested aqueous species found in DEW database.")


def aqueous_species_by_elements(dew_db, elements):
    """Return aqueous species names composed only of the given elements."""
    allowed = set(elements)
    pattern = re.compile(r"[A-Z][a-z]?")
    names = []
    for species in dew_db.species():
        formula = str(species.formula())
        elems = set(pattern.findall(formula))
        if elems and elems.issubset(allowed):
            names.append(species.name())
    return sorted(set(names))


def element_species_coeffs(dew_db, species_names, element):
    """Return (species, coefficient) pairs for the given element."""
    entries = []
    for name in species_names:
        try:
            species = dew_db.species(name)
            coeff = float(species.elements().coefficient(element))
        except Exception:
            continue
        if coeff != 0.0:
            entries.append((name, coeff))
    return entries


def build_system(dew_db, supcrt_db, mineral_config, water_config=None):
    """Build ChemicalSystem combining DEW aqueous species with mineral from SUPCRTBL.

    Args:
        dew_db: DEW database
        supcrt_db: SUPCRT database
        mineral_config: Dictionary with mineral-specific settings
        water_config: Optional water model configuration

    Returns:
        ChemicalSystem configured for the specified mineral
    """
    if water_config is None:
        water_config = DEW_CONFIG

    # Get mineral species from database
    mineral_name = mineral_config["mineral_name"]
    mineral_species = supcrt_db.species(mineral_name)
    combined_db = Database(dew_db.species())
    combined_db.addSpecies(mineral_species)

    # Build aqueous phase species list
    base_species = "H2O(aq) H+(aq) OH-(aq)"
    solute = mineral_config["solute_species"]
    additional = mineral_config.get("aqueous_species", "")

    include_elements = mineral_config.get("include_elements")
    if include_elements:
        element_species = aqueous_species_by_elements(dew_db, include_elements)
        additional = " ".join(element_species)

    if additional:
        aqueous_species_str = f"{base_species} {solute} {additional}"
    else:
        aqueous_species_str = f"{base_species} {solute}"

    # Cache element-bearing species for explicit element molality sums
    species_names = aqueous_species_str.split()
    target_element = mineral_config.get("target_element")
    if target_element:
        mineral_config["_element_species_map"] = {
            target_element: element_species_coeffs(
                dew_db, species_names, target_element
            )
        }

    validate_aqueous_species(dew_db, aqueous_species_str)
    aqueous = AqueousPhase(aqueous_species_str)

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

    mineral = MineralPhase(mineral_name)
    system = ChemicalSystem(combined_db, aqueous, mineral)

    print(f"✓ System built for {mineral_name} solubility")
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
    mineral_name = MINERAL_CONFIG["mineral_name"]
    print(f"{mineral_name} Solubility Analysis - Reaktoro DEW2024")
    print(f"Mineral: {MINERAL_CONFIG['mineral_formula']}")
    print(f"Solute species: {MINERAL_CONFIG['solute_species']}")
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
        dew_db = DEWDatabase("dew2024-aqueous")
        print("    Successfully loaded DEW2024 aqueous database")
    except Exception as e:
        print(f"    ERROR: Failed to load DEW database: {e}")
        raise

    try:
        if not os.path.exists(SUPCRTBL_DATABASE_FILE):
            raise FileNotFoundError(
                f"SUPCRTBL database not found: {SUPCRTBL_DATABASE_FILE}"
            )
        supcrt_db = Database.fromFile(SUPCRTBL_DATABASE_FILE)
        print(
            "    Successfully loaded SUPCRTBL mineral database from "
            f"{SUPCRTBL_DATABASE_FILE}"
        )
    except Exception as e:
        print(f"    ERROR: Failed to load SUPCRTBL database: {e}")
        raise

    system = build_system(dew_db, supcrt_db, MINERAL_CONFIG)

    # Calculate solubility curves for each experimental pressure (drops NaN)
    mineral_name = MINERAL_CONFIG["mineral_name"]
    solute_species = MINERAL_CONFIG["solute_species"]
    solute_species_list = get_solute_species_list(MINERAL_CONFIG)
    print(f"\n[3] Calculating {mineral_name.lower()} solubility curves...")
    solubility_curves = {}
    element_curves = {"Si": {}}

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

        # Modern EquilibriumSpecs pattern (matches official tutorial)
        specs = EquilibriumSpecs(system)
        specs.temperature()
        specs.pressure()

        solver = EquilibriumSolver(specs)
        conditions = EquilibriumConditions(specs)

        state = ChemicalState(system)
        state.set("H2O(aq)", 1.0, "kg")
        state.set("H+(aq)", 1e-8, "mol")
        state.set("OH-(aq)", 1e-8, "mol")
        state.set(solute_species, 1e-6, "mol")
        state.set(mineral_name, 10.0, "mol")

        molalities = []
        element_m = {"Si": []}
        for T_C in T_range:
            conditions.temperature(float(T_C), "celsius")
            conditions.pressure(float(P_bar), "bar")
            result = solver.solve(state, conditions)

            if result.succeeded():
                try:
                    aqprops = AqueousProps(state)
                    molality = total_element_molality(
                        aqprops, MINERAL_CONFIG, solute_species_list
                    )
                    for elem in ("Si",):
                        try:
                            element_m[elem].append(float(aqprops.elementMolality(elem)))
                        except Exception:
                            element_m[elem].append(np.nan)
                except Exception:
                    molality = np.nan
                    for elem in ("Si",):
                        element_m[elem].append(np.nan)
            else:
                molality = np.nan
                for elem in ("Si",):
                    element_m[elem].append(np.nan)

            molalities.append(molality)

        solubility_curves[P_kbar] = {"T_C": T_range, "molality": np.array(molalities)}
        for elem in ("Si",):
            element_curves[elem][P_kbar] = {
                "T_C": T_range,
                "molality": np.array(element_m[elem]),
            }

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

    # Modern EquilibriumSpecs pattern for Psat calculations
    specs_psat = EquilibriumSpecs(system)
    specs_psat.temperature()
    specs_psat.pressure()

    solver_psat = EquilibriumSolver(specs_psat)
    conditions_psat = EquilibriumConditions(specs_psat)

    state_psat = ChemicalState(system)
    state_psat.set("H2O(aq)", 1.0, "kg")
    state_psat.set("H+(aq)", 1e-8, "mol")
    state_psat.set("OH-(aq)", 1e-8, "mol")
    state_psat.set(solute_species, 1e-6, "mol")
    state_psat.set(mineral_name, 10.0, "mol")

    psat_molalities = []
    psat_element_m = {"Si": []}
    for i, T_C in enumerate(T_psat_range):
        if not valid_temps[i]:
            psat_molalities.append(np.nan)
            for elem in ("Si",):
                psat_element_m[elem].append(np.nan)
            continue

        P_bar = P_psat_values[i] * 1000.0

        conditions_psat.temperature(float(T_C), "celsius")
        conditions_psat.pressure(float(P_bar), "bar")

        result = solver_psat.solve(state_psat, conditions_psat)
        if result.succeeded():
            try:
                aqprops = AqueousProps(state_psat)
                molality = total_element_molality(
                    aqprops, MINERAL_CONFIG, solute_species_list
                )
                for elem in ("Si",):
                    try:
                        psat_element_m[elem].append(
                            float(aqprops.elementMolality(elem))
                        )
                    except Exception:
                        psat_element_m[elem].append(np.nan)
            except Exception:
                molality = np.nan
                for elem in ("Si",):
                    psat_element_m[elem].append(np.nan)
        else:
            molality = np.nan
            for elem in ("Si",):
                psat_element_m[elem].append(np.nan)

        psat_molalities.append(molality)

    solubility_curves["Psat"] = {
        "T_C": T_psat_range,
        "P_kbar": P_psat_values,
        "molality": np.array(psat_molalities),
    }
    for elem in ("Si",):
        element_curves[elem]["Psat"] = {
            "T_C": T_psat_range,
            "P_kbar": P_psat_values,
            "molality": np.array(psat_element_m[elem]),
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

    author_markers = (
        build_author_markers(exp_data["reference"]) if len(exp_data) > 0 else {}
    )

    # =========================================================================
    # PLOT 1: High Pressure (>=1 kbar)
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

    # Plot calculated curves for high pressures (Si)
    for P_kbar in high_P_all_pressures:
        if P_kbar not in element_curves["Si"]:
            continue
        si_curve = element_curves["Si"][P_kbar]

        si_mol = np.array(si_curve["molality"], dtype=float)
        si_valid = np.isfinite(si_mol)

        if np.any(si_valid):
            ax2.plot(
                si_curve["T_C"][si_valid],
                si_mol[si_valid],
                color=P_to_color_high[P_kbar],
                linewidth=2.0,
                linestyle="-",
                label=f"Calc Si P={P_kbar:.2f} kbar",
                zorder=5,
            )

    ax2.set_yscale("log")
    ax2.set_xlabel("Temperature (°C)", fontsize=14, fontweight="bold")
    ax2.set_ylabel(MINERAL_CONFIG["y_label"], fontsize=14, fontweight="bold")
    ax2.set_title(
        f"{MINERAL_CONFIG['plot_title']}: High Pressure (>=1 kbar)",
        fontsize=16,
        fontweight="bold",
        pad=20,
    )
    ax2.grid(True, which="both", alpha=0.3, linestyle="--")

    # Add database/model info annotation
    info_text = (
        f"DEW24 (species) + SUPCRTBL ({mineral_name}) + Zhang-Duan 2005 EOS (H₂O)"
    )
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
    # PLOT 2: Absolute and Relative Differences (measured - predicted)
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
        state.set("H2O(aq)", 1.0, "kg")
        state.set("H+(aq)", 1e-8, "mol")
        state.set("OH-(aq)", 1e-8, "mol")
        state.set(solute_species, 1e-6, "mol")
        state.set(mineral_name, 10.0, "mol")
        for extra_mineral in MINERAL_CONFIG.get("additional_minerals", []):
            state.set(extra_mineral, 10.0, "mol")
        state.pressure(P_bar, "bar")
        state.temperature(float(T_C), "celsius")

        result = solver_resid.solve(state)
        if result.succeeded():
            try:
                aqprops = AqueousProps(state)
                molality = total_element_molality(
                    aqprops, MINERAL_CONFIG, solute_species_list
                )
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

