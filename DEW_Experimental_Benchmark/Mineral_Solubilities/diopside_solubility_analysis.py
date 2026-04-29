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
import matplotlib.pyplot as plt
import os
import sys
import re

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BENCHMARK_DIR = os.path.dirname(SCRIPT_DIR)
ROOT_DIR = os.path.dirname(BENCHMARK_DIR)
PYD_DIR = os.path.join(ROOT_DIR, "build", "Reaktoro", "Release")
if os.path.isdir(PYD_DIR):
    if PYD_DIR not in sys.path:
        sys.path.insert(0, PYD_DIR)
    try:
        os.add_dll_directory(PYD_DIR)
    except Exception:
        pass

# Try to import Reaktoro; fall back to local extension module if not installed
try:
    from reaktoro import *  # noqa: F401,F403
except ModuleNotFoundError:
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

# ============================================================
# MINERAL CONFIGURATION - Change these for different minerals
# ============================================================
MINERAL_CONFIG = {
    # Mineral identification
    "mineral_name": "Diopside",  # Name in database
    "mineral_formula": "CaMgSi2O6",  # Chemical formula
    "target_element": "Ca",  # Element to report total dissolved molality
    "solute_species": "Ca+2",  # Primary aqueous species
    "include_elements": ["H", "O", "Ca", "Mg", "Si"],
    "exclude_organics": True,
    "excluded_species": ["CaO_aq", "MgO_aq"],
    "additional_minerals": ["Forsterite"],
    # Aqueous species to include (besides water, H+, OH-)
    "aqueous_species": "",
    # File paths
    "csv_file": "diopside_DEW_testset.csv",
    "output_prefix": "diopside",  # Prefix for output files
    # Plot settings
    "plot_title": "Diopside Solubility",
    "y_label": "Element totals (m_Ca, m_Mg, 0.5 m_Si; mol/kg-Hâ‚‚O)",
}
# ============================================================

# Generate file paths from mineral config
CSV_FILE = os.path.join(SCRIPT_DIR, MINERAL_CONFIG["csv_file"])
OUTPUT_PLOT_LOW = os.path.join(
    SCRIPT_DIR,
    f"{MINERAL_CONFIG['output_prefix']}_solubility_comparison_low_P_dew24.png",
)
OUTPUT_PLOT_HIGH = os.path.join(
    SCRIPT_DIR,
    f"{MINERAL_CONFIG['output_prefix']}_solubility_comparison_high_P_dew24.png",
)
RESIDUALS_PLOT = os.path.join(
    SCRIPT_DIR, f"{MINERAL_CONFIG['output_prefix']}_solubility_residuals_dew24.png"
)
SPECIES_PLOT = os.path.join(
    SCRIPT_DIR, f"{MINERAL_CONFIG['output_prefix']}_speciation_dew24.png"
)

# Calculation settings
T_MIN, T_MAX = 650, 900
N_POINTS = 100
DEFAULT_PRESSURES = [7.0, 10.0, 15.0]

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


def build_reversal_markers(values):
    """Assign marker shapes for undersaturated/supersaturated reversal points."""
    marker_map = {
        "U": "^",  # undersaturated
        "S": "v",  # supersaturated
        "?": "o",  # unknown/other
    }
    return {key: marker_map.get(key, "^") for key in set(values)}


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

    excluded = set(MINERAL_CONFIG.get("excluded_species", []))

    if MINERAL_CONFIG.get("exclude_organics", True):
        organic_tokens = (
            "ACET",
            "FORM",
            "METH",
            "ETH",
            "PROP",
            "BUT",
            "PENT",
            "HEX",
            "HEPT",
            "OCT",
            "BENZ",
            "TOLU",
            "LACT",
            "GLYCOL",
            "SUCCIN",
            "GLUTAR",
            "ISOBUT",
        )
        for name in names:
            if any(tok in name for tok in organic_tokens):
                excluded.add(name)

    return sorted(set(names) - excluded)


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
        return float(aqprops.elementMolality(element))
    return total_solute_molality(aqprops, solute_species_list)


def diopside_element_totals(aqprops):
    """Return element totals for Ca, Mg, and 0.5*Si."""
    try:
        m_ca = float(aqprops.elementMolality("Ca"))
    except Exception:
        m_ca = np.nan
    try:
        m_mg = float(aqprops.elementMolality("Mg"))
    except Exception:
        m_mg = np.nan
    try:
        m_si = float(aqprops.elementMolality("Si"))
    except Exception:
        m_si = np.nan

    return m_ca, m_mg, 0.5 * m_si


def validate_aqueous_species(dew_db, aqueous_species_str):
    """Check that all aqueous species exist in the DEW database."""
    missing = []
    if isinstance(aqueous_species_str, (list, tuple)):
        names = aqueous_species_str
    else:
        names = aqueous_species_str.split()
    for name in names:
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

    additional_minerals = mineral_config.get("additional_minerals", [])
    additional_mineral_phases = []
    for name in additional_minerals:
        try:
            combined_db.addSpecies(supcrt_db.species(name))
            additional_mineral_phases.append(MineralPhase(name))
        except Exception as e:
            print(f"    WARNING: Could not add mineral '{name}': {e}")

    # Build aqueous phase species list
    base_species = "WATER,AQ H+ OH-"
    solute = mineral_config["solute_species"]
    additional = mineral_config.get("aqueous_species", "")

    include_elements = mineral_config.get("include_elements")
    if include_elements:
        aqueous_species_list = aqueous_species_by_elements(dew_db, include_elements)
        print(
            f"    Included {len(aqueous_species_list)} aqueous species after filtering."
        )
        validate_aqueous_species(dew_db, aqueous_species_list)
        aqueous = AqueousPhase(" ".join(aqueous_species_list))
    else:
        if additional:
            aqueous_species_str = f"{base_species} {solute} {additional}"
        else:
            aqueous_species_str = f"{base_species} {solute}"

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
        print(f" DEW configured: EOS={eos_name}; phase activity=ActivityModelDEW")

    except Exception as e:
        print(f"Warning: Could not configure DEW: {e}")
        print("  Falling back to default ActivityModelDEW()")
        try:
            aqueous.setActivityModel(ActivityModelDEW())
        except NameError:
            aqueous.setActivityModel(ActivityModelHKF())

    mineral = MineralPhase(mineral_name)
    system = ChemicalSystem(combined_db, aqueous, mineral, *additional_mineral_phases)

    print(f" System built for {mineral_name} solubility")
    return system


def load_experimental_data(csv_file):
    """Load and organize experimental data from CSV."""
    df = pd.read_csv(csv_file, encoding="cp1252")
    if {"T_C", "P_kbar", "molality_m"}.issubset(df.columns):
        df = df[["T_C", "P_kbar", "molality_m", "reference", "experiment_type"]].copy()
        df["reversal"] = df.get("reversal", "?")
    else:
        t_col = "T (Â°C)" if "T (Â°C)" in df.columns else "T (ï¿½C)"
        p_col = "P (bar)"
        m_col = "Molality (mol/kg Hâ‚‚O)"
        if m_col not in df.columns:
            m_col = "Molality (mol/kg H?O)"

        df = df[
            [t_col, p_col, m_col, "Reversal", "reference", "experiment_type"]
        ].copy()
        df = df.rename(
            columns={
                t_col: "T_C",
                p_col: "P_bar",
                m_col: "molality_m",
                "Reversal": "reversal",
            }
        )
        df["P_kbar"] = df["P_bar"] / 1000.0

    df["molality_m"] = pd.to_numeric(df["molality_m"], errors="coerce")
    df["P_bar"] = df["P_kbar"] * 1000.0
    # Keep NaN pressures for Hemley and other saturation curve experiments
    df["experiment_id"] = df["reference"] + " (" + df["experiment_type"] + ")"

    df["reversal"] = df["reversal"].astype(str).str.strip().str.upper()
    df["reversal"] = df["reversal"].where(df["reversal"].isin(["U", "S"]), "?")

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
    print("    WARNING: CaO_aq and MgO_aq excluded (not real aqueous species).")
    if not os.path.exists(CSV_FILE):
        print(f"    WARNING: Experimental data file not found: {CSV_FILE}")
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

    if len(exp_data) > 0:
        experiments = exp_data["experiment_id"].unique()
        print(f"    Experiments: {len(experiments)}")
    else:
        experiments = []
        pressures_kbar = DEFAULT_PRESSURES

    # Initialize databases
    print("\n[2] Initializing Reaktoro databases...")
    try:
        dew_db = DEWDatabase("dew2024-aqueous")
        print("    Successfully loaded DEW2024 aqueous database")
    except Exception as e:
        print(f"    ERROR: Failed to load DEW database: {e}")
        raise

    try:
        supcrt_db = SupcrtDatabase("supcrtbl")
        print("    Successfully loaded SUPCRTBL mineral database")
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
    element_curves = {"Ca": {}, "Mg": {}, "Si": {}}

    pressures_for_curves = sorted(exp_data["P_kbar"].dropna().unique())

    for P_kbar in pressures_for_curves:
        P_bar = P_kbar * 1000.0
        print(f"    P = {P_kbar:.3f} kbar ({P_bar:.0f} bar)...")

        # Determine T range from experiments at this pressure (Â±5%)
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
        state.set("WATER,AQ", 1.0, "kg")
        state.set("H+", 1e-8, "mol")
        state.set("OH-", 1e-8, "mol")
        state.set(solute_species, 1e-6, "mol")
        state.set(mineral_name, 10.0, "mol")
        for extra_mineral in MINERAL_CONFIG.get("additional_minerals", []):
            state.set(extra_mineral, 10.0, "mol")

        molalities = []
        element_m = {"Ca": [], "Mg": [], "Si": []}
        for T_C in T_range:
            conditions.temperature(float(T_C), "celsius")
            conditions.pressure(float(P_bar), "bar")
            result = solver.solve(state, conditions)

            if result.succeeded():
                try:
                    aqprops = AqueousProps(state)
                    m_ca, m_mg, m_si_half = diopside_element_totals(aqprops)
                    molality = np.nan
                    for elem in ("Ca", "Mg", "Si"):
                        try:
                            if elem == "Ca":
                                element_m[elem].append(m_ca)
                            elif elem == "Mg":
                                element_m[elem].append(m_mg)
                            else:
                                element_m[elem].append(m_si_half)
                        except Exception:
                            element_m[elem].append(np.nan)
                except Exception:
                    molality = np.nan
                    for elem in ("Ca", "Mg", "Si"):
                        element_m[elem].append(np.nan)
            else:
                molality = np.nan
                for elem in ("Ca", "Mg", "Si"):
                    element_m[elem].append(np.nan)

            molalities.append(molality)

        solubility_curves[P_kbar] = {"T_C": T_range, "molality": np.array(molalities)}
        for elem in ("Ca", "Mg", "Si"):
            element_curves[elem][P_kbar] = {
                "T_C": T_range,
                "molality": np.array(element_m[elem]),
            }

        valid_points = np.sum(~np.isnan(molalities))
        if valid_points > 0:
            valid_idx = np.where(~np.isnan(molalities))[0]
            first_T = T_range[valid_idx[0]]
            last_T = T_range[valid_idx[-1]]
            first_m = molalities[valid_idx[0]]
            last_m = molalities[valid_idx[-1]]
            print(
                f"       Calculated {valid_points}/{N_POINTS} points (T: {first_T:.0f}-{last_T:.0f}Â°C)"
            )

    # Plotting
    print("\n[4] Creating plots...")

    reversal_markers = (
        build_reversal_markers(exp_data["reversal"]) if len(exp_data) > 0 else {}
    )

    # =========================================================================
    # PLOT 2: High Pressure (>=1 kbar)
    # =========================================================================
    print("    Creating high-pressure plot (>=1 kbar)...")
    fig2, ax2 = plt.subplots(figsize=(14, 8))

    # Collect all experiments >= 1.0 kbar
    high_P_all_data = exp_data[exp_data["P_kbar"] >= 1.0].copy()

    # Get all unique pressures for these high-P experiments
    high_P_all_pressures = sorted(high_P_all_data["P_kbar"].unique())

    # Generate colors for high-pressure experiments
    n_high = max(len(high_P_all_pressures), 1)
    colors_high = plt.cm.plasma(np.linspace(0, 0.9, n_high))
    P_to_color_high = {
        P: colors_high[i % len(colors_high)] for i, P in enumerate(high_P_all_pressures)
    }

    # Plot high-pressure experimental data by actual pressure
    exp_labels_used = set()
    for P_kbar in high_P_all_pressures:
        P_tol = 0.05 * P_kbar
        subset = high_P_all_data[
            (high_P_all_data["P_kbar"] >= P_kbar - P_tol)
            & (high_P_all_data["P_kbar"] <= P_kbar + P_tol)
        ]
        if len(subset) == 0:
            continue

        # Plot by author within this pressure group
        for reversal in subset["reversal"].unique():
            rev_subset = subset[subset["reversal"] == reversal]
            marker = reversal_markers.get(reversal, "^")
            author = str(rev_subset["reference"].iloc[0])
            exp_label = f"{author} Exp P={P_kbar:.2f} kbar ({reversal})"
            if exp_label in exp_labels_used:
                exp_label = "_nolegend_"
            else:
                exp_labels_used.add(exp_label)
            ax2.scatter(
                rev_subset["T_C"],
                rev_subset["molality_m"],
                c=[P_to_color_high[P_kbar]],
                marker=marker,
                s=70,
                alpha=0.7,
                edgecolors="black",
                linewidths=0.4,
                label=exp_label,
                zorder=10,
            )

    # Plot calculated curves for high pressures (Ca, Mg, 0.5*Si)
    y_positive = []
    for P_kbar in high_P_all_pressures:
        for elem, linestyle, suffix in (
            ("Ca", "-", "Ca"),
            ("Mg", "--", "Mg"),
            ("Si", ":", "0.5Si"),
        ):
            if P_kbar not in element_curves[elem]:
                continue
            curve = element_curves[elem][P_kbar]
            molality = np.array(curve["molality"], dtype=float)
            valid = np.isfinite(molality) & (molality > 0)
            if np.any(valid):
                ax2.plot(
                    curve["T_C"][valid],
                    molality[valid],
                    color=P_to_color_high[P_kbar],
                    linewidth=2.5,
                    linestyle=linestyle,
                    alpha=0.9,
                    label=f"Calc {suffix} P={P_kbar:.2f} kbar",
                    zorder=15,
                )
                y_positive.extend(molality[valid].tolist())

    exp_positive = high_P_all_data.loc[
        high_P_all_data["molality_m"] > 0, "molality_m"
    ].to_numpy()
    if len(exp_positive) > 0:
        y_positive.extend(exp_positive.tolist())

    ax2.set_yscale("log")
    if y_positive:
        y_min = min(y_positive)
        y_max = max(y_positive)
        if y_min > 0 and y_max > 0:
            ax2.set_ylim(y_min * 0.7, y_max * 1.3)
    ax2.set_xlabel("Temperature (Â°C)", fontsize=14, fontweight="bold")
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
        f"DEW24 (species) + SUPCRTBL ({mineral_name}) + Zhang-Duan 2005 EOS (Hâ‚‚O)"
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
    # PLOT 2B: Element totals (Ca, Mg, Si) to assess incongruent dissolution
    # =========================================================================
    print("    Creating element total plots (Ca, Mg, Si)...")
    fig_el, axes_el = plt.subplots(3, 1, figsize=(12, 14), sharex=True)
    element_labels = {
        "Ca": "m_Ca_tot",
        "Mg": "m_Mg_tot",
        "Si": "0.5 m_Si_tot",
    }
    for ax_el, elem in zip(axes_el, ("Ca", "Mg", "Si")):
        for P_kbar in high_P_all_pressures:
            if P_kbar not in element_curves[elem]:
                continue
            curve = element_curves[elem][P_kbar]
            molality = np.array(curve["molality"], dtype=float)
            valid = np.isfinite(molality) & (molality > 0)
            if np.any(valid):
                ax_el.plot(
                    curve["T_C"][valid],
                    molality[valid],
                    color=P_to_color_high.get(P_kbar, "black"),
                    linewidth=2.5,
                    label=f"P={P_kbar:.2f} kbar",
                )

        ax_el.set_yscale("log")
        ax_el.set_ylabel(
            f"{element_labels[elem]} (mol/kg-Hâ‚‚O)", fontsize=12, fontweight="bold"
        )
        ax_el.grid(True, which="both", alpha=0.3, linestyle="--")

    axes_el[-1].set_xlabel("Temperature (Â°C)", fontsize=12, fontweight="bold")
    axes_el[0].set_title(
        "Diopside element totals: Ca, Mg, and 0.5Â·Si",
        fontsize=13,
        fontweight="bold",
        pad=12,
    )
    axes_el[0].legend(loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=9)
    plt.tight_layout()
    element_plot = os.path.join(
        SCRIPT_DIR, f"{MINERAL_CONFIG['output_prefix']}_element_totals_dew24.png"
    )
    plt.savefig(element_plot, dpi=300, bbox_inches="tight")
    print(f"    Element totals plot saved to: {element_plot}")
    plt.close(fig_el)

    # =========================================================================
    # PLOT 3: Aqueous speciation vs temperature (representative pressure)
    # =========================================================================
    if len(pressures_for_curves) > 0:
        spec_P_kbar = pressures_for_curves[len(pressures_for_curves) // 2]
        if spec_P_kbar in solubility_curves:
            print(f"    Creating speciation plot at P={spec_P_kbar:.3f} kbar...")
            curve = solubility_curves[spec_P_kbar]
            T_spec = curve["T_C"]
            spec_names = get_solute_species_list(MINERAL_CONFIG)
            spec_data = {name: [] for name in spec_names}

            specs_spec = EquilibriumSpecs(system)
            specs_spec.temperature()
            specs_spec.pressure()
            solver_spec = EquilibriumSolver(specs_spec)
            conditions_spec = EquilibriumConditions(specs_spec)

            state_spec = ChemicalState(system)
            state_spec.set("WATER,AQ", 1.0, "kg")
            state_spec.set("H+", 1e-8, "mol")
            state_spec.set("OH-", 1e-8, "mol")
            state_spec.set(solute_species, 1e-6, "mol")
            state_spec.set(mineral_name, 10.0, "mol")
            for extra_mineral in MINERAL_CONFIG.get("additional_minerals", []):
                state_spec.set(extra_mineral, 10.0, "mol")

            for T_C in T_spec:
                conditions_spec.temperature(float(T_C), "celsius")
                conditions_spec.pressure(float(spec_P_kbar * 1000.0), "bar")
                result = solver_spec.solve(state_spec, conditions_spec)
                if result.succeeded():
                    aqprops = AqueousProps(state_spec)
                    for name in spec_names:
                        try:
                            spec_data[name].append(float(aqprops.speciesMolality(name)))
                        except Exception:
                            spec_data[name].append(np.nan)
                else:
                    for name in spec_names:
                        spec_data[name].append(np.nan)

            fig3, ax3 = plt.subplots(figsize=(12, 7))
            for name, values in spec_data.items():
                values = np.array(values)
                valid = ~np.isnan(values)
                if np.any(valid):
                    ax3.plot(T_spec[valid], values[valid], label=name, linewidth=2.0)

            ax3.set_yscale("log")
            ax3.set_xlabel("Temperature (Â°C)", fontsize=12, fontweight="bold")
            ax3.set_ylabel(
                "Species molality (mol/kg-Hâ‚‚O)", fontsize=12, fontweight="bold"
            )
            ax3.set_title(
                f"{MINERAL_CONFIG['plot_title']} Speciation at {spec_P_kbar:.3f} kbar",
                fontsize=14,
                fontweight="bold",
                pad=16,
            )
            ax3.grid(True, which="both", alpha=0.3, linestyle="--")
            ax3.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=9)
            plt.tight_layout()
            plt.savefig(SPECIES_PLOT, dpi=300, bbox_inches="tight")
            print(f"    Speciation plot saved to: {SPECIES_PLOT}")
            plt.close(fig3)

    # =========================================================================
    # PLOT 4: Absolute and Relative Differences (measured - predicted)
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
                m_ca, m_mg, m_si_half = diopside_element_totals(aqprops)
                molality = min(m_ca, m_mg, m_si_half)
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
        reversal = row.get("reversal", "?")
        marker = reversal_markers.get(reversal, "^")

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
    ax_rel.set_xlabel("Temperature (Â°C)", fontsize=12, fontweight="bold")
    ax_rel.grid(True, alpha=0.3, linestyle="--")

    plt.tight_layout()
    plt.savefig(RESIDUALS_PLOT, dpi=300, bbox_inches="tight")
    print(f"    Residual plots saved to: {RESIDUALS_PLOT}")

    print("\n" + "=" * 80)
    print("Analysis complete!")
    print("=" * 80)


if __name__ == "__main__":
    main()

