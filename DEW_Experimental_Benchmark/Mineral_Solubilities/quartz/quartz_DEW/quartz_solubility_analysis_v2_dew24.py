"""
Mineral Solubility Analysis using Reaktoro with DEW2024
Generic framework for comparing calculated solubilities with experimental data
Includes per-kbar-category temperature ranges, Psat curve, and uncertainty analysis

Easily adaptable for different minerals by changing MINERAL_CONFIG section
"""

import os
import sys
import argparse
import importlib

# Keep DLL resolution focused on the active conda env on Windows.
# This avoids accidental runtime conflicts with other env/toolchain paths.
if os.name == "nt":
    env_prefix = sys.prefix
    env_paths = [
        env_prefix,
        os.path.join(env_prefix, "Library", "mingw-w64", "bin"),
        os.path.join(env_prefix, "Library", "usr", "bin"),
        os.path.join(env_prefix, "Library", "bin"),
        os.path.join(env_prefix, "Scripts"),
        os.path.join(env_prefix, "bin"),
    ]
    system_root = os.environ.get("SystemRoot", r"C:\Windows")
    system_paths = [
        os.path.join(system_root, "System32"),
        system_root,
        os.path.join(system_root, "System32", "Wbem"),
    ]
    os.environ["PATH"] = ";".join(
        [p for p in env_paths + system_paths if os.path.isdir(p)]
    )

try:
    import autodiff
except ModuleNotFoundError:

    class _AutoDiffShim:
        @staticmethod
        def real(value):
            return value

    autodiff = _AutoDiffShim()
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Try to import Reaktoro; fall back to local extension module if not installed
try:
    from reaktoro import *  # noqa: F401,F403
except ModuleNotFoundError:
    # Add local build path for reaktoro4py if available
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    ROOT_DIR = os.path.dirname(SCRIPT_DIR)
    pyd_candidates = [
        os.path.join(ROOT_DIR, "build", "Reaktoro", "Release"),
        os.path.join(ROOT_DIR, "build", "Reaktoro", "Release"),
    ]
    loaded_from = None
    for pyd_dir in pyd_candidates:
        if not os.path.isdir(pyd_dir):
            continue

        if pyd_dir in sys.path:
            sys.path.remove(pyd_dir)
        sys.path.insert(0, pyd_dir)
        sys.modules.pop("reaktoro4py", None)

        try:
            local_mod = importlib.import_module("reaktoro4py")
        except ModuleNotFoundError:
            continue

        # Mirror wildcard import behavior while letting us control load location.
        globals().update(
            {
                name: getattr(local_mod, name)
                for name in dir(local_mod)
                if not name.startswith("_")
            }
        )
        loaded_from = os.path.dirname(getattr(local_mod, "__file__", pyd_dir))
        break

    if loaded_from is None:
        raise ModuleNotFoundError(
            "Could not import Reaktoro. Install the 'reaktoro' package or ensure reaktoro4py is on PYTHONPATH."
        )

    print(f"Using local reaktoro4py extension from {loaded_from}.")

# Silence repeated non-convergence warnings
try:
    Warnings.disable(906)
except Exception:
    pass

# =============================================================================
# Configuration
# =============================================================================

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# ============================================================
# MINERAL CONFIGURATION - Change these for different minerals
# ============================================================
MINERAL_CONFIG = {
    # Mineral identification
    "mineral_name": "Quartz",  # Name in database
    "mineral_formula": "SiO2",  # Chemical formula
    "target_element": "Si",  # Element to report total dissolved molality
    "solute_species": "SiO2(aq)",  # Primary aqueous species
    # Aqueous species to include (besides water, H+, OH-)
    "aqueous_species": "H2(aq) O2(aq) HO2-(aq) HSiO3-(aq) Si2O4(aq) Si3O6(aq)",
    # File paths
    "csv_file": "quartz_DEW_testset.csv",
    "output_prefix": "quartz",  # Prefix for output files
    # Plot settings
    "plot_title": "Quartz Solubility",
    "y_label": "Quartz Solubility (mol/kg-Hâ‚‚O)",
}
# ============================================================

CSV_FILE = os.path.join(SCRIPT_DIR, MINERAL_CONFIG["csv_file"])

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

# Switch between the original DEW backend and Perple_X-linked DEW-style backend.
# Accepted values: "DEW", "PerplexDEW"
MODEL_BACKEND = "PerplexDEW"
PERPLEXDEW_REQUIRED_SYMBOLS = (
    "ActivityModelPerplexDEW",
    "ActivityDHModel",
    "StandardThermoModelParamsPerplexDEW",
    "StandardThermoModelPerplexDEW",
)


def backend_tag(name):
    """Normalize backend name for file names and labels."""
    return str(name or "DEW").strip()


def output_paths(
    mineral_config,
    backend_name,
    dh_model=None,
    perplex_activity_model="PerplexDEW",
):
    """Return output file paths scoped by backend to avoid overwrite."""
    tag = backend_tag(backend_name)
    if tag.lower() == "perplexdew" and dh_model:
        tag = f"{tag}_{dh_model}"
    if (
        tag.lower().startswith("perplexdew")
        and str(perplex_activity_model or "PerplexDEW").strip().lower() == "dew"
    ):
        tag = f"{tag}_DEWActivity"
    prefix = mineral_config["output_prefix"]
    return {
        "low_plot": os.path.join(
            SCRIPT_DIR,
            f"{prefix}_solubility_comparison_low_P_dew24_{tag}.png",
        ),
        "high_plot": os.path.join(
            SCRIPT_DIR,
            f"{prefix}_solubility_comparison_high_P_dew24_{tag}.png",
        ),
        "residuals_plot": os.path.join(
            SCRIPT_DIR,
            f"{prefix}_solubility_residuals_dew24_{tag}.png",
        ),
        "residuals_csv": os.path.join(
            SCRIPT_DIR,
            f"{prefix}_residuals_dew24_{tag}.csv",
        ),
        "curves_csv": os.path.join(
            SCRIPT_DIR,
            f"{prefix}_curves_dew24_{tag}.csv",
        ),
    }


def parse_args():
    parser = argparse.ArgumentParser(
        description="Quartz solubility benchmark with selectable backend."
    )
    parser.add_argument(
        "--backend",
        default=MODEL_BACKEND,
        choices=["DEW", "PerplexDEW"],
        help="Aqueous backend to use.",
    )
    parser.add_argument(
        "--dh-model",
        default="Davies",
        choices=["Davies", "ExtendedDH"],
        help="Debye-HÃ¼ckel variant for PerplexDEW backend (Davies or ExtendedDH).",
    )
    parser.add_argument(
        "--perplex-activity-model",
        default="PerplexDEW",
        choices=["PerplexDEW", "DEW"],
        help=(
            "Activity model path for PerplexDEW backend: native PerplexDEW or "
            "DEW-compatible activity handling."
        ),
    )
    return parser.parse_args()


def _local_pyd_candidates():
    root_dir = os.path.dirname(SCRIPT_DIR)
    return [
        os.path.join(root_dir, "build", "Reaktoro", "Release"),
        os.path.join(root_dir, "build", "Reaktoro", "Release"),
        os.path.join(root_dir, "build", "Reaktoro", "Release"),
    ]


def ensure_perplexdew_symbols():
    """Ensure PerplexDEW symbols are loaded, trying local reaktoro4py builds if needed."""
    missing = [name for name in PERPLEXDEW_REQUIRED_SYMBOLS if name not in globals()]
    if not missing:
        return

    searched_paths = []
    for pyd_dir in _local_pyd_candidates():
        searched_paths.append(pyd_dir)
        if not os.path.isdir(pyd_dir):
            continue

        if pyd_dir in sys.path:
            sys.path.remove(pyd_dir)
        sys.path.insert(0, pyd_dir)

        # Force import from this specific candidate path instead of reusing
        # a previously imported reaktoro4py from another build folder.
        sys.modules.pop("reaktoro4py", None)
        importlib.invalidate_caches()

        try:
            local_mod = importlib.import_module("reaktoro4py")
        except ModuleNotFoundError:
            continue

        for name in PERPLEXDEW_REQUIRED_SYMBOLS:
            if hasattr(local_mod, name):
                globals()[name] = getattr(local_mod, name)

        missing = [
            name for name in PERPLEXDEW_REQUIRED_SYMBOLS if name not in globals()
        ]
        if not missing:
            print(f"Using local reaktoro4py extension from {pyd_dir}.")
            return

    missing_str = ", ".join(missing)
    searched_str = "\n  - ".join(searched_paths)
    raise RuntimeError(
        "PerplexDEW backend requested but required symbols are unavailable: "
        f"{missing_str}.\nSearched local reaktoro4py build folders:\n  - {searched_str}"
    )


def to_real(value):
    """Convert Python scalars to autodiff real when required by pybind signatures."""
    try:
        return autodiff.real(value)
    except Exception:
        return value


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


def build_system(
    dew_db,
    supcrt_db,
    mineral_config,
    water_config=None,
    model_backend="DEW",
    dh_model="Davies",
    perplex_activity_model="PerplexDEW",
):
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

    if additional:
        aqueous_species_str = f"{base_species} {solute} {additional}"
    else:
        aqueous_species_str = f"{base_species} {solute}"

    validate_aqueous_species(dew_db, aqueous_species_str)
    aqueous = AqueousPhase(aqueous_species_str)

    backend = str(model_backend or "DEW").strip().lower()

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

        eos_name = water_config.get("eos_model", "ZhangDuan2005")

        if backend == "perplexdew":
            ensure_perplexdew_symbols()

            _dh = (
                ActivityDHModel.ExtendedDH
                if dh_model == "ExtendedDH"
                else ActivityDHModel.Davies
            )
            activity_mode = str(perplex_activity_model or "PerplexDEW").strip().lower()
            if activity_mode == "dew":
                aqueous.setActivityModel(ActivityModelDEW())
                print(
                    f"âœ“ PerplexDEW configured: EOS={eos_name}; DH model={dh_model}; phase activity=ActivityModelDEW (compatibility mode)"
                )
            else:
                aqueous.setActivityModel(ActivityModelPerplexDEW(_dh))
                print(
                    f"âœ“ PerplexDEW configured: EOS={eos_name}; DH model={dh_model}; phase activity=ActivityModelPerplexDEW"
                )
        else:
            _ = StandardThermoModelDEW(params)
            aqueous.setActivityModel(ActivityModelDEW())
            print(
                f"âœ“ DEW configured: EOS={eos_name}; phase activity=ActivityModelDEW"
            )

    except Exception as e:
        if backend == "perplexdew":
            raise RuntimeError(
                f"Could not configure requested aqueous backend ({model_backend}): {e}"
            ) from e

        print(
            f"Warning: Could not configure requested aqueous backend ({model_backend}): {e}"
        )
        print("  Falling back to DEW/HKF default activity model")
        try:
            aqueous.setActivityModel(ActivityModelDEW())
        except NameError:
            aqueous.setActivityModel(ActivityModelHKF())

    mineral = MineralPhase(mineral_name)
    system = ChemicalSystem(combined_db, aqueous, mineral)

    print(f"âœ“ System built for {mineral_name} solubility (backend={model_backend})")
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
    args = parse_args()
    model_backend = backend_tag(args.backend)
    dh_model = args.dh_model  # "Davies" or "ExtendedDH"
    perplex_activity_model = args.perplex_activity_model
    outputs = output_paths(
        MINERAL_CONFIG,
        args.backend,
        dh_model,
        perplex_activity_model=perplex_activity_model,
    )

    print("=" * 80)
    mineral_name = MINERAL_CONFIG["mineral_name"]
    print(f"{mineral_name} Solubility Analysis - Reaktoro DEW2024")
    print(f"Mineral: {MINERAL_CONFIG['mineral_formula']}")
    print(f"Solute species: {MINERAL_CONFIG['solute_species']}")
    backend_display = model_backend + (
        f" ({dh_model}, activity={perplex_activity_model})"
        if "perplexdew" in model_backend.lower()
        else ""
    )
    print(f"Aqueous backend: {backend_display}")
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
                f"    Temperature range: {exp_data['T_C'].min():.0f} - {exp_data['T_C'].max():.0f} Â°C"
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
        supcrt_db = SupcrtDatabase("supcrtbl")
        print("    Successfully loaded SUPCRTBL mineral database")
    except Exception as e:
        print(f"    ERROR: Failed to load SUPCRTBL database: {e}")
        raise

    system = build_system(
        dew_db,
        supcrt_db,
        MINERAL_CONFIG,
        water_config=DEW_CONFIG,
        model_backend=model_backend,
        dh_model=dh_model,
        perplex_activity_model=perplex_activity_model,
    )

    # Calculate solubility curves for each experimental pressure (drops NaN)
    mineral_name = MINERAL_CONFIG["mineral_name"]
    solute_species = MINERAL_CONFIG["solute_species"]
    solute_species_list = get_solute_species_list(MINERAL_CONFIG)
    print(f"\n[3] Calculating {mineral_name.lower()} solubility curves...")
    solubility_curves = {}

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
        state.set("H2O(aq)", to_real(1.0), "kg")
        state.set("H+(aq)", to_real(1e-8), "mol")
        state.set("OH-(aq)", to_real(1e-8), "mol")
        state.set(solute_species, to_real(1e-6), "mol")
        state.set(mineral_name, to_real(10.0), "mol")

        molalities = []
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
            first_m = molalities[valid_idx[0]]
            last_m = molalities[valid_idx[-1]]
            print(
                f"       Calculated {valid_points}/{N_POINTS} points (T: {first_T:.0f}-{last_T:.0f}Â°C)"
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
    for i, T_C in enumerate(T_psat_range):
        if not valid_temps[i]:
            psat_molalities.append(np.nan)
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

    # Separate data into low (<1 kbar) and high (â‰¥1 kbar) pressure ranges
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
            marker = author_markers.get(author, "o")
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
    ax1.set_xlabel("Temperature (Â°C)", fontsize=14, fontweight="bold")
    ax1.set_ylabel(MINERAL_CONFIG["y_label"], fontsize=14, fontweight="bold")
    ax1.set_title(
        f"{MINERAL_CONFIG['plot_title']}: Low Pressure (<1 kbar)",
        fontsize=16,
        fontweight="bold",
        pad=20,
    )
    ax1.grid(True, which="both", alpha=0.3, linestyle="--")

    # Add database/model info annotation
    info_text = (
        f"DEW24 (species) + SUPCRTBL ({mineral_name}) + Zhang-Duan 2005 EOS (Hâ‚‚O)"
    )
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
    plt.savefig(outputs["low_plot"], dpi=300, bbox_inches="tight")
    print(f"    Low-P plot saved to: {outputs['low_plot']}")
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
    plt.savefig(outputs["high_plot"], dpi=300, bbox_inches="tight")
    print(f"    High-P plot saved to: {outputs['high_plot']}")
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
        state.set("H2O(aq)", to_real(1.0), "kg")
        state.set("H+(aq)", to_real(1e-8), "mol")
        state.set("OH-(aq)", to_real(1e-8), "mol")
        state.set(solute_species, to_real(1e-6), "mol")
        state.set(mineral_name, to_real(10.0), "mol")
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
    ax_rel.set_xlabel("Temperature (Â°C)", fontsize=12, fontweight="bold")
    ax_rel.grid(True, alpha=0.3, linestyle="--")

    plt.tight_layout()
    plt.savefig(outputs["residuals_plot"], dpi=300, bbox_inches="tight")
    print(f"    Residual plots saved to: {outputs['residuals_plot']}")

    exp_resid["backend"] = backend_display
    exp_resid.to_csv(outputs["residuals_csv"], index=False)
    print(f"    Residual data saved to: {outputs['residuals_csv']}")

    curve_rows = []
    for pressure_key, curve in solubility_curves.items():
        if pressure_key == "Psat":
            curve_type = "psat"
            pressures = curve["P_kbar"]
        else:
            curve_type = "isobar"
            pressures = np.full_like(curve["T_C"], float(pressure_key), dtype=float)

        for t_c, p_kbar, molality in zip(curve["T_C"], pressures, curve["molality"]):
            curve_rows.append(
                {
                    "backend": backend_display,
                    "curve_type": curve_type,
                    "P_kbar": float(p_kbar),
                    "T_C": float(t_c),
                    "molality": float(molality) if not np.isnan(molality) else np.nan,
                }
            )

    curves_df = pd.DataFrame(curve_rows)
    curves_df.to_csv(outputs["curves_csv"], index=False)
    print(f"    Curve data saved to: {outputs['curves_csv']}")

    print("\n" + "=" * 80)
    print("Analysis complete!")
    print("=" * 80)


if __name__ == "__main__":
    main()
