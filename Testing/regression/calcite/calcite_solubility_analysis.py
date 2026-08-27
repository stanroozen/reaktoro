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
import re

# Keep DLL resolution focused on the active conda env on Windows.
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
    import autodiff  # noqa: F401
except ModuleNotFoundError:
    autodiff = None
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BENCHMARK_DIR = os.path.dirname(SCRIPT_DIR)
ROOT_DIR = os.path.dirname(BENCHMARK_DIR)

# Try to import Reaktoro; fall back to local extension modules if not installed
try:
    from reaktoro import *  # noqa: F401,F403
except ModuleNotFoundError:
    # First-resort: try reaktoro4py directly from sys.path (covers PYTHONPATH-set builds).
    _loaded_from = None
    try:
        _local_mod = importlib.import_module("reaktoro4py")
        globals().update(
            {
                k: getattr(_local_mod, k)
                for k in dir(_local_mod)
                if not k.startswith("_")
            }
        )
        _loaded_from = os.path.dirname(getattr(_local_mod, "__file__", ""))
    except ModuleNotFoundError:
        pass

    if _loaded_from is None:
        # Fallback: search for reaktoro4py in known relative build locations.
        _pyd_candidates = [
            os.path.join(ROOT_DIR, "build", "Reaktoro", "Release"),
        ]
        for _pyd_dir in _pyd_candidates:
            if not os.path.isdir(_pyd_dir):
                continue
            if _pyd_dir in sys.path:
                sys.path.remove(_pyd_dir)
            sys.path.insert(0, _pyd_dir)
            sys.modules.pop("reaktoro4py", None)
            try:
                _local_mod = importlib.import_module("reaktoro4py")
            except ModuleNotFoundError:
                continue
            globals().update(
                {
                    k: getattr(_local_mod, k)
                    for k in dir(_local_mod)
                    if not k.startswith("_")
                }
            )
            _loaded_from = _pyd_dir
            break
    if _loaded_from is None:
        raise ModuleNotFoundError(
            "Could not import Reaktoro. Install the 'reaktoro' package or "
            "ensure reaktoro4py is on PYTHONPATH."
        )
    print(f"Using local reaktoro4py extension from {_loaded_from}.")

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
    "mineral_name": "Calcite",  # Name in database
    "mineral_formula": "CaCO3",  # Chemical formula
    "target_element": "Ca",  # Element to report total dissolved molality
    "solute_species": "Ca+2(aq)",  # Primary aqueous species
    "include_elements": ["H", "O", "Ca", "C"],
    "exclude_organics": True,
    # CaO(aq) is not a real aqueous species; CO(aq)/H2(aq)/O2(aq) are dissolved redox
    # gases irrelevant to carbonate chemistry and excluded to keep the system clean.
    "excluded_species": ["CaO(aq)", "CO(aq)", "H2(aq)", "O2(aq)"],
    # Full Ca-C-H-O species list (used when include_elements is None; otherwise the
    # element filter is used and returns: WATER,AQ H+ OH- CO2_aq HCO3- CO3-2
    # H2CO3_aq CaCO3_aq Ca(HCO3)+ Ca(OH)+ Ca+2)
    "aqueous_species": "CO2_aq HCO3- CO3-2 H2CO3_aq CaCO3_aq Ca(HCO3)+ Ca(OH)+",
    # File paths
    "csv_file": "calcite_DEW_testset.csv",
    "output_prefix": "calcite",  # Prefix for output files
    # Plot settings
    "plot_title": "Calcite Solubility",
    "y_label": "Calcite Solubility (mol/kg-H₂O)",
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
MODEL_BACKEND = "DEW"
PERPLEXDEW_REQUIRED_SYMBOLS = (
    "ActivityModelPerplexDEW",
    "ActivityDHModel",
    "StandardThermoModelParamsPerplexDEW",
    "StandardThermoModelPerplexDEW",
    "ActivityModelPerplexGFSM",
    "ActivityModelParamsPerplexGFSM",
    "PerpleXHybridEosOptions",
    "PerpleXCO2Eos",
)


def backend_tag(name):
    """Normalize backend name for file names and labels."""
    return str(name or "DEW").strip()


def output_paths(mineral_config, backend_name, dh_model=None, fluid=None):
    """Return output file paths scoped by backend and fluid composition."""
    tag = backend_tag(backend_name)
    if tag.lower() == "perplexdew" and dh_model:
        tag = f"{tag}_{dh_model}"
    fluid_tag = str(fluid or "H2O").upper().replace("-", "")
    if fluid_tag != "H2O":
        tag = f"{tag}_{fluid_tag}"
    prefix = mineral_config["output_prefix"]
    return {
        "high_plot": os.path.join(
            SCRIPT_DIR,
            f"{prefix}_solubility_comparison_high_P_dew24_{tag}.png",
        ),
        "residuals_plot": os.path.join(
            SCRIPT_DIR,
            f"{prefix}_solubility_residuals_dew24_{tag}.png",
        ),
        "species_plot": os.path.join(
            SCRIPT_DIR,
            f"{prefix}_speciation_dew24_{tag}.png",
        ),
    }


def parse_args():
    parser = argparse.ArgumentParser(
        description="Calcite solubility benchmark with selectable backend."
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
        help="Debye-Hückel variant for PerplexDEW backend (Davies or ExtendedDH).",
    )
    parser.add_argument(
        "--fluid",
        default="H2O",
        choices=["H2O", "H2OCO2", "H2OCO2aq"],
        help="Bulk fluid composition: pure water (H2O), water-CO2 gas-phase mixture (H2OCO2), or water with dissolved CO2_aq solute (H2OCO2aq).",
    )
    parser.add_argument(
        "--xco2",
        type=float,
        default=0.1,
        help="CO2 mole fraction in the bulk H2O-CO2 fluid (only used when --fluid H2OCO2).",
    )
    parser.add_argument(
        "--co2aq-molality",
        type=float,
        default=6.168,
        dest="co2aq_molality",
        help="Initial CO2_aq molality (mol/kg-H2O) added to the pure-H2O system (only used when --fluid H2OCO2aq).",
    )
    return parser.parse_args()


def _local_pyd_candidates():
    return [
        os.path.join(ROOT_DIR, "build", "Reaktoro", "Release"),
        os.path.join(ROOT_DIR, "build", "Reaktoro", "Release"),
        os.path.join(ROOT_DIR, "build", "Reaktoro", "Release"),
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
    excluded.update(
        [
            "Ca(CH3COO)+",
            "Ca(CH3COO)2_aq",
            "Ca(CH3COO)3-",
            "Ca(HCOO)+",
            "Isobutane_aq",
            "HO2-",
        ]
    )

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


def build_system(
    dew_db,
    supcrt_db,
    mineral_config,
    water_config=None,
    model_backend="DEW",
    dh_model="Davies",
    fluid="H2O",
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

    _fluid_tag = str(fluid or "H2O").upper().replace("-", "")

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
        aqueous_species_list = aqueous_species_by_elements(dew_db, include_elements)
        if _fluid_tag == "H2OCO2" and "CO2(aq)" in aqueous_species_list:
            aqueous_species_list = [s for s in aqueous_species_list if s != "CO2(aq)"]
            print("    (CO2(aq) excluded: CO2 treated as gas-phase solvent only)")
        print(
            f"    Included {len(aqueous_species_list)} aqueous species after filtering."
        )
        target_element = mineral_config.get("target_element")
        if target_element:
            mineral_config["_element_species_map"] = {
                target_element: element_species_coeffs(
                    dew_db, aqueous_species_list, target_element
                )
            }
        validate_aqueous_species(dew_db, aqueous_species_list)
        aqueous = AqueousPhase(" ".join(aqueous_species_list))
    else:
        if additional:
            aqueous_species_str = f"{base_species} {solute} {additional}"
        else:
            aqueous_species_str = f"{base_species} {solute}"

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

        eos_name = water_config.get("eos_model", "ZhangDuan2005")
        backend = str(model_backend or "DEW").strip().lower()
        if backend == "perplexdew":
            ensure_perplexdew_symbols()
            _dh = (
                ActivityDHModel.ExtendedDH
                if dh_model == "ExtendedDH"
                else ActivityDHModel.Davies
            )
            aqueous.setActivityModel(ActivityModelPerplexDEW(_dh))
            print(
                f"✓ PerplexDEW configured: EOS={eos_name}; DH model={dh_model}; phase activity=ActivityModelPerplexDEW"
            )
        else:
            _ = StandardThermoModelDEW(params)
            aqueous.setActivityModel(ActivityModelDEW())
            print(f"✓ DEW configured: EOS={eos_name}; phase activity=ActivityModelDEW")

    except Exception as e:
        if str(model_backend or "").strip().lower() == "perplexdew":
            raise RuntimeError(
                f"Could not configure requested aqueous backend ({model_backend}): {e}"
            ) from e
        print(f"Warning: Could not configure DEW: {e}")
        print("  Falling back to default ActivityModelDEW()")
        try:
            aqueous.setActivityModel(ActivityModelDEW())
        except NameError:
            aqueous.setActivityModel(ActivityModelHKF())

    mineral = MineralPhase(mineral_name)

    # Optionally add supercritical CO2 gas phase for H2O-CO2 mixed fluid
    if _fluid_tag == "H2OCO2":
        try:
            co2_gas_species = supcrt_db.species("CO2(g)")
            combined_db.addSpecies(co2_gas_species)
            co2_phase = GaseousPhase("CO2(g)")
            gfsm_params = ActivityModelParamsPerplexGFSM()
            hybrid_opts = PerpleXHybridEosOptions()
            hybrid_opts.co2 = PerpleXCO2Eos.ZhangDuan09
            gfsm_params.hybridEosOptions = hybrid_opts
            co2_phase.setActivityModel(ActivityModelPerplexGFSM(gfsm_params))
            system = ChemicalSystem(combined_db, aqueous, co2_phase, mineral)
            print("✓ CO2(g) gas phase added (Zhang-Duan 2009 EOS) for H2O-CO2 fluid")
        except Exception as e:
            raise RuntimeError(f"Failed to add CO2(g) gaseous phase: {e}") from e
    else:
        system = ChemicalSystem(combined_db, aqueous, mineral)

    print(f"✓ System built for {mineral_name} solubility")
    return system


def load_experimental_data(csv_file):
    """Load and organize experimental data from CSV."""
    df = pd.read_csv(csv_file, encoding="cp1252")
    if {"T_C", "P_kbar", "molality_m"}.issubset(df.columns):
        df = df[["T_C", "P_kbar", "molality_m", "reference", "experiment_type"]].copy()
        df["reversal"] = df.get("reversal", "?")
    else:
        t_col = "T (°C)" if "T (°C)" in df.columns else "T (�C)"
        p_col = "P (bar)"
        m_col = "Molality (mol/kg H₂O)"
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
    args = parse_args()
    dh_model = args.dh_model
    model_backend = backend_tag(args.backend)
    fluid = args.fluid
    xco2 = args.xco2
    outputs = output_paths(MINERAL_CONFIG, args.backend, dh_model, fluid)

    # Moles of CO2 to add per 1 kg H2O (n_co2 = 0 for pure H2O)
    n_co2 = xco2 / (1.0 - xco2) * 55.508 if fluid == "H2OCO2" else 0.0
    # Moles of dissolved CO2_aq to add as solute in pure H2O (H2OCO2aq mode)
    n_co2_aq = args.co2aq_molality if fluid == "H2OCO2aq" else 0.0

    print("=" * 80)
    mineral_name = MINERAL_CONFIG["mineral_name"]
    print(f"{mineral_name} Solubility Analysis - Reaktoro DEW2024")
    print(f"Mineral: {MINERAL_CONFIG['mineral_formula']}")
    print(f"Solute species: {MINERAL_CONFIG['solute_species']}")
    backend_display = model_backend + (
        f" ({dh_model})" if "perplexdew" in model_backend.lower() else ""
    )
    print(f"Aqueous backend: {backend_display}")
    if fluid == "H2OCO2":
        print(
            f"Fluid: H2O-CO2 mixture (XCO2 = {xco2:.2f}, n_CO2 = {n_co2:.3f} mol/kg-H2O)"
        )
    elif fluid == "H2OCO2aq":
        print(
            f"Fluid: Pure H2O solvent + CO2_aq solute (initial CO2_aq = {n_co2_aq:.3f} mol/kg-H2O)"
        )
    print("=" * 80)

    # Load experimental data
    print("\n[1] Loading experimental data...")
    print("    WARNING: CaO_aq excluded (not a real aqueous species).")
    if not os.path.exists(CSV_FILE):
        print(f"    WARNING: Experimental data file not found: {CSV_FILE}")
        exp_data = pd.DataFrame(
            columns=["T_C", "P_kbar", "molality_m", "reference", "experiment_id"]
        )
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

    system = build_system(
        dew_db,
        supcrt_db,
        MINERAL_CONFIG,
        model_backend=model_backend,
        dh_model=dh_model,
        fluid=fluid,
    )

    # Calculate solubility curves for each experimental pressure (drops NaN)
    mineral_name = MINERAL_CONFIG["mineral_name"]
    solute_species = MINERAL_CONFIG["solute_species"]
    solute_species_list = get_solute_species_list(MINERAL_CONFIG)
    print(f"\n[3] Calculating {mineral_name.lower()} solubility curves...")
    solubility_curves = {}

    pressures_for_curves = (
        sorted(exp_data["P_kbar"].dropna().unique())
        if len(exp_data) > 0
        else DEFAULT_PRESSURES
    )

    for P_kbar in pressures_for_curves:
        P_bar = P_kbar * 1000.0
        print(f"    P = {P_kbar:.3f} kbar ({P_bar:.0f} bar)...")

        # Determine T range from experiments at this pressure (�5%)
        P_tol = 0.05 * P_kbar
        if len(exp_data) > 0:
            exp_at_P = exp_data[
                (exp_data["P_kbar"] >= P_kbar - P_tol)
                & (exp_data["P_kbar"] <= P_kbar + P_tol)
            ]
        else:
            exp_at_P = pd.DataFrame()

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
        if n_co2 > 0.0:
            state.set("CO2(g)", n_co2, "mol")
        if n_co2_aq > 0.0:
            state.set("CO2(aq)", n_co2_aq, "mol")

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
                f"       Calculated {valid_points}/{N_POINTS} points (T: {first_T:.0f}-{last_T:.0f}°C)"
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

    # Plot calculated curves for high pressures
    y_positive = []
    for P_kbar in high_P_all_pressures:
        if P_kbar not in solubility_curves:
            continue
        curve = solubility_curves[P_kbar]
        molality = np.array(curve["molality"], dtype=float)
        valid = np.isfinite(molality) & (molality > 0)
        if np.any(valid):
            ax2.plot(
                curve["T_C"][valid],
                molality[valid],
                color=P_to_color_high[P_kbar],
                linewidth=3.0,
                linestyle="-",
                alpha=0.9,
                label=f"Calc P={P_kbar:.2f} kbar",
                zorder=15,
            )
            y_positive.extend(molality[valid].tolist())
        else:
            print(
                f"    WARNING: No positive molalities to plot at P={P_kbar:.2f} kbar."
            )

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
    plt.savefig(outputs["high_plot"], dpi=300, bbox_inches="tight")
    print(f"    High-P plot saved to: {outputs['high_plot']}")
    plt.close(fig2)

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
            state_spec.set("H2O(aq)", 1.0, "kg")
            state_spec.set("H+(aq)", 1e-8, "mol")
            state_spec.set("OH-(aq)", 1e-8, "mol")
            state_spec.set(solute_species, 1e-6, "mol")
            state_spec.set(mineral_name, 10.0, "mol")
            if n_co2 > 0.0:
                state_spec.set("CO2(g)", n_co2, "mol")

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
            ax3.set_xlabel("Temperature (°C)", fontsize=12, fontweight="bold")
            ax3.set_ylabel(
                "Species molality (mol/kg-H₂O)", fontsize=12, fontweight="bold"
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
            plt.savefig(outputs["species_plot"], dpi=300, bbox_inches="tight")
            print(f"    Speciation plot saved to: {outputs['species_plot']}")
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
        state.set("H2O(aq)", 1.0, "kg")
        state.set("H+(aq)", 1e-8, "mol")
        state.set("OH-(aq)", 1e-8, "mol")
        state.set(solute_species, 1e-6, "mol")
        state.set(mineral_name, 10.0, "mol")
        if n_co2 > 0.0:
            state.set("CO2(g)", n_co2, "mol")
        if n_co2_aq > 0.0:
            state.set("CO2(aq)", n_co2_aq, "mol")
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
    ax_rel.set_xlabel("Temperature (°C)", fontsize=12, fontweight="bold")
    ax_rel.grid(True, alpha=0.3, linestyle="--")

    plt.tight_layout()
    plt.savefig(outputs["residuals_plot"], dpi=300, bbox_inches="tight")
    print(f"    Residual plots saved to: {outputs['residuals_plot']}")

    print("\n" + "=" * 80)
    print("Analysis complete!")
    print("=" * 80)


if __name__ == "__main__":
    main()
