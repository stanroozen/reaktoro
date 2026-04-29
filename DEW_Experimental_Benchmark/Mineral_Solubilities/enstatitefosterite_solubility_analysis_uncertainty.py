"""
Mineral Solubility Analysis using Reaktoro with DEW2024
Generic framework for comparing calculated solubilities with experimental data
Includes per-kbar-category temperature ranges, Psat curve, and uncertainty analysis

Easily adaptable for different minerals by changing MINERAL_CONFIG section
"""

import os
import sys
import re
import json
import copy
import tempfile
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
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# Try to import Reaktoro; fall back to local extension module if not installed
try:
    from reaktoro import *  # noqa: F401,F403
except ModuleNotFoundError:
    # Add local build path for reaktoro4py if available
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    BENCHMARK_DIR = os.path.dirname(SCRIPT_DIR)
    ROOT_DIR = os.path.dirname(BENCHMARK_DIR)
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
    "mineral_name": "Enstatite",  # Name in database
    "mineral_formula": "MgSiO3",  # Chemical formula
    "target_element": "Si",  # Element to report total dissolved molality
    "solute_species": "SiO2_aq",  # Primary aqueous species
    "include_elements": ["H", "O", "Mg", "Si"],
    "exclude_organics": True,
    "excluded_species": ["MgO_aq"],
    "additional_minerals": ["Forsterite"],
    # Aqueous species to include (besides water, H+, OH-)
    "aqueous_species": "",
    # File paths
    "csv_file": "enstatite_DEW_testset.csv",
    "output_prefix": "enstatitefosterite",  # Prefix for output files
    # Plot settings
    "plot_title": "Enstatite+Forsterite Apparent Solubility",
    "y_label": "Apparent Enstatite+Forsterite Solubility (m_Si, mol/kg-Hâ‚‚O)",
}
# ============================================================

CSV_FILE = os.path.join(SCRIPT_DIR, MINERAL_CONFIG["csv_file"])
ROOT_DIR = os.path.dirname(os.path.dirname(SCRIPT_DIR))
DEFAULT_TC_DS62_COVARIANCE_JSON = os.path.join(
    ROOT_DIR,
    "embedded",
    "databases",
    "hollandpowell",
    "tc-ds62-covariance.json",
)
DEFAULT_TC_DS62_MINERAL_JSON = os.path.join(
    ROOT_DIR,
    "embedded",
    "databases",
    "hollandpowell",
    "tc-ds62-reaktoro.json",
)

# Calculation settings
T_MIN, T_MAX = 800, 900
N_POINTS = 100
DEFAULT_PRESSURES = [10.0]

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


def output_paths(mineral_config, backend_name, dh_model=None, uncertainty_enabled=True):
    """Return output file paths scoped by backend to avoid overwrite."""
    tag = backend_tag(backend_name)
    if tag.lower() == "perplexdew" and dh_model:
        tag = f"{tag}_{dh_model}"
    mode_tag = "uncertainty" if uncertainty_enabled else "deterministic"
    tag = f"{tag}_{mode_tag}"
    prefix = mineral_config["output_prefix"]
    return {
        "low_plot": os.path.join(
            SCRIPT_DIR,
            f"{prefix}_solubility_comparison_all_P_dew24_{tag}.png",
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
        "uncertainty_csv": os.path.join(
            SCRIPT_DIR,
            f"{prefix}_uncertainty_summary_dew24_{tag}.csv",
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
        "--covariance-json",
        default=DEFAULT_TC_DS62_COVARIANCE_JSON,
        help="Path to covariance JSON (packed upper triangle format).",
    )
    parser.add_argument(
        "--uncertainty-entity",
        default="en",
        help="Entity code in covariance JSON used for uncertainty propagation (default: en for Enstatite).",
    )
    parser.add_argument(
        "--uncertainty-samples",
        type=int,
        default=300,
        help="Number of Monte Carlo samples for uncertainty envelopes.",
    )
    parser.add_argument(
        "--uncertainty-ci",
        type=float,
        default=95.0,
        help="Central confidence interval in percent (e.g., 95 gives 2.5-97.5%).",
    )
    parser.add_argument(
        "--uncertainty-seed",
        type=int,
        default=42,
        help="Random seed for Monte Carlo sampling.",
    )
    parser.add_argument(
        "--disable-uncertainty",
        action="store_true",
        help="Disable uncertainty propagation and plot only deterministic curves.",
    )
    parser.add_argument(
        "--mineral-db-json",
        default=DEFAULT_TC_DS62_MINERAL_JSON,
        help="Path to mineral Reaktoro JSON used for uncertainty sampling (e.g., tc-ds62-reaktoro.json).",
    )
    parser.add_argument(
        "--mineral-species-code",
        default="en",
        help="Mineral species code in mineral JSON used in equilibrium system (default: en for Enstatite).",
    )
    parser.add_argument(
        "--quick-eval",
        action="store_true",
        help="Use only a small representative subset of low/high/Psat experimental points for quick runs.",
    )
    parser.add_argument(
        "--quick-n-per-group",
        type=int,
        default=2,
        help="Number of points per group (low-P, high-P, Psat) when --quick-eval is enabled.",
    )
    parser.add_argument(
        "--quick-npoints",
        type=int,
        default=25,
        help="Number of temperature points per curve when --quick-eval is enabled.",
    )
    return parser.parse_args()


def _sample_evenly_by_temperature(df, n):
    """Pick up to n points spread across temperature range."""
    if len(df) <= n:
        return df.copy()
    ordered = df.sort_values("T_C").reset_index(drop=True)
    idx = np.linspace(0, len(ordered) - 1, num=n, dtype=int)
    return ordered.iloc[np.unique(idx)].copy()


def select_quick_experimental_subset(exp_data, n_per_group=2):
    """Return a compact subset with representative low-P, high-P, and Psat points."""
    n = max(1, int(n_per_group))

    psat = exp_data[exp_data["is_psat"]]
    non_psat = exp_data[~exp_data["is_psat"] & exp_data["P_kbar"].notna()]
    low = non_psat[non_psat["P_kbar"] < 1.0]
    high = non_psat[non_psat["P_kbar"] >= 1.0]

    parts = [
        _sample_evenly_by_temperature(low, n),
        _sample_evenly_by_temperature(high, n),
        _sample_evenly_by_temperature(psat, n),
    ]
    subset = pd.concat(parts, ignore_index=True)
    if len(subset) == 0:
        return exp_data.copy()

    subset = (
        subset.drop_duplicates(subset=["T_C", "P_kbar", "reference", "experiment_type"])
        .sort_values(["is_psat", "P_kbar", "T_C"], na_position="last")
        .reset_index(drop=True)
    )
    return subset


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


def load_covariance_matrix(covariance_json_path):
    """Load covariance matrix from packed upper-triangle JSON export.

    NOTE: Source covariance from THERMOCALC/Perple_X is in (kJ/mol)Â².
    PackedUpperTriangle is stored as source / 1000, so we scale by 1000
    to recover the original (kJ/mol)Â² units.
    """
    with open(covariance_json_path, "r", encoding="utf-8") as f:
        data = json.load(f)

    entities = data["Entities"]
    packed = data["PackedUpperTriangle"]
    n = len(entities)
    expected = n * (n + 1) // 2
    if len(packed) != expected:
        raise ValueError(
            f"Invalid packed covariance length ({len(packed)}). Expected {expected} for n={n}."
        )

    sigma = np.zeros((n, n), dtype=float)
    k = 0
    for i in range(n):
        for j in range(i, n):
            v = float(packed[k]) * 1000.0  # Scale from (J/mol)Â² to (kJ/mol)Â²
            sigma[i, j] = v
            sigma[j, i] = v
            k += 1

    return entities, sigma


def sample_from_covariance_cholesky(sigma, nsamples, rng):
    """Draw multivariate normal samples with numerically robust Cholesky factorization."""
    n = sigma.shape[0]
    eye = np.eye(n)
    jitter = 0.0
    last_err = None
    for _ in range(8):
        try:
            L = np.linalg.cholesky(sigma + jitter * eye)
            z = rng.standard_normal(size=(nsamples, n))
            return z @ L.T
        except np.linalg.LinAlgError as err:
            last_err = err
            jitter = 1e-20 if jitter == 0.0 else jitter * 10.0

    raise np.linalg.LinAlgError(
        f"Could not Cholesky-factor covariance matrix after jitter regularization: {last_err}"
    )


def covariance_diagnostics(sigma, rel_tol=1e-12, abs_tol=1e-20):
    """Return PSD/eigenvalue diagnostics for covariance quality assessment."""
    sym = 0.5 * (sigma + sigma.T)
    eigvals = np.linalg.eigvalsh(sym)
    maxeig = float(np.max(eigvals)) if eigvals.size else np.nan
    mineig = float(np.min(eigvals)) if eigvals.size else np.nan
    tol_psd = max(abs_tol, rel_tol * max(maxeig, 1.0))
    is_psd = bool(mineig >= -tol_psd)

    positive = eigvals[eigvals > max(abs_tol, rel_tol * max(maxeig, 1.0))]
    effective_rank = int(len(positive))
    cond_number = float(np.inf)
    if len(positive) > 0:
        cond_number = float(np.max(positive) / np.min(positive))

    return {
        "is_psd": is_psd,
        "min_eigenvalue": mineig,
        "max_eigenvalue": maxeig,
        "effective_rank": effective_rank,
        "condition_number": cond_number,
    }


def relevant_mineral_entities(cov_entities, mineral_json_data):
    """Select covariance entities that map to mineral species with Holland-Powell Gf entries."""
    species_map = mineral_json_data.get("Species", {})
    relevant = []
    for code in cov_entities:
        sp = species_map.get(code)
        if not isinstance(sp, dict):
            continue
        hp = sp.get("StandardThermoModel", {}).get("HollandPowell", {})
        if isinstance(hp, dict) and "Gf" in hp:
            relevant.append(code)
    return relevant


def sampled_mineral_db_file(
    base_mineral_json_data, shifts_by_entity, tmp_dir, sample_id
):
    """Create a sampled mineral DB JSON file by perturbing all relevant species Gf values."""
    sampled = copy.deepcopy(base_mineral_json_data)
    species_map = sampled.get("Species", {})
    for code, shift in shifts_by_entity.items():
        sp = species_map.get(code)
        if not isinstance(sp, dict):
            continue
        hp = sp.get("StandardThermoModel", {}).get("HollandPowell", {})
        if not isinstance(hp, dict) or "Gf" not in hp:
            continue
        hp["Gf"] = float(hp["Gf"]) + float(shift)

    path = os.path.join(tmp_dir, f"sampled_mineral_db_{sample_id:05d}.json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(sampled, f)
    return path


def evaluate_curves_on_fixed_grids(
    system,
    pressure_grids,
    psat_grid,
    mineral_species_name,
    solute_species,
    solute_species_list,
):
    """Forward-evaluate molality curves on fixed T/P grids for one sampled system."""
    curves = {}

    for p_kbar, t_range in pressure_grids.items():
        specs = EquilibriumSpecs(system)
        specs.temperature()
        specs.pressure()
        solver = EquilibriumSolver(specs)
        conds = EquilibriumConditions(specs)

        state = ChemicalState(system)
        state.set("WATER,AQ", to_real(1.0), "kg")
        state.set("H+", to_real(1e-8), "mol")
        state.set("OH-", to_real(1e-8), "mol")
        state.set(solute_species, to_real(1e-6), "mol")
        state.set(mineral_species_name, to_real(10.0), "mol")

        vals = []
        p_bar = float(p_kbar) * 1000.0
        for t_c in t_range:
            conds.temperature(float(t_c), "celsius")
            conds.pressure(float(p_bar), "bar")
            result = solver.solve(state, conds)
            if result.succeeded():
                try:
                    aqprops = AqueousProps(state)
                    vals.append(
                        total_element_molality(
                            aqprops,
                            MINERAL_CONFIG,
                            solute_species_list,
                        )
                    )
                except Exception:
                    vals.append(np.nan)
            else:
                vals.append(np.nan)

        curves[p_kbar] = np.array(vals, dtype=float)

    if psat_grid is not None:
        t_range = psat_grid["T_C"]
        p_kbar = psat_grid["P_kbar"]

        specs = EquilibriumSpecs(system)
        specs.temperature()
        specs.pressure()
        solver = EquilibriumSolver(specs)
        conds = EquilibriumConditions(specs)

        state = ChemicalState(system)
        state.set("WATER,AQ", to_real(1.0), "kg")
        state.set("H+", to_real(1e-8), "mol")
        state.set("OH-", to_real(1e-8), "mol")
        state.set(solute_species, to_real(1e-6), "mol")
        state.set(mineral_species_name, to_real(10.0), "mol")

        vals = []
        for t_c, p in zip(t_range, p_kbar):
            if np.isnan(p):
                vals.append(np.nan)
                continue
            conds.temperature(float(t_c), "celsius")
            conds.pressure(float(p) * 1000.0, "bar")
            result = solver.solve(state, conds)
            if result.succeeded():
                try:
                    aqprops = AqueousProps(state)
                    vals.append(
                        total_element_molality(
                            aqprops,
                            MINERAL_CONFIG,
                            solute_species_list,
                        )
                    )
                except Exception:
                    vals.append(np.nan)
            else:
                vals.append(np.nan)

        curves["Psat"] = np.array(vals, dtype=float)

    return curves


def compute_uncertainty_envelope(
    base_molality, temperatures_c, gibbs_samples_j_per_mol, ci
):
    """Propagate Gibbs uncertainty to molality envelopes using log-K sensitivity."""
    base = np.asarray(base_molality, dtype=float)
    t_c = np.asarray(temperatures_c, dtype=float)
    t_k = t_c + 273.15
    gas_r = 8.31446261815324

    lo_pct = 50.0 - ci / 2.0
    hi_pct = 50.0 + ci / 2.0

    lower = np.full_like(base, np.nan)
    median = np.full_like(base, np.nan)
    upper = np.full_like(base, np.nan)

    valid = np.isfinite(base) & (base > 0.0) & np.isfinite(t_k) & (t_k > 0.0)
    if not np.any(valid):
        return lower, median, upper

    factors = np.exp(np.outer(gibbs_samples_j_per_mol, 1.0 / (gas_r * t_k[valid])))
    sampled = base[valid][None, :] * factors

    lower[valid] = np.nanpercentile(sampled, lo_pct, axis=0)
    median[valid] = np.nanpercentile(sampled, 50.0, axis=0)
    upper[valid] = np.nanpercentile(sampled, hi_pct, axis=0)
    return lower, median, upper


def interval_bounds_from_samples(arr, interval_percent):
    """Return lower/upper bounds for a central interval from sample array.

    Parameters
    ----------
    arr : np.ndarray
        Sample matrix with shape (nsamples, npoints).
    interval_percent : float
        Central interval width in percent (e.g., 50, 80, 95).
    """
    alpha = max(0.0, min(1.0, 0.5 * (1.0 - float(interval_percent) / 100.0)))
    return np.nanquantile(arr, alpha, axis=0), np.nanquantile(arr, 1.0 - alpha, axis=0)


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


def aqueous_species_by_elements(dew_db, elements, mineral_config):
    """Return aqueous species names composed only of the given elements."""
    allowed = set(elements)
    pattern = re.compile(r"[A-Z][a-z]?")
    names = []
    for species in dew_db.species():
        formula = str(species.formula())
        elems = set(pattern.findall(formula))
        if elems and elems.issubset(allowed):
            names.append(species.name())

    excluded = set(mineral_config.get("excluded_species", []))

    if mineral_config.get("exclude_organics", True):
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

    # Ensure water and acid/base species are always present for charge balance.
    for base in ("WATER,AQ", "H+", "OH-"):
        if base not in excluded:
            names.append(base)

    return sorted(set(names) - excluded)


def validate_aqueous_species(dew_db, aqueous_species):
    """Check that all aqueous species exist in the DEW database."""
    if isinstance(aqueous_species, str):
        names = [n for n in aqueous_species.split() if n]
    else:
        names = [str(n) for n in aqueous_species if str(n)]

    missing = []
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
    mineral_species_override=None,
    mineral_phase_name_override=None,
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
    if mineral_species_override is not None:
        mineral_species = mineral_species_override
    else:
        mineral_species = supcrt_db.species(mineral_name)

    mineral_phase_name = (
        str(mineral_phase_name_override)
        if mineral_phase_name_override is not None
        else str(mineral_species.name())
    )

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
        aqueous_species_list = aqueous_species_by_elements(
            dew_db, include_elements, mineral_config
        )
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
            aqueous.setActivityModel(ActivityModelPerplexDEW(_dh))
            print(
                f"[OK] PerplexDEW configured: EOS={eos_name}; DH model={dh_model}; phase activity=ActivityModelPerplexDEW"
            )
        else:
            _ = StandardThermoModelDEW(params)
            aqueous.setActivityModel(ActivityModelDEW())
            print(
                f"[OK] DEW configured: EOS={eos_name}; phase activity=ActivityModelDEW"
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

    mineral = MineralPhase(mineral_phase_name)
    system = ChemicalSystem(combined_db, aqueous, mineral, *additional_mineral_phases)

    print(
        f"[OK] System built for {mineral_phase_name} solubility "
        f"(display={mineral_name}, backend={model_backend})"
    )
    return system


def load_experimental_data(csv_file):
    """Load and organize experimental data from CSV."""
    df = pd.read_csv(csv_file)
    df = df[["T_C", "P_kbar", "molality_m", "reference", "experiment_type"]].copy()
    df["P_bar"] = df["P_kbar"] * 1000.0
    # Keep NaN pressures for Hemley and other saturation curve experiments
    df["experiment_id"] = df["reference"] + " (" + df["experiment_type"] + ")"

    # Mark experiments on Psat only when pressure is unspecified.
    df["is_psat"] = df["P_kbar"].isna()

    df = df.sort_values(["is_psat", "P_kbar", "T_C"]).reset_index(drop=True)

    return df


# =============================================================================
# Main Script
# =============================================================================


def main():
    args = parse_args()
    model_backend = backend_tag(args.backend)
    dh_model = args.dh_model  # "Davies" or "ExtendedDH"
    quick_eval = bool(args.quick_eval)
    quick_n_per_group = max(1, int(args.quick_n_per_group))
    quick_npoints = max(5, int(args.quick_npoints))
    uncertainty_enabled = not args.disable_uncertainty
    uncertainty_entity = str(args.uncertainty_entity)
    uncertainty_samples = max(1, int(args.uncertainty_samples))
    uncertainty_ci = float(args.uncertainty_ci)
    uncertainty_seed = int(args.uncertainty_seed)
    covariance_json = args.covariance_json
    mineral_db_json = args.mineral_db_json
    mineral_species_code = str(args.mineral_species_code)
    outputs = output_paths(
        MINERAL_CONFIG,
        args.backend,
        dh_model,
        uncertainty_enabled=uncertainty_enabled,
    )

    print("=" * 80)
    mineral_name = MINERAL_CONFIG["mineral_name"]
    print(f"{mineral_name} Solubility Analysis - Reaktoro DEW2024")
    print(f"Mineral: {MINERAL_CONFIG['mineral_formula']}")
    print(f"Solute species: {MINERAL_CONFIG['solute_species']}")
    backend_display = model_backend + (
        f" ({dh_model})" if "perplexdew" in model_backend.lower() else ""
    )
    print(f"Aqueous backend: {backend_display}")
    if quick_eval:
        print(
            f"Quick evaluation mode: enabled (n/group={quick_n_per_group}, points/curve={quick_npoints})"
        )
    else:
        print("Quick evaluation mode: disabled")
    if uncertainty_enabled:
        print(
            "Uncertainty: "
            f"entity={uncertainty_entity}, samples={uncertainty_samples}, CI={uncertainty_ci:.1f}%"
        )
        print(f"Covariance JSON: {covariance_json}")
        print(f"Mineral DB JSON: {mineral_db_json}")
        print(f"Mineral species code: {mineral_species_code}")
        print(
            "Note: uncertainty currently propagates mineral database covariance only (aqueous covariance unavailable)."
        )
    else:
        print("Uncertainty: disabled")
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

    if quick_eval and len(exp_data) > 0:
        exp_data = select_quick_experimental_subset(exp_data, quick_n_per_group)
        print(
            f"    Quick subset selected: {len(exp_data)} points "
            "(representative low-P, high-P, and Psat)."
        )

    if len(exp_data) > 0:
        experiments = exp_data["experiment_id"].unique()
        # Separate Psat and non-Psat experiments
        psat_data = exp_data[exp_data["is_psat"]]
        non_psat_data = exp_data[~exp_data["is_psat"]]

        # Report actual non-Psat pressure levels present in the loaded dataset.
        pressures_kbar = sorted(non_psat_data["P_kbar"].dropna().unique())

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

    mineral_species_name = MINERAL_CONFIG["mineral_name"]
    mineral_base_json_data = None
    sampled_theta = None
    sampled_entities = None
    sampled_entity_indices = None
    cov_diag = None
    samples_evaluated = 0

    if uncertainty_enabled:
        if not os.path.exists(mineral_db_json):
            raise FileNotFoundError(f"Mineral DB JSON not found: {mineral_db_json}")
        with open(mineral_db_json, "r", encoding="utf-8") as f:
            mineral_base_json_data = json.load(f)

        base_mineral_db = Database.fromFile(mineral_db_json)
        base_mineral_species = base_mineral_db.species(mineral_species_code)
        mineral_species_name = mineral_species_code
        system = build_system(
            dew_db,
            supcrt_db,
            MINERAL_CONFIG,
            water_config=DEW_CONFIG,
            model_backend=model_backend,
            dh_model=dh_model,
            mineral_species_override=base_mineral_species,
            mineral_phase_name_override=mineral_species_name,
        )
    else:
        system = build_system(
            dew_db,
            supcrt_db,
            MINERAL_CONFIG,
            water_config=DEW_CONFIG,
            model_backend=model_backend,
            dh_model=dh_model,
        )

    gibbs_samples = None
    q_sigma = np.nan
    if uncertainty_enabled:
        print("\n[2b] Loading covariance and drawing Monte Carlo samples...")
        entities, sigma = load_covariance_matrix(covariance_json)
        if uncertainty_entity not in entities:
            raise ValueError(
                f"Entity '{uncertainty_entity}' not found in covariance file. "
                f"Available entities include: {', '.join(entities[:12])}..."
            )

        idx = entities.index(uncertainty_entity)
        q_sigma = float(np.sqrt(max(sigma[idx, idx], 0.0)))
        cov_diag = covariance_diagnostics(sigma)
        sampled_entities = relevant_mineral_entities(entities, mineral_base_json_data)
        if not sampled_entities:
            raise RuntimeError(
                "No overlap between covariance entities and mineral JSON species with HollandPowell Gf."
            )
        sampled_entity_indices = [entities.index(code) for code in sampled_entities]

        rng = np.random.default_rng(uncertainty_seed)
        theta_samples = sample_from_covariance_cholesky(sigma, uncertainty_samples, rng)
        sampled_theta = theta_samples[:, sampled_entity_indices]
        gibbs_samples = theta_samples[:, idx]

        print(
            f"    Loaded covariance for {len(entities)} entities. "
            f"{uncertainty_entity} sigma = {q_sigma:.3f} J/mol"
        )
        print(
            "    Covariance diagnostics: "
            f"PSD={cov_diag['is_psd']}, "
            f"rank={cov_diag['effective_rank']}/{sigma.shape[0]}, "
            f"cond={cov_diag['condition_number']:.3e}"
        )
        print(
            f"    Relevant mineral entities with uncertainty: {len(sampled_entities)} "
            f"(example: {', '.join(sampled_entities[:10])})"
        )
        print(f"    Generated {len(gibbs_samples)} covariance-consistent samples.")

    # Calculate solubility curves for each experimental pressure (drops NaN)
    mineral_name = MINERAL_CONFIG["mineral_name"]
    solute_species = MINERAL_CONFIG["solute_species"]
    solute_species_list = get_solute_species_list(MINERAL_CONFIG)
    print(f"\n[3] Calculating {mineral_name.lower()} solubility curves...")
    solubility_curves = {}

    pressures_for_curves = sorted(exp_data["P_kbar"].dropna().unique())

    n_curve_points = quick_npoints if quick_eval else N_POINTS

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
            if quick_eval:
                # In quick mode, evaluate exactly at the representative experimental temperatures.
                T_range = np.array(
                    sorted(exp_at_P["T_C"].dropna().unique()), dtype=float
                )
            else:
                T_min = (
                    max(25, T_min_cat - 0.05 * T_span) if T_span > 0 else T_min_cat - 50
                )
                T_max = (
                    min(1000, T_max_cat + 0.05 * T_span)
                    if T_span > 0
                    else T_max_cat + 50
                )
                T_range = np.linspace(T_min, T_max, n_curve_points)
        else:
            if quick_eval:
                # Defensive fallback in quick mode if pressure selection ever yields no matches.
                T_range = np.array([float(T_MIN), float(T_MAX)], dtype=float)
            else:
                T_range = np.linspace(T_MIN, T_MAX, n_curve_points)

        # Modern EquilibriumSpecs pattern (matches official tutorial)
        specs = EquilibriumSpecs(system)
        specs.temperature()
        specs.pressure()

        solver = EquilibriumSolver(specs)
        conditions = EquilibriumConditions(specs)

        state = ChemicalState(system)
        state.set("WATER,AQ", to_real(1.0), "kg")
        state.set("H+", to_real(1e-8), "mol")
        state.set("OH-", to_real(1e-8), "mol")
        state.set(solute_species, to_real(1e-6), "mol")
        state.set(mineral_species_name, to_real(10.0), "mol")

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

        molality_array = np.array(molalities, dtype=float)
        curve_data = {"T_C": T_range, "molality": molality_array}

        solubility_curves[P_kbar] = curve_data

        valid_points = np.sum(~np.isnan(molalities))
        if valid_points > 0:
            valid_idx = np.where(~np.isnan(molalities))[0]
            first_T = T_range[valid_idx[0]]
            last_T = T_range[valid_idx[-1]]
            first_m = molalities[valid_idx[0]]
            last_m = molalities[valid_idx[-1]]
            print(
                f"       Calculated {valid_points}/{n_curve_points} points (T: {first_T:.0f}-{last_T:.0f}Â°C)"
            )

    has_psat_experiments = len(psat_data) > 0
    if has_psat_experiments:
        # Calculate Psat solubility curve (following saturation pressure)
        print("    P = Psat curve...")
        if quick_eval:
            T_psat_range = np.array(
                sorted(psat_data["T_C"].dropna().unique()), dtype=float
            )
        else:
            T_psat_min = max(25, psat_data["T_C"].min() - 25)
            T_psat_max = min(374, psat_data["T_C"].max() + 25)
            T_psat_range = np.linspace(T_psat_min, T_psat_max, n_curve_points)

        P_psat_values = np.array([psat_kbar(T) for T in T_psat_range])
        valid_temps = ~np.isnan(P_psat_values)

        # Modern EquilibriumSpecs pattern for Psat calculations
        specs_psat = EquilibriumSpecs(system)
        specs_psat.temperature()
        specs_psat.pressure()

        solver_psat = EquilibriumSolver(specs_psat)
        conditions_psat = EquilibriumConditions(specs_psat)

        state_psat = ChemicalState(system)
        state_psat.set("WATER,AQ", 1.0, "kg")
        state_psat.set("H+", 1e-8, "mol")
        state_psat.set("OH-", 1e-8, "mol")
        state_psat.set(solute_species, 1e-6, "mol")
        state_psat.set(mineral_species_name, to_real(10.0), "mol")

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

        curve_psat = {
            "T_C": T_psat_range,
            "P_kbar": P_psat_values,
            "molality": np.array(psat_molalities),
        }

        solubility_curves["Psat"] = curve_psat

        valid_psat_points = np.sum(~np.isnan(psat_molalities))
        print(
            f"       Calculated {valid_psat_points}/{n_curve_points} points along Psat curve"
        )
    else:
        print("    P = Psat curve skipped (no Psat experimental points in dataset).")

    if (
        uncertainty_enabled
        and sampled_theta is not None
        and sampled_entities is not None
    ):
        print(
            "\n[3b] Full-forward uncertainty propagation (sampled mineral databases)..."
        )
        print(
            f"    Running {uncertainty_samples} forward solves over "
            f"{len(pressures_for_curves)} pressure curves"
            + (" + Psat..." if has_psat_experiments else "...")
        )
        pressures_for_sampling = [float(p) for p in pressures_for_curves]
        with tempfile.TemporaryDirectory(prefix="qz_uncert_") as temp_dir:
            all_sampled_curves = []
            for i in range(uncertainty_samples):
                if (i + 1) % max(1, uncertainty_samples // 10) == 0 or i == 0:
                    print(f"    Sample {i + 1}/{uncertainty_samples}...")
                try:
                    shifts_by_entity = {
                        code: float(shift)
                        * 1000.0  # Convert kJ/mol to J/mol (Reaktoro uses J/mol)
                        for code, shift in zip(sampled_entities, sampled_theta[i, :])
                    }
                    sampled_file = sampled_mineral_db_file(
                        mineral_base_json_data,
                        shifts_by_entity,
                        temp_dir,
                        i,
                    )
                    sampled_db = Database.fromFile(sampled_file)
                    sampled_species = sampled_db.species(mineral_species_code)
                    sampled_system = build_system(
                        dew_db,
                        supcrt_db,
                        MINERAL_CONFIG,
                        water_config=DEW_CONFIG,
                        model_backend=model_backend,
                        dh_model=dh_model,
                        mineral_species_override=sampled_species,
                        mineral_phase_name_override=mineral_species_name,
                    )
                    psat_grid = None
                    if has_psat_experiments and "Psat" in solubility_curves:
                        psat_grid = {
                            "T_C": np.asarray(
                                solubility_curves["Psat"]["T_C"], dtype=float
                            ),
                            "P_kbar": np.asarray(
                                solubility_curves["Psat"]["P_kbar"], dtype=float
                            ),
                        }

                    sampled_curves = evaluate_curves_on_fixed_grids(
                        sampled_system,
                        {
                            float(p): np.asarray(
                                solubility_curves[p]["T_C"], dtype=float
                            )
                            for p in pressures_for_sampling
                        },
                        psat_grid,
                        mineral_species_name,
                        solute_species,
                        solute_species_list,
                    )
                    all_sampled_curves.append(sampled_curves)
                    samples_evaluated += 1
                except Exception as e:
                    print(f"    Warning: sample {i + 1} failed and was skipped: {e}")

        if samples_evaluated == 0:
            raise RuntimeError(
                "All full-forward uncertainty samples failed. Cannot compute uncertainty bands."
            )

        alpha = max(0.0, min(1.0, 0.5 * (1.0 - uncertainty_ci / 100.0)))
        q_lo = alpha
        q_med = 0.5
        q_hi = 1.0 - alpha
        fan_intervals = (50, 80, 95)

        for p in pressures_for_sampling:
            arr = np.array([s[p] for s in all_sampled_curves], dtype=float)
            base_curve = np.asarray(solubility_curves[p]["molality"], dtype=float)
            arr_with_base = np.vstack([arr, base_curve[None, :]])
            # Use pure sample quantiles (including the nominal/base case) without clamping.
            lo = np.nanquantile(arr_with_base, q_lo, axis=0)
            med = np.nanquantile(arr_with_base, q_med, axis=0)
            hi = np.nanquantile(arr_with_base, q_hi, axis=0)
            solubility_curves[p]["molality_lo"] = lo
            solubility_curves[p]["molality_med"] = med
            solubility_curves[p]["molality_hi"] = hi
            for interval in fan_intervals:
                lo_i, hi_i = interval_bounds_from_samples(arr_with_base, interval)
                solubility_curves[p][f"molality_ci{interval}_lo"] = lo_i
                solubility_curves[p][f"molality_ci{interval}_hi"] = hi_i

        if has_psat_experiments and "Psat" in solubility_curves:
            psat_samples = [s.get("Psat") for s in all_sampled_curves]
            psat_samples = [x for x in psat_samples if x is not None]
            if len(psat_samples) > 0:
                arr_psat = np.array(psat_samples, dtype=float)
                base_psat = np.asarray(
                    solubility_curves["Psat"]["molality"], dtype=float
                )
                arr_psat_with_base = np.vstack([arr_psat, base_psat[None, :]])
                lo_psat = np.nanquantile(arr_psat_with_base, q_lo, axis=0)
                med_psat = np.nanquantile(arr_psat_with_base, q_med, axis=0)
                hi_psat = np.nanquantile(arr_psat_with_base, q_hi, axis=0)
                solubility_curves["Psat"]["molality_lo"] = lo_psat
                solubility_curves["Psat"]["molality_med"] = med_psat
                solubility_curves["Psat"]["molality_hi"] = hi_psat
                for interval in fan_intervals:
                    lo_i, hi_i = interval_bounds_from_samples(
                        arr_psat_with_base, interval
                    )
                    solubility_curves["Psat"][f"molality_ci{interval}_lo"] = lo_i
                    solubility_curves["Psat"][f"molality_ci{interval}_hi"] = hi_i

        print(
            f"    Full-forward uncertainty complete: {samples_evaluated}/{uncertainty_samples} "
            "samples evaluated successfully."
        )

    # Plotting
    print("\n[4] Creating plots...")

    # Single combined panel across all pressures (no low/high split).
    all_P_data = non_psat_data
    all_P_pressures = (
        sorted(all_P_data["P_kbar"].unique()) if len(all_P_data) > 0 else []
    )

    author_markers = (
        build_author_markers(exp_data["reference"]) if len(exp_data) > 0 else {}
    )

    # =========================================================================
    # PLOT 1: Combined pressure panel with Psat
    # =========================================================================
    print("    Creating combined-pressure plot...")
    fig1, ax1 = plt.subplots(figsize=(14, 8))

    def interp_curve_value(curve, t_query, key):
        x = np.asarray(curve.get("T_C", []), dtype=float)
        y = np.asarray(curve.get(key, []), dtype=float)
        valid = np.isfinite(x) & np.isfinite(y)
        n_valid = int(np.sum(valid))
        if n_valid == 0:
            return np.nan
        xv = x[valid]
        yv = y[valid]
        t = float(t_query)

        # Quick-eval often has one computed point per curve; allow exact/near match.
        close = np.isclose(xv, t, rtol=0.0, atol=1e-8)
        if np.any(close):
            return float(yv[np.argmax(close)])

        if n_valid == 1:
            # If only one point exists and temperature does not match, do not extrapolate.
            return np.nan

        order = np.argsort(xv)
        xv = xv[order]
        yv = yv[order]
        if t < xv[0] or t > xv[-1]:
            return np.nan
        return float(np.interp(t, xv, yv))

    def curve_for_experimental_row(row, pressure_curve_keys):
        p_row = row.get("P_kbar", np.nan)
        keys = [float(p) for p in pressure_curve_keys]

        # Prefer mapping to an isobar whenever a finite pressure is provided.
        # This prevents high-P points marked as psat by metadata heuristics from
        # being forced onto the Psat curve (which may be undefined at high T).
        if not pd.isna(p_row) and len(keys) > 0:
            p_row = float(p_row)
            nearest = min(keys, key=lambda p: abs(p - p_row))
            tol = max(0.01, 0.05 * max(abs(p_row), 1.0))
            if abs(nearest - p_row) <= tol:
                return solubility_curves.get(nearest), nearest

        if bool(row.get("is_psat", False)) and "Psat" in solubility_curves:
            return solubility_curves["Psat"], "Psat"

        return None, None

    def plotting_yerr_with_visibility_floor(m_lo, m_med, m_hi, min_rel=None):
        """Return asymmetric yerr with a small plotting-only floor for log-scale visibility."""
        low_err = float(m_med - m_lo)
        high_err = float(m_hi - m_med)

        if min_rel is None:
            # Quick mode can produce nearly collapsed CIs at single-point evaluations.
            min_rel = 0.20 if quick_eval else 0.03

        # Keep error bars visible even when Monte Carlo spread collapses numerically.
        floor = float(max(min_rel * float(m_med), 1e-16))
        if not np.isfinite(low_err) or low_err <= 0.0:
            low_err = floor
        if not np.isfinite(high_err) or high_err <= 0.0:
            high_err = floor

        # Prevent lower whisker from crossing <= 0 on log axis.
        low_err = min(low_err, 0.99 * float(m_med))
        return np.array([[low_err], [high_err]], dtype=float)

    def _collect_positive_values(*arrays):
        vals = []
        for arr in arrays:
            if arr is None:
                continue
            a = np.asarray(arr, dtype=float)
            good = a[np.isfinite(a) & (a > 0.0)]
            if good.size:
                vals.append(good)
        if not vals:
            return np.array([], dtype=float)
        return np.concatenate(vals)

    def dynamic_log_ylim(
        values, fallback=(1e-6, 1e1), lower_pad_decades=0.12, upper_pad_decades=0.15
    ):
        """Compute log-scale limits that include all positive values with small margins."""
        vals = _collect_positive_values(values)
        if vals.size == 0:
            return fallback

        vmin = float(np.min(vals))
        vmax = float(np.max(vals))
        if not np.isfinite(vmin) or not np.isfinite(vmax) or vmin <= 0.0:
            return fallback
        if vmax <= vmin:
            vmax = vmin * 1.1

        ylo = 10 ** (np.log10(vmin) - lower_pad_decades)
        yhi = 10 ** (np.log10(vmax) + upper_pad_decades)
        ylo = max(ylo, 1e-20)
        if yhi <= ylo:
            yhi = ylo * 10.0
        return ylo, yhi

    def draw_uncertainty_fan(ax, x_values, curve, color, zbase=4):
        """Draw nested central interval bands (95/80/50%) for visual density cue."""
        specs = [
            # (interval, alpha, facecolor)
            (95, 0.18, "#9ecae1"),
            (80, 0.28, "#fd8d3c"),
            (50, 0.40, "#e31a1c"),
        ]
        x = np.asarray(x_values, dtype=float)
        for interval, alpha_band, facecolor in specs:
            k_lo = f"molality_ci{interval}_lo"
            k_hi = f"molality_ci{interval}_hi"
            if k_lo not in curve or k_hi not in curve:
                continue
            lo = np.asarray(curve[k_lo], dtype=float)
            hi = np.asarray(curve[k_hi], dtype=float)
            valid_band = (
                np.isfinite(x)
                & np.isfinite(lo)
                & np.isfinite(hi)
                & (lo > 0.0)
                & (hi > 0.0)
                & (hi >= lo)
            )
            if np.any(valid_band):
                ax.fill_between(
                    x[valid_band],
                    lo[valid_band],
                    hi[valid_band],
                    color=facecolor,
                    alpha=alpha_band,
                    linewidth=0,
                    zorder=zbase,
                )

    # Generate colors for low-pressure experiments
    n_low = max(len(all_P_pressures), 1)
    colors_low = plt.cm.viridis(np.linspace(0, 0.9, n_low))
    P_to_color_low = {
        P: colors_low[i % len(colors_low)] for i, P in enumerate(all_P_pressures)
    }
    low_plot_positive_values = []

    # Plot all non-Psat experimental data
    for P_kbar in all_P_pressures:
        P_tol = 0.05 * P_kbar if P_kbar > 0.1 else 0.01
        subset = all_P_data[
            (all_P_data["P_kbar"] >= P_kbar - P_tol)
            & (all_P_data["P_kbar"] <= P_kbar + P_tol)
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
            low_plot_positive_values.extend(
                _collect_positive_values(
                    author_subset["molality_m"].to_numpy(dtype=float)
                )
            )

    # Plot Psat experiments (if present)
    psat_plot_data = pd.DataFrame()
    if len(psat_data) > 0:
        psat_plot_data = psat_data
        for author in psat_plot_data["reference"].unique():
            author_psat = psat_plot_data[psat_plot_data["reference"] == author]
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
            low_plot_positive_values.extend(
                _collect_positive_values(
                    author_psat["molality_m"].to_numpy(dtype=float)
                )
            )

    # Plot calculated curves for all non-Psat pressures
    for P_kbar in all_P_pressures:
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
        low_plot_positive_values.extend(
            _collect_positive_values(curve["molality"][valid])
        )
        if "molality_lo" in curve and "molality_hi" in curve:
            valid_band = (
                np.isfinite(curve["molality_lo"])
                & np.isfinite(curve["molality_hi"])
                & (curve["molality_lo"] > 0.0)
                & (curve["molality_hi"] > 0.0)
            )
            draw_uncertainty_fan(
                ax1, curve["T_C"], curve, P_to_color_low[P_kbar], zbase=4
            )
            low_plot_positive_values.extend(
                _collect_positive_values(
                    curve["molality_lo"][valid_band],
                    curve.get("molality_med", None),
                    curve["molality_hi"][valid_band],
                )
            )

    if len(all_P_pressures) == 0 and len(psat_plot_data) == 0:
        ax1.text(
            0.5,
            0.5,
            "No experimental points in this dataset",
            transform=ax1.transAxes,
            ha="center",
            va="center",
            fontsize=11,
            color="dimgray",
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.7),
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
        low_plot_positive_values.extend(
            _collect_positive_values(curve_psat["molality"][valid_psat_m])
        )
        if "molality_lo" in curve_psat and "molality_hi" in curve_psat:
            valid_psat_band = (
                np.isfinite(curve_psat["molality_lo"])
                & np.isfinite(curve_psat["molality_hi"])
                & (curve_psat["molality_lo"] > 0.0)
                & (curve_psat["molality_hi"] > 0.0)
            )
            draw_uncertainty_fan(ax1, curve_psat["T_C"], curve_psat, "purple", zbase=5)
            low_plot_positive_values.extend(
                _collect_positive_values(
                    curve_psat["molality_lo"][valid_psat_band],
                    curve_psat.get("molality_med", None),
                    curve_psat["molality_hi"][valid_psat_band],
                )
            )

    # Optional predicted-point overlays (disabled to keep only CI fan bands).
    show_predicted_point_overlay = False
    if uncertainty_enabled and show_predicted_point_overlay:
        predicted_marker_x_offset_c = 1.0
        pred_label_low = (
            "Predicted at experimental T,P "
            f"(+{predicted_marker_x_offset_c:.1f}Â°C offset; nested 50/80/95% CI error bars)"
        )
        drew_label_low = False
        low_overlay_data = pd.concat([all_P_data, psat_plot_data], ignore_index=True)
        for _, row in low_overlay_data.iterrows():
            curve, curve_key = curve_for_experimental_row(row, all_P_pressures)
            if curve is None:
                continue
            t_c = float(row["T_C"])
            m_med = interp_curve_value(curve, t_c, "molality_med")

            if not (np.isfinite(m_med) and m_med > 0.0):
                continue

            if curve_key == "Psat":
                c = "purple"
            else:
                c = P_to_color_low.get(float(curve_key), "black")

            ci_bar_specs = [
                (95, "#9ecae1", 2.1, 6.0),
                (80, "#fd8d3c", 2.3, 5.0),
                (50, "#e31a1c", 2.5, 4.0),
            ]
            drew_any_bar = False
            x_pt = t_c + predicted_marker_x_offset_c
            for interval, ecolor, lw, cap in ci_bar_specs:
                m_lo_i = interp_curve_value(curve, t_c, f"molality_ci{interval}_lo")
                m_hi_i = interp_curve_value(curve, t_c, f"molality_ci{interval}_hi")
                if not (
                    np.isfinite(m_lo_i)
                    and np.isfinite(m_hi_i)
                    and m_lo_i > 0.0
                    and m_hi_i > 0.0
                    and m_hi_i >= m_med >= m_lo_i
                ):
                    continue
                yerr = plotting_yerr_with_visibility_floor(m_lo_i, m_med, m_hi_i)
                low_plot_positive_values.extend(
                    _collect_positive_values(
                        np.array(
                            [
                                m_med - float(yerr[0, 0]),
                                m_med,
                                m_med + float(yerr[1, 0]),
                            ],
                            dtype=float,
                        )
                    )
                )
                ax1.errorbar(
                    [x_pt],
                    [m_med],
                    yerr=yerr,
                    fmt="none",
                    ecolor=ecolor,
                    elinewidth=lw,
                    capsize=cap,
                    capthick=lw,
                    alpha=0.95,
                    zorder=13,
                )
                drew_any_bar = True

            if not drew_any_bar:
                m_lo = interp_curve_value(curve, t_c, "molality_lo")
                m_hi = interp_curve_value(curve, t_c, "molality_hi")
                if (
                    np.isfinite(m_lo)
                    and np.isfinite(m_hi)
                    and m_lo > 0.0
                    and m_hi > 0.0
                    and m_hi >= m_med >= m_lo
                ):
                    yerr = plotting_yerr_with_visibility_floor(m_lo, m_med, m_hi)
                    low_plot_positive_values.extend(
                        _collect_positive_values(
                            np.array(
                                [
                                    m_med - float(yerr[0, 0]),
                                    m_med,
                                    m_med + float(yerr[1, 0]),
                                ],
                                dtype=float,
                            )
                        )
                    )
                    ax1.errorbar(
                        [x_pt],
                        [m_med],
                        yerr=yerr,
                        fmt="none",
                        ecolor="black",
                        elinewidth=2.2,
                        capsize=6.0,
                        capthick=2.2,
                        alpha=0.95,
                        zorder=13,
                    )

            ax1.plot(
                [x_pt],
                [m_med],
                marker="o",
                linestyle="none",
                markerfacecolor="white",
                markeredgecolor="black",
                markeredgewidth=1.1,
                markersize=8,
                alpha=1.0,
                zorder=14,
                label=(pred_label_low if not drew_label_low else None),
            )
            drew_label_low = True

    ax1.set_yscale("log")
    low_ymin, low_ymax = dynamic_log_ylim(
        np.asarray(low_plot_positive_values, dtype=float),
        fallback=(1e-4, 1e-1),
    )
    ax1.set_ylim(low_ymin, low_ymax)
    ax1.set_xlabel("Temperature (Â°C)", fontsize=14, fontweight="bold")
    ax1.set_ylabel(MINERAL_CONFIG["y_label"], fontsize=14, fontweight="bold")
    ax1.set_title(
        f"{MINERAL_CONFIG['plot_title']}: All Experimental Pressures",
        fontsize=16,
        fontweight="bold",
        pad=20,
    )
    ax1.grid(True, which="both", alpha=0.3, linestyle="--")

    # Add database/model info annotation
    info_text = (
        f"DEW24 (species) + SUPCRTBL ({mineral_name}) + Zhang-Duan 2005 EOS (Hâ‚‚O)"
    )
    if gibbs_samples is not None:
        info_text += (
            "\n"
            f"Shaded fan = nested 50/80/95% central CI from covariance({uncertainty_entity}), "
            f"N={len(gibbs_samples)}"
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

    fan_handles = []
    if uncertainty_enabled:
        fan_handles = [
            Patch(
                facecolor="#e31a1c", alpha=0.40, edgecolor="none", label="50% CI band"
            ),
            Patch(
                facecolor="#fd8d3c", alpha=0.28, edgecolor="none", label="80% CI band"
            ),
            Patch(
                facecolor="#9ecae1", alpha=0.18, edgecolor="none", label="95% CI band"
            ),
        ]
    handles1, labels1 = ax1.get_legend_handles_labels()
    ax1.legend(
        handles1 + fan_handles,
        labels1 + [h.get_label() for h in fan_handles],
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        fontsize=9,
        framealpha=0.9,
        ncol=1,
    )
    fig1.subplots_adjust(left=0.10, right=0.76, top=0.90, bottom=0.12)
    plt.savefig(
        outputs["low_plot"],
        dpi=180,
        bbox_inches="tight",
        pil_kwargs={"compress_level": 1},
    )
    print(f"    Combined plot saved to: {outputs['low_plot']}")
    print("    High-P split plot generation disabled (single combined plot mode).")
    plt.close(fig1)

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
        state.set("WATER,AQ", to_real(1.0), "kg")
        state.set("H+", to_real(1e-8), "mol")
        state.set("OH-", to_real(1e-8), "mol")
        state.set(solute_species, to_real(1e-6), "mol")
        state.set(mineral_species_name, to_real(10.0), "mol")
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

    fig_res.subplots_adjust(left=0.10, right=0.96, top=0.90, bottom=0.12)
    plt.savefig(
        outputs["residuals_plot"],
        dpi=180,
        bbox_inches="tight",
        pil_kwargs={"compress_level": 1},
    )
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

        lo_values = curve.get("molality_lo")
        med_values = curve.get("molality_med")
        hi_values = curve.get("molality_hi")

        if lo_values is None:
            lo_values = np.full_like(curve["molality"], np.nan, dtype=float)
        if med_values is None:
            med_values = np.full_like(curve["molality"], np.nan, dtype=float)
        if hi_values is None:
            hi_values = np.full_like(curve["molality"], np.nan, dtype=float)

        for t_c, p_kbar, molality, m_lo, m_med, m_hi in zip(
            curve["T_C"],
            pressures,
            curve["molality"],
            lo_values,
            med_values,
            hi_values,
        ):
            curve_rows.append(
                {
                    "backend": backend_display,
                    "curve_type": curve_type,
                    "P_kbar": float(p_kbar),
                    "T_C": float(t_c),
                    "molality": float(molality) if not np.isnan(molality) else np.nan,
                    "molality_lo": float(m_lo) if not np.isnan(m_lo) else np.nan,
                    "molality_med": float(m_med) if not np.isnan(m_med) else np.nan,
                    "molality_hi": float(m_hi) if not np.isnan(m_hi) else np.nan,
                }
            )

    curves_df = pd.DataFrame(curve_rows)
    curves_df.to_csv(outputs["curves_csv"], index=False)
    print(f"    Curve data saved to: {outputs['curves_csv']}")

    if gibbs_samples is not None:
        summary_df = pd.DataFrame(
            [
                {
                    "backend": backend_display,
                    "covariance_file": covariance_json,
                    "mineral_db_file": mineral_db_json,
                    "entity": uncertainty_entity,
                    "samples": len(gibbs_samples),
                    "samples_forward_success": int(samples_evaluated),
                    "relevant_entity_count": int(len(sampled_entities or [])),
                    "seed": uncertainty_seed,
                    "ci_percent": uncertainty_ci,
                    "sigma_j_per_mol": q_sigma,
                    "gibbs_shift_mean": float(np.mean(gibbs_samples)),
                    "gibbs_shift_std": float(np.std(gibbs_samples, ddof=1)),
                    "gibbs_shift_p2p5": float(np.percentile(gibbs_samples, 2.5)),
                    "gibbs_shift_p97p5": float(np.percentile(gibbs_samples, 97.5)),
                    "cov_is_psd": bool(cov_diag["is_psd"]) if cov_diag else np.nan,
                    "cov_min_eig": float(cov_diag["min_eigenvalue"])
                    if cov_diag
                    else np.nan,
                    "cov_effective_rank": float(cov_diag["effective_rank"])
                    if cov_diag
                    else np.nan,
                    "cov_condition_number": float(cov_diag["condition_number"])
                    if cov_diag
                    else np.nan,
                }
            ]
        )
        summary_df.to_csv(outputs["uncertainty_csv"], index=False)
        print(f"    Uncertainty summary saved to: {outputs['uncertainty_csv']}")

    print("\n" + "=" * 80)
    print("Analysis complete!")
    print("=" * 80)


if __name__ == "__main__":
    main()

