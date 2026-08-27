"""
Brucite Solubility Analysis using Reaktoro with DEW2024
Holland-Powell tc-ds62 mineral database + DEW2024 aqueous species
Includes uncertainty propagation from tc-ds62 thermodynamic covariance matrix

Brucite: Mg(OH)2  entity code "br" in Holland-Powell tc-ds62
Dissolution: Mg(OH)2 <-> Mg2+ + 2 OH-
Solubility metric: total dissolved Mg molality (mol/kg H2O)
"""

import os
import sys
import json
import copy
import tempfile
import argparse
import importlib

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
    SCRIPT_DIR_TMP = os.path.dirname(os.path.abspath(__file__))
    _d1 = os.path.dirname(SCRIPT_DIR_TMP)  # Mineral_Solubilities/
    _d2 = os.path.dirname(_d1)  # DEW_Experimental_Benchmark/
    _ROOT_TMP = os.path.dirname(_d2)  # reaktoro/ root
    _seen = set()

    def _dedupe(paths):
        out = []
        for p in paths:
            pp = os.path.normpath(p)
            if pp in _seen:
                continue
            _seen.add(pp)
            out.append(pp)
        return out

    pyd_candidates = [
        os.path.join(
            _ROOT_TMP, "build", "python", "package", "build", "lib", "reaktoro"
        ),
        os.path.join(_ROOT_TMP, "build", "Reaktoro", "Release"),
        os.path.join(_ROOT_TMP, "build-msvc", "Reaktoro", "Release"),
        os.path.join(_ROOT_TMP, "build-dew", "Reaktoro", "Release"),
        os.path.join(_ROOT_TMP, "build", "python", "Release"),
        os.path.join(_ROOT_TMP, "build-msvc", "python", "Release"),
        os.path.join(_ROOT_TMP, "install", "lib"),
    ]
    pyd_candidates = _dedupe(pyd_candidates)
    loaded_from = None
    for pyd_dir in pyd_candidates:
        if not os.path.isdir(pyd_dir):
            continue
        # Ensure DLLs in this directory can be found on Windows
        if os.name == "nt":
            os.environ["PATH"] = pyd_dir + os.pathsep + os.environ.get("PATH", "")
            if hasattr(os, "add_dll_directory"):
                try:
                    os.add_dll_directory(pyd_dir)
                except OSError:
                    pass
        if pyd_dir in sys.path:
            sys.path.remove(pyd_dir)
        sys.path.insert(0, pyd_dir)
        sys.modules.pop("reaktoro4py", None)
        try:
            local_mod = importlib.import_module("reaktoro4py")
        except (ModuleNotFoundError, ImportError):
            continue
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
        searched = "\n  - ".join(pyd_candidates)
        raise ModuleNotFoundError(
            "Could not import Reaktoro. Install the 'reaktoro' package or ensure "
            "reaktoro4py is on PYTHONPATH.\n"
            f"Searched local extension paths:\n  - {searched}"
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

# ---- Directory layout (script is 3 levels deep in the repo) ----
_MINERAL_SOL_DIR = os.path.dirname(SCRIPT_DIR)  # Mineral_Solubilities/
_DEW_BENCHMARK_DIR = os.path.dirname(_MINERAL_SOL_DIR)  # DEW_Experimental_Benchmark/
ROOT_DIR = os.path.dirname(_DEW_BENCHMARK_DIR)  # reaktoro/ root

# ============================================================
# MINERAL CONFIGURATION
# ============================================================
MINERAL_CONFIG = {
    "mineral_name": "br",  # Entity code in Holland-Powell tc-ds62
    "mineral_formula": "Mg(OH)2",  # Human-readable formula
    "target_element": "Mg",  # Element for total dissolved molality
    "solute_species": "Mg+2(aq)",  # Primary aqueous solute species
    "aqueous_species": "MgOH+(aq)",  # Additional aqueous species (space-separated)
    "csv_file": "brucite_DEW_testset.csv",
    "output_prefix": "brucite",
    "plot_title": "Brucite Solubility (Holland-Powell tc-ds62 + DEW2024)",
    "y_label": "Brucite Solubility (mol/kg H\u2082O)",
}
# ============================================================

# Default file paths
DEFAULT_TC_DS62_COVARIANCE_JSON = os.path.join(
    ROOT_DIR, "embedded", "databases", "hollandpowell", "tc-ds62-covariance.json"
)
DEFAULT_TC_DS62_MINERAL_JSON = os.path.join(
    ROOT_DIR, "embedded", "databases", "hollandpowell", "tc-ds62-reaktoro.json"
)
CSV_FILE = os.path.join(SCRIPT_DIR, MINERAL_CONFIG["csv_file"])

# Calculation settings
T_MIN, T_MAX = 350, 700
N_POINTS = 60
DEFAULT_PRESSURES_KBAR = [1.0, 2.0, 3.0]

# DEW water model configuration
DEW_CONFIG = {
    "eos_model": "ZhangDuan2005",
    "dielectric_model": "PowerFunction",
    "gibbs_model": "DewIntegral",
    "born_model": "Shock92Dew",
}

MODEL_BACKEND = "PerplexDEW"

PERPLEXDEW_REQUIRED_SYMBOLS = (
    "ActivityModelPerplexDEW",
    "ActivityDHModel",
    "StandardThermoModelParamsPerplexDEW",
    "StandardThermoModelPerplexDEW",
)


# =============================================================================
# Argument Parsing
# =============================================================================


def parse_args():
    parser = argparse.ArgumentParser(
        description="Brucite solubility benchmark with DEW2024 + Holland-Powell tc-ds62."
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
        help="Debye-Hückel variant for PerplexDEW backend.",
    )
    parser.add_argument(
        "--covariance-json",
        default=DEFAULT_TC_DS62_COVARIANCE_JSON,
        help="Path to tc-ds62 covariance JSON.",
    )
    parser.add_argument(
        "--mineral-db-json",
        default=DEFAULT_TC_DS62_MINERAL_JSON,
        help="Path to tc-ds62-reaktoro.json for mineral database.",
    )
    parser.add_argument(
        "--mineral-species-code",
        default="br",
        help="Entity code for brucite in mineral JSON (default: br).",
    )
    parser.add_argument(
        "--uncertainty-entity",
        default="br",
        help="Entity code in covariance JSON used for uncertainty propagation.",
    )
    parser.add_argument(
        "--uncertainty-samples",
        type=int,
        default=200,
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
        help="Disable uncertainty propagation.",
    )
    parser.add_argument(
        "--quick-eval",
        action="store_true",
        help="Use reduced point counts for fast evaluation.",
    )
    parser.add_argument(
        "--quick-npoints",
        type=int,
        default=20,
        help="Number of temperature points per curve in quick mode.",
    )
    return parser.parse_args()


# =============================================================================
# PerplexDEW symbol loading
# =============================================================================


def _local_pyd_candidates():
    return [
        os.path.join(ROOT_DIR, "build", "Reaktoro", "Release"),
        os.path.join(ROOT_DIR, "build-msvc", "Reaktoro", "Release"),
        os.path.join(ROOT_DIR, "build-dew", "Reaktoro", "Release"),
        os.path.join(ROOT_DIR, "build", "python", "Release"),
        os.path.join(ROOT_DIR, "build-msvc", "python", "Release"),
        os.path.join(ROOT_DIR, "install", "lib"),
    ]


def ensure_perplexdew_symbols():
    """Ensure PerplexDEW symbols are available, trying local reaktoro4py builds if needed."""
    missing = [name for name in PERPLEXDEW_REQUIRED_SYMBOLS if name not in globals()]
    if not missing:
        return

    for pyd_dir in _local_pyd_candidates():
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
            return

    missing_str = ", ".join(missing)
    raise RuntimeError(
        f"PerplexDEW backend requested but required symbols unavailable: {missing_str}."
    )


def to_real(value):
    try:
        return autodiff.real(value)
    except Exception:
        return value


# =============================================================================
# Covariance / Uncertainty Helpers
# =============================================================================


def load_covariance_matrix(covariance_json_path):
    """Load covariance matrix from packed upper-triangle JSON export.

    The PackedUpperTriangle is stored scaled (divide by 1000 vs source),
    so we multiply by 1000 to recover (J/mol)^2 units for sampling.
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
            v = float(packed[k]) * 1000.0  # scale to (J/mol)^2
            sigma[i, j] = v
            sigma[j, i] = v
            k += 1

    return entities, sigma


def sample_from_covariance_cholesky(sigma, nsamples, rng):
    """Draw multivariate normal samples via numerically robust Cholesky factorization."""
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
        f"Could not Cholesky-factor covariance after jitter regularization: {last_err}"
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
    """Select covariance entities that map to Holland-Powell species in the mineral JSON."""
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
    """Write a perturbed mineral DB JSON file for one Monte Carlo sample."""
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


def interval_bounds_from_samples(arr, interval_percent):
    """Return (lower, upper) quantile bounds for a central interval."""
    alpha = max(0.0, min(1.0, 0.5 * (1.0 - float(interval_percent) / 100.0)))
    return np.nanquantile(arr, alpha, axis=0), np.nanquantile(arr, 1.0 - alpha, axis=0)


# =============================================================================
# System Building
# =============================================================================


def build_system(
    dew_db,
    mineral_species,
    mineral_phase_name,
    model_backend="PerplexDEW",
    dh_model="Davies",
):
    """Build ChemicalSystem combining DEW2024 aqueous + HP brucite mineral.

    Parameters
    ----------
    dew_db : DEWDatabase
        Loaded DEW2024 aqueous database.
    mineral_species : Species
        Brucite species object from Holland-Powell database.
    mineral_phase_name : str
        Name to use for the MineralPhase (must match species name in combined db).
    model_backend : str
        "PerplexDEW" or "DEW".
    dh_model : str
        Debye-Hückel variant: "Davies" or "ExtendedDH".
    """
    combined_db = Database(dew_db.species())
    combined_db.addSpecies(mineral_species)

    base_species = "H2O(aq) H+(aq) OH-(aq)"
    solute = MINERAL_CONFIG["solute_species"]
    additional = MINERAL_CONFIG.get("aqueous_species", "")
    aqueous_species_str = f"{base_species} {solute} {additional}".strip()

    aqueous = AqueousPhase(aqueous_species_str)

    backend = str(model_backend or "DEW").strip().lower()

    try:
        if backend == "perplexdew":
            ensure_perplexdew_symbols()
            _dh = (
                ActivityDHModel.ExtendedDH
                if dh_model == "ExtendedDH"
                else ActivityDHModel.Davies
            )
            aqueous.setActivityModel(ActivityModelPerplexDEW(_dh))
        else:
            params = StandardThermoModelParamsDEW()
            eos_map = {
                "ZhangDuan2005": WaterEosModel.ZhangDuan2005,
                "ZhangDuan2009": WaterEosModel.ZhangDuan2009,
                "WagnerPruss": WaterEosModel.WagnerPruss,
                "HGK": WaterEosModel.HGK,
            }
            params.waterOptions.eosModel = eos_map.get(
                DEW_CONFIG.get("eos_model", "ZhangDuan2005"),
                WaterEosModel.ZhangDuan2005,
            )
            aqueous.setActivityModel(ActivityModelDEW())
    except Exception as e:
        if backend == "perplexdew":
            raise RuntimeError(f"Could not configure PerplexDEW backend: {e}") from e
        print(f"  Warning: aqueous backend config failed ({e}); using HKF fallback.")
        try:
            aqueous.setActivityModel(ActivityModelDEW())
        except NameError:
            aqueous.setActivityModel(ActivityModelHKF())

    mineral = MineralPhase(mineral_phase_name)
    system = ChemicalSystem(combined_db, aqueous, mineral)
    return system


# =============================================================================
# Experimental Data Loading
# =============================================================================


def load_experimental_data(csv_file):
    """Load brucite experimental data from CSV.

    Brucite CSV columns:
        Sample, P (bar), T (°C), Molality (mol/kg H2O), Reversal, reference, experiment_type

    Returns a DataFrame with standardized column names:
        T_C, P_bar, P_kbar, molality_m, reference, experiment_type
    """
    df = pd.read_csv(csv_file, encoding="latin-1")

    # Rename / standardize columns
    rename_map = {}
    for col in df.columns:
        col_lower = col.lower().strip()
        if "p (bar)" in col_lower or col_lower == "p (bar)":
            rename_map[col] = "P_bar"
        elif "t (" in col_lower or col_lower.startswith("t ("):
            rename_map[col] = "T_C"
        elif "molality" in col_lower:
            rename_map[col] = "molality_m"

    df = df.rename(columns=rename_map)

    # Ensure numeric types; coerce non-numeric (e.g., "−", "?") to NaN
    for col in ("P_bar", "T_C", "molality_m"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    df["P_kbar"] = df["P_bar"] / 1000.0

    if "reference" not in df.columns:
        df["reference"] = "unknown"
    if "experiment_type" not in df.columns:
        df["experiment_type"] = "unknown"

    # Drop rows without valid T, P, or molality
    df = df.dropna(subset=["T_C", "P_bar", "molality_m"]).copy()
    df = df.sort_values(["P_kbar", "T_C"]).reset_index(drop=True)
    return df


# =============================================================================
# Solubility Calculation Helpers
# =============================================================================


def get_solute_species_list(mineral_config):
    """Return the unique list of aqueous solute species to sum for total molality."""
    species = []
    for key in ("solute_species", "aqueous_species"):
        value = mineral_config.get(key, "")
        for name in value.split():
            if name and name not in species:
                species.append(name)
    return species


def total_element_molality(aqprops, mineral_config):
    """Return total dissolved Mg molality using element stoichiometry."""
    element = mineral_config.get("target_element")
    if element:
        return float(aqprops.elementMolality(element))
    # fallback: sum Mg-bearing solute species
    total = 0.0
    for name in get_solute_species_list(mineral_config):
        total += float(aqprops.speciesMolality(name))
    return total


def evaluate_solubility_curve(system, t_range_c, p_bar, mineral_species_name):
    """Compute brucite solubility over a T range at fixed P.

    Returns (T_C_array, molality_array).
    """
    specs = EquilibriumSpecs(system)
    specs.temperature()
    specs.pressure()
    solver = EquilibriumSolver(specs)
    conds = EquilibriumConditions(specs)

    state = ChemicalState(system)
    state.set("H2O(aq)", to_real(1.0), "kg")
    state.set("H+(aq)", to_real(1e-8), "mol")
    state.set("OH-(aq)", to_real(1e-8), "mol")
    state.set(MINERAL_CONFIG["solute_species"], to_real(1e-7), "mol")
    state.set(mineral_species_name, to_real(10.0), "mol")

    molalities = []
    for t_c in t_range_c:
        conds.temperature(float(t_c), "celsius")
        conds.pressure(float(p_bar), "bar")
        result = solver.solve(state, conds)
        if result.succeeded():
            try:
                aqprops = AqueousProps(state)
                molalities.append(total_element_molality(aqprops, MINERAL_CONFIG))
            except Exception:
                molalities.append(np.nan)
        else:
            molalities.append(np.nan)

    return np.asarray(molalities, dtype=float)


def evaluate_curves_on_fixed_grids(system, pressure_grids, mineral_species_name):
    """Evaluate molality curves on fixed T grids for one sampled system.

    Parameters
    ----------
    pressure_grids : dict
        {p_kbar: T_array_celsius}

    Returns
    -------
    dict
        {p_kbar: np.ndarray of molalities}
    """
    curves = {}
    for p_kbar, t_range in pressure_grids.items():
        p_bar = float(p_kbar) * 1000.0
        curves[p_kbar] = evaluate_solubility_curve(
            system, t_range, p_bar, mineral_species_name
        )
    return curves


# =============================================================================
# Plotting Helpers
# =============================================================================


def build_author_markers(references):
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
        a: marker_cycle[i % len(marker_cycle)] for i, a in enumerate(unique_authors)
    }


def _collect_positive(arr):
    a = np.asarray(arr, dtype=float)
    return a[np.isfinite(a) & (a > 0.0)]


def dynamic_log_ylim(values, fallback=(1e-8, 1e0), pad=0.15):
    vals = _collect_positive(values)
    if vals.size == 0:
        return fallback
    vmin, vmax = float(np.min(vals)), float(np.max(vals))
    if not np.isfinite(vmin) or vmin <= 0.0:
        return fallback
    if vmax <= vmin:
        vmax = vmin * 10.0
    return 10 ** (np.log10(vmin) - pad), 10 ** (np.log10(vmax) + pad)


def draw_uncertainty_fan(ax, x_values, curve, color, zbase=4):
    """Draw nested 50/80/95% central interval bands."""
    specs = [(95, 0.18, "#9ecae1"), (80, 0.28, "#fd8d3c"), (50, 0.40, "#e31a1c")]
    x = np.asarray(x_values, dtype=float)
    for interval, alpha_band, facecolor in specs:
        lo = curve.get(f"molality_ci{interval}_lo")
        hi = curve.get(f"molality_ci{interval}_hi")
        if lo is None or hi is None:
            continue
        lo = np.asarray(lo, dtype=float)
        hi = np.asarray(hi, dtype=float)
        valid = (
            np.isfinite(x) & np.isfinite(lo) & np.isfinite(hi) & (lo > 0) & (hi >= lo)
        )
        if np.any(valid):
            ax.fill_between(
                x[valid],
                lo[valid],
                hi[valid],
                color=facecolor,
                alpha=alpha_band,
                linewidth=0,
                zorder=zbase,
            )


def output_paths(backend_name, dh_model=None):
    prefix = MINERAL_CONFIG["output_prefix"]
    tag = str(backend_name)
    if tag.lower() == "perplexdew" and dh_model:
        tag = f"{tag}_{dh_model}"
    tag += "_uncertainty"
    return {
        "main_plot": os.path.join(
            SCRIPT_DIR, f"{prefix}_solubility_hp_dew24_{tag}.png"
        ),
        "residuals_plot": os.path.join(
            SCRIPT_DIR, f"{prefix}_residuals_hp_dew24_{tag}.png"
        ),
        "residuals_csv": os.path.join(
            SCRIPT_DIR, f"{prefix}_residuals_hp_dew24_{tag}.csv"
        ),
        "curves_csv": os.path.join(SCRIPT_DIR, f"{prefix}_curves_hp_dew24_{tag}.csv"),
        "uncertainty_csv": os.path.join(
            SCRIPT_DIR, f"{prefix}_uncertainty_hp_dew24_{tag}.csv"
        ),
    }


# =============================================================================
# Main
# =============================================================================


def main():
    args = parse_args()
    model_backend = str(args.backend)
    dh_model = args.dh_model
    quick_eval = bool(args.quick_eval)
    quick_npoints = max(5, int(args.quick_npoints))
    uncertainty_enabled = not args.disable_uncertainty
    uncertainty_entity = str(args.uncertainty_entity)
    uncertainty_samples = max(1, int(args.uncertainty_samples))
    uncertainty_ci = float(args.uncertainty_ci)
    uncertainty_seed = int(args.uncertainty_seed)
    covariance_json = args.covariance_json
    mineral_db_json = args.mineral_db_json
    mineral_species_code = str(args.mineral_species_code)
    outputs = output_paths(model_backend, dh_model)
    n_curve_points = quick_npoints if quick_eval else N_POINTS

    print("=" * 80)
    print("Brucite Solubility Analysis — Holland-Powell tc-ds62 + DEW2024")
    print(
        f"Mineral: {MINERAL_CONFIG['mineral_formula']} (entity: {mineral_species_code})"
    )
    print(f"Target element: {MINERAL_CONFIG['target_element']}")
    backend_display = model_backend + (
        f" ({dh_model})" if "perplexdew" in model_backend.lower() else ""
    )
    print(f"Aqueous backend: {backend_display}")
    print(f"Quick eval: {quick_eval}")
    if uncertainty_enabled:
        print(
            f"Uncertainty: entity={uncertainty_entity}, samples={uncertainty_samples}, CI={uncertainty_ci:.1f}%"
        )
        print(f"Covariance JSON : {covariance_json}")
        print(f"Mineral DB JSON : {mineral_db_json}")
    else:
        print("Uncertainty: disabled")
    print("=" * 80)

    # ------------------------------------------------------------------
    # [1] Load experimental data
    # ------------------------------------------------------------------
    print("\n[1] Loading experimental data...")
    if not os.path.exists(CSV_FILE):
        print(f"    WARNING: CSV not found: {CSV_FILE}")
        exp_data = pd.DataFrame(
            columns=[
                "T_C",
                "P_bar",
                "P_kbar",
                "molality_m",
                "reference",
                "experiment_type",
            ]
        )
    else:
        exp_data = load_experimental_data(CSV_FILE)
        print(f"    Loaded {len(exp_data)} experimental data points")
        if len(exp_data) > 0:
            print(
                f"    T range : {exp_data['T_C'].min():.0f} – {exp_data['T_C'].max():.0f} °C"
            )
            print(
                f"    P range : {exp_data['P_kbar'].min():.3f} – {exp_data['P_kbar'].max():.3f} kbar"
            )

    # Pressures to compute curves at: union of unique experimental pressures + defaults
    if len(exp_data) > 0:
        exp_pressures = sorted(exp_data["P_kbar"].dropna().unique())
        # Round to nearest 0.05 kbar to merge near-identical pressures
        rounded = sorted({round(float(p) / 0.05) * 0.05 for p in exp_pressures})
        pressures_for_curves = rounded
    else:
        pressures_for_curves = DEFAULT_PRESSURES_KBAR
        exp_pressures = []

    print(f"    Curve pressures: {pressures_for_curves} kbar")

    # ------------------------------------------------------------------
    # [2] Initialize databases
    # ------------------------------------------------------------------
    print("\n[2] Initializing databases...")
    try:
        dew_db = DEWDatabase("dew2024-aqueous")
        print("    DEW2024 aqueous: OK")
    except Exception as e:
        print(f"    ERROR: DEW2024 failed: {e}")
        raise

    if not os.path.exists(mineral_db_json):
        raise FileNotFoundError(f"Mineral DB JSON not found: {mineral_db_json}")

    with open(mineral_db_json, "r", encoding="utf-8") as f:
        mineral_base_json_data = json.load(f)

    base_mineral_db = Database.fromFile(mineral_db_json)
    base_mineral_species = base_mineral_db.species(mineral_species_code)
    mineral_species_name = mineral_species_code
    print(f"    HP mineral '{mineral_species_code}' loaded from tc-ds62: OK")
    print(f"    Formula: {base_mineral_species.formula()}")

    system = build_system(
        dew_db,
        base_mineral_species,
        mineral_species_name,
        model_backend=model_backend,
        dh_model=dh_model,
    )
    print(
        f"    ChemicalSystem: OK ({system.species().size()} species, {system.phases().size()} phases)"
    )

    # ------------------------------------------------------------------
    # [2b] Load covariance and draw MC samples
    # ------------------------------------------------------------------
    sampled_theta = None
    sampled_entities = None
    cov_diag = None
    br_sigma = np.nan
    samples_evaluated = 0

    if uncertainty_enabled:
        print("\n[2b] Loading covariance and drawing Monte Carlo samples...")
        if not os.path.exists(covariance_json):
            raise FileNotFoundError(f"Covariance JSON not found: {covariance_json}")
        entities, sigma = load_covariance_matrix(covariance_json)

        if uncertainty_entity not in entities:
            print(f"    WARNING: entity '{uncertainty_entity}' not in covariance file.")
            print(f"    Available (first 15): {', '.join(entities[:15])}")
            uncertainty_enabled = False
        else:
            idx = entities.index(uncertainty_entity)
            br_sigma = float(np.sqrt(max(sigma[idx, idx], 0.0)))
            cov_diag = covariance_diagnostics(sigma)
            sampled_entities = relevant_mineral_entities(
                entities, mineral_base_json_data
            )
            if not sampled_entities:
                raise RuntimeError(
                    "No overlap between covariance entities and HP mineral species in JSON."
                )
            sampled_entity_indices = [entities.index(c) for c in sampled_entities]

            rng = np.random.default_rng(uncertainty_seed)
            theta_samples = sample_from_covariance_cholesky(
                sigma, uncertainty_samples, rng
            )
            sampled_theta = theta_samples[:, sampled_entity_indices]

            print(
                f"    Covariance: {len(entities)} entities; '{uncertainty_entity}' σ = {br_sigma:.2f} J/mol"
            )
            print(
                f"    Diagnostics: PSD={cov_diag['is_psd']}, "
                f"rank={cov_diag['effective_rank']}/{sigma.shape[0]}, "
                f"cond={cov_diag['condition_number']:.2e}"
            )
            print(
                f"    Relevant entities with HP Gf: {len(sampled_entities)} "
                f"(e.g. {', '.join(sampled_entities[:8])})"
            )
            print(f"    Generated {uncertainty_samples} MC samples.")

    # ------------------------------------------------------------------
    # [3] Calculate deterministic solubility curves
    # ------------------------------------------------------------------
    print(
        f"\n[3] Calculating brucite solubility curves ({n_curve_points} T points / curve)..."
    )
    solubility_curves = {}

    for p_kbar in pressures_for_curves:
        p_bar = float(p_kbar) * 1000.0
        print(f"    P = {p_kbar:.3f} kbar ({p_bar:.0f} bar)...")

        # T range: span experimental data at this P ± 5%, or default range
        P_tol = 0.05 * max(p_kbar, 0.1)
        exp_at_P = (
            exp_data[
                (exp_data["P_kbar"] >= p_kbar - P_tol)
                & (exp_data["P_kbar"] <= p_kbar + P_tol)
            ]
            if len(exp_data) > 0
            else pd.DataFrame()
        )

        if len(exp_at_P) > 0:
            T_min_cat = float(exp_at_P["T_C"].min())
            T_max_cat = float(exp_at_P["T_C"].max())
            T_span = T_max_cat - T_min_cat
            T_min = max(25.0, T_min_cat - max(0.1 * T_span, 25.0))
            T_max = min(900.0, T_max_cat + max(0.1 * T_span, 25.0))
        else:
            T_min, T_max = float(T_MIN), float(T_MAX)

        T_range = np.linspace(T_min, T_max, n_curve_points)
        molality = evaluate_solubility_curve(
            system, T_range, p_bar, mineral_species_name
        )

        valid = int(np.sum(~np.isnan(molality)))
        print(
            f"       {valid}/{n_curve_points} valid points "
            f"(T: {T_range[0]:.0f}–{T_range[-1]:.0f} °C)"
        )
        solubility_curves[p_kbar] = {"T_C": T_range, "molality": molality}

    # ------------------------------------------------------------------
    # [3b] Full-forward uncertainty propagation
    # ------------------------------------------------------------------
    if (
        uncertainty_enabled
        and sampled_theta is not None
        and sampled_entities is not None
    ):
        print(
            f"\n[3b] Full-forward uncertainty propagation ({uncertainty_samples} samples)..."
        )
        pressure_grids = {
            float(p): np.asarray(solubility_curves[p]["T_C"], dtype=float)
            for p in pressures_for_curves
        }

        all_sampled_curves = []
        with tempfile.TemporaryDirectory(prefix="br_uncert_") as tmp_dir:
            for i in range(uncertainty_samples):
                if (i + 1) % max(1, uncertainty_samples // 10) == 0 or i == 0:
                    print(f"    Sample {i + 1}/{uncertainty_samples}...")
                try:
                    # Shifts are in kJ/mol from the covariance → convert to J/mol for Reaktoro
                    shifts_by_entity = {
                        code: float(shift) * 1000.0
                        for code, shift in zip(sampled_entities, sampled_theta[i, :])
                    }
                    sampled_file = sampled_mineral_db_file(
                        mineral_base_json_data, shifts_by_entity, tmp_dir, i
                    )
                    sampled_db = Database.fromFile(sampled_file)
                    sampled_species = sampled_db.species(mineral_species_code)
                    sampled_system = build_system(
                        dew_db,
                        sampled_species,
                        mineral_species_name,
                        model_backend=model_backend,
                        dh_model=dh_model,
                    )
                    sampled_curves = evaluate_curves_on_fixed_grids(
                        sampled_system, pressure_grids, mineral_species_name
                    )
                    all_sampled_curves.append(sampled_curves)
                    samples_evaluated += 1
                except Exception as e:
                    print(f"    Warning: sample {i + 1} failed and was skipped: {e}")

        if samples_evaluated == 0:
            raise RuntimeError("All MC uncertainty samples failed.")

        alpha = max(0.0, min(1.0, 0.5 * (1.0 - uncertainty_ci / 100.0)))
        q_lo, q_med, q_hi = alpha, 0.5, 1.0 - alpha
        fan_intervals = (50, 80, 95)

        for p in pressures_for_curves:
            arr = np.array([s[float(p)] for s in all_sampled_curves], dtype=float)
            solubility_curves[p]["molality_lo"] = np.nanquantile(arr, q_lo, axis=0)
            solubility_curves[p]["molality_med"] = np.nanquantile(arr, q_med, axis=0)
            solubility_curves[p]["molality_hi"] = np.nanquantile(arr, q_hi, axis=0)
            for ivl in fan_intervals:
                lo_i, hi_i = interval_bounds_from_samples(arr, ivl)
                solubility_curves[p][f"molality_ci{ivl}_lo"] = lo_i
                solubility_curves[p][f"molality_ci{ivl}_hi"] = hi_i

        print(
            f"    Uncertainty complete: {samples_evaluated}/{uncertainty_samples} samples succeeded."
        )

    # ------------------------------------------------------------------
    # [4] Plotting
    # ------------------------------------------------------------------
    print("\n[4] Creating main solubility plot...")

    n_pressures = len(pressures_for_curves)
    colors = plt.cm.viridis(np.linspace(0.05, 0.92, max(n_pressures, 1)))
    P_to_color = {
        p: colors[i % len(colors)] for i, p in enumerate(pressures_for_curves)
    }

    author_markers = (
        build_author_markers(exp_data["reference"]) if len(exp_data) > 0 else {}
    )

    fig, ax = plt.subplots(figsize=(14, 8))

    # Collect all positive molality values for log-scale y limits
    all_positive = []

    # --- Experimental data ---
    for _, row in exp_data.iterrows():
        p_row = float(row["P_kbar"])
        m_row = float(row["molality_m"])
        if not np.isfinite(m_row) or m_row <= 0:
            continue
        all_positive.append(m_row)
        # Match to nearest computed pressure curve for color
        if pressures_for_curves:
            nearest_p = min(pressures_for_curves, key=lambda p: abs(p - p_row))
            pt_color = P_to_color.get(nearest_p, "gray")
        else:
            pt_color = "gray"
        author = str(row.get("reference", "?"))
        marker = author_markers.get(author, "o")
        ax.scatter(
            row["T_C"],
            m_row,
            color=pt_color,
            marker=marker,
            s=70,
            alpha=0.75,
            edgecolors="black",
            linewidths=0.4,
            zorder=10,
        )

    # --- Calculated curves + uncertainty bands ---
    for p_kbar in pressures_for_curves:
        curve = solubility_curves[p_kbar]
        color = P_to_color[p_kbar]
        T = curve["T_C"]
        m = curve["molality"]
        valid = ~np.isnan(m)

        if np.any(valid):
            all_positive.extend(_collect_positive(m[valid]))
            ax.plot(
                T[valid],
                m[valid],
                color=color,
                linewidth=2.0,
                linestyle="-",
                label=f"Calc P={p_kbar:.2f} kbar",
                zorder=5,
            )

        if "molality_lo" in curve:
            lo = np.asarray(curve["molality_lo"], dtype=float)
            hi = np.asarray(curve["molality_hi"], dtype=float)
            valid_band = np.isfinite(lo) & np.isfinite(hi) & (lo > 0) & (hi > 0)
            if np.any(valid_band):
                all_positive.extend(_collect_positive(lo[valid_band]))
                all_positive.extend(_collect_positive(hi[valid_band]))
                ax.fill_between(
                    T[valid_band],
                    lo[valid_band],
                    hi[valid_band],
                    color=color,
                    alpha=0.18,
                    linewidth=0,
                    zorder=4,
                    label=f"{uncertainty_ci:.0f}% CI P={p_kbar:.2f} kbar",
                )
            draw_uncertainty_fan(ax, T, curve, color, zbase=4)

    ax.set_yscale("log")
    if all_positive:
        ylo, yhi = dynamic_log_ylim(np.array(all_positive))
        ax.set_ylim(ylo, yhi)

    # Add a simple color-coded pressure legend for experimental data
    if len(exp_data) > 0:
        import matplotlib.patches as mpatches

        exp_legend_patches = []
        for p in pressures_for_curves:
            P_tol = 0.05 * max(p, 0.1)
            subset = exp_data[
                (exp_data["P_kbar"] >= p - P_tol) & (exp_data["P_kbar"] <= p + P_tol)
            ]
            if len(subset) > 0:
                exp_legend_patches.append(
                    mpatches.Patch(color=P_to_color[p], label=f"Exp P≈{p:.2f} kbar")
                )
        if exp_legend_patches:
            exp_legend = ax.legend(
                handles=exp_legend_patches,
                loc="upper right",
                fontsize=8,
                title="Experimental data",
            )
            ax.add_artist(exp_legend)

    ax.set_xlabel("Temperature (°C)", fontsize=12)
    ax.set_ylabel(MINERAL_CONFIG["y_label"], fontsize=12)
    title = MINERAL_CONFIG["plot_title"]
    if uncertainty_enabled and samples_evaluated > 0:
        title += f"\n(uncertainty: {uncertainty_ci:.0f}% CI, n={samples_evaluated} MC samples)"
    ax.set_title(title, fontsize=13)
    ax.legend(loc="upper left", fontsize=8)
    ax.grid(True, which="both", alpha=0.3)

    plt.tight_layout()
    plt.savefig(outputs["main_plot"], dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"    Saved: {outputs['main_plot']}")

    # ------------------------------------------------------------------
    # [5] Residuals
    # ------------------------------------------------------------------
    if len(exp_data) > 0:
        print("\n[5] Computing residuals (calculated vs experimental)...")
        residuals = []
        for _, row in exp_data.iterrows():
            p_row = float(row["P_kbar"])
            t_row = float(row["T_C"])
            m_exp = float(row["molality_m"])
            if not np.isfinite(m_exp) or m_exp <= 0:
                continue

            # Find nearest pressure curve
            if not pressures_for_curves:
                continue
            nearest_p = min(pressures_for_curves, key=lambda p: abs(p - p_row))
            if abs(nearest_p - p_row) > 0.15:
                continue

            curve = solubility_curves.get(nearest_p)
            if curve is None:
                continue

            # Interpolate calculated molality at this T
            xv = np.asarray(curve["T_C"], dtype=float)
            yv = np.asarray(curve["molality"], dtype=float)
            valid = np.isfinite(xv) & np.isfinite(yv)
            if not np.any(valid):
                continue
            xv, yv = xv[valid], yv[valid]
            if t_row < xv.min() or t_row > xv.max():
                continue
            m_calc = float(np.interp(t_row, xv, yv))
            if not np.isfinite(m_calc) or m_calc <= 0:
                continue

            log10_ratio = np.log10(m_calc / m_exp)
            residuals.append(
                {
                    "reference": row.get("reference", "?"),
                    "experiment_type": row.get("experiment_type", "?"),
                    "T_C": t_row,
                    "P_kbar": p_row,
                    "P_kbar_nearest": nearest_p,
                    "molality_exp": m_exp,
                    "molality_calc": m_calc,
                    "log10_ratio": log10_ratio,
                }
            )

        if residuals:
            resid_df = pd.DataFrame(residuals)
            resid_df.to_csv(outputs["residuals_csv"], index=False)
            print(f"    {len(resid_df)} residuals computed.")
            rmse = float(np.sqrt(np.mean(resid_df["log10_ratio"] ** 2)))
            bias = float(np.mean(resid_df["log10_ratio"]))
            print(f"    RMSE (log10 ratio): {rmse:.4f}   |   Bias: {bias:.4f}")
            print(f"    Saved residuals CSV: {outputs['residuals_csv']}")

            # Residuals plot
            fig2, ax2 = plt.subplots(figsize=(10, 6))
            author_markers2 = build_author_markers(resid_df["reference"])
            for author in resid_df["reference"].unique():
                sub = resid_df[resid_df["reference"] == author]
                marker = author_markers2.get(author, "o")
                ax2.scatter(
                    sub["T_C"],
                    sub["log10_ratio"],
                    marker=marker,
                    s=70,
                    alpha=0.75,
                    edgecolors="black",
                    linewidths=0.4,
                    label=author,
                    zorder=5,
                )
            ax2.axhline(0, color="black", linewidth=1.2, linestyle="--", alpha=0.7)
            ax2.axhline(0.3, color="gray", linewidth=0.8, linestyle=":", alpha=0.6)
            ax2.axhline(-0.3, color="gray", linewidth=0.8, linestyle=":", alpha=0.6)
            ax2.set_xlabel("Temperature (°C)", fontsize=12)
            ax2.set_ylabel("log10(calc / exp)", fontsize=12)
            ax2.set_title(
                f"Brucite Residuals — HP tc-ds62 + DEW2024 ({backend_display})\n"
                f"RMSE={rmse:.4f}  Bias={bias:.4f}",
                fontsize=12,
            )
            ax2.legend(fontsize=8, loc="best")
            ax2.grid(True, alpha=0.3)
            plt.tight_layout()
            plt.savefig(outputs["residuals_plot"], dpi=150, bbox_inches="tight")
            plt.close(fig2)
            print(f"    Saved residuals plot: {outputs['residuals_plot']}")

    # ------------------------------------------------------------------
    # [6] Save curves CSV
    # ------------------------------------------------------------------
    print("\n[6] Saving curves data...")
    rows = []
    for p_kbar, curve in solubility_curves.items():
        T = curve["T_C"]
        m = curve["molality"]
        m_lo = curve.get("molality_lo", np.full_like(m, np.nan))
        m_med = curve.get("molality_med", np.full_like(m, np.nan))
        m_hi = curve.get("molality_hi", np.full_like(m, np.nan))
        for j in range(len(T)):
            rows.append(
                {
                    "P_kbar": p_kbar,
                    "T_C": T[j],
                    "molality_calc": m[j],
                    "molality_lo": m_lo[j],
                    "molality_med": m_med[j],
                    "molality_hi": m_hi[j],
                }
            )
    curves_df = pd.DataFrame(rows)
    curves_df.to_csv(outputs["curves_csv"], index=False)
    print(f"    Saved: {outputs['curves_csv']}")

    # ------------------------------------------------------------------
    # [7] Summary
    # ------------------------------------------------------------------
    print("\n" + "=" * 80)
    print("Summary")
    print(
        f"  Mineral      : Brucite Mg(OH)2 | HP tc-ds62 entity '{mineral_species_code}'"
    )
    print(f"  Aqueous DB   : DEW2024 | backend={backend_display}")
    print(f"  Pressures    : {pressures_for_curves} kbar")
    if uncertainty_enabled and samples_evaluated > 0:
        print(
            f"  Uncertainty  : {uncertainty_ci:.0f}% CI from {samples_evaluated} MC samples "
            f"| br σ = {br_sigma:.2f} J/mol"
        )
        print(
            f"  Note: Uncertainty reflects Holland-Powell Gf covariance only "
            f"(aqueous species covariance not available)."
        )
    print(f"  Output dir   : {SCRIPT_DIR}")
    print("=" * 80)


if __name__ == "__main__":
    main()
