"""
quartz_H2OCO2_solubility.py

Quartz (SiO2) solubility in H2O-CO2 mixtures as a function of XCO2.

Fixed conditions : T = 800 Â°C,  P = 10 kbar
Fluid            : H2O-CO2 mixture, XCO2 swept from 0 to ~0.85
Aqueous model    : PerplexDEW (Davies or ExtendedDH)
CO2 EOS          : Zhang-Duan 2009 (ActivityModelPerplexGFSM)
Experimental data: Newton & Manning (2000), Table 2 H2O-CO2 series
"""

import os
import sys
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
    import autodiff  # noqa: F401
except ModuleNotFoundError:
    autodiff = None
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BENCHMARK_DIR = os.path.dirname(SCRIPT_DIR)
ROOT_DIR = os.path.dirname(BENCHMARK_DIR)

# Try to import Reaktoro; fall back to local extension modules if not installed.
try:
    from reaktoro import *  # noqa: F401,F403
except ModuleNotFoundError:
    _pyd_candidates = [
        os.path.join(ROOT_DIR, "build", "Reaktoro", "Release"),
        os.path.join(ROOT_DIR, "build", "Reaktoro", "Release"),
        os.path.join(ROOT_DIR, "build", "Reaktoro", "Release"),
    ]
    _loaded_from = None
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
            "ensure reaktoro4py is available in a build/build/build folder."
        )
    print(f"Using local reaktoro4py extension from {_loaded_from}.")

try:
    Warnings.disable(906)
except Exception:
    pass

# =============================================================================
# Configuration
# =============================================================================

MINERAL_CONFIG = {
    "mineral_name": "Quartz",
    "mineral_formula": "SiO2",
    "target_element": "Si",
    "solute_species": "SiO2(aq)",
    # H2_aq / O2_aq / HO2- are dissolved redox gases irrelevant to silicate
    # chemistry and are excluded to keep the system clean.
    "aqueous_species": "HSiO3-(aq) Si2O4(aq) Si3O6(aq)",
    "csv_file": "quartz_H2OCO2_DEW_testset.csv",
    "output_prefix": "quartz_H2OCO2",
    "plot_title": "Quartz Solubility in H\u2082O-CO\u2082 (T\u2009=\u2009800\u00b0C)",
    "y_label": "Quartz Solubility (mol/kg-H\u2082O)",
}

CSV_FILE = os.path.join(SCRIPT_DIR, MINERAL_CONFIG["csv_file"])

# Fixed T and model pressures matching both experimental datasets
T_FIXED_C = 800.0
MODEL_PRESSURES_KBAR = [
    9.0,
    10.0,
]  # 9 kbar = Shmulovich+01; 10 kbar = Newton&Manning 2000

# XCO2 sweep for the model curve
XCO2_MODEL = np.concatenate([[0.0], np.linspace(0.01, 0.85, 50)])

# Water EOS configuration (shared DEW settings)
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
    "ActivityModelPerplexGFSM",
    "ActivityModelParamsPerplexGFSM",
    "PerpleXHybridEosOptions",
    "PerpleXCO2Eos",
    "PerpleXWaterEos",
)


# =============================================================================
# CLI
# =============================================================================


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Quartz solubility in H2O-CO2 vs XCO2. "
            "Newton & Manning (2000) P=10 kbar + Shmulovich et al. (2001) P=9 kbar, T=800Â°C."
        )
    )
    parser.add_argument(
        "--backend",
        default=MODEL_BACKEND,
        choices=["DEW", "PerplexDEW"],
        help="Aqueous backend to use (default: PerplexDEW).",
    )
    parser.add_argument(
        "--dh-model",
        default="Davies",
        choices=["Davies", "ExtendedDH"],
        help="Debye-HÃ¼ckel variant for PerplexDEW backend (default: Davies).",
    )
    return parser.parse_args()


def backend_tag(name):
    return str(name or "DEW").strip()


def output_paths(mineral_config, backend_name, dh_model=None):
    tag = backend_tag(backend_name)
    if tag.lower() == "perplexdew" and dh_model:
        tag = f"{tag}_{dh_model}"
    prefix = mineral_config["output_prefix"]
    return {
        "solubility_plot": os.path.join(
            SCRIPT_DIR, f"{prefix}_solubility_vs_xco2_{tag}.png"
        ),
        "residuals_plot": os.path.join(SCRIPT_DIR, f"{prefix}_residuals_{tag}.png"),
    }


# =============================================================================
# PerplexDEW Symbol Loading
# =============================================================================


def _local_pyd_candidates():
    return [
        os.path.join(ROOT_DIR, "build", "Reaktoro", "Release"),
        os.path.join(ROOT_DIR, "build", "Reaktoro", "Release"),
        os.path.join(ROOT_DIR, "build", "Reaktoro", "Release"),
    ]


def ensure_perplexdew_symbols():
    """Ensure PerplexDEW symbols are accessible in globals()."""
    missing = [n for n in PERPLEXDEW_REQUIRED_SYMBOLS if n not in globals()]
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
        missing = [n for n in PERPLEXDEW_REQUIRED_SYMBOLS if n not in globals()]
        if not missing:
            print(f"Using local reaktoro4py extension from {pyd_dir}.")
            return
    raise RuntimeError(
        "PerplexDEW backend requested but required symbols are unavailable: "
        + ", ".join(missing)
    )


# =============================================================================
# Helpers
# =============================================================================


def element_species_coeffs(dew_db, species_names, element):
    """Return (species_name, stoichiometric_coeff) pairs for the given element."""
    entries = []
    for name in species_names:
        try:
            sp = dew_db.species(name)
            coeff = float(sp.elements().coefficient(element))
        except Exception:
            continue
        if coeff != 0.0:
            entries.append((name, coeff))
    return entries


def total_element_molality(aqprops, mineral_config):
    """Sum molality of all Si-bearing aqueous species, weighted by stoichiometry."""
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
    # Fallback: primary solute only
    return float(aqprops.speciesMolality(mineral_config["solute_species"]))


def validate_aqueous_species(dew_db, aqueous_species_str):
    names = (
        aqueous_species_str.split()
        if isinstance(aqueous_species_str, str)
        else list(aqueous_species_str)
    )
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


# =============================================================================
# System Builder
# =============================================================================


def build_system(
    dew_db,
    supcrt_db,
    mineral_config,
    water_config=None,
    model_backend="PerplexDEW",
    dh_model="Davies",
):
    """Build ChemicalSystem: quartz + Si-H-O aqueous (DEW) + CO2(g) (ZhangDuan09)."""
    if water_config is None:
        water_config = DEW_CONFIG

    mineral_name = mineral_config["mineral_name"]
    mineral_species = supcrt_db.species(mineral_name)
    combined_db = Database(dew_db.species())
    combined_db.addSpecies(mineral_species)

    # Aqueous phase: base + SiO2_aq + HSiO3- + Si2O4_aq + Si3O6_aq
    base_species = "H2O(aq) H+(aq) OH-(aq)"
    solute = mineral_config["solute_species"]
    additional = mineral_config.get("aqueous_species", "")
    if additional:
        aqueous_species_str = f"{base_species} {solute} {additional}"
    else:
        aqueous_species_str = f"{base_species} {solute}"

    # Register element stoichiometry for reporting total dissolved Si
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

    # --- Configure aqueous activity / EOS model ---
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
        born_map = {"Shock92Dew": WaterBornModel.Shock92Dew}

        params.waterOptions.eosModel = eos_map.get(
            water_config.get("eos_model", "ZhangDuan2005"),
            WaterEosModel.ZhangDuan2005,
        )
        params.waterOptions.dielectricModel = dielectric_map.get(
            water_config.get("dielectric_model", "PowerFunction"),
            WaterDielectricModel.PowerFunction,
        )
        params.waterOptions.gibbsModel = gibbs_map.get(
            water_config.get("gibbs_model", "DewIntegral"),
            WaterGibbsModel.DewIntegral,
        )
        params.waterOptions.bornModel = born_map.get(
            water_config.get("born_model", "Shock92Dew"),
            WaterBornModel.Shock92Dew,
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
            print(f"âœ“ PerplexDEW configured: EOS={eos_name}; DH model={dh_model}")
        else:
            aqueous.setActivityModel(ActivityModelDEW())
            print(f"âœ“ DEW configured: EOS={eos_name}")

    except Exception as e:
        if str(model_backend or "").strip().lower() == "perplexdew":
            raise RuntimeError(f"Could not configure PerplexDEW: {e}") from e
        try:
            aqueous.setActivityModel(ActivityModelDEW())
        except NameError:
            aqueous.setActivityModel(ActivityModelHKF())

    mineral = MineralPhase(mineral_name)

    # No gas phase: the DEW model computes species standard states from pure-water
    # EOS at T,P. CO2-induced water activity changes are applied analytically via
    # a separate GFSM gas computation (see compute_aH2O_gfsm / build_gfsm_system).
    system = ChemicalSystem(combined_db, aqueous, mineral)
    print(f"âœ“ System built for {mineral_name} solubility (pure-H2O reference)")
    return system


# =============================================================================
# Experimental Data
# =============================================================================


def load_experimental_data(csv_file):
    """Load Newton & Manning (2000) quartz H2O-CO2 data."""
    df = pd.read_csv(csv_file)
    for col in ("molality_m", "xco2", "P_kbar", "T_C"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df.dropna(subset=["xco2", "molality_m"]).reset_index(drop=True)


# =============================================================================
# XCO2 Sweep Solver
# =============================================================================


def solve_xco2_sweep(system, mineral_config, xco2_array, T_C, P_bar):
    """Compute quartz solubility (pure-H2O DEW reference) at fixed T and P.

    Returns a flat curve: the DEW model gives the pure-water solubility
    regardless of XCO2 because Si species have anhydrous formulas (SiO2, no H)
    in the database.  CO2-induced reduction is applied separately via the
    activity-corrected curve (see compute_activity_corrected_curve).
    """
    solver = EquilibriumSolver(system)
    conditions = EquilibriumConditions(system)
    conditions.temperature(T_C, "celsius")
    conditions.pressure(P_bar, "bar")

    results_xco2 = []
    results_m = []

    for xco2 in xco2_array:
        state = ChemicalState(system)
        state.set("H2O(aq)", 1.0, "kg")
        state.set("SiO2(aq)", 1e-6, "mol")
        state.set("Quartz", 10.0, "mol")

        try:
            result = solver.solve(state, conditions)
            if not result.succeeded():
                print(f"    XCO2 = {xco2:.4f}: solver did not converge â€” skipped.")
                continue
        except Exception as exc:
            print(f"    XCO2 = {xco2:.4f}: exception â†’ {exc}")
            continue

        aqprops = AqueousProps(state)
        m_si = total_element_molality(aqprops, mineral_config)
        if np.isfinite(m_si) and m_si > 0.0:
            results_xco2.append(xco2)
            results_m.append(m_si)

    return np.array(results_xco2), np.array(results_m)


def get_baseline_species_molalities(system, mineral_config, T_C, P_bar):
    """Solve DEW at pure H2O (XCO2=0) and return per-species molalities.

    Returns dict: {species_name: molality_mol_per_kg}
    Also returns total dissolved Si (mol/kg).
    """
    solver = EquilibriumSolver(system)
    conditions = EquilibriumConditions(system)
    conditions.temperature(T_C, "celsius")
    conditions.pressure(P_bar, "bar")

    state = ChemicalState(system)
    state.set("H2O(aq)", 1.0, "kg")
    state.set("SiO2(aq)", 1e-6, "mol")
    state.set("Quartz", 10.0, "mol")
    result = solver.solve(state, conditions)
    if not result.succeeded():
        raise RuntimeError(
            f"Baseline DEW solve failed at T={T_C}Â°C P={P_bar / 1000:.1f} kbar"
        )

    aqprops = AqueousProps(state)
    si_species = ["SiO2(aq)", "HSiO3-(aq)", "Si2O4(aq)", "Si3O6(aq)"]
    m0 = {}
    for sp in si_species:
        try:
            m0[sp] = float(aqprops.speciesMolality(sp))
        except Exception:
            m0[sp] = 0.0
    m0["total_Si"] = total_element_molality(aqprops, mineral_config)
    return m0


def build_gfsm_system(supcrt_db):
    """Build a gas-only ChemicalSystem (CO2(g)+H2O(g)) with GFSM Zhang-Duan 2009.

    Used exclusively to compute the H2O-CO2 water activity ratio a_H2O = f_H2O_mix/f_H2O_pure
    at given T, P, XCO2.  No equilibrium solving is performed â€” ChemicalProps is
    evaluated directly at the prescribed fluid composition.
    """
    ensure_perplexdew_symbols()
    co2_sp = supcrt_db.species("CO2(g)")
    h2o_sp = supcrt_db.species("H2O(g)")
    gas_db = Database([co2_sp, h2o_sp])
    gas_phase = GaseousPhase("CO2(g) H2O(g)")
    gfsm_params = ActivityModelParamsPerplexGFSM()
    hybrid_opts = PerpleXHybridEosOptions()
    hybrid_opts.co2 = PerpleXCO2Eos.ZhangDuan09
    hybrid_opts.water = PerpleXWaterEos.ZhangDuan09
    gfsm_params.hybridEosOptions = hybrid_opts
    gas_phase.setActivityModel(ActivityModelPerplexGFSM(gfsm_params))
    return ChemicalSystem(gas_db, gas_phase)


def compute_aH2O_gfsm(gfsm_system, xco2_array, T_C, P_bar):
    """Compute a_H2O = f_H2O_mix / f_H2O_pure for each XCO2 using the GFSM EOS.

    The fugacity ratio is the thermodynamic water activity in the H2O-CO2 fluid
    relative to pure water at the same T and P.  This ratio determines how much
    the mixed fluid lowers the effective 'availability' of water for dissolution
    reactions.

    Parameters
    ----------
    gfsm_system : ChemicalSystem (gas CO2(g)+H2O(g) with GFSM ZhangDuan09)
    xco2_array  : array of CO2 mole fractions
    T_C         : temperature (Â°C)
    P_bar       : pressure (bar)

    Returns
    -------
    numpy array of a_H2O values (1.0 at XCO2=0, decreasing toward 0 at XCO2â†’1)
    """
    _N = 1000.0  # total moles for state (arbitrary; only mole fraction matters)

    # Pure water baseline
    state_pure = ChemicalState(gfsm_system)
    state_pure.set("H2O(g)", _N, "mol")
    state_pure.set("CO2(g)", 1e-10, "mol")
    state_pure.setTemperature(T_C, "celsius")
    state_pure.setPressure(P_bar, "bar")
    lna_pure = float(ChemicalProps(state_pure).speciesActivityLn("H2O(g)"))

    result = []
    for xco2 in xco2_array:
        if xco2 <= 0.0:
            result.append(1.0)
            continue
        if xco2 >= 1.0:
            result.append(0.0)
            continue
        state = ChemicalState(gfsm_system)
        state.set("H2O(g)", (1.0 - xco2) * _N, "mol")
        state.set("CO2(g)", xco2 * _N, "mol")
        state.setTemperature(T_C, "celsius")
        state.setPressure(P_bar, "bar")
        lna_mix = float(ChemicalProps(state).speciesActivityLn("H2O(g)"))
        result.append(float(np.exp(lna_mix - lna_pure)))
    return np.array(result)


def compute_activity_corrected_curve(m0_dict, aH2O_array, xco2_array):
    """Apply the a_H2O^n correction to the DEW baseline molalities.

    Physical basis
    --------------
    Each Si species dissolved in water corresponds to a hydrated form:
      SiO2_aq  = Hâ‚„SiOâ‚„ = SiOâ‚‚ + 2 Hâ‚‚O   â†’ n_Hâ‚‚O = 2
      HSiOâ‚ƒâ»  dissolves as SiOâ‚‚ + Hâ‚‚O     â†’ n_Hâ‚‚O = 1
      Siâ‚‚Oâ‚„aq = 2Ã—SiOâ‚‚ + 4 Hâ‚‚O            â†’ n_Hâ‚‚O = 4
      Siâ‚ƒOâ‚†aq = 3Ã—SiOâ‚‚ + 6 Hâ‚‚O            â†’ n_Hâ‚‚O = 6

    Although the DEW database stores anhydrous Si formulas (SiOâ‚‚), the
    Gibbs energies were tabulated using the pure-water DEW EOS which
    implicitly encodes the 2-Hâ‚‚O hydration of orthosilicic acid.  In a
    mixed Hâ‚‚O-COâ‚‚ fluid the effective equilibrium constant is reduced by
    a factor  a_Hâ‚‚O^n_Hâ‚‚O  relative to pure water.

    So the corrected total dissolved Si (mol Si / kg Hâ‚‚O) is:
        m_Si = 1Â·mâ‚€(SiOâ‚‚_aq)Â·aáµ¥Â²  +  1Â·mâ‚€(HSiOâ‚ƒâ»)Â·aáµ¥Â¹
             + 2Â·mâ‚€(Siâ‚‚Oâ‚„_aq)Â·aáµ¥â´  +  3Â·mâ‚€(Siâ‚ƒOâ‚†_aq)Â·aáµ¥â¶

    Parameters
    ----------
    m0_dict    : dict with keys 'SiO2(aq)', 'HSiO3-(aq)', 'Si2O4(aq)', 'Si3O6(aq)'
                 giving pure-Hâ‚‚O DEW molalities
    aH2O_array : array of a_Hâ‚‚O values from compute_aH2O_gfsm (same length as xco2_array)
    xco2_array : XCOâ‚‚ values (for filtering valid range)

    Returns
    -------
    (xco2_out, m_corrected): arrays of valid XCOâ‚‚ and corrected molalities
    """
    m0_sio2 = float(m0_dict.get("SiO2(aq)", 0.0))
    m0_hsi = float(m0_dict.get("HSiO3-(aq)", 0.0))
    m0_si2 = float(m0_dict.get("Si2O4(aq)", 0.0))
    m0_si3 = float(m0_dict.get("Si3O6(aq)", 0.0))

    a = np.asarray(aH2O_array, dtype=float)
    m_corr = (
        1 * m0_sio2 * a**2 + 1 * m0_hsi * a**1 + 2 * m0_si2 * a**4 + 3 * m0_si3 * a**6
    )
    valid = np.isfinite(m_corr) & (m_corr > 0.0)
    return np.asarray(xco2_array)[valid], m_corr[valid]


# =============================================================================
# Plotting
# =============================================================================


def plot_solubility_vs_xco2(
    model_results_flat,
    model_results_corr,
    exp_df,
    mineral_config,
    output_path,
    dh_model,
):
    """
    model_results_flat : dict {P_kbar: (xco2_array, m_array)} -- pure-H2O DEW (flat, dashed)
    model_results_corr : dict {P_kbar: (xco2_array, m_array)} -- a(H2O) corrected (solid)
    exp_df             : DataFrame with columns xco2, molality_m, P_kbar, reference, Name
    """
    pressure_palette = {
        10.0: {"color": "#d62728", "label": "10 kbar"},
        9.0: {"color": "#1f77b4", "label": "9 kbar"},
    }
    ref_markers = {
        "Newton_Manning_2000": ("o", "Newton & Manning (2000)  10 kbar"),
        "Shmulovich_Graham_Yardley_2001": ("s", "Shmulovich et al. (2001)  9 kbar"),
    }

    fig, ax = plt.subplots(figsize=(9, 6))

    # --- Flat DEW curves (pure-H2O reference, dashed) ---
    for P_kbar in sorted(model_results_flat.keys()):
        xco2_arr, m_arr = model_results_flat[P_kbar]
        if len(xco2_arr) == 0:
            continue
        style = pressure_palette.get(
            P_kbar, {"color": "gray", "label": f"{P_kbar} kbar"}
        )
        ax.plot(
            xco2_arr,
            m_arr,
            color=style["color"],
            linewidth=1.8,
            linestyle="--",
            alpha=0.65,
            label=f"DEW ({dh_model}) pure H\u2082O \u2014 {style['label']}",
            zorder=3,
        )

    # --- Activity-corrected curves (solid) ---
    for P_kbar in sorted(model_results_corr.keys()):
        xco2_arr, m_arr = model_results_corr[P_kbar]
        if len(xco2_arr) == 0:
            continue
        style = pressure_palette.get(
            P_kbar, {"color": "gray", "label": f"{P_kbar} kbar"}
        )
        ax.plot(
            xco2_arr,
            m_arr,
            color=style["color"],
            linewidth=2,
            linestyle="-",
            label=f"DEW + a(H\u2082O) corr. \u2014 {style['label']}",
            zorder=4,
        )

    # --- Experimental data ---
    for ref, (marker, label) in ref_markers.items():
        subset = exp_df[exp_df["reference"] == ref]
        if subset.empty:
            continue
        p_vals = subset["P_kbar"].unique()
        color = (
            pressure_palette.get(float(p_vals[0]), {"color": "gray"})["color"]
            if len(p_vals) == 1
            else "gray"
        )
        ax.scatter(
            subset["xco2"],
            subset["molality_m"],
            color=color,
            marker=marker,
            s=80,
            zorder=5,
            label=label,
            edgecolors="black",
            linewidths=0.7,
        )
        if "Name" in subset.columns:
            for _, row in subset.iterrows():
                ax.annotate(
                    str(row["Name"]),
                    (row["xco2"], row["molality_m"]),
                    textcoords="offset points",
                    xytext=(5, 3),
                    fontsize=7.5,
                    color=color,
                )

    ax.set_yscale("log")
    ax.set_xlabel(r"$X_{\mathrm{CO_2}}$ (mole fraction)", fontsize=13)
    ax.set_ylabel(mineral_config["y_label"], fontsize=13)
    ax.set_title(
        mineral_config["plot_title"]
        + "\n(dashed = pure-H\u2082O DEW; solid = DEW + GFSM a(H\u2082O) correction)",
        fontsize=11,
        fontweight="bold",
    )
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(True, which="both", alpha=0.3)
    ax.set_xlim(left=0.0)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"    Saved: {output_path}")
    plt.close()


def plot_residuals(
    model_results_corr,
    exp_df,
    mineral_config,
    output_path,
    dh_model,
):
    """Bar chart of percent residuals for the activity-corrected curve vs experiment.

    model_results_corr: dict {P_kbar: (xco2_array, m_array)}
    """
    if not model_results_corr or len(exp_df) == 0:
        return
    model_results = model_results_corr  # alias for readability below

    pressure_palette = {
        10.0: {"color": "#d62728", "label": "10 kbar"},
        9.0: {"color": "#1f77b4", "label": "9 kbar"},
    }

    all_rows = []
    print(f"\n    Residuals (PerplexDEW {dh_model}):")
    for _, row in exp_df.iterrows():
        P = float(row["P_kbar"])
        if P not in model_results:
            continue
        xco2_arr, m_arr = model_results[P]
        if len(xco2_arr) == 0:
            continue
        m_mod = float(
            np.interp(row["xco2"], xco2_arr, m_arr, left=np.nan, right=np.nan)
        )
        if not np.isfinite(m_mod) or row["molality_m"] <= 0:
            continue
        res_pct = (m_mod - row["molality_m"]) / row["molality_m"] * 100.0
        name = row.get("Name", str(_)) if "Name" in exp_df.columns else str(_)
        print(
            f"      {name} ({P:.0f} kbar): XCO2={row['xco2']:.3f}  "
            f"exp={row['molality_m']:.4f}  model={m_mod:.4f}  res={res_pct:+.1f}%"
        )
        style = pressure_palette.get(P, {"color": "gray", "label": f"{P:.0f} kbar"})
        all_rows.append(
            {
                "xco2": row["xco2"],
                "res": res_pct,
                "name": name,
                "color": style["color"],
                "P_label": style["label"],
            }
        )

    if not all_rows:
        print("    No valid residuals to plot.")
        return

    rdf = pd.DataFrame(all_rows).sort_values("xco2")
    colors = ["steelblue" if r >= 0 else "tomato" for r in rdf["res"]]

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.bar(
        rdf["xco2"],
        rdf["res"],
        width=0.018,
        color=colors,
        edgecolor="black",
        linewidth=0.5,
    )
    ax.axhline(0, color="black", linewidth=1)
    for level in (50, -50, 100, -100):
        ax.axhline(level, color="gray", linestyle="--", linewidth=0.7, alpha=0.6)
    for _, r in rdf.iterrows():
        va = "bottom" if r["res"] >= 0 else "top"
        ax.text(r["xco2"], r["res"], f"  {r['name']}", ha="center", va=va, fontsize=7.5)
    ax.set_xlabel(r"$X_{\mathrm{CO_2}}$", fontsize=12)
    ax.set_ylabel("Residual (%)", fontsize=12)
    ax.set_title(
        f"Quartz Solubility â€” Model vs Experiment Residuals\n"
        f"PerplexDEW ({dh_model}), T = 800Â°C  "
        f"(\u25cf\u202f=\u202f10 kbar, \u25a0\u202f=\u202f9 kbar)",
        fontsize=11,
    )
    ax.set_xlim(left=0.0)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"\n    Saved: {output_path}")
    plt.close()


# =============================================================================
# Main
# =============================================================================


def main():
    args = parse_args()
    dh_model = args.dh_model
    model_backend = backend_tag(args.backend)
    outputs = output_paths(MINERAL_CONFIG, model_backend, dh_model)

    print("=" * 80)
    print("Quartz Solubility in H2O-CO2  â€”  Reaktoro / PerplexDEW")
    print(f"Mineral   : {MINERAL_CONFIG['mineral_formula']}")
    print(f"Solute    : {MINERAL_CONFIG['solute_species']}")
    print(f"Backend   : {model_backend} ({dh_model})")
    print(f"Gas EOS   : H2O-CO2 water activity â€” Zhang-Duan 2009 (GFSM, separate pass)")
    print(f"Conditions: T = {T_FIXED_C} Â°C,  P = {MODEL_PRESSURES_KBAR} kbar")
    print(
        f"XCO2 grid : {len(XCO2_MODEL)} values from {XCO2_MODEL[0]:.3f} to {XCO2_MODEL[-1]:.3f}"
    )
    print(
        f"Si species: SiO2_aq, HSiO3-, Si2O4_aq, Si3O6_aq (no Si-C species in DEW2024)"
    )
    print("=" * 80)

    # --- Databases ---
    print("\n[1] Loading databases...")
    dew_db = DEWDatabase("dew2024-aqueous")
    supcrt_db = SupcrtDatabase("supcrtbl")
    print("    DEW2024 and SUPCRTBL databases loaded.")

    # --- Build DEW system (aqueous + mineral, no gas phase) ---
    print("\n[2] Building chemical system...")
    if model_backend.lower() == "perplexdew":
        ensure_perplexdew_symbols()
    system = build_system(
        dew_db,
        supcrt_db,
        MINERAL_CONFIG,
        model_backend=model_backend,
        dh_model=dh_model,
    )

    # --- Load experimental data ---
    print("\n[3] Loading experimental data...")
    exp_df = load_experimental_data(CSV_FILE)
    print(f"    {len(exp_df)} experimental data points loaded.")
    for ref, grp in exp_df.groupby("reference"):
        print(
            f"      {ref}: {len(grp)} pts, P={sorted(grp['P_kbar'].unique().tolist())} kbar, "
            f"XCO2=[{grp['xco2'].min():.3f}\u2013{grp['xco2'].max():.3f}], "
            f"m=[{grp['molality_m'].min():.4f}\u2013{grp['molality_m'].max():.4f}] mol/kg"
        )

    # --- Build GFSM gas system for a(H2O) computation ---
    print("\n[4] Building GFSM gas system for H2O activity (Zhang-Duan 2009)...")
    gfsm_system = build_gfsm_system(supcrt_db)
    print("    GFSM CO2(g)+H2O(g) system ready.")

    # --- Pure-H2O DEW baseline + activity-corrected curves per pressure ---
    print("\n[5] Computing quartz solubility vs XCO2...")
    model_results_flat = {}  # {P_kbar: (xco2_array, m_flat)}  -- horizontal line
    model_results_corr = {}  # {P_kbar: (xco2_array, m_corrected)}
    for P_kbar in MODEL_PRESSURES_KBAR:
        P_bar = P_kbar * 1000.0
        print(f"\n    P = {P_kbar} kbar ({P_bar:.0f} bar):")

        # Baseline: solve at XCO2 = 0 (pure H2O) â€” gives individual species mâ‚€
        print("      Solving DEW baseline (pure H2O)...")
        m0 = get_baseline_species_molalities(system, MINERAL_CONFIG, T_FIXED_C, P_bar)
        print(
            f"      mâ‚€(SiO2_aq)={m0['SiO2(aq)']:.4f}  mâ‚€(Si2O4_aq)={m0['Si2O4(aq)']:.4f}  "
            f"mâ‚€(Si3O6_aq)={m0['Si3O6(aq)']:.4e}  mâ‚€(HSiO3-)={m0['HSiO3-(aq)']:.4e}  "
            f"total_Si={m0['total_Si']:.4f} mol/kg"
        )

        # Flat curve: same mâ‚€ for all XCO2 (pure-H2O DEW reference)
        m_flat = np.full(len(XCO2_MODEL), m0["total_Si"])
        model_results_flat[P_kbar] = (XCO2_MODEL.copy(), m_flat)

        # Water activity from GFSM at each XCO2
        print("      Computing a(H2O) from GFSM...")
        aH2O = compute_aH2O_gfsm(gfsm_system, XCO2_MODEL, T_FIXED_C, P_bar)
        print(
            f"      a_H2O range: {aH2O.min():.4f} (XCO2={XCO2_MODEL[np.argmin(aH2O)]:.3f}) "
            f"to {aH2O.max():.4f} (XCO2={XCO2_MODEL[np.argmax(aH2O)]:.3f})"
        )

        # Activity-corrected curve
        xco2_corr, m_corr = compute_activity_corrected_curve(m0, aH2O, XCO2_MODEL)
        model_results_corr[P_kbar] = (xco2_corr, m_corr)
        print(
            f"      Corrected m range: {m_corr.min():.4f} â€“ {m_corr.max():.4f} mol/kg"
        )

    # --- Generate plots ---
    print("\n[6] Generating plots...")
    plot_solubility_vs_xco2(
        model_results_flat,
        model_results_corr,
        exp_df,
        MINERAL_CONFIG,
        outputs["solubility_plot"],
        dh_model,
    )
    plot_residuals(
        model_results_corr,
        exp_df,
        MINERAL_CONFIG,
        outputs["residuals_plot"],
        dh_model,
    )

    print("\nâœ“ Done.")


if __name__ == "__main__":
    main()

