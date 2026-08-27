"""
Beginner Tutorial: Willemite Solubility with DEW17HP622_Zn (PerpleX JSON)

What this script does
- builds one aqueous phase with either:
    - Willemite-only mineral mode, or
    - multi-mineral Zn competition mode for true stability-limited solubility
- computes Zn solubility versus temperature at 3 pressures
- computes Zn solubility sensitivity versus pH, silica chemical potential,
    sulfur fugacity, and fO2
- saves temperature, sensitivity, and true-fO2 plots in the same folder as this script

Important note on equilibrium controls
- pH-based sweep is imposed with equilibrium constraints (specs.pH + conditions.pH).
- Optional pH+charge mode is available (specs.pH + specs.charge + specs.openTo).
- The SiO2,aq chemical-potential sweep uses native specs.chemicalPotential when
    available, and otherwise falls back to a thermodynamically equivalent
    μ↔activity mapping at the chosen T,P reference point.
- SiO2,aq can also be constrained via native specs.chemicalPotential when the
  local binding supports species names with punctuation.
- fO2 and fH2S are imposed with real gas phases and specs.fugacity(...).
- For DEW17HP622_Zn_2025, H2O exists as aqueous/liquid species only (no H2O gas).
    Therefore, this tutorial does not build a full GFSM mixed-fluid gas phase.
"""

import json
import os
import re
import sys

import matplotlib.pyplot as plt
import numpy as np


# -----------------------------------------------------------------------------
# 1) Import Reaktoro from local build when available
# -----------------------------------------------------------------------------

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(os.path.dirname(THIS_DIR))
LOCAL_BUILD_PYD_DIR = os.path.join(
    REPO_ROOT, "build", "python", "package", "build", "lib", "reaktoro"
)
LOCAL_REAKTORO_RELEASE_DIR = os.path.join(REPO_ROOT, "build", "Reaktoro", "Release")
LOCAL_REAKTORO_DEBUG_DIR = os.path.join(REPO_ROOT, "build", "Reaktoro", "Debug")

_default_use_local_build = "0" if sys.platform.startswith("win") else "1"
USE_LOCAL_REAKTORO_BUILD = (
    os.environ.get("REAKTORO_USE_LOCAL_BUILD", _default_use_local_build) != "0"
)


def _normpath(path):
    return os.path.normcase(os.path.normpath(path))


if sys.platform.startswith("win") and os.path.isdir(LOCAL_REAKTORO_RELEASE_DIR):
    # On Windows, loading the Debug extension into a standard Python process can
    # trigger CRT heap assertions. Always prefer Release when available.
    if LOCAL_REAKTORO_RELEASE_DIR not in sys.path:
        sys.path.insert(0, LOCAL_REAKTORO_RELEASE_DIR)

    debug_dir_norm = _normpath(LOCAL_REAKTORO_DEBUG_DIR)
    sys.path = [p for p in sys.path if _normpath(p) != debug_dir_norm]

if (
    USE_LOCAL_REAKTORO_BUILD
    and os.path.isdir(LOCAL_BUILD_PYD_DIR)
    and LOCAL_BUILD_PYD_DIR not in sys.path
):
    # Keep package-build path lower priority than explicit Release folder.
    sys.path.append(LOCAL_BUILD_PYD_DIR)

try:
    import autodiff  # noqa: E402
except ModuleNotFoundError:

    class _AutodiffShim:
        @staticmethod
        def real(value):
            return value

    autodiff = _AutodiffShim()

try:
    import reaktoro4py as _reaktoro4py  # noqa: E402
    from reaktoro4py import *  # noqa: F401,F403

    print(
        "Using 'reaktoro4py' extension from: "
        f"{getattr(_reaktoro4py, '__file__', '<unknown>')} "
        f"(REAKTORO_USE_LOCAL_BUILD={USE_LOCAL_REAKTORO_BUILD})"
    )
except ModuleNotFoundError:
    from reaktoro import *  # noqa: F401,F403

    print("Using installed 'reaktoro' package (local build disabled or not found).")


# -----------------------------------------------------------------------------
# 2) User input section
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# 2A) Mineral and aqueous chemistry
# -----------------------------------------------------------------------------

MINERAL_NAME = "Wlm"
MINERAL_FORMULA = "Zn2SiO4"
TARGET_DISSOLVED_ELEMENT = "Zn"

# Mineral competition mode:
# - False: Willemite-only (metastable, constrained to Wlm as the sole Zn solid)
# - True: include competing Zn minerals to compute true stability-limited solubility
USE_COMPETING_ZN_MINERALS = True

# If enabled, auto-include Zn-free gangue solids from the database that are
# compositionally compatible with the initialized chemical system elements.
INCLUDE_COMPATIBLE_GANGUE_MINERALS = True

# Zn-bearing minerals available in DEW17HP622_Zn_2025 (from zn_species_summary.md).
# Keep Willemite first as the reference mineral in reporting.
COMPETING_ZN_MINERALS = [
    "Wlm",  # Willemite
    "Sph",  # Sphalerite
    "Smth",  # Smithsonite
    "Znc",  # Zincite
    "Znks",  # Zinkosite
    "Ghn",  # Gahnite
    "Frk",  # Franklinite
    "Hrds",  # Hardystonite
    "ZnSp",  # Zn-spinel
    "HZnc",  # Hydrozincite
    "Zn",  # Native zinc
    "Wrt",  # Wurtzite
    "Zn-St",  # Zn-stilpnomelane
]

# Populated automatically from the database when enabled.
AUTO_COMPATIBLE_GANGUE_MINERALS = []

# Regex for extracting chemical element symbols from formulas.
ELEMENT_SYMBOL_REGEX = re.compile(r"[A-Z][a-z]?")

AQUEOUS_SPECIES = [
    "H2O",
    "H+",
    "OH-",
    "Na+",
    # Zinc aqueous speciation (DEW17HP622_Zn_2025)
    "Zn2+",
    "ZnOH+",
    "ZnO",
    "ZnO2-2",
    "HZnO2-",
    "ZnCl+",
    "ZnCl2",
    "ZnCl3-",
    "ZnCl4-2",
    "ZnF+",
    "ZnHCO3+",
    "Zn(HS)2",
    "Zn(HS)2OH-",
    "Zn(HS)3-",
    "Zn(HS)4-2",
    # Silicon species
    "SiO2,aq",
    "HSiO3-",
    "Si2O4,aq",
    # Sulfur species
    "HS-",
    "SO3-2",
    "HSO4-",
    "SO4-2",
    # Carbonate (for ZnHCO3+ coupling)
    "HCO3-",
    "CO3-2",
    # Chloride and fluoride (for ZnCl+, ZnF+ coupling)
    "Cl-",
    "F-",
]

SOLVENT_SPECIES_NAME = "H2O"


# -----------------------------------------------------------------------------
# 2B) Thermodynamic database and PerpleX-DEW model settings
# -----------------------------------------------------------------------------

PERPLEX_DATABASE_FILENAME = "DEW17HP622_Zn_2025-reaktoro.json"
AQUEOUS_ACTIVITY_MODEL = ActivityModelPerplexDEW

DEW_WATER_EOS_MODEL = "ZhangDuan2005"
DEW_WATER_DIELECTRIC_MODEL = "PowerFunction"
DEW_WATER_GIBBS_MODEL = "DewIntegral"
DEW_WATER_BORN_MODEL = "Shock92Dew"
DEW_USE_PSAT_POLYNOMIALS = True
DEW_PSAT_RELATIVE_TOLERANCE = 1.0e-3


# -----------------------------------------------------------------------------
# 2C) Temperature-pressure sweep
# -----------------------------------------------------------------------------

TEMPERATURE_MIN_C = 50.0
TEMPERATURE_MAX_C = 500.0
NUMBER_OF_TEMPERATURE_POINTS = 90
PRESSURES_KBAR = [1.0, 2.0, 5.0]


# -----------------------------------------------------------------------------
# 2D) Initial chemical state
# -----------------------------------------------------------------------------

TRACE_ZN_AQUEOUS_SEED_MOL = 1.0e-20

INITIAL_SPECIES_AMOUNTS_MOL = {
    "H2O": 55.5,
    "H+": 1.0e-7,
    "OH-": 1.0e-7,
    "Na+": 1.0e-10,
    # Zinc
    "Zn2+": TRACE_ZN_AQUEOUS_SEED_MOL,
    "ZnOH+": TRACE_ZN_AQUEOUS_SEED_MOL,
    "ZnO": TRACE_ZN_AQUEOUS_SEED_MOL,
    "ZnO2-2": TRACE_ZN_AQUEOUS_SEED_MOL,
    "HZnO2-": TRACE_ZN_AQUEOUS_SEED_MOL,
    "ZnCl+": TRACE_ZN_AQUEOUS_SEED_MOL,
    "ZnCl2": TRACE_ZN_AQUEOUS_SEED_MOL,
    "ZnCl3-": TRACE_ZN_AQUEOUS_SEED_MOL,
    "ZnCl4-2": TRACE_ZN_AQUEOUS_SEED_MOL,
    "ZnF+": TRACE_ZN_AQUEOUS_SEED_MOL,
    "ZnHCO3+": TRACE_ZN_AQUEOUS_SEED_MOL,
    "Zn(HS)2": TRACE_ZN_AQUEOUS_SEED_MOL,
    "Zn(HS)2OH-": TRACE_ZN_AQUEOUS_SEED_MOL,
    "Zn(HS)3-": TRACE_ZN_AQUEOUS_SEED_MOL,
    "Zn(HS)4-2": TRACE_ZN_AQUEOUS_SEED_MOL,
    # Silicon
    "SiO2,aq": 1.0e-6,
    "HSiO3-": 1.0e-9,
    "Si2O4,aq": 1.0e-12,
    # Sulfur
    "HS-": 1.0e-10,
    "SO3-2": 1.0e-10,
    "HSO4-": 1.0e-10,
    "SO4-2": 1.0e-10,
    # Carbonate
    "HCO3-": 1.0e-10,
    "CO3-2": 1.0e-12,
    # Chloride and fluoride
    "Cl-": 1.0e-10,
    "F-": 1.0e-12,
    MINERAL_NAME: 10.0,
}

# Optional NaCl brine background used by all sweeps.
# With 55.5 mol H2O (~1 kg), setting 0.5 gives approximately 0.5 mol/kg NaCl.
USE_NACL_BRINE_BACKGROUND = True
NACL_BRINE_MOL_PER_KG_H2O = 0.5


# -----------------------------------------------------------------------------
# 2E) Sensitivity sweeps (single T,P reference point)
# -----------------------------------------------------------------------------

SENS_TEMPERATURE_C = 300.0
SENS_PRESSURE_KBAR = 2.0
SENS_POINTS = 31

PH_RANGE = np.linspace(3.0, 10.0, SENS_POINTS)
LOG_FO2_RANGE = np.linspace(-40.0, -10.0, SENS_POINTS)
LOG_FH2S_RANGE = np.linspace(-20.0, -6.0, SENS_POINTS)

# Optional pH+charge mode for aqueous electroneutral sensitivity calculations.
USE_PH_CHARGE_MODE = False
PH_TARGET_CHARGE_MOL = 0.0
PH_CHARGE_TITRANT_SPECIES = "Cl-"

# Chemical potential ranges (J/mol) for thermodynamic sensitivity sweeps.
# These are calibrated at the 300°C, 2 kbar reference point so that the
# converted activity windows remain in a numerically tractable range.
MU_SIO2_MIN_J_PER_MOL = -902079.0  # ≈ lg(a[SiO2,aq]) = -5.0 at 300°C, 2 kbar
MU_SIO2_MAX_J_PER_MOL = -874647.0  # ≈ lg(a[SiO2,aq]) = -2.5 at 300°C, 2 kbar

MU_SIO2_RANGE = np.linspace(MU_SIO2_MIN_J_PER_MOL, MU_SIO2_MAX_J_PER_MOL, SENS_POINTS)

# Species names used in sensitivity controls.
SIO2_PROXY_SPECIES = "SiO2,aq"
SULFUR_FUGACITY_SPECIES = "H2S"


# -----------------------------------------------------------------------------
# 2F) Plot and output settings
# -----------------------------------------------------------------------------

RUN_TITLE = "Willemite solubility tutorial (DEW17HP622_Zn PerpleX JSON)"

OUTPUT_T_CURVE_FILENAME = "willemite_solubility_vs_temperature_dew17hp622_zn.png"
OUTPUT_SENS_FILENAME = "willemite_solubility_sensitivity_dew17hp622_zn.png"
OUTPUT_TRUE_FO2_FILENAME = "willemite_solubility_true_fo2_dew17hp622_zn.png"

PLOT_DPI = 250
PLOT_Y_SCALE = "log"
WATER_MOLAR_MASS_KG_PER_MOL = 0.01801528


# -----------------------------------------------------------------------------
# 3) Derived paths
# -----------------------------------------------------------------------------

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PERPLEX_DATABASE_FILE = os.path.join(
    os.path.dirname(os.path.dirname(SCRIPT_DIR)),
    "embedded",
    "databases",
    "perplex",
    PERPLEX_DATABASE_FILENAME,
)

PERPLEX_DATABASE_SUBSET_FILE = os.path.join(
    SCRIPT_DIR,
    "DEW17HP622_Zn_2025-reaktoro-willemite-subset.json",
)

OUTPUT_T_CURVE = os.path.join(SCRIPT_DIR, OUTPUT_T_CURVE_FILENAME)
OUTPUT_SENS = os.path.join(SCRIPT_DIR, OUTPUT_SENS_FILENAME)
OUTPUT_TRUE_FO2 = os.path.join(SCRIPT_DIR, OUTPUT_TRUE_FO2_FILENAME)


# -----------------------------------------------------------------------------
# 4) Helpers
# -----------------------------------------------------------------------------


def infer_solvent_species_name(aqueous_species, configured_name, initial_amounts):
    if configured_name is not None:
        return configured_name

    for candidate_name in ["H2O", "H2O(aq)", "H2O(l)"]:
        if candidate_name in aqueous_species:
            return candidate_name

    largest_amount = -1.0
    best_name = None
    for species_name in aqueous_species:
        amount = float(initial_amounts.get(species_name, 0.0))
        if amount > largest_amount:
            largest_amount = amount
            best_name = species_name

    if best_name is None:
        raise ValueError(
            "Could not infer solvent species name. Set SOLVENT_SPECIES_NAME explicitly."
        )

    return best_name


def formula_element_symbols(formula):
    return set(ELEMENT_SYMBOL_REGEX.findall(str(formula)))


def discover_compatible_gangue_minerals(database_data):
    species_data = database_data.get("Species", {})

    # Determine elements effectively represented in the initialized system.
    initialized_elements = set()
    for species_name, amount in INITIAL_SPECIES_AMOUNTS_MOL.items():
        if float(amount) <= 0.0:
            continue
        entry = species_data.get(species_name)
        if entry is None:
            continue
        initialized_elements.update(formula_element_symbols(entry.get("Formula", "")))

    allowed_non_zn_elements = initialized_elements - {"Zn"}
    if not allowed_non_zn_elements:
        return []

    gangue = []
    for species_name, entry in species_data.items():
        aggregate_state = str(entry.get("AggregateState", "")).lower()
        if aggregate_state != "solid":
            continue

        element_symbols = formula_element_symbols(entry.get("Formula", ""))
        if not element_symbols or "Zn" in element_symbols:
            continue

        if element_symbols.issubset(allowed_non_zn_elements):
            gangue.append(species_name)

    return sorted(dict.fromkeys(gangue))


def validate_user_inputs():
    if not AQUEOUS_SPECIES:
        raise ValueError("AQUEOUS_SPECIES must contain at least one species.")

    if NUMBER_OF_TEMPERATURE_POINTS < 2:
        raise ValueError("NUMBER_OF_TEMPERATURE_POINTS must be at least 2.")

    if TEMPERATURE_MAX_C <= TEMPERATURE_MIN_C:
        raise ValueError("TEMPERATURE_MAX_C must be greater than TEMPERATURE_MIN_C.")

    if MINERAL_NAME not in INITIAL_SPECIES_AMOUNTS_MOL:
        raise ValueError(
            "INITIAL_SPECIES_AMOUNTS_MOL must include the mineral species."
        )

    if USE_COMPETING_ZN_MINERALS:
        if not COMPETING_ZN_MINERALS:
            raise ValueError(
                "COMPETING_ZN_MINERALS must contain at least one mineral when "
                "USE_COMPETING_ZN_MINERALS is enabled."
            )
        if MINERAL_NAME not in COMPETING_ZN_MINERALS:
            raise ValueError(
                "COMPETING_ZN_MINERALS must include MINERAL_NAME so the primary "
                "Willemite reference mineral remains represented."
            )

    if INCLUDE_COMPATIBLE_GANGUE_MINERALS and not os.path.isfile(PERPLEX_DATABASE_FILE):
        raise ValueError(f"Database file not found: {PERPLEX_DATABASE_FILE}")


def selected_mineral_names():
    global AUTO_COMPATIBLE_GANGUE_MINERALS

    if USE_COMPETING_ZN_MINERALS:
        # Deduplicate while preserving order.
        minerals = list(dict.fromkeys(COMPETING_ZN_MINERALS))
    else:
        minerals = [MINERAL_NAME]

    if INCLUDE_COMPATIBLE_GANGUE_MINERALS:
        if not AUTO_COMPATIBLE_GANGUE_MINERALS:
            with open(PERPLEX_DATABASE_FILE, "r", encoding="utf-8") as file:
                database_data = json.load(file)
            AUTO_COMPATIBLE_GANGUE_MINERALS = discover_compatible_gangue_minerals(
                database_data
            )
        minerals.extend(AUTO_COMPATIBLE_GANGUE_MINERALS)

    return list(dict.fromkeys(minerals))


def make_mineral_phases():
    minerals = selected_mineral_names()
    if len(minerals) == 1:
        return MineralPhase(minerals[0])

    return MineralPhases(StringList(minerals))


def print_run_configuration(solvent_species_name):
    print("=" * 78)
    print(RUN_TITLE)
    print("=" * 78)
    print("[Mineral and aqueous chemistry]")
    print(f"Primary mineral: {MINERAL_NAME} ({MINERAL_FORMULA})")
    minerals = selected_mineral_names()
    if len(minerals) == 1:
        print("Mineral competition mode: Willemite-only")
    else:
        print(
            "Mineral competition mode: multi-mineral Zn competition "
            f"({len(minerals)} minerals)"
        )
        print("Competing minerals: " + ", ".join(minerals))
    if INCLUDE_COMPATIBLE_GANGUE_MINERALS:
        print(
            "Auto-included compatible gangue minerals: "
            + (
                ", ".join(AUTO_COMPATIBLE_GANGUE_MINERALS)
                if AUTO_COMPATIBLE_GANGUE_MINERALS
                else "none"
            )
        )
    print("Aqueous species included: " + ", ".join(AQUEOUS_SPECIES))
    print(f"Solvent species used for molality: {solvent_species_name}")
    print(f"Target dissolved element: {TARGET_DISSOLVED_ELEMENT}")
    print("Gas-phase fugacity controls: O2 and H2S fugacity sweeps enabled")
    print()
    print("[Database and model settings]")
    print(f"Database file: {PERPLEX_DATABASE_FILENAME}")
    print(
        "Tutorial subset database file: "
        f"{os.path.basename(PERPLEX_DATABASE_SUBSET_FILE)}"
    )
    print("PerplexDEW water convention metadata:")
    print(f"  Water EOS model: {DEW_WATER_EOS_MODEL}")
    print(f"  Water dielectric model: {DEW_WATER_DIELECTRIC_MODEL}")
    print(f"  Water Gibbs model: {DEW_WATER_GIBBS_MODEL}")
    print(f"  Water Born model: {DEW_WATER_BORN_MODEL}")
    print(f"  Use Psat polynomials: {DEW_USE_PSAT_POLYNOMIALS}")
    print(f"  Psat relative tolerance: {DEW_PSAT_RELATIVE_TOLERANCE}")
    print()
    print("[Temperature-pressure sweep]")
    print(
        f"Temperature range: {TEMPERATURE_MIN_C} to {TEMPERATURE_MAX_C} C "
        f"with {NUMBER_OF_TEMPERATURE_POINTS} points"
    )
    print("Pressures (kbar): " + ", ".join(str(value) for value in PRESSURES_KBAR))
    print()
    print("[Sensitivity sweep reference point]")
    print(f"T = {SENS_TEMPERATURE_C} C, P = {SENS_PRESSURE_KBAR} kbar")
    if USE_PH_CHARGE_MODE:
        print(
            "pH mode: pH + charge balance "
            f"(target charge={PH_TARGET_CHARGE_MOL}, titrant={PH_CHARGE_TITRANT_SPECIES})"
        )
    else:
        print("pH mode: pH-only")
    if USE_NACL_BRINE_BACKGROUND:
        print(f"NaCl brine background: {NACL_BRINE_MOL_PER_KG_H2O} mol/kg-H2O")
    else:
        print("NaCl brine background: disabled")
    print("=" * 78)


def prepare_tutorial_database_subset():
    with open(PERPLEX_DATABASE_FILE, "r", encoding="utf-8") as file:
        database_data = json.load(file)

    required_species_names = set(AQUEOUS_SPECIES)
    required_species_names.update(selected_mineral_names())

    species_data = database_data.get("Species", {})

    missing_species = sorted(
        species_name
        for species_name in required_species_names
        if species_name not in species_data
    )
    if missing_species:
        raise ValueError(
            "Missing required species in source database: " + ", ".join(missing_species)
        )

    reduced_species_data = {}
    for species_name in sorted(required_species_names):
        species_entry = dict(species_data[species_name])

        # Keep proton formula explicit for AqueousProps compatibility.
        if species_name == "H+":
            species_entry["Formula"] = "H+"

        # In PerpleX mixed DEW databases, H2O can be tagged as Gas/GFSM.
        # For aqueous modeling in this tutorial, force this entry to Aqueous.
        if species_name == "H2O":
            species_entry["AggregateState"] = "Aqueous"

        has_formation_reaction = "FormationReaction" in species_entry
        has_standard_model = "StandardThermoModel" in species_entry
        if not has_formation_reaction and not has_standard_model:

            def _as_float(value, default=0.0):
                try:
                    return float(value)
                except (TypeError, ValueError):
                    return default

            thermo_reference = species_entry.get("ThermoReference", {})
            metadata = species_entry.get("Metadata", {})
            perplex_params = {}
            if isinstance(metadata, dict):
                params = metadata.get("PerpleX_Params", {})
                if isinstance(params, dict):
                    perplex_params = params

            reference_gibbs = _as_float(
                thermo_reference.get("Gf", 0.0)
                if isinstance(thermo_reference, dict)
                else 0.0
            )
            reference_entropy = _as_float(
                thermo_reference.get("S0", perplex_params.get("S0", 0.0))
                if isinstance(thermo_reference, dict)
                else perplex_params.get("S0", 0.0)
            )
            # Perple_X V0 is in J/bar; convert to m3/mol with factor 1e-5.
            reference_v0_jbar = _as_float(
                thermo_reference.get("V0", perplex_params.get("V0", 0.0))
                if isinstance(thermo_reference, dict)
                else perplex_params.get("V0", 0.0)
            )

            reference_h0 = reference_gibbs + 298.15 * reference_entropy
            reference_v0 = reference_v0_jbar * 1.0e-5
            species_entry["StandardThermoModel"] = {
                "Constant": {
                    "G0": reference_gibbs,
                    "H0": reference_h0,
                    "V0": reference_v0,
                }
            }

        reduced_species_data[species_name] = species_entry

    reduced_database_data = dict(database_data)
    reduced_database_data["Species"] = reduced_species_data

    with open(PERPLEX_DATABASE_SUBSET_FILE, "w", encoding="utf-8") as file:
        json.dump(reduced_database_data, file, indent=2)


def build_tutorial_system():
    prepare_tutorial_database_subset()
    database = Database.fromFile(PERPLEX_DATABASE_SUBSET_FILE)

    aqueous_phase = AqueousPhase(" ".join(AQUEOUS_SPECIES))
    aqueous_phase.setActivityModel(AQUEOUS_ACTIVITY_MODEL())

    mineral_phase = make_mineral_phases()

    return ChemicalSystem(database, aqueous_phase, mineral_phase)


def build_system_with_o2_gas_phase():
    """Build a chemical system with the Zn aqueous phase, willemite mineral,
    and a single-species O2 gas phase using ActivityModelIdealGas.

    This enables a thermodynamically rigorous oxygen-fugacity constraint via
    specs.fugacity("O2") without requiring the full GFSM mixed-fluid setup.
    The O2 species in the DEW17HP622_Zn database is a proper gas-phase entry
    (PerplexGFSM index 7, AggregateState=Gas).  Using ActivityModelIdealGas
    for the single-species gas phase is consistent with fO2 as an intensive
    variable: the activity of pure O2(g) equals its fugacity in bar at the
    standard-state reference pressure of 1 bar.
    """
    database = Database.fromFile(PERPLEX_DATABASE_FILE)

    aqueous_phase = AqueousPhase(" ".join(AQUEOUS_SPECIES))
    aqueous_phase.setActivityModel(AQUEOUS_ACTIVITY_MODEL())

    mineral_phase = make_mineral_phases()

    # Single-species gas phase for O2: ideal-gas activity avoids GFSM H2O
    # requirement while correctly representing the O2 chemical potential.
    o2_gas_phase = GaseousPhase("O2")
    o2_gas_phase.setActivityModel(ActivityModelIdealGas())

    return ChemicalSystem(database, aqueous_phase, mineral_phase, o2_gas_phase)


def build_system_with_h2s_gas_phase():
    """Build a chemical system with Zn aqueous phase, willemite mineral,
    and a single-species H2S gas phase for native sulfur fugacity constraints.
    """
    database = Database.fromFile(PERPLEX_DATABASE_FILE)

    aqueous_phase = AqueousPhase(" ".join(AQUEOUS_SPECIES))
    aqueous_phase.setActivityModel(AQUEOUS_ACTIVITY_MODEL())

    mineral_phase = make_mineral_phases()

    h2s_gas_phase = GaseousPhase("H2S")
    h2s_gas_phase.setActivityModel(ActivityModelIdealGas())

    return ChemicalSystem(database, aqueous_phase, mineral_phase, h2s_gas_phase)


def species_index_map(system):
    return {system.species(i).name(): i for i in range(system.species().size())}


def apply_species_amount_overrides(state, amount_overrides_mol):
    index_map = species_index_map(state.system())
    # Use vector-based assignment, which accepts NumPy arrays in this binding.
    amounts = np.full(state.system().species().size(), 1.0e-16, dtype=float)

    for species_name, amount in amount_overrides_mol.items():
        idx = index_map.get(species_name)
        if idx is not None:
            amounts[idx] = float(amount)

    state.setSpeciesAmounts(amounts)


def make_base_state(system, include_gas_seed=False):
    state = ChemicalState(system)

    amounts_mol = dict(INITIAL_SPECIES_AMOUNTS_MOL)

    if USE_NACL_BRINE_BACKGROUND:
        water_moles = float(INITIAL_SPECIES_AMOUNTS_MOL.get("H2O", 0.0))
        water_mass_kg = water_moles * WATER_MOLAR_MASS_KG_PER_MOL
        nacl_moles = NACL_BRINE_MOL_PER_KG_H2O * water_mass_kg
        amounts_mol["Na+"] = float(amounts_mol.get("Na+", 0.0)) + nacl_moles
        amounts_mol["Cl-"] = float(amounts_mol.get("Cl-", 0.0)) + nacl_moles

    apply_species_amount_overrides(state, amounts_mol)

    return state


def dissolved_element_molality(state, solvent_species_name):
    props = ChemicalProps(state)
    dissolved_element_moles = float(
        props.elementAmountInPhase(TARGET_DISSOLVED_ELEMENT, "AqueousPhase")
    )

    solvent_moles = float(state.speciesAmount(solvent_species_name))
    solvent_mass_kg = solvent_moles * WATER_MOLAR_MASS_KG_PER_MOL
    if solvent_mass_kg <= 0.0:
        return np.nan

    return dissolved_element_moles / solvent_mass_kg


def make_equilibrium_solver(system, specs):
    try:
        return EquilibriumSolver(specs)
    except Exception as exc:
        # Some local bindings/builds can fail constructing solver from specs.
        # Fallback keeps the same solve(state, conditions) workflow.
        print(
            f"Warning: EquilibriumSolver(specs) failed ({exc}); falling back to EquilibriumSolver(system)."
        )
        return EquilibriumSolver(system)


def compute_temperature_curves(system, solvent_species_name):
    specs = EquilibriumSpecs(system)
    specs.temperature()
    specs.pressure()

    solver = make_equilibrium_solver(system, specs)
    conditions = EquilibriumConditions(specs)

    temperatures_c = np.linspace(
        TEMPERATURE_MIN_C,
        TEMPERATURE_MAX_C,
        NUMBER_OF_TEMPERATURE_POINTS,
    )

    curves = {}
    for pressure_kbar in PRESSURES_KBAR:
        pressure_bar = pressure_kbar * 1000.0
        values = []
        state = make_base_state(system)

        for temperature_c in temperatures_c:
            conditions.temperature(float(temperature_c), "celsius")
            conditions.pressure(float(pressure_bar), "bar")

            result = solver.solve(state, conditions)
            if result.succeeded():
                values.append(dissolved_element_molality(state, solvent_species_name))
            else:
                values.append(np.nan)

        curves[pressure_kbar] = np.array(values)
        valid_count = np.isfinite(curves[pressure_kbar]).sum()
        print(
            f"Temperature sweep: {valid_count}/{len(temperatures_c)} points at {pressure_kbar:.1f} kbar"
        )

    return temperatures_c, curves


def compute_pH_sensitivity(system, solvent_species_name):
    specs = EquilibriumSpecs(system)
    specs.temperature()
    specs.pressure()
    specs.pH()
    if USE_PH_CHARGE_MODE:
        specs.charge()
        specs.openTo(PH_CHARGE_TITRANT_SPECIES)

    solver = make_equilibrium_solver(system, specs)
    conditions = EquilibriumConditions(specs)

    pressure_bar = SENS_PRESSURE_KBAR * 1000.0
    values = []

    for pH_value in PH_RANGE:
        state = make_base_state(system)

        conditions.temperature(float(SENS_TEMPERATURE_C), "celsius")
        conditions.pressure(float(pressure_bar), "bar")
        conditions.pH(float(pH_value))
        if USE_PH_CHARGE_MODE:
            conditions.charge(float(PH_TARGET_CHARGE_MOL))

        result = solver.solve(state, conditions)
        if result.succeeded():
            values.append(dissolved_element_molality(state, solvent_species_name))
        else:
            values.append(np.nan)

    return PH_RANGE, np.array(values)


def compute_chemical_potential_sensitivity(
    system, solvent_species_name, constrained_species_name, mu_range_j_per_mol
):
    """Compute solubility vs. constrained chemical potential (thermodynamically equivalent).

    Reaktoro API limitation: chemicalPotential() doesn't accept species names with special
    characters (e.g., commas in "SiO2,aq"). Instead, we use the thermodynamic equivalence:

        μ(species) = μ°(species) + RT·ln(a[species])

    By constraining activity (lgActivity), we are indirectly constraining chemical potential.
    However, for rigorous control over absolute μ, we map the μ values to activity values:

        a[species] = exp((μ - μ°) / RT)

    Reference state (μ°) is computed at the equilibrium conditions (T, P).

    Args:
        system: ChemicalSystem
        solvent_species_name: e.g., "H2O"
        constrained_species_name: e.g., "SiO2,aq" or "HS-"
        mu_range_j_per_mol: numpy array of chemical potential values (J/mol)

    Returns:
        (mu_range, solubility_values): tuple of arrays
    """
    # Convert μ range to activity range via the relationship μ = μ° + RT·ln(a).
    # Only μ° at the chosen T,P reference point is needed for this mapping.
    R = universalGasConstant  # J/(mol·K)
    T_sens = float(SENS_TEMPERATURE_C) + 273.15  # Kelvin

    # Find the species and get its standard chemical potential
    species_idx = None
    for i in range(system.species().size()):
        if system.species(i).name() == constrained_species_name:
            species_idx = i
            break

    if species_idx is None:
        raise ValueError(f"Species '{constrained_species_name}' not found in system")

    # Get standard Gibbs energy of formation at the sensitivity reference point.
    reference_pressure_pa = float(SENS_PRESSURE_KBAR) * 1.0e8
    species_obj = system.species(species_idx)
    G0 = species_obj.props(T_sens, reference_pressure_pa).G0

    # Convert chemical potential values to lgActivity values: a = exp((μ - G0) / RT)
    # where μ is the chemical potential of the current state at equilibrium
    # and G0 is the standard state chemical potential at that T, P
    # So: lg(a) = (μ - G0) / (RT·ln(10))

    lg_activity_range = (mu_range_j_per_mol - G0) / (R * T_sens * np.log(10.0))

    # Now use standard lgActivity constraint with these converted values
    specs = EquilibriumSpecs(system)
    specs.temperature()
    specs.pressure()
    specs.lgActivity(constrained_species_name)

    solver = make_equilibrium_solver(system, specs)
    conditions = EquilibriumConditions(specs)

    pressure_bar = SENS_PRESSURE_KBAR * 1000.0
    values = []

    for lg_activity in lg_activity_range:
        state = make_base_state(system)

        conditions.temperature(float(SENS_TEMPERATURE_C), "celsius")
        conditions.pressure(float(pressure_bar), "bar")
        conditions.lgActivity(constrained_species_name, float(lg_activity))

        result = solver.solve(state, conditions)
        if result.succeeded():
            values.append(dissolved_element_molality(state, solvent_species_name))
        else:
            values.append(np.nan)

    return (mu_range_j_per_mol, np.array(values))


def compute_silica_potential_sensitivity(
    system, solvent_species_name, mu_range_j_per_mol
):
    """Compute SiO2,aq potential sensitivity.

    The function first tries native Reaktoro chemicalPotential constraints.
    If the binding cannot parse species names with punctuation (e.g., "SiO2,aq"),
    it falls back to the thermodynamically equivalent activity-mapped approach.
    """
    specs = EquilibriumSpecs(system)
    specs.temperature()
    specs.pressure()
    try:
        specs.chemicalPotential(SIO2_PROXY_SPECIES)
    except Exception as exc:
        print(
            "Native silica chemicalPotential constraint unavailable; "
            f"falling back to μ↔activity mapping ({exc})."
        )
        return compute_chemical_potential_sensitivity(
            system,
            solvent_species_name,
            SIO2_PROXY_SPECIES,
            mu_range_j_per_mol,
        )

    solver = make_equilibrium_solver(system, specs)
    conditions = EquilibriumConditions(specs)

    pressure_bar = SENS_PRESSURE_KBAR * 1000.0
    values = []

    for mu_value in mu_range_j_per_mol:
        state = make_base_state(system)

        conditions.temperature(float(SENS_TEMPERATURE_C), "celsius")
        conditions.pressure(float(pressure_bar), "bar")
        conditions.chemicalPotential(SIO2_PROXY_SPECIES, float(mu_value))

        result = solver.solve(state, conditions)
        if result.succeeded():
            values.append(dissolved_element_molality(state, solvent_species_name))
        else:
            values.append(np.nan)

    print("Using native specs.chemicalPotential for SiO2,aq.")
    return (mu_range_j_per_mol, np.array(values))


def compute_true_fo2_sensitivity(system_with_o2_gas, solvent_species_name):
    """Rigorous fO2 sweep using specs.fugacity('O2') with a real O2 gas phase.

    This is thermodynamically equivalent to the Reaktoro tutorial at
    https://reaktoro.org/tutorials/equilibrium/equilibrium-with-fixed-fugacity.html
    The system must have been built with build_system_with_o2_gas_phase().
    """
    specs = EquilibriumSpecs(system_with_o2_gas)
    specs.temperature()
    specs.pressure()
    specs.fugacity("O2")

    solver = make_equilibrium_solver(system_with_o2_gas, specs)
    conditions = EquilibriumConditions(specs)

    pressure_bar = SENS_PRESSURE_KBAR * 1000.0
    values = []

    for logf in LOG_FO2_RANGE:
        state = make_base_state(system_with_o2_gas)
        # Seed a tiny O2 gas amount so the gas phase is non-empty initially.
        try:
            state.set("O2", 1.0e-20, "mol")
        except Exception:
            pass

        conditions.temperature(float(SENS_TEMPERATURE_C), "celsius")
        conditions.pressure(float(pressure_bar), "bar")
        conditions.fugacity("O2", float(10.0**logf), "bar")

        result = solver.solve(state, conditions)
        if result.succeeded():
            values.append(dissolved_element_molality(state, solvent_species_name))
        else:
            values.append(np.nan)

    return (LOG_FO2_RANGE, np.array(values))


def compute_true_fh2s_sensitivity(system_with_h2s_gas, solvent_species_name):
    """Rigorous sulfur-fugacity sweep using specs.fugacity('H2S')."""
    specs = EquilibriumSpecs(system_with_h2s_gas)
    specs.temperature()
    specs.pressure()
    specs.fugacity(SULFUR_FUGACITY_SPECIES)

    solver = make_equilibrium_solver(system_with_h2s_gas, specs)
    conditions = EquilibriumConditions(specs)

    pressure_bar = SENS_PRESSURE_KBAR * 1000.0
    values = []

    for logf in LOG_FH2S_RANGE:
        state = make_base_state(system_with_h2s_gas)
        try:
            state.set(SULFUR_FUGACITY_SPECIES, 1.0e-20, "mol")
        except Exception:
            pass

        conditions.temperature(float(SENS_TEMPERATURE_C), "celsius")
        conditions.pressure(float(pressure_bar), "bar")
        conditions.fugacity(SULFUR_FUGACITY_SPECIES, float(10.0**logf), "bar")

        result = solver.solve(state, conditions)
        if result.succeeded():
            values.append(dissolved_element_molality(state, solvent_species_name))
        else:
            values.append(np.nan)

    return (LOG_FH2S_RANGE, np.array(values))


def save_true_fo2_plot(logfO2_x, fo2_y):
    fig, ax = plt.subplots(figsize=(8, 5.5))

    valid = np.isfinite(fo2_y) & (fo2_y > 0.0)
    valid_count = int(valid.sum())
    ax.plot(logfO2_x[valid], fo2_y[valid], color="tab:red", linewidth=2.0)
    ax.set_yscale("log")
    ax.set_xlabel("log10(fO2 / bar)")
    ax.set_ylabel(f"{TARGET_DISSOLVED_ELEMENT} molality (mol/kg-H2O)")
    ax.set_title(
        f"Zn Solubility vs. True O2 Fugacity\n"
        f"T={SENS_TEMPERATURE_C:.0f} °C, P={SENS_PRESSURE_KBAR:.1f} kbar  "
        f"({valid_count}/{len(logfO2_x)} points converged)"
    )
    ax.grid(True, which="both", alpha=0.3, linestyle="--")

    fig.text(
        0.01,
        0.01,
        "Constraint: specs.fugacity('O2') with GaseousPhase('O2') + ActivityModelIdealGas()\n"
        "Thermodynamically rigorous: µ(O2,aq) = µ°(O2,g) + RT·ln(fO2)",
        fontsize=8,
        color="gray",
    )

    fig.tight_layout(rect=(0.0, 0.06, 1.0, 1.0))
    fig.savefig(OUTPUT_TRUE_FO2, dpi=PLOT_DPI)
    plt.close(fig)

    print(
        f"Saved true fO2 plot ({valid_count}/{len(logfO2_x)} points): {OUTPUT_TRUE_FO2}"
    )


def save_temperature_plot(temperatures_c, curves):
    plt.figure(figsize=(9, 6))

    for pressure_kbar in PRESSURES_KBAR:
        y = curves[pressure_kbar]
        valid = np.isfinite(y) & (y > 0.0)
        plt.plot(
            temperatures_c[valid],
            y[valid],
            linewidth=2.0,
            label=f"{pressure_kbar:.1f} kbar",
        )

    plt.yscale(PLOT_Y_SCALE)
    plt.xlabel("Temperature (C)")
    plt.ylabel(f"Total dissolved {TARGET_DISSOLVED_ELEMENT} molality (mol/kg-H2O)")
    if len(selected_mineral_names()) == 1:
        title_prefix = "Willemite"
    else:
        title_prefix = "Zn (with competing minerals)"
    plt.title(f"{title_prefix} Solubility vs Temperature (DEW17HP622_Zn)")
    plt.grid(True, which="both", alpha=0.3, linestyle="--")
    plt.legend(title="Pressure")
    plt.tight_layout()
    plt.savefig(OUTPUT_T_CURVE, dpi=PLOT_DPI)
    plt.close()

    print(f"Saved temperature plot: {OUTPUT_T_CURVE}")


def save_sensitivity_plot(
    pH_x,
    pH_y,
    mu_sio2_x,
    mu_sio2_y,
    fH2S_x,
    fH2S_y,
    fO2_x,
    fO2_y,
):
    """Save sensitivity plots for thermodynamically rigorous constraints.

    Args:
        pH_x, pH_y: pH sweep data
        mu_sio2_x, mu_sio2_y: Chemical potential of SiO2,aq (J/mol) vs solubility
        fH2S_x, fH2S_y: H2S fugacity sweep data
        fO2_x, fO2_y: O2 fugacity sweep data
    """
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    axes = axes.ravel()

    series = [
        (
            axes[0],
            pH_x,
            pH_y,
            "pH",
            "linear",
            "Sensitivity to pH",
        ),
        (
            axes[1],
            mu_sio2_x / 1000.0,  # Convert J/mol to kJ/mol for readability
            mu_sio2_y,
            "μ(SiO2,aq) / kJ·mol⁻¹",
            "log",
            "Sensitivity to SiO2,aq chemical potential (thermodynamic)",
        ),
        (
            axes[2],
            fH2S_x,
            fH2S_y,
            "log10(fH2S / bar)",
            "log",
            "Sensitivity to sulfur fugacity (H2S master variable)",
        ),
        (
            axes[3],
            fO2_x,
            fO2_y,
            "log10(fO2 / bar)",
            "log",
            "Sensitivity to O2 fugacity (redox state)",
        ),
    ]

    for axis, x, y, xlabel, yscale, title in series:
        valid = np.isfinite(y) & (y > 0.0)
        axis.plot(x[valid], y[valid], color="tab:blue", linewidth=2.0)
        if yscale == "log":
            axis.set_yscale("log")
        axis.set_xlabel(xlabel)
        axis.set_ylabel(f"{TARGET_DISSOLVED_ELEMENT} molality (mol/kg-H2O)")
        axis.set_title(title)
        axis.grid(True, which="both", alpha=0.3, linestyle="--")

    fig.suptitle(
        "Zn Solubility Sensitivity (Thermodynamically Rigorous)\n"
        "DEW17HP622_Zn Database",
        fontsize=14,
    )
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.96))
    fig.savefig(OUTPUT_SENS, dpi=PLOT_DPI)
    plt.close(fig)

    print(f"Saved sensitivity plot: {OUTPUT_SENS}")


# -----------------------------------------------------------------------------
# 5) Main
# -----------------------------------------------------------------------------


def main():
    validate_user_inputs()

    try:
        Warnings.disable(906)
    except Exception:
        pass

    solvent_species_name = infer_solvent_species_name(
        AQUEOUS_SPECIES,
        SOLVENT_SPECIES_NAME,
        INITIAL_SPECIES_AMOUNTS_MOL,
    )

    print_run_configuration(solvent_species_name)

    system_no_gas = build_tutorial_system()
    temperatures_c, curves = compute_temperature_curves(
        system_no_gas, solvent_species_name
    )
    save_temperature_plot(temperatures_c, curves)

    pH_x, pH_y = compute_pH_sensitivity(system_no_gas, solvent_species_name)

    # Thermodynamically rigorous chemical potential constraints
    print("\nRunning SiO2,aq chemical potential sweep...")
    mu_sio2_x, mu_sio2_y = compute_silica_potential_sensitivity(
        system_no_gas,
        solvent_species_name,
        MU_SIO2_RANGE,
    )
    print(
        f"SiO2 sweep: {np.isfinite(mu_sio2_y).sum()}/{len(MU_SIO2_RANGE)} points converged"
    )

    print("Running rigorous sulfur fugacity sweep (H2S master variable)...")
    system_with_h2s_gas = build_system_with_h2s_gas_phase()
    fH2S_x, fH2S_y = compute_true_fh2s_sensitivity(
        system_with_h2s_gas,
        solvent_species_name,
    )
    print(
        f"H2S fugacity sweep: {np.isfinite(fH2S_y).sum()}/{len(LOG_FH2S_RANGE)} points converged"
    )

    print("Running rigorous fO2 sweep (true fugacity constraint)...")
    system_with_o2_gas = build_system_with_o2_gas_phase()
    fO2_x, fO2_y = compute_true_fo2_sensitivity(
        system_with_o2_gas,
        solvent_species_name,
    )

    save_sensitivity_plot(
        pH_x,
        pH_y,
        mu_sio2_x,
        mu_sio2_y,
        fH2S_x,
        fH2S_y,
        fO2_x,
        fO2_y,
    )

    try:
        save_true_fo2_plot(fO2_x, fO2_y)
    except Exception as exc:
        print(f"Skipping true fO2 plot due to: {exc}")


if __name__ == "__main__":
    main()
