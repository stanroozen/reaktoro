"""
Beginner Tutorial: Generic Mineral Solubility with Reaktoro

What this script does
- builds one aqueous phase plus one mineral phase
- equilibrates that system over a temperature range at several pressures
- reports the dissolved molality of one target element
- saves one plot in the same folder as this script

How to run this script
  Option 1 — from a terminal:
      cd <path-to-this-folder>
      python brucite_solubility_tutorial.py

  Option 2 — from VS Code:
      Open this file and press F5, or right-click and choose "Run Python File in Terminal".

  Requirements:
      - Python environment with Reaktoro installed (e.g. the reaktoro conda environment), OR
      - a local Reaktoro build at <repo-root>/build/Reaktoro/Release (detected automatically).

  Output:
      A PNG plot is saved in the same folder as this script.
      The filename is set by OUTPUT_PLOT_FILENAME in the user settings section below.

"""

# Import os so we can build file paths in a platform-safe way.
import os

# Import json to create a reduced tutorial-specific database subset.
import json

# Import sys so we can add the local Reaktoro build folder to Python's import path.
import sys

# Import NumPy for arrays and numeric helpers.
import numpy as np

# Import Matplotlib for plotting the final solubility curves.
import matplotlib.pyplot as plt


# -----------------------------------------------------------------------------
# 1) Import Reaktoro from the local build when available
# -----------------------------------------------------------------------------

# Store the absolute folder of this script.
THIS_DIR = os.path.dirname(os.path.abspath(__file__))

# Move up to the repository root from the Tutorial folder.
REPO_ROOT = os.path.dirname(os.path.dirname(THIS_DIR))

# Point to the local compiled Python extension folder in this repository.
LOCAL_BUILD_PYD_DIR = os.path.join(REPO_ROOT, "build", "Reaktoro", "Release")

# If the local compiled extension folder exists and is not already on sys.path,
# add it at the front so Python prefers the local build.
if os.path.isdir(LOCAL_BUILD_PYD_DIR) and LOCAL_BUILD_PYD_DIR not in sys.path:
    sys.path.insert(0, LOCAL_BUILD_PYD_DIR)

# Import autodiff after the local build path is on sys.path.
# The local reaktoro4py build may ship autodiff alongside the extension module.
try:
    import autodiff  # noqa: E402
except ModuleNotFoundError:
    # Fallback for environments where autodiff is not exposed as a module.
    class _AutodiffShim:
        @staticmethod
        def real(value):
            return value

    autodiff = _AutodiffShim()

try:
    # First try the local compiled extension from this repository.
    from reaktoro4py import *  # noqa: F401,F403

    # Print a message so the user knows the local build is being used.
    print(f"Using local Reaktoro build from: {LOCAL_BUILD_PYD_DIR}")
except ModuleNotFoundError:
    # If the local compiled extension is not available, fall back to an installed package.
    from reaktoro import *  # noqa: F401,F403

    # Print a message so the user knows the installed package is being used instead.
    print("Using installed 'reaktoro' package (local build not found).")


# -----------------------------------------------------------------------------
# 2) User input section
# Edit this block first. Most users only need to change values here.
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# 2A) Mineral and aqueous chemistry
# -----------------------------------------------------------------------------

# Choose the mineral name exactly as it appears in the DEW24HP622 PerpleX JSON database.
MINERAL_NAME = "br"

# Give the mineral a human-readable formula string for printed output only.
MINERAL_FORMULA = "Mg(OH)2"

# Choose which dissolved element should be plotted as the solubility signal.
TARGET_DISSOLVED_ELEMENT = "Mg"

# List the aqueous species that will be included in the aqueous phase.
# This is the main list to change when adapting the tutorial to another system.
AQUEOUS_SPECIES = [
    "H2O",
    "CO2",
    "CO2,aq",
    "H+",
    "OH-",
    "Mg+2",
    "MgOH+",
    "HCO3-",
    "CO3-2",
]

# Optionally set the solvent species explicitly.
# Leave this as None if you want the script to infer a likely solvent automatically.
SOLVENT_SPECIES_NAME = "H2O"


# -----------------------------------------------------------------------------
# 2B) Thermodynamic database and PerpleX-DEW model settings
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# PerplexDEW settable API explained (Python)
#
# 1) StandardThermoModelParamsPerplexDEW (species-level HKF-style parameters)
#    - Gf, Hf, Sr, a1, a2, a3, a4, c1, c2, wref, charge, Tmax
#    - Note: this tutorial uses species thermo directly from the selected
#      PerpleX-derived database, so these fields are not manually assigned here.
#
# 2) ActivityDHModel enum
#    - ActivityDHModel.Davies
#    - ActivityDHModel.ExtendedDH
#
# 3) ActivityModelParamsPerplexDEW
#    - dhModel
#    - errorOnConflictingStandardState
#    - warnOnUnmappedGFSMCoupling
#
# 4) Model factories used in this tutorial
#    - StandardThermoModelPerplexDEW(params)
#    - ActivityModelPerplexDEW(params) or ActivityModelPerplexDEW(model)
# -----------------------------------------------------------------------------

# Choose the PerpleX-derived Reaktoro JSON database file.
PERPLEX_DATABASE_FILENAME = "DEW24HP622ver_elements-reaktoro.json"

# Choose the Debye-Hückel variant used by ActivityModelPerplexDEW.
# Physical meaning: both models compute the electrostatic contribution to the mean
# activity coefficient of ions arising from long-range Coulomb screening by the
# surrounding ion cloud; the choice affects how strongly ionic strength depresses
# the activity of charged species and therefore mineral solubility.
# Options:
#   "Davies"     — Davies equation; no ionic-size parameter required; good to ~0.5 mol/kg.
#   "ExtendedDH" — extended Debye-Hückel with ion-size parameters; good to ~1 mol/kg.
PERPLEX_ACTIVITY_DH_MODEL = "Davies"

# Validate GFSM/HKF standard-state conflicts during model construction.
# Physical meaning: prevents combining incompatible reference states for the same
# species (mole-fraction GFSM solvent and molal HKF solute), which would double-count
# excess chemical potential contributions.
# Options: False (warn only, default), True (raise error and stop).
PERPLEX_ERROR_ON_CONFLICTING_STANDARD_STATE = False

# Emit a warning when a GFSM phase species cannot be matched to any Reaktoro species
# in the aqueous phase during PerplexDEW evaluation.
# Physical meaning: an unmapped GFSM coupling means the Perple_X fluid activity
# for that species is silently ignored; the species uses purely Reaktoro's internal
# activity model instead. This is expected when the PerpleX database contains species
# not listed in AQUEOUS_SPECIES, and benign if those species have negligible amounts.
# Options: True (print a one-time stderr warning per unmapped species, default),
#          False (suppress all unmapped-coupling warnings).
PERPLEX_WARN_ON_UNMAPPED_GFSM_COUPLING = True

# The settings below document the DEW water conventions used when generating
# this PerpleX-derived database.

# Water EOS used when the PerpleX-derived database was generated with the DEW spreadsheet.
# Physical meaning: sets the P-T-density relationship of water; determines volumetric
# and caloric properties that underpin all HKF standard-state calculations in the database.
# Options: "ZhangDuan2005" (DEW default), "ZhangDuan2009", "WagnerPruss", "HGK"
# NOTE: changing these values does NOT recompute the PerpleX database; they are
# documentation-only to record which convention the database was generated with.
DEW_WATER_EOS_MODEL = "ZhangDuan2005"

# Dielectric-constant model used when the database was generated.
# Physical meaning: the relative permittivity of water controls how strongly the
# solvent screens electrostatic interactions; it directly sets the Debye-Huckel A
# and B parameters and shifts the Born solvation energies of ions.
# Options: "PowerFunction" (DEW default), "JohnsonNorton1991", "Franck1990", "Fernandez1997"
DEW_WATER_DIELECTRIC_MODEL = "PowerFunction"

# Water Gibbs energy model used when the database was generated.
# Physical meaning: the Gibbs energy of liquid water at P-T enters the calculation
# of the Born coefficient omega, which corrects standard partial molar volumes
# and heat capacities of charged species for electrostatic solvation.
# Options: "DewIntegral" (DEW default), "DelaneyHelgeson1978"
DEW_WATER_GIBBS_MODEL = "DewIntegral"

# Born model used when the database was generated.
# Physical meaning: the Born model computes the electrostatic free energy of
# transferring an ion from vacuum into the high-P dielectric medium; omitting
# it ("None") removes all pressure dependence of ionic solvation energies.
# Options: "Shock92Dew" (DEW default), "None" (Born correction disabled)
DEW_WATER_BORN_MODEL = "Shock92Dew"

# Whether Psat polynomials were used when the database was generated.
# Physical meaning: Psat(T) locates the liquid-vapour boundary used to decide
# which branch of the water EOS applies at a given T; only relevant near saturation.
# Options: True (DEW default), False (iterative solver)
DEW_USE_PSAT_POLYNOMIALS = True

# Relative Psat solver tolerance used when the database was generated.
# Physical meaning: controls precision of the liquid-vapour boundary location;
# only significant for calculations very close to the saturation curve.
# Typical values: 1e-3 (loose, DEW default), 1e-6 (tight).
DEW_PSAT_RELATIVE_TOLERANCE = 1.0e-3


# -----------------------------------------------------------------------------
# 2C) Temperature-pressure sweep
# -----------------------------------------------------------------------------

# Choose the temperature range in degrees Celsius.
TEMPERATURE_MIN_C = 50.0
TEMPERATURE_MAX_C = 450.0

# Choose how many temperature points to compute.
NUMBER_OF_TEMPERATURE_POINTS = 80

# Choose the pressures to calculate, in kbar.
PRESSURES_KBAR = [1.0, 2.0, 5.0]


# -----------------------------------------------------------------------------
# 2D) Initial chemical state
# -----------------------------------------------------------------------------

# Set the initial amounts for the starting chemical state, in mol.
# Keep the mineral amount comfortably above saturation so the mineral can remain present.
# Keep the solvent amount large enough to represent the bulk fluid.
# Add or remove species here when adapting the script to another chemistry setup.
INITIAL_SPECIES_AMOUNTS_MOL = {
    "H2O": 27.75,
    "CO2": 27.75,
    "CO2,aq": 1e-8,
    "H+": 1e-8,
    "OH-": 1e-8,
    "Mg+2": 1e-8,
    "MgOH+": 1e-10,
    "HCO3-": 1e-10,
    "CO3-2": 1e-12,
    MINERAL_NAME: 10.0,
}


# -----------------------------------------------------------------------------
# 2E) Graph and output settings
# -----------------------------------------------------------------------------

# Give the tutorial run a short human-readable title.
RUN_TITLE = "Brucite solubility tutorial (DEW24HP622 PerpleX JSON, 50/50 H2O-CO2 fluid)"

# Choose the output figure name.
OUTPUT_PLOT_FILENAME = "brucite_solubility_tutorial_dew24hp622.png"

# Choose the plot title.
PLOT_TITLE = (
    "Brucite Solubility Tutorial (DEW24HP622 PerpleX JSON, 50/50 H2O-CO2 fluid)"
)

# Choose the x-axis label.
PLOT_X_LABEL = "Temperature (C)"

# Choose the y-axis label.
# This is placed after TARGET_DISSOLVED_ELEMENT so the label updates automatically when you change the element.
PLOT_Y_LABEL = f"Total dissolved {TARGET_DISSOLVED_ELEMENT} molality (mol/kg-H2O)"

# Set water molar mass (kg/mol) for manual molality conversion.
WATER_MOLAR_MASS_KG_PER_MOL = 0.01801528

# Choose the plot size in inches.
PLOT_FIGURE_SIZE = (9, 6)

# Choose the saved PNG resolution.
PLOT_DPI = 250

# Choose the y-axis scaling.
PLOT_Y_SCALE = "log"


# -----------------------------------------------------------------------------
# 3) Derived paths and small helper choices
# You usually do not need to edit anything below this line.
# -----------------------------------------------------------------------------

# Store the script folder again with a beginner-friendly name.
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# Build the absolute path to the output plot file.
OUTPUT_PLOT = os.path.join(SCRIPT_DIR, OUTPUT_PLOT_FILENAME)

# Build the absolute path to the PerpleX-derived Reaktoro database file in this repository.
PERPLEX_DATABASE_FILE = os.path.join(
    os.path.dirname(os.path.dirname(SCRIPT_DIR)),
    "embedded",
    "databases",
    "perplex",
    PERPLEX_DATABASE_FILENAME,
)

# Build a tutorial-local subset database path that contains only the selected species.
PERPLEX_DATABASE_SUBSET_FILE = os.path.join(
    SCRIPT_DIR,
    "DEW24HP622ver_elements-reaktoro-brucite-subset.json",
)


# -----------------------------------------------------------------------------
# 4) Helper functions
# These functions keep the main calculation logic easier to read.
# -----------------------------------------------------------------------------


# Define a helper that tries to determine the solvent species name automatically.
def infer_solvent_species_name(aqueous_species, configured_name, initial_amounts):
    """Return the solvent species name using explicit choice first, then safe guesses."""

    # If the user provided a solvent species explicitly, use it.
    if configured_name is not None:
        return configured_name

    # Prefer common water-like solvent names if they appear in the aqueous species list.
    preferred_names = ["H2O", "H2O(aq)", "H2O(l)"]

    # Check each preferred solvent name in order.
    for candidate_name in preferred_names:
        if candidate_name in aqueous_species:
            return candidate_name

    # If no common water name was found, try to use the aqueous species with the largest initial amount.
    # This makes the script more general for custom solvent naming.
    largest_amount = -1.0
    best_name = None

    # Loop through the explicit aqueous species only.
    for species_name in aqueous_species:
        amount = float(initial_amounts.get(species_name, 0.0))
        if amount > largest_amount:
            largest_amount = amount
            best_name = species_name

    # If we found a species with the largest amount, use it.
    if best_name is not None:
        return best_name

    # If we get here, the script cannot infer a solvent safely.
    raise ValueError(
        "Could not infer the solvent species name. "
        "Set SOLVENT_SPECIES_NAME explicitly or add the solvent to INITIAL_SPECIES_AMOUNTS_MOL."
    )


# Define a helper that checks the most important user inputs early.
def validate_user_inputs():
    """Raise a clear error if the beginner-facing configuration is inconsistent."""

    # Make sure the aqueous species list is not empty.
    if not AQUEOUS_SPECIES:
        raise ValueError("AQUEOUS_SPECIES must contain at least one aqueous species.")

    # Make sure at least two temperatures are requested so a curve can be drawn.
    if NUMBER_OF_TEMPERATURE_POINTS < 2:
        raise ValueError("NUMBER_OF_TEMPERATURE_POINTS must be at least 2.")

    # Make sure the temperature range is ordered correctly.
    if TEMPERATURE_MAX_C <= TEMPERATURE_MIN_C:
        raise ValueError("TEMPERATURE_MAX_C must be greater than TEMPERATURE_MIN_C.")

    # Make sure the mineral has an initial amount so it is present in the starting state.
    if MINERAL_NAME not in INITIAL_SPECIES_AMOUNTS_MOL:
        raise ValueError(
            "INITIAL_SPECIES_AMOUNTS_MOL must include the mineral phase name so the mineral is present initially."
        )


# Define a helper that prints the grouped user settings in a clear beginner-friendly way.
def print_run_configuration(solvent_species_name):
    """Print the most important tutorial settings in the same grouped order used at the top."""

    print("=" * 78)
    print(RUN_TITLE)
    print("=" * 78)
    print("[Mineral and aqueous chemistry]")
    print(f"Mineral phase: {MINERAL_NAME} ({MINERAL_FORMULA})")
    print("Aqueous species included: " + ", ".join(AQUEOUS_SPECIES))
    print(f"Solvent species used for interpretation: {solvent_species_name}")
    print(f"Target dissolved element: {TARGET_DISSOLVED_ELEMENT}")
    print()
    print("[Database and model settings]")
    print(f"PerpleX-derived database file: {PERPLEX_DATABASE_FILENAME}")
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
    print(f"  Activity DH mode: {PERPLEX_ACTIVITY_DH_MODEL}")
    print(
        "  Note: ActivityModelPerplexDEW is built from ActivityModelParamsPerplexDEW in this tutorial."
    )
    print()
    print("[Temperature-pressure sweep]")
    print(
        f"Temperature range: {TEMPERATURE_MIN_C} to {TEMPERATURE_MAX_C} C "
        f"with {NUMBER_OF_TEMPERATURE_POINTS} points"
    )
    print("Pressures (kbar): " + ", ".join(str(value) for value in PRESSURES_KBAR))
    print()
    print("[Initial chemical state in mol]")
    for species_name, amount in INITIAL_SPECIES_AMOUNTS_MOL.items():
        print(f"  {species_name}: {amount}")
    print("  Fluid composition target: 50 mol% H2O and 50 mol% CO2 (by initial moles)")
    print("=" * 78)


# Define a helper that builds the chemical system used for the equilibrium calculations.
def prepare_tutorial_database_subset():
    """Create a reduced database with only the species needed by this tutorial."""

    # Load the full converted PerpleX database.
    with open(PERPLEX_DATABASE_FILE, "r", encoding="utf-8") as file:
        database_data = json.load(file)

    # Collect the exact species required for the selected chemistry setup.
    required_species_names = set(AQUEOUS_SPECIES)
    required_species_names.add(MINERAL_NAME)

    # Access the species dictionary in the converted database.
    species_data = database_data.get("Species", {})

    # Validate that all required species exist in the source database.
    missing_species = sorted(
        species_name
        for species_name in required_species_names
        if species_name not in species_data
    )
    if missing_species:
        raise ValueError(
            "The following required species are missing in the source database: "
            + ", ".join(missing_species)
        )

    # Build a reduced species dictionary keeping only the required entries.
    reduced_species_data = {}
    for species_name in sorted(required_species_names):
        species_entry = dict(species_data[species_name])

        # AqueousProps expects proton formula to be H+ (or H3O+) explicitly.
        if species_name == "H+":
            species_entry["Formula"] = "H+"

        # Some converted DEW entries intentionally carry only ThermoReference values.
        # Add a constant fallback model with G0/H0/V0 so derived reporting
        # properties are not forced to zero.
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

    # Copy metadata and replace only the species collection.
    reduced_database_data = dict(database_data)
    reduced_database_data["Species"] = reduced_species_data

    # Write the reduced database next to this tutorial script.
    with open(PERPLEX_DATABASE_SUBSET_FILE, "w", encoding="utf-8") as file:
        json.dump(reduced_database_data, file, indent=2)


def build_tutorial_system():
    """Create the ChemicalSystem for the current tutorial configuration."""

    # Create a reduced tutorial database to avoid unrelated malformed entries.
    prepare_tutorial_database_subset()

    # Load the reduced PerpleX-derived database subset.
    database = Database.fromFile(PERPLEX_DATABASE_SUBSET_FILE)

    # Join the aqueous species names into the single space-separated format Reaktoro expects.
    aqueous_species_string = " ".join(AQUEOUS_SPECIES)

    # Create the aqueous phase from the chosen aqueous species.
    aqueous_phase = AqueousPhase(aqueous_species_string)

    # Configure PerplexDEW activity-model parameters explicitly.
    activity_params = ActivityModelParamsPerplexDEW()
    activity_params.dhModel = getattr(ActivityDHModel, PERPLEX_ACTIVITY_DH_MODEL)
    activity_params.errorOnConflictingStandardState = (
        PERPLEX_ERROR_ON_CONFLICTING_STANDARD_STATE
    )
    activity_params.warnOnUnmappedGFSMCoupling = PERPLEX_WARN_ON_UNMAPPED_GFSM_COUPLING

    # Set the activity model for the aqueous phase using explicit PerplexDEW options.
    aqueous_phase.setActivityModel(ActivityModelPerplexDEW(activity_params))

    # Create the mineral phase for the chosen mineral.
    mineral_phase = MineralPhase(MINERAL_NAME)

    # Return the full chemical system.
    return ChemicalSystem(database, aqueous_phase, mineral_phase)


# Define a helper that sets a species amount in a way that works with both local and installed builds.
def set_state_amount(state, species_name, amount, unit):
    """Set a species amount with compatibility for typed local reaktoro4py builds."""

    # First try the simple Python float path.
    try:
        state.set(species_name, float(amount), unit)

    # If the local build expects autodiff.real, convert and retry.
    except TypeError:
        state.set(species_name, autodiff.real(float(amount)), unit)


# Define a helper that creates the initial chemical state for one pressure sweep.
def make_initial_state(system):
    """Create and populate a ChemicalState from INITIAL_SPECIES_AMOUNTS_MOL."""

    # Start with an empty chemical state for the current system.
    state = ChemicalState(system)

    # Loop over every user-provided initial amount.
    for species_name, amount in INITIAL_SPECIES_AMOUNTS_MOL.items():
        # Set that amount in mol units.
        set_state_amount(state, species_name, amount, "mol")

    # Return the populated starting state.
    return state


# -----------------------------------------------------------------------------
# 5) Main calculation
# -----------------------------------------------------------------------------


# Define the main function so the file can be run as a script.
def main():
    """Run the full mineral-solubility tutorial workflow."""

    # Check the beginner-facing settings before doing any thermodynamics work.
    validate_user_inputs()

    # Determine the solvent species name once.
    solvent_species_name = infer_solvent_species_name(
        AQUEOUS_SPECIES,
        SOLVENT_SPECIES_NAME,
        INITIAL_SPECIES_AMOUNTS_MOL,
    )

    # Print the grouped run configuration so beginners can confirm everything before solving.
    print_run_configuration(solvent_species_name)

    # Build the chemical system from the selected user inputs.
    system = build_tutorial_system()

    # Create an equilibrium specification object for this system.
    specs = EquilibriumSpecs(system)

    # Declare temperature as an input variable.
    specs.temperature()

    # Declare pressure as an input variable.
    specs.pressure()

    # Create the equilibrium solver from those specs.
    solver = EquilibriumSolver(specs)

    # Create a conditions object that will receive temperature and pressure values.
    conditions = EquilibriumConditions(specs)

    # Create a uniformly spaced temperature array in degrees Celsius.
    temperatures_c = np.linspace(
        TEMPERATURE_MIN_C,
        TEMPERATURE_MAX_C,
        NUMBER_OF_TEMPERATURE_POINTS,
    )

    # Create a dictionary that will store one solubility curve for each pressure.
    curves = {}

    # Loop over each requested pressure in kbar.
    for pressure_kbar in PRESSURES_KBAR:
        # Convert kbar to bar because the equilibrium conditions use bar units.
        pressure_bar = pressure_kbar * 1000.0

        # Create an empty Python list to collect solubility values for this pressure.
        dissolved_element_molality = []

        # Create the starting chemical state for this pressure branch.
        state = make_initial_state(system)

        # Loop over every temperature on the temperature grid.
        for temperature_c in temperatures_c:
            # Set the current temperature in Celsius.
            conditions.temperature(float(temperature_c), "celsius")

            # Set the current pressure in bar.
            conditions.pressure(float(pressure_bar), "bar")

            # Solve equilibrium for the current state and conditions.
            result = solver.solve(state, conditions)

            # If the solver succeeded, read the dissolved element molality.
            if result.succeeded():
                # Read aqueous dissolved element moles directly from the phase properties.
                props = ChemicalProps(state)
                dissolved_element_moles = float(
                    props.elementAmountInPhase(TARGET_DISSOLVED_ELEMENT, "AqueousPhase")
                )

                # Convert solvent species amount to kilograms of H2O for molality.
                solvent_moles = float(state.speciesAmount(solvent_species_name))
                solvent_mass_kg = solvent_moles * WATER_MOLAR_MASS_KG_PER_MOL

                # Guard against zero or negative solvent mass.
                if solvent_mass_kg > 0.0:
                    element_molality = dissolved_element_moles / solvent_mass_kg
                else:
                    element_molality = np.nan

            # If the solver failed, store NaN so the plot skips that point.
            else:
                element_molality = np.nan

            # Add this temperature-point result to the current pressure curve.
            dissolved_element_molality.append(element_molality)

        # Convert the current pressure curve to a NumPy array and store it.
        curves[pressure_kbar] = np.array(dissolved_element_molality)

        # Count how many finite values were produced as a simple progress check.
        valid_count = np.isfinite(curves[pressure_kbar]).sum()

        # Print a short progress line for this pressure.
        print(
            f"Computed {valid_count}/{len(temperatures_c)} points at {pressure_kbar:.1f} kbar"
        )

    # Create the Matplotlib figure.
    plt.figure(figsize=PLOT_FIGURE_SIZE)

    # Loop over each pressure again to draw one curve per pressure.
    for pressure_kbar in PRESSURES_KBAR:
        # Load the stored y-values for this pressure.
        y = curves[pressure_kbar]

        # Build a mask that is True only where the calculation succeeded.
        valid = np.isfinite(y)

        # Plot temperature against dissolved-element molality for the valid points.
        plt.plot(
            temperatures_c[valid],
            y[valid],
            linewidth=2.0,
            label=f"{pressure_kbar:.1f} kbar",
        )

    # Apply the requested y-axis scaling.
    plt.yscale(PLOT_Y_SCALE)

    # Set the x-axis label.
    plt.xlabel(PLOT_X_LABEL)

    # Set the y-axis label.
    plt.ylabel(PLOT_Y_LABEL)

    # Set the plot title.
    plt.title(PLOT_TITLE)

    # Draw a light grid to help beginners read the figure.
    plt.grid(True, which="both", alpha=0.3, linestyle="--")

    # Draw the legend with pressure labels.
    plt.legend(title="Pressure")

    # Tighten layout so labels are not clipped.
    plt.tight_layout()

    # Save the figure as a PNG file.
    plt.savefig(OUTPUT_PLOT, dpi=PLOT_DPI)

    # Print the file path of the saved output plot.
    print(f"Saved plot: {OUTPUT_PLOT}")


# Run the main function only when the script is executed directly.
if __name__ == "__main__":
    # Start the tutorial run.
    main()
