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
import autodiff  # noqa: E402

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

# Choose the mineral name exactly as it appears in the mineral database file.
MINERAL_NAME = "Brucite"

# Give the mineral a human-readable formula string for printed output only.
MINERAL_FORMULA = "Mg(OH)2"

# Choose which dissolved element should be plotted as the solubility signal.
TARGET_DISSOLVED_ELEMENT = "Mg"

# List the aqueous species that will be included in the aqueous phase.
# This is the main list to change when adapting the tutorial to another system.
AQUEOUS_SPECIES = [
    "H2O(aq)",
    "H+(aq)",
    "OH-(aq)",
    "Mg+2(aq)",
    "MgOH+(aq)",
]

# Optionally set the solvent species explicitly.
# Leave this as None if you want the script to infer a likely solvent automatically.
SOLVENT_SPECIES_NAME = None


# -----------------------------------------------------------------------------
# 2B) Thermodynamic databases and DEW model settings
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# DEW settable API explained (Python)
#
# 1) DEWDatabase
#    - DEWDatabase(name)
#    - DEWDatabase.withName(name), fromFile(path), fromContents(text), load(...)
#
# 2) WaterModelOptions (settable fields)
#    - eosModel
#    - dielectricModel
#    - gibbsModel
#    - bornModel
#    - usePsatPolynomials
#    - psatRelTol
#    - densityTolerance
#
# 3) StandardThermoModelParamsDEW (species-level HKF/DEW parameters)
#    - Gf, Hf, Sr, a1, a2, a3, a4, c1, c2, wref, charge, Tmax, waterOptions
#    - Note: this tutorial uses species thermo directly from DEW database entries,
#      so these fields are not manually assigned here.
#
# 4) ActivityModelParamsDEW (aqueous activity-model options)
#    - waterOptions
#    - bExtended
#
# 5) Model factories used in this tutorial
#    - StandardThermoModelDEW(params)
#    - ActivityModelDEW(params)
# -----------------------------------------------------------------------------

# Choose the name of the DEW aqueous database.
DEW_DATABASE_NAME = "dew2024-aqueous"

# Choose the mineral database file name.
SUPCRTBL_DATABASE_FILENAME = "supcrtbl.yaml"

# Choose the activity model for the aqueous phase.
# Available options in Reaktoro: ActivityModelDEW (this tutorial), ActivityModelPitzer,
# ActivityModelHKF, ActivityModelDebyeHuckel, ActivityModelIdealAqueous
# Change this to switch to a different model without editing the rest of the script.
# The DEW (Deep Earth Water) activity model is used for the aqueous phase.
# We configure it with ActivityModelParamsDEW so all DEW activity options are explicit.

# Water equation-of-state used for density/Cp calculations.
# Physical meaning: determines water density and volumetric properties at high P-T,
# which feed directly into the Debye-Huckel A and B parameters and ion solvation energies.
# Options: "ZhangDuan2005" (DEW default), "ZhangDuan2009", "WagnerPruss", "HGK"
DEW_WATER_EOS_MODEL = "ZhangDuan2005"

# Dielectric-constant model for the relative permittivity (epsilon) of water.
# Physical meaning: the dielectric constant measures how strongly water screens
# electrostatic interactions between ions; a higher epsilon means charge is felt
# less strongly, reducing activity coefficients at a given ionic strength.
# Options: "PowerFunction" (DEW default), "JohnsonNorton1991", "Franck1990", "Fernandez1997"
DEW_WATER_DIELECTRIC_MODEL = "PowerFunction"

# Water Gibbs energy model used to compute the solvent contribution to omega.
# Physical meaning: the Gibbs energy of water at P-T enters the Born solvation
# energy calculation; different models give slightly different solvent G curves
# at deep-Earth conditions and thus shift the solvation free energies of ions.
# Options: "DewIntegral" (DEW default), "DelaneyHelgeson1978"
DEW_WATER_GIBBS_MODEL = "DewIntegral"

# Born model for the conventional Born coefficient omega of solvent water.
# Physical meaning: omega quantifies the electrostatic work required to transfer
# an ion from vacuum into the dielectric medium; this correction shifts standard
# partial molar volumes and heat capacities of aqueous species at high P.
# Options: "Shock92Dew" (DEW default), "None" (disables Born correction entirely)
DEW_WATER_BORN_MODEL = "Shock92Dew"

# Whether to use DEW polynomial fits for the saturation pressure Psat(T).
# Physical meaning: Psat sets the liquid-vapour boundary; it is used internally
# to distinguish compressed-liquid from steam conditions when evaluating water
# thermodynamic properties above the critical isochore.
# Options: True (faster polynomial fit, DEW default), False (full iterative Psat solver)
DEW_USE_PSAT_POLYNOMIALS = True

# Relative tolerance for the iterative Psat solver (only used when DEW_USE_PSAT_POLYNOMIALS=False).
# Physical meaning: tighter tolerances give a more precise liquid-vapour boundary
# at the cost of extra Newton iterations; only relevant for near-saturation conditions.
# Typical values: 1e-3 (loose), 1e-6 (tight). Has no effect when DEW_USE_PSAT_POLYNOMIALS=True.
DEW_PSAT_RELATIVE_TOLERANCE = 1.0e-3

# Density tolerance in bar used in the DEW water-state solver convergence check.
# Physical meaning: the solver iterates on water density to match the target P-T;
# looser tolerances converge faster but may produce slightly inconsistent densities
# that propagate into the dielectric constant and hence into activity coefficients.
# Typical values: 1e-3 (DEW default, loose), 1e-6 (tighter, slower).
DEW_DENSITY_TOLERANCE_BAR = 1.0e-3

# Extended Debye-Huckel correction term b_c,k used in ActivityModelDEW.
# Physical meaning: adds a linear ionic-strength term (b * I) to log(gamma) that
# accounts for short-range ion-solvent interactions beyond the limiting Debye-Huckel
# law; positive b raises activity coefficients at high ionic strength.
# Options: 0.0 (standard DEW default — no extended correction),
#          any positive float (adds a linear ionic-strength term b*I to log gamma).
DEW_ACTIVITY_B_EXTENDED = 0.0


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
    "H2O(aq)": 55.5,
    "H+(aq)": 1e-8,
    "OH-(aq)": 1e-8,
    "Mg+2(aq)": 1e-8,
    MINERAL_NAME: 10.0,
}


# -----------------------------------------------------------------------------
# 2E) Graph and output settings
# -----------------------------------------------------------------------------

# Give the tutorial run a short human-readable title.
RUN_TITLE = "Brucite solubility tutorial"

# Choose the output figure name.
OUTPUT_PLOT_FILENAME = "brucite_solubility_tutorial.png"

# Choose the plot title.
PLOT_TITLE = "Brucite Solubility Tutorial (DEW2024 + SUPCRTBL)"

# Choose the x-axis label.
PLOT_X_LABEL = "Temperature (C)"

# Choose the y-axis label.
# This is placed after TARGET_DISSOLVED_ELEMENT so the label updates automatically when you change the element.
PLOT_Y_LABEL = f"Total dissolved {TARGET_DISSOLVED_ELEMENT} molality (mol/kg-H2O)"

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

# Build the absolute path to the SUPCRTBL mineral database file in this repository.
SUPCRTBL_DATABASE_FILE = os.path.join(
    os.path.dirname(os.path.dirname(SCRIPT_DIR)),
    "embedded",
    "databases",
    "reaktoro",
    SUPCRTBL_DATABASE_FILENAME,
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
    preferred_names = ["H2O(aq)", "H2O", "H2O(l)"]

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
    print("[Databases and DEW model settings]")
    print(f"Aqueous database: {DEW_DATABASE_NAME}")
    print(f"Mineral database file: {SUPCRTBL_DATABASE_FILENAME}")
    print("DEW activity model choices used by this tutorial:")
    print(f"  Water EOS model: {DEW_WATER_EOS_MODEL}")
    print(f"  Water dielectric model: {DEW_WATER_DIELECTRIC_MODEL}")
    print(f"  Water Gibbs model: {DEW_WATER_GIBBS_MODEL}")
    print(f"  Water Born model: {DEW_WATER_BORN_MODEL}")
    print(f"  Use Psat polynomials: {DEW_USE_PSAT_POLYNOMIALS}")
    print(f"  Psat relative tolerance: {DEW_PSAT_RELATIVE_TOLERANCE}")
    print(f"  Density tolerance (bar): {DEW_DENSITY_TOLERANCE_BAR}")
    print(f"  Activity bExtended: {DEW_ACTIVITY_B_EXTENDED}")
    print("  ActivityModelDEW is built from ActivityModelParamsDEW in this tutorial.")
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
    print("=" * 78)


# Define a helper that builds the chemical system used for the equilibrium calculations.
def build_tutorial_system():
    """Create the ChemicalSystem for the current tutorial configuration."""

    # Load the DEW aqueous database.
    dew_db = DEWDatabase(DEW_DATABASE_NAME)

    # Load the mineral database file.
    supcrt_db = Database.fromFile(SUPCRTBL_DATABASE_FILE)

    # Start a combined database with all aqueous species from the DEW database.
    combined_db = Database(dew_db.species())

    # Add the chosen mineral species from the mineral database into the combined database.
    combined_db.addSpecies(supcrt_db.species(MINERAL_NAME))

    # Join the aqueous species names into the single space-separated format Reaktoro expects.
    aqueous_species_string = " ".join(AQUEOUS_SPECIES)

    # Create the aqueous phase from the chosen aqueous species.
    aqueous_phase = AqueousPhase(aqueous_species_string)

    # Configure explicit DEW water options used by the activity model.
    water_options = makeWaterModelOptionsDEW()
    water_options.eosModel = getattr(WaterEosModel, DEW_WATER_EOS_MODEL)
    water_options.dielectricModel = getattr(
        WaterDielectricModel, DEW_WATER_DIELECTRIC_MODEL
    )
    water_options.gibbsModel = getattr(WaterGibbsModel, DEW_WATER_GIBBS_MODEL)
    water_options.bornModel = getattr(WaterBornModel, DEW_WATER_BORN_MODEL)
    water_options.usePsatPolynomials = bool(DEW_USE_PSAT_POLYNOMIALS)
    water_options.psatRelTol = float(DEW_PSAT_RELATIVE_TOLERANCE)
    water_options.densityTolerance = float(DEW_DENSITY_TOLERANCE_BAR)

    # Configure DEW activity-model parameters.
    activity_params = ActivityModelParamsDEW()
    activity_params.waterOptions = water_options
    activity_params.bExtended = float(DEW_ACTIVITY_B_EXTENDED)

    # Set the activity model for the aqueous phase using explicit DEW parameters.
    aqueous_phase.setActivityModel(ActivityModelDEW(activity_params))

    # -------------------------------------------------------------------------
    # Model-swap reference: use PerplexDEW instead of DEW
    #
    # To switch from ActivityModelDEW to ActivityModelPerplexDEW, replace the
    # three lines above (water_options build + setActivityModel) with:
    #
    #   from reaktoro4py import ActivityModelParamsPerplexDEW, ActivityModelPerplexDEW
    #   perp_params = ActivityModelParamsPerplexDEW()
    #   perp_params.dhModel = ActivityDHModel.Davies      # or ActivityDHModel.ExtendedDH
    #   perp_params.errorOnConflictingStandardState = False
    #   aqueous_phase.setActivityModel(ActivityModelPerplexDEW(perp_params))
    #
    # PerplexDEW uses Perple_X solvent internals (ZD05 + Looyenga mixing) instead
    # of the DEW spreadsheet water sub-models, but both backends share the same
    # ActivityDHModel enum and the same ChemicalSystem/EquilibriumSolver interface.
    # The AQUEOUS_SPECIES list and solver loop below are unchanged for either model.
    # -------------------------------------------------------------------------

    # Create the mineral phase for the chosen mineral.
    mineral_phase = MineralPhase(MINERAL_NAME)

    # Return the full chemical system.
    return ChemicalSystem(combined_db, aqueous_phase, mineral_phase)


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
                # Build the aqueous-properties helper from the updated equilibrium state.
                aqprops = AqueousProps(state)

                # Read the molality of the chosen dissolved element.
                element_molality = float(
                    aqprops.elementMolality(TARGET_DISSOLVED_ELEMENT)
                )

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
