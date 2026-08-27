"""
Beginner Tutorial: Anhydrite Solubility with an Oxygen Fugacity Buffer

Summary
- This script calculates anhydrite solubility in water versus temperature.
- It applies a mineral oxygen fugacity buffer (for example NNO, HM, or Mn2O3-MnO).
- It saves one plot with solubility curves at fixed pressures.

Goal
- Keep the code minimal and beginner-friendly.
- Show clearly where to choose mineral, aqueous species, and fO2 buffer.


"""

import os
import sys
import numpy as np
import matplotlib

matplotlib.use("Agg")  # Non-interactive backend; required for savefig without a display
import matplotlib.pyplot as plt

# Prefer the local build in build/Reaktoro/Release when available.
THIS_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(os.path.dirname(THIS_DIR))
LOCAL_BUILD_PYD_DIR = os.path.join(REPO_ROOT, "build", "Reaktoro", "Release")

if os.path.isdir(LOCAL_BUILD_PYD_DIR) and LOCAL_BUILD_PYD_DIR not in sys.path:
    sys.path.insert(0, LOCAL_BUILD_PYD_DIR)

# autodiff must be imported AFTER the local build path is on sys.path,
# because it is bundled inside build/Reaktoro/Release alongside reaktoro4py.
import autodiff  # noqa: E402

try:
    # Use the local compiled extension from this repository build.
    from reaktoro4py import *  # noqa: F401,F403

    print(f"Using local Reaktoro build from: {LOCAL_BUILD_PYD_DIR}")
except ModuleNotFoundError:
    # Fallback only if local build is not available in the environment.
    from reaktoro import *  # noqa: F401,F403

    print("Using installed 'reaktoro' package (local build not found).")


# -----------------------------------------------------------------------------
# 1) User input section (edit this first)
# -----------------------------------------------------------------------------

# Mineral phase in SUPCRTBL.
MINERAL_NAME = "Anhydrite"
MINERAL_FORMULA = "CaSO4"

# Explicit aqueous species list (DEW database names).
# Use explicit names confirmed in this build's DEW species list.
# Include reduced and intermediate sulfur species plus dissolved oxygen so
# sulfur speciation can respond to the imposed oxygen fugacity.
AQUEOUS_SPECIES = [
    "H2O(aq)",
    "H+(aq)",
    "OH-(aq)",
    "O2(aq)",
    "Ca+2(aq)",
    "HS-(aq)",
    "H2S(aq)",
    "HSO3-(aq)",
    "SO3-2(aq)",
    "S2O3-2(aq)",
    "S2-2(aq)",
    "SO4-2(aq)",
    "HSO4-(aq)",
]

# Species/element to report as solubility.
TARGET_DISSOLVED_ELEMENT = "Ca"

# Diagnostic species to plot alongside dissolved Ca so the redox effect is
# visible directly in the output figure.
DIAGNOSTIC_SPECIES = [
    ("Ca", "element", "Total dissolved Ca molality (mol/kg-H2O)", "log"),
    ("HS-(aq)", "species", "HS-(aq) molality (mol/kg-H2O)", "log"),
    ("H2S(aq)", "species", "H2S(aq) molality (mol/kg-H2O)", "log"),
    ("SO4-2(aq)", "species", "SO4-2(aq) molality (mol/kg-H2O)", "log"),
    (
        "log10((HS-+H2S)/SO4-2)",
        "ratio",
        "log10((HS-(aq)+H2S(aq))/SO4-2(aq))",
        "linear",
    ),
]

# Oxygen fugacity buffers to compare.
# Common options from the helper: NNO, HM, Mn2O3-MnO, IW, FMQ, CCO.
FUGACITY_BUFFER_NAMES = ["NNO", "FMQ"]

# Fugacity species constrained in the equilibrium problem.
FUGACITY_SPECIES_NAME = "O2(g)"

# Include sulfur-bearing gases supported by SUPCRTBL so gas/aqueous sulfur
# speciation can shift with the imposed oxygen fugacity.
GAS_SPECIES = ["O2(g)", "S2(g)", "H2S(g)", "SO2(g)"]

# Pressure and temperature grid.
PRESSURES_KBAR = [1.0, 2.0, 5.0]
TEMPERATURE_MIN_C = 100.0
TEMPERATURE_MAX_C = 500.0
NUMBER_OF_TEMPERATURE_POINTS = 80

# Output plot location.
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_PLOT = os.path.join(SCRIPT_DIR, "anhydrite_solubility_buffer_tutorial.png")


# -----------------------------------------------------------------------------
# 2) Buffer fugacity helper
# -----------------------------------------------------------------------------

# The buffer helper lives in Mineral_Solubilities/anhydrite.
BUFFER_HELPER_DIR = os.path.join(
    os.path.dirname(SCRIPT_DIR), "Mineral_Solubilities", "anhydrite"
)
if BUFFER_HELPER_DIR not in sys.path:
    sys.path.insert(0, BUFFER_HELPER_DIR)

import buffer_fO2_from_supcrtbl as buffer_supcrtbl


# -----------------------------------------------------------------------------
# 3) Build chemical system
# -----------------------------------------------------------------------------


# Path to the SUPCRTBL database file embedded in this repository.
SUPCRTBL_DATABASE_FILE = os.path.join(
    REPO_ROOT, "embedded", "databases", "reaktoro", "supcrtbl.yaml"
)


def build_anhydrite_system():
    """Create system with redox-sensitive sulfur chemistry and gas species."""

    dew_db = DEWDatabase("dew2024-aqueous")
    supcrt_db = Database.fromFile(SUPCRTBL_DATABASE_FILE)

    combined_db = Database(dew_db.species())
    combined_db.addSpecies(supcrt_db.species(MINERAL_NAME))
    for gas_species_name in GAS_SPECIES:
        combined_db.addSpecies(supcrt_db.species(gas_species_name))

    aqueous_phase = AqueousPhase(" ".join(AQUEOUS_SPECIES))
    aqueous_phase.setActivityModel(ActivityModelDEW())

    mineral_phase = MineralPhase(MINERAL_NAME)
    gas_phase = GaseousPhase(" ".join(GAS_SPECIES))

    return ChemicalSystem(combined_db, aqueous_phase, mineral_phase, gas_phase)


def set_state_amount(state, species_name, amount, unit):
    """Set species amount with compatibility for local reaktoro4py typed API."""
    try:
        state.set(species_name, float(amount), unit)
    except TypeError:
        state.set(species_name, autodiff.real(float(amount)), unit)


def get_diagnostic_value(aqprops, name, value_kind):
    """Read one diagnostic quantity from aqueous properties."""
    if value_kind == "element":
        return float(aqprops.elementMolality(name))
    if value_kind == "ratio":
        hs = float(aqprops.speciesMolality("HS-(aq)"))
        h2s = float(aqprops.speciesMolality("H2S(aq)"))
        so4 = float(aqprops.speciesMolality("SO4-2(aq)"))
        numerator = hs + h2s
        if numerator <= 0.0 or so4 <= 0.0:
            return np.nan
        return float(np.log10(numerator / so4))
    return float(aqprops.speciesMolality(name))


# -----------------------------------------------------------------------------
# 4) Main calculation
# -----------------------------------------------------------------------------


def main():
    buffer_label = " vs ".join(FUGACITY_BUFFER_NAMES)
    print("=" * 78)
    print("Anhydrite solubility tutorial with oxygen fugacity buffer")
    print(f"Mineral: {MINERAL_NAME} ({MINERAL_FORMULA})")
    print("Aqueous species: " + ", ".join(AQUEOUS_SPECIES))
    print(f"fO2 buffers: {buffer_label}")
    print("=" * 78)

    # Suppress the non-convergence warning spam (warning 906) — these are
    # expected for some T/P points and are handled gracefully as NaN.
    Warnings.disable(906)

    system = build_anhydrite_system()

    # This solver uses temperature, pressure, and a fixed fugacity condition.
    specs = EquilibriumSpecs(system)
    specs.temperature()
    specs.pressure()
    specs.fugacity(FUGACITY_SPECIES_NAME)

    solver = EquilibriumSolver(specs)
    conditions = EquilibriumConditions(specs)

    temperatures_c = np.linspace(
        TEMPERATURE_MIN_C,
        TEMPERATURE_MAX_C,
        NUMBER_OF_TEMPERATURE_POINTS,
    )

    # curves[buffer_name][diagnostic_name][pressure_kbar] = diagnostic value array
    all_curves = {}

    for buffer_name in FUGACITY_BUFFER_NAMES:
        curves = {name: {} for name, _, _, _ in DIAGNOSTIC_SPECIES}
        for pressure_kbar in PRESSURES_KBAR:
            pressure_bar = pressure_kbar * 1000.0
            diagnostic_values = {name: [] for name, _, _, _ in DIAGNOSTIC_SPECIES}

            state = ChemicalState(system)
            set_state_amount(state, "H2O(aq)", 55.5, "mol")
            set_state_amount(state, "H+(aq)", 1e-8, "mol")
            set_state_amount(state, "OH-(aq)", 1e-8, "mol")
            set_state_amount(state, "Ca+2(aq)", 1e-6, "mol")
            set_state_amount(state, MINERAL_NAME, 10.0, "mol")
            set_state_amount(state, FUGACITY_SPECIES_NAME, 1e-12, "mol")

            for temperature_c in temperatures_c:
                f_o2_bar = float(
                    buffer_supcrtbl.buffer_fugacity_bar(
                        "fO2",
                        buffer_name,
                        float(temperature_c),
                        float(pressure_bar),
                    )
                )

                conditions.temperature(float(temperature_c), "celsius")
                conditions.pressure(float(pressure_bar), "bar")
                conditions.fugacity(FUGACITY_SPECIES_NAME, f_o2_bar, "bar")

                result = solver.solve(state, conditions)

                if result.succeeded():
                    aqprops = AqueousProps(state)
                    for name, value_kind, _, _ in DIAGNOSTIC_SPECIES:
                        diagnostic_values[name].append(
                            get_diagnostic_value(aqprops, name, value_kind)
                        )
                else:
                    for name, _, _, _ in DIAGNOSTIC_SPECIES:
                        diagnostic_values[name].append(np.nan)

            for name, _, _, _ in DIAGNOSTIC_SPECIES:
                curves[name][pressure_kbar] = np.array(diagnostic_values[name])

            valid = np.isfinite(curves[TARGET_DISSOLVED_ELEMENT][pressure_kbar]).sum()
            print(
                f"[{buffer_name}] Computed {valid}/{len(temperatures_c)} points at {pressure_kbar:.1f} kbar",
                flush=True,
            )

        all_curves[buffer_name] = curves

    # One linestyle per buffer, one colour per pressure.
    linestyles = ["-", "--", "-.", ":"]
    colors = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple"]

    # Plot Ca plus diagnostic aqueous species in one figure.
    fig, axes = plt.subplots(3, 2, figsize=(13, 12), sharex=True)
    axes = axes.ravel()

    for plot_index, (diagnostic_name, _, ylabel, yscale) in enumerate(
        DIAGNOSTIC_SPECIES
    ):
        axis = axes[plot_index]
        for bi, buffer_name in enumerate(FUGACITY_BUFFER_NAMES):
            ls = linestyles[bi % len(linestyles)]
            for ci, pressure_kbar in enumerate(PRESSURES_KBAR):
                y = all_curves[buffer_name][diagnostic_name][pressure_kbar]
                if yscale == "log":
                    valid = np.isfinite(y) & (y > 0.0)
                else:
                    valid = np.isfinite(y)
                axis.plot(
                    temperatures_c[valid],
                    y[valid],
                    linestyle=ls,
                    color=colors[ci % len(colors)],
                    linewidth=2.0,
                    label=f"{buffer_name} {pressure_kbar:.1f} kbar",
                )

        if yscale == "log":
            axis.set_yscale("log")
        axis.set_ylabel(ylabel)
        axis.grid(True, which="both", alpha=0.3, linestyle="--")
        axis.set_title(diagnostic_name)

    axes[4].set_xlabel("Temperature (°C)")
    axes[5].set_xlabel("Temperature (°C)")
    axes[5].axis("off")
    handles, labels = axes[0].get_legend_handles_labels()
    axes[5].legend(handles, labels, title="Buffer / Pressure", loc="center")

    fig.suptitle(
        f"Anhydrite Solubility and Redox Diagnostics: {buffer_label} Buffer Comparison"
    )
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.97))
    fig.savefig(OUTPUT_PLOT, dpi=250)
    plt.close(fig)

    print(f"Saved plot: {OUTPUT_PLOT}", flush=True)


if __name__ == "__main__":
    main()
