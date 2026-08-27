"""
Brucite Database-Agnostic Solubility Tutorial / Test Harness

Goal
- Run the same Brucite solubility workflow against different database setups.
- Switch aqueous activity model with one variable: DEW or PerplexDEW.
- Optionally run a full matrix of database combinations and save results.

Supported database modes
1) PerpleX mixed database (aqueous + mineral in one JSON)
2) DEW aqueous database + mineral species from a PerpleX mineral database JSON

Notes
- The Brucite mineral in PerpleX files is typically named "br".
- DEW aqueous species usually use suffixes like H2O(aq), H+(aq), Mg+2(aq).
- PerpleX DEW-style aqueous species typically use names like H2O, H+, Mg+2.
"""

import csv
import json
import os
import argparse
import subprocess
import sys
import traceback
import copy
from dataclasses import dataclass
from dataclasses import asdict
from typing import Dict
from typing import List
from typing import Optional

import numpy as np

# Optional plotting, disabled by default for matrix runs.
import matplotlib.pyplot as plt


# -----------------------------------------------------------------------------
# 1) Import Reaktoro from local build when available
# -----------------------------------------------------------------------------

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(THIS_DIR)))
LOCAL_BUILD_PYD_DIR = os.path.join(REPO_ROOT, "build", "Reaktoro", "Release")

if os.path.isdir(LOCAL_BUILD_PYD_DIR) and LOCAL_BUILD_PYD_DIR not in sys.path:
    sys.path.insert(0, LOCAL_BUILD_PYD_DIR)

try:
    import autodiff  # noqa: E402
except ModuleNotFoundError:

    class _AutodiffShim:
        @staticmethod
        def real(value):
            return value

    autodiff = _AutodiffShim()

try:
    from reaktoro4py import *  # noqa: F401,F403

    print(f"Using local Reaktoro build from: {LOCAL_BUILD_PYD_DIR}")
except ModuleNotFoundError:
    from reaktoro import *  # noqa: F401,F403

    print("Using installed 'reaktoro' package (local build not found).")


# -----------------------------------------------------------------------------
# 2) User controls (simple switches)
# -----------------------------------------------------------------------------

# Main switch: "DEW" or "PerplexDEW"
AQUEOUS_ACTIVITY_MODEL_NAME = "PerplexDEW"

# Run either one selected case, or all matrix cases below.
RUN_ALL_CASES = True
SELECTED_CASE_NAME = "perplex_mixed_dew24hp622"

# Numerical settings
TEMPERATURE_MIN_C = 50.0
TEMPERATURE_MAX_C = 450.0
NUMBER_OF_TEMPERATURE_POINTS = 80
PRESSURES_KBAR = [1.0, 2.0, 5.0]

# Initial amounts are species-name dependent and set per case type below.
INITIAL_MINERAL_AMOUNT_MOL = 10.0
INITIAL_WATER_AMOUNT_MOL = 55.5
INITIAL_TRACE_AMOUNT_MOL = 1.0e-8
INITIAL_MGOH_TRACE_AMOUNT_MOL = 1.0e-10

# Plot control
SAVE_PLOT_FOR_SINGLE_CASE = True

# Output files
RESULTS_JSON_FILE = os.path.join(THIS_DIR, "brucite_database_test_results.json")
RESULTS_CSV_FILE = os.path.join(THIS_DIR, "brucite_database_test_results.csv")
PLOT_FILE = os.path.join(THIS_DIR, "brucite_database_test_plot.png")
GENERATED_DB_DIR = os.path.join(THIS_DIR, "_generated")
WATER_MOLAR_MASS_KG_PER_MOL = 0.01801528


# -----------------------------------------------------------------------------
# 3) Standardized test-case definitions
# -----------------------------------------------------------------------------

PERPLEX_DB_DIR = os.path.join(REPO_ROOT, "embedded", "databases", "perplex")


@dataclass
class BruciteTestCase:
    name: str
    source_mode: str  # "perplex_mixed" or "dew_plus_perplex_mineral"
    aqueous_db_name: Optional[str]
    database_file: Optional[str]
    mineral_db_file: Optional[str]
    mineral_name: str
    aqueous_species: List[str]
    solvent_species: str


def build_test_cases() -> List[BruciteTestCase]:
    """Return standardized database combinations requested for brucite tests."""

    return [
        # PerpleX mixed: aqueous + mineral in same file
        BruciteTestCase(
            name="perplex_mixed_dew24hp622",
            source_mode="perplex_mixed",
            aqueous_db_name=None,
            database_file=os.path.join(
                PERPLEX_DB_DIR, "DEW24HP622ver_elements-reaktoro.json"
            ),
            mineral_db_file=None,
            mineral_name="br",
            aqueous_species=["H2O", "H+", "OH-", "Mg+2", "MgOH+"],
            solvent_species="H2O",
        ),
        BruciteTestCase(
            name="perplex_mixed_dew17hp622_zn_2025",
            source_mode="perplex_mixed",
            aqueous_db_name=None,
            database_file=os.path.join(
                PERPLEX_DB_DIR, "DEW17HP622_Zn_2025-reaktoro.json"
            ),
            mineral_db_file=None,
            mineral_name="br",
            aqueous_species=["H2O", "H+", "OH-", "Mg+2", "MgOH+"],
            solvent_species="H2O",
        ),
        # DEW aqueous + HP62 mineral database
        BruciteTestCase(
            name="dew2024_aqueous_plus_hp62_mineral",
            source_mode="dew_plus_perplex_mineral",
            aqueous_db_name="dew2024-aqueous",
            database_file=None,
            mineral_db_file=os.path.join(PERPLEX_DB_DIR, "hp62ver-reaktoro.json"),
            mineral_name="br",
            aqueous_species=["H2O(aq)", "H+(aq)", "OH-(aq)", "Mg+2(aq)", "MgOH+(aq)"],
            solvent_species="H2O(aq)",
        ),
        # "Holland & Powell 636" was requested; hp633/hp634 are available in this workspace.
        BruciteTestCase(
            name="dew2024_aqueous_plus_hp633_mineral",
            source_mode="dew_plus_perplex_mineral",
            aqueous_db_name="dew2024-aqueous",
            database_file=None,
            mineral_db_file=os.path.join(PERPLEX_DB_DIR, "hp633ver-reaktoro.json"),
            mineral_name="br",
            aqueous_species=["H2O(aq)", "H+(aq)", "OH-(aq)", "Mg+2(aq)", "MgOH+(aq)"],
            solvent_species="H2O(aq)",
        ),
        BruciteTestCase(
            name="dew2024_aqueous_plus_hp634_mineral",
            source_mode="dew_plus_perplex_mineral",
            aqueous_db_name="dew2024-aqueous",
            database_file=None,
            mineral_db_file=os.path.join(PERPLEX_DB_DIR, "hp634ver-reaktoro.json"),
            mineral_name="br",
            aqueous_species=["H2O(aq)", "H+(aq)", "OH-(aq)", "Mg+2(aq)", "MgOH+(aq)"],
            solvent_species="H2O(aq)",
        ),
        # DEW aqueous + PerpleX mineral-only source (taking only br from mixed file)
        BruciteTestCase(
            name="dew2024_aqueous_plus_perplex_dew24_mineral",
            source_mode="dew_plus_perplex_mineral",
            aqueous_db_name="dew2024-aqueous",
            database_file=None,
            mineral_db_file=os.path.join(
                PERPLEX_DB_DIR, "DEW24HP622ver_elements-reaktoro.json"
            ),
            mineral_name="br",
            aqueous_species=["H2O(aq)", "H+(aq)", "OH-(aq)", "Mg+2(aq)", "MgOH+(aq)"],
            solvent_species="H2O(aq)",
        ),
        BruciteTestCase(
            name="dew2024_aqueous_plus_perplex_dew17zn2025_mineral",
            source_mode="dew_plus_perplex_mineral",
            aqueous_db_name="dew2024-aqueous",
            database_file=None,
            mineral_db_file=os.path.join(
                PERPLEX_DB_DIR, "DEW17HP622_Zn_2025-reaktoro.json"
            ),
            mineral_name="br",
            aqueous_species=["H2O(aq)", "H+(aq)", "OH-(aq)", "Mg+2(aq)", "MgOH+(aq)"],
            solvent_species="H2O(aq)",
        ),
    ]


# -----------------------------------------------------------------------------
# 4) Core helpers
# -----------------------------------------------------------------------------


def validate_paths(cases: List[BruciteTestCase]) -> None:
    for case in cases:
        if case.database_file and not os.path.isfile(case.database_file):
            raise FileNotFoundError(
                f"Missing database file for case '{case.name}': {case.database_file}"
            )
        if case.mineral_db_file and not os.path.isfile(case.mineral_db_file):
            raise FileNotFoundError(
                f"Missing mineral DB file for case '{case.name}': {case.mineral_db_file}"
            )


def set_state_amount(state, species_name: str, amount: float, unit: str) -> None:
    try:
        state.set(species_name, float(amount), unit)
    except TypeError:
        state.set(species_name, autodiff.real(float(amount)), unit)


def sanitize_species_entry(species_entry: Dict[str, object]) -> Dict[str, object]:
    """Apply minimal fixes so Database.fromFile can instantiate the species."""

    out = copy.deepcopy(species_entry)

    # Keep proton formula explicit for AqueousProps compatibility.
    if out.get("Name") == "H+":
        out["Formula"] = "H+"

    # In PerpleX mixed DEW databases, water may appear as GFSM H2O with Gas aggregate.
    # For aqueous-phase brucite tests we want H2O usable as aqueous solvent species.
    if out.get("Name") == "H2O":
        out["AggregateState"] = "Aqueous"

    has_formation_reaction = "FormationReaction" in out
    has_standard_model = "StandardThermoModel" in out

    # Some PerpleX conversions only carry ThermoReference; add a constant fallback
    # with G0/H0/V0 so derived reporting properties are not forced to zero.
    if not has_formation_reaction and not has_standard_model:

        def _as_float(value: object, default: float = 0.0) -> float:
            try:
                return float(value)
            except (TypeError, ValueError):
                return default

        thermo_ref = out.get("ThermoReference", {})
        metadata = out.get("Metadata", {})
        perplex_params = {}
        if isinstance(metadata, dict):
            params = metadata.get("PerpleX_Params", {})
            if isinstance(params, dict):
                perplex_params = params

        g0 = _as_float(
            thermo_ref.get("Gf", 0.0) if isinstance(thermo_ref, dict) else 0.0
        )
        s0 = _as_float(
            thermo_ref.get("S0", perplex_params.get("S0", 0.0))
            if isinstance(thermo_ref, dict)
            else perplex_params.get("S0", 0.0)
        )
        # Perple_X V0 is in J/bar; convert to m3/mol with factor 1e-5.
        v0_j_per_bar = _as_float(
            thermo_ref.get("V0", perplex_params.get("V0", 0.0))
            if isinstance(thermo_ref, dict)
            else perplex_params.get("V0", 0.0)
        )

        h0 = g0 + 298.15 * s0
        v0 = v0_j_per_bar * 1.0e-5

        out["StandardThermoModel"] = {
            "Constant": {
                "G0": g0,
                "H0": h0,
                "V0": v0,
            }
        }

    return out


def create_subset_database_file(
    source_json_path: str,
    required_species_names: List[str],
    output_basename: str,
) -> str:
    """Create a sanitized subset JSON containing only required species."""

    with open(source_json_path, "r", encoding="utf-8") as file:
        data = json.load(file)

    species_data = data.get("Species", {})

    missing = [name for name in required_species_names if name not in species_data]
    if missing:
        raise ValueError(
            f"Missing required species in '{source_json_path}': {', '.join(missing)}"
        )

    reduced_species = {}
    for name in required_species_names:
        reduced_species[name] = sanitize_species_entry(species_data[name])

    reduced_data = dict(data)
    reduced_data["Species"] = reduced_species

    os.makedirs(GENERATED_DB_DIR, exist_ok=True)
    out_path = os.path.join(GENERATED_DB_DIR, output_basename)

    with open(out_path, "w", encoding="utf-8") as file:
        json.dump(reduced_data, file, indent=2, allow_nan=False)

    return out_path


def create_initial_amounts(case: BruciteTestCase) -> Dict[str, float]:
    amounts = {
        case.solvent_species: INITIAL_WATER_AMOUNT_MOL,
        case.mineral_name: INITIAL_MINERAL_AMOUNT_MOL,
    }

    # Fill common solutes if they are present in the selected aqueous species list.
    for species_name in case.aqueous_species:
        if species_name in (case.solvent_species,):
            continue
        if species_name in ("MgOH+", "MgOH+(aq)"):
            amounts[species_name] = INITIAL_MGOH_TRACE_AMOUNT_MOL
        else:
            amounts[species_name] = INITIAL_TRACE_AMOUNT_MOL

    return amounts


def build_database_and_system(case: BruciteTestCase, model_name: str):
    if case.source_mode == "perplex_mixed":
        required = list(case.aqueous_species) + [case.mineral_name]
        subset_file = create_subset_database_file(
            case.database_file,
            required,
            f"{case.name}-subset.json",
        )
        db = Database.fromFile(subset_file)

    elif case.source_mode == "dew_plus_perplex_mineral":
        dew_db = DEWDatabase(case.aqueous_db_name)
        db = Database(dew_db.species())
        mineral_subset_file = create_subset_database_file(
            case.mineral_db_file,
            [case.mineral_name],
            f"{case.name}-mineral-subset.json",
        )
        mineral_db = Database.fromFile(mineral_subset_file)
        db.addSpecies(mineral_db.species(case.mineral_name))

    else:
        raise ValueError(f"Unknown source mode '{case.source_mode}'")

    aqueous_phase = AqueousPhase(" ".join(case.aqueous_species))

    if model_name == "DEW":
        aqueous_phase.setActivityModel(ActivityModelDEW())
    elif model_name == "PerplexDEW":
        aqueous_phase.setActivityModel(ActivityModelPerplexDEW())
    else:
        raise ValueError(
            f"Unknown model '{model_name}', expected 'DEW' or 'PerplexDEW'"
        )

    mineral_phase = MineralPhase(case.mineral_name)
    system = ChemicalSystem(db, aqueous_phase, mineral_phase)
    return system


def estimate_total_mg_molality_from_state(state, case: BruciteTestCase) -> float:
    """Fallback when AqueousProps cannot be constructed due proton-formula constraints."""

    water_name = case.solvent_species
    n_water = float(state.speciesAmount(water_name))
    if n_water <= 0.0:
        return float("nan")

    # Sum dissolved Mg-bearing species for this tutorial chemistry.
    mg_species_coeffs = {
        "Mg+2": 1.0,
        "Mg+2(aq)": 1.0,
        "MgOH+": 1.0,
        "MgOH+(aq)": 1.0,
    }

    n_mg_total = 0.0
    for species_name in case.aqueous_species:
        coeff = mg_species_coeffs.get(species_name)
        if coeff is None:
            continue
        n_mg_total += coeff * float(state.speciesAmount(species_name))

    kg_water = n_water * WATER_MOLAR_MASS_KG_PER_MOL
    if kg_water <= 0.0:
        return float("nan")

    return n_mg_total / kg_water


def run_case(case: BruciteTestCase, model_name: str) -> Dict[str, object]:
    result = {
        "case": case.name,
        "model": model_name,
        "status": "failed",
        "message": "",
        "valid_points": 0,
        "total_points": len(PRESSURES_KBAR) * NUMBER_OF_TEMPERATURE_POINTS,
        "pressures_kbar": PRESSURES_KBAR,
    }

    try:
        system = build_database_and_system(case, model_name)

        specs = EquilibriumSpecs(system)
        specs.temperature()
        specs.pressure()

        solver = EquilibriumSolver(specs)
        conditions = EquilibriumConditions(specs)

        temperatures_c = np.linspace(
            TEMPERATURE_MIN_C,
            TEMPERATURE_MAX_C,
            NUMBER_OF_TEMPERATURE_POINTS,
        )

        valid_points = 0
        all_curves = {}

        for pressure_kbar in PRESSURES_KBAR:
            pressure_bar = pressure_kbar * 1000.0
            state = ChemicalState(system)

            for species_name, amount in create_initial_amounts(case).items():
                set_state_amount(state, species_name, amount, "mol")

            y = []
            for temperature_c in temperatures_c:
                conditions.temperature(float(temperature_c), "celsius")
                conditions.pressure(float(pressure_bar), "bar")

                solve_result = solver.solve(state, conditions)
                if solve_result.succeeded():
                    try:
                        aqprops = AqueousProps(state)
                        m = float(aqprops.elementMolality("Mg"))
                    except Exception:
                        m = estimate_total_mg_molality_from_state(state, case)
                    y.append(m)
                    if np.isfinite(m):
                        valid_points += 1
                else:
                    y.append(np.nan)

            all_curves[pressure_kbar] = np.array(y)

        result["valid_points"] = int(valid_points)
        result["status"] = "ok" if valid_points > 0 else "failed"
        result["message"] = "completed" if valid_points > 0 else "no valid points"

        # Optional plot only when running a single case to keep matrix runs light.
        if SAVE_PLOT_FOR_SINGLE_CASE and not RUN_ALL_CASES:
            plt.figure(figsize=(9, 6))
            for pressure_kbar in PRESSURES_KBAR:
                y = all_curves[pressure_kbar]
                valid = np.isfinite(y)
                plt.plot(
                    temperatures_c[valid],
                    y[valid],
                    linewidth=2.0,
                    label=f"{pressure_kbar:.1f} kbar",
                )
            plt.yscale("log")
            plt.xlabel("Temperature (C)")
            plt.ylabel("Total dissolved Mg molality (mol/kg-H2O)")
            plt.title(f"Brucite test: {case.name} | {model_name}")
            plt.grid(True, which="both", alpha=0.3, linestyle="--")
            plt.legend(title="Pressure")
            plt.tight_layout()
            plt.savefig(PLOT_FILE, dpi=220)
            plt.close()

    except Exception as exc:
        result["status"] = "failed"
        result["message"] = f"{type(exc).__name__}: {exc}"
        result["traceback"] = traceback.format_exc()

    return result


def run_case_in_subprocess(case_name: str, model_name: str) -> Dict[str, object]:
    """Run one case/model in a child process to isolate potential native crashes."""

    cmd = [
        sys.executable,
        os.path.abspath(__file__),
        "--worker",
        "--case",
        case_name,
        "--model",
        model_name,
    ]

    try:
        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=900,
            cwd=REPO_ROOT,
        )

        payload = {
            "case": case_name,
            "model": model_name,
            "status": "failed",
            "valid_points": 0,
            "total_points": len(PRESSURES_KBAR) * NUMBER_OF_TEMPERATURE_POINTS,
            "message": f"worker exit code {proc.returncode}",
            "worker_stdout": proc.stdout,
            "worker_stderr": proc.stderr,
        }

        if proc.returncode == 0:
            for line in reversed(proc.stdout.splitlines()):
                line = line.strip()
                if not line:
                    continue
                try:
                    row = json.loads(line)
                    if isinstance(row, dict) and "case" in row and "model" in row:
                        return row
                except json.JSONDecodeError:
                    continue

        return payload

    except subprocess.TimeoutExpired:
        return {
            "case": case_name,
            "model": model_name,
            "status": "failed",
            "valid_points": 0,
            "total_points": len(PRESSURES_KBAR) * NUMBER_OF_TEMPERATURE_POINTS,
            "message": "worker timeout",
        }


def parse_args():
    parser = argparse.ArgumentParser(description="Brucite database test harness")
    parser.add_argument(
        "--worker", action="store_true", help="Run one case/model and emit one JSON row"
    )
    parser.add_argument("--case", default="", help="Case name for worker mode")
    parser.add_argument("--model", default="", help="Model name for worker mode")
    return parser.parse_args()


# -----------------------------------------------------------------------------
# 5) Main runner
# -----------------------------------------------------------------------------


def main() -> None:
    args = parse_args()

    cases = build_test_cases()
    validate_paths(cases)

    if args.worker:
        target_case = [case for case in cases if case.name == args.case]
        if not target_case:
            raise ValueError(f"Unknown worker case '{args.case}'")
        if args.model not in ("DEW", "PerplexDEW"):
            raise ValueError(f"Unknown worker model '{args.model}'")
        row = run_case(target_case[0], args.model)
        print(json.dumps(row))
        return

    if RUN_ALL_CASES:
        selected_cases = cases
        selected_models = ["DEW", "PerplexDEW"]
    else:
        selected_cases = [case for case in cases if case.name == SELECTED_CASE_NAME]
        if not selected_cases:
            available = ", ".join(case.name for case in cases)
            raise ValueError(
                f"SELECTED_CASE_NAME '{SELECTED_CASE_NAME}' not found. Available: {available}"
            )
        selected_models = [AQUEOUS_ACTIVITY_MODEL_NAME]

    print("=" * 90)
    print("Brucite Database-Agnostic Test Matrix")
    print("=" * 90)
    print(f"Cases: {len(selected_cases)}")
    print("Models:", ", ".join(selected_models))
    print(
        "Temperatures:",
        f"{TEMPERATURE_MIN_C}..{TEMPERATURE_MAX_C} C ({NUMBER_OF_TEMPERATURE_POINTS} points)",
    )
    print("Pressures (kbar):", ", ".join(str(x) for x in PRESSURES_KBAR))
    print("-" * 90)

    results = []

    for case in selected_cases:
        print(f"Case: {case.name}")
        for model_name in selected_models:
            print(f"  Running model={model_name} ...")
            row = run_case_in_subprocess(case.name, model_name)
            results.append(row)
            print(
                "    ->",
                row["status"],
                f"valid={row['valid_points']}/{row['total_points']}",
                f"message={row.get('message', '')}",
            )

    payload = {
        "script": os.path.abspath(__file__),
        "settings": {
            "run_all_cases": RUN_ALL_CASES,
            "selected_case_name": SELECTED_CASE_NAME,
            "selected_model": AQUEOUS_ACTIVITY_MODEL_NAME,
            "temperature_min_c": TEMPERATURE_MIN_C,
            "temperature_max_c": TEMPERATURE_MAX_C,
            "n_temperature_points": NUMBER_OF_TEMPERATURE_POINTS,
            "pressures_kbar": PRESSURES_KBAR,
        },
        "cases": [asdict(case) for case in selected_cases],
        "results": results,
    }

    with open(RESULTS_JSON_FILE, "w", encoding="utf-8") as file:
        json.dump(payload, file, indent=2)

    with open(RESULTS_CSV_FILE, "w", newline="", encoding="utf-8") as file:
        fieldnames = [
            "case",
            "model",
            "status",
            "valid_points",
            "total_points",
            "message",
        ]
        writer = csv.DictWriter(file, fieldnames=fieldnames)
        writer.writeheader()
        for row in results:
            writer.writerow({key: row.get(key, "") for key in fieldnames})

    print("-" * 90)
    print(f"Saved JSON results: {RESULTS_JSON_FILE}")
    print(f"Saved CSV results:  {RESULTS_CSV_FILE}")


if __name__ == "__main__":
    main()
