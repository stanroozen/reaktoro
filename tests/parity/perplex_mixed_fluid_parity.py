"""
Direct Perple_X vs PerplexDEW mixed-fluid parity checks for quartz and calcite.

This script runs a small regression-style case matrix at selected H2O-CO2 mixed-fluid
conditions. For each case it:

1. Solves a PerplexDEW aqueous + CO2(g) mixed-fluid system in Reaktoro.
2. Generates a temporary Perple_X meemum project with the same mineral/fluid bulk basis.
3. Parses the Perple_X .prn solute molalities.
4. Compares total dissolved mineral molality and associated aqueous species.

The output is a CSV summary plus a concise terminal report. The script exits non-zero
only if --fail-threshold-log10 is provided and a case exceeds that tolerance.
"""

from __future__ import annotations

import argparse
import csv
import importlib
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path


if os.name == "nt":
    ep = sys.prefix
    env_paths = [
        ep,
        os.path.join(ep, "Library", "mingw-w64", "bin"),
        os.path.join(ep, "Library", "usr", "bin"),
        os.path.join(ep, "Library", "bin"),
        os.path.join(ep, "Scripts"),
        os.path.join(ep, "bin"),
    ]
    sr = os.environ.get("SystemRoot", r"C:\Windows")
    os.environ["PATH"] = ";".join(
        [
            p
            for p in env_paths
            + [os.path.join(sr, "System32"), sr, os.path.join(sr, "System32", "Wbem")]
            if os.path.isdir(p)
        ]
    )


SCRIPT_DIR = Path(__file__).resolve().parent
BENCHMARK_DIR = SCRIPT_DIR.parent
ROOT_DIR = BENCHMARK_DIR.parent

for _build_pkg in [
    ROOT_DIR / "build" / "python" / "package",
    ROOT_DIR / "build" / "python" / "package",
    ROOT_DIR / "build" / "python" / "package",
]:
    _rkt_inner = _build_pkg / "reaktoro"
    if _rkt_inner.is_dir():
        if str(_build_pkg) not in sys.path:
            sys.path.insert(0, str(_build_pkg))
        os.environ["PATH"] = str(_rkt_inner) + os.pathsep + os.environ.get("PATH", "")
        break

try:
    from reaktoro import *
except (ModuleNotFoundError, ImportError):
    _pyd_dir = None
    for _d in [
        ROOT_DIR / "build" / "Reaktoro" / "Release",
        ROOT_DIR / "build" / "Reaktoro" / "Release",
        ROOT_DIR / "build" / "Reaktoro" / "Debug",
        ROOT_DIR / "build" / "Reaktoro" / "Release",
    ]:
        if not _d.is_dir():
            continue
        sys.path.insert(0, str(_d))
        sys.modules.pop("reaktoro4py", None)
        try:
            _m = importlib.import_module("reaktoro4py")
            globals().update(
                {k: getattr(_m, k) for k in dir(_m) if not k.startswith("_")}
            )
            _pyd_dir = _d
            break
        except (ModuleNotFoundError, ImportError):
            continue
    if _pyd_dir is None:
        raise ModuleNotFoundError(
            "Reaktoro import failed for all local build candidates."
        )

try:
    Warnings.disable(906)
except Exception:
    pass


H2O_MOL_PER_KG = 55.508
MINERAL_MOLES = 10.0
N_GFSM_GAS_MOLES = 1000.0
OUTPUT_CSV = SCRIPT_DIR / "perplex_mixed_fluid_parity_results.csv"
OUTPUT_CSV_NO_HYDROUS = (
    SCRIPT_DIR / "perplex_mixed_fluid_parity_results_no_hydrous_water_activity.csv"
)
DEFAULT_MEEMUM = Path(r"C:\Program Files (x86)\Perplex\meemum.exe")
TEMPLATE_DAT = Path(
    r"C:\Users\stanroozen\Documents\Projects\Perplex\Perple_X\test\weigang\gfsm_fluid_probe.dat"
)
PERPLEX_DATAFILES_DIR = Path(
    r"C:\Users\stanroozen\Documents\Projects\Perplex\Perple_X\datafiles"
)
PERPLEX_DATABASE_FILE = "DEW24HP62ver_elements.dat"
# Perple_X mixed-fluid solution model used in this parity harness.
PERPLEX_FLUID_PHASE = "COH-Fluid+"
PERPLEX_HYBRID_EOS_H2O = "6"
PERPLEX_HYBRID_EOS_CO2 = "7"


@dataclass(frozen=True)
class Case:
    mineral: str
    name: str
    t_c: float
    p_kbar: float
    xco2: float


CASE_MATRIX = [
    Case("Quartz", "quartz_q37", 800.0, 9.0, 0.157),
    Case("Quartz", "quartz_q34", 800.0, 9.0, 0.350),
    Case("Quartz", "quartz_q33", 800.0, 9.0, 0.500),
    Case("Quartz", "quartz_c4", 800.0, 10.0, 0.337),
    Case("Quartz", "quartz_c3", 800.0, 10.0, 0.676),
    Case("Calcite", "calcite_500C_5kbar_x010", 500.0, 5.0, 0.10),
    Case("Calcite", "calcite_600C_8kbar_x010", 600.0, 8.0, 0.10),
    Case("Calcite", "calcite_700C_10kbar_x010", 700.0, 10.0, 0.10),
    Case("Calcite", "calcite_700C_10kbar_x250", 700.0, 10.0, 0.25),
]


MINERAL_CONFIGS = {
    "Quartz": {
        "mineral_phase": "Quartz",
        "seed_species": "SiO2_aq",
        "aq_species": [
            "WATER,AQ",
            "H+",
            "OH-",
            "SiO2_aq",
            "HSiO3-",
            "Si2O4_aq",
            "Si3O6_aq",
        ],
        "species_map": {
            "SiO2_aq": "SiO2,aq",
            "HSiO3-": "HSiO3-",
            "Si2O4_aq": "Si2O4,aq",
            "Si3O6_aq": "Si3O6,aq",
        },
        "element_coeffs": {
            "SiO2_aq": 1.0,
            "HSiO3-": 1.0,
            "Si2O4_aq": 2.0,
            "Si3O6_aq": 3.0,
        },
        "bulk_components": ("O2", "H2", "C", "Si"),
        "mineral_component": "Si",
        "mineral_component_add": {"O2": 1.0, "Si": 1.0},
    },
    "Calcite": {
        "mineral_phase": "Calcite",
        "seed_species": "Ca+2",
        "aq_species": [
            "WATER,AQ",
            "H+",
            "OH-",
            "Ca+2",
            "HCO3-",
            "CO3-2",
            "H2CO3_aq",
            "CaCO3_aq",
            "Ca(HCO3)+",
            "Ca(OH)+",
        ],
        "species_map": {
            "Ca+2": "Ca+2",
            "HCO3-": "HCO3-",
            "CO3-2": "CO3-2",
            "H2CO3_aq": "H2CO3,aq",
            "CaCO3_aq": "CaCO3,aq",
            "Ca(HCO3)+": "Ca(HCO3)+",
            "Ca(OH)+": "Ca(OH)+",
        },
        "element_coeffs": {
            "Ca+2": 1.0,
            "CaCO3_aq": 1.0,
            "Ca(HCO3)+": 1.0,
            "Ca(OH)+": 1.0,
        },
        "bulk_components": ("O2", "H2", "C", "Ca"),
        "mineral_component": "Ca",
        "mineral_component_add": {"O2": 1.5, "C": 1.0, "Ca": 1.0},
    },
}


def normalize_species_name(name: str) -> str:
    return name.replace(",aq", "_aq")


def safe_log10_delta(a: float, b: float, floor: float = 1e-30) -> float:
    import math

    return math.log10(max(a, floor)) - math.log10(max(b, floor))


def parse_fortran_float(token: str) -> float:
    token = token.strip()
    if not token:
        return 0.0
    try:
        return float(token)
    except ValueError:
        pass

    match = re.match(r"^([+-]?(?:\d+(?:\.\d*)?|\.\d+))([+-]\d+)$", token)
    if match:
        return float(f"{match.group(1)}E{match.group(2)}")
    raise ValueError(f"Could not parse numeric token: {token}")


def n_co2_for_xco2(xco2: float) -> float:
    if xco2 <= 0.0:
        return 0.0
    if xco2 >= 1.0:
        raise ValueError("xCO2 must be < 1.0")
    return xco2 / (1.0 - xco2) * H2O_MOL_PER_KG


def gfsm_gas_amounts(
    xco2: float, total_moles: float = N_GFSM_GAS_MOLES
) -> tuple[float, float]:
    if not 0.0 <= xco2 < 1.0:
        raise ValueError("xCO2 must be in [0, 1)")
    n_co2 = max(xco2, 1e-12) * total_moles
    n_h2o = max(1.0 - xco2, 1e-12) * total_moles
    return n_h2o, n_co2


def build_reaktoro_system(
    mineral: str,
    dew_db,
    supcrt_db,
    *,
    dh_model: str,
    enable_hydrous_species_correction: bool,
):
    cfg = MINERAL_CONFIGS[mineral]
    mineral_species = supcrt_db.species(cfg["mineral_phase"])
    h2o_species = supcrt_db.species("H2O(g)")
    co2_species = supcrt_db.species("CO2(g)")

    db = Database(dew_db.species())
    db.addSpecies(mineral_species)
    db.addSpecies(h2o_species)
    db.addSpecies(co2_species)

    # Use the exact same target species set for both PerplexDEW calculation and
    # Perple_X parsing/comparison. This keeps the parity comparison apples-to-apples.
    # Mixed-fluid H2O and CO2 are still supplied through the GFSM gas phase.
    target_species = ["WATER,AQ", "H+", "OH-"] + list(cfg["species_map"].keys())
    aq = AqueousPhase(" ".join(target_species))
    os.environ["REAKTORO_PERPLEXDEW_ENABLE_HYDROUS_SPECIES_CORRECTION"] = (
        "1" if enable_hydrous_species_correction else "0"
    )
    dh = ActivityDHModel.Davies
    if str(dh_model).strip().lower() == "extendeddh":
        dh = ActivityDHModel.ExtendedDH
    aq.setActivityModel(ActivityModelPerplexDEW(dh))

    gas = GaseousPhase("H2O(g) CO2(g)")
    gfsm_params = ActivityModelParamsPerplexGFSM()
    hybrid_opts = PerpleXHybridEosOptions()
    hybrid_opts.water = PerpleXWaterEos.ZhangDuan05
    hybrid_opts.co2 = PerpleXCO2Eos.ZhangDuan09
    gfsm_params.hybridEosOptions = hybrid_opts
    gas.setActivityModel(ActivityModelPerplexGFSM(gfsm_params))

    mineral_phase = MineralPhase(cfg["mineral_phase"])
    return ChemicalSystem(db, aq, gas, mineral_phase)


def solve_reaktoro_case(system, case: Case) -> dict[str, float]:
    cfg = MINERAL_CONFIGS[case.mineral]
    solver = EquilibriumSolver(system)
    conditions = EquilibriumConditions(system)
    conditions.temperature(case.t_c, "celsius")
    conditions.pressure(case.p_kbar * 1000.0, "bar")

    state = ChemicalState(system)
    state.set("WATER,AQ", 1.0, "kg")
    state.set("H+", 1e-8, "mol")
    state.set("OH-", 1e-8, "mol")
    state.set(cfg["seed_species"], 1e-6, "mol")
    state.set(cfg["mineral_phase"], MINERAL_MOLES, "mol")

    # Mixed-fluid composition comes directly from the GFSM EOS control variable
    # xCO2. No aqueous CO2 source species are seeded here.
    n_h2o, n_co2 = gfsm_gas_amounts(case.xco2)
    state.set("H2O(g)", n_h2o, "mol")
    state.set("CO2(g)", n_co2, "mol")

    result = solver.solve(state, conditions)
    if not result.succeeded():
        raise RuntimeError(f"PerplexDEW solve failed for {case.name}")

    aqp = AqueousProps(state)
    values = {}
    for species in cfg["species_map"]:
        try:
            values[species] = float(aqp.speciesMolality(species))
        except Exception:
            values[species] = 0.0

    total = 0.0
    for species, coeff in cfg["element_coeffs"].items():
        total += coeff * values.get(species, 0.0)
    values["total"] = total
    return values


def rewrite_dat_for_case(
    template: Path, target: Path, option_name: str, case: Case
) -> None:
    cfg = MINERAL_CONFIGS[case.mineral]
    n_co2 = n_co2_for_xco2(case.xco2)

    components = {
        "O2": 0.5 * H2O_MOL_PER_KG + n_co2,
        "H2": H2O_MOL_PER_KG,
        "C": n_co2,
    }
    for name, amount in cfg["mineral_component_add"].items():
        components[name] = components.get(name, 0.0) + amount * MINERAL_MOLES

    ordered_lines = []
    for component in cfg["bulk_components"]:
        ordered_lines.append(
            f"{component:<5} 1  {components.get(component, 0.0):.6f}      0.00000      0.00000     molar amount"
        )

    out = []
    lines = template.read_text(encoding="utf-8").splitlines()
    in_components = False
    for line in lines:
        stripped = line.strip()
        if stripped.endswith("thermodynamic data file"):
            out.append(f"{PERPLEX_DATABASE_FILE}     thermodynamic data file")
            continue
        if stripped == "gfsm_fluid_probe_perplex_option.dat     | Perple_X option file":
            out.append(f"{option_name}     | Perple_X option file")
            continue
        if stripped == "GFSM fluid probe":
            out.append(f"{case.mineral} mixed-fluid parity: {case.name}")
            continue
        if stripped == "COH-Fluid" or stripped == PERPLEX_FLUID_PHASE:
            out.append(PERPLEX_FLUID_PHASE)
            continue
        if stripped == "begin thermodynamic component list":
            in_components = True
            out.append(line)
            out.extend(ordered_lines)
            continue
        if stripped == "end thermodynamic component list":
            in_components = False
            out.append(line)
            continue
        if in_components:
            continue
        out.append(line)

    target.write_text("\n".join(out) + "\n", encoding="utf-8")


def configure_option_file(path: Path) -> None:
    text = path.read_text(encoding="utf-8")
    text = text.replace("warn_interactive T", "warn_interactive F")
    text = text.replace("pause_on_error T", "pause_on_error F")
    text = re.sub(
        r"(?m)^hybrid_EoS_H2O\s+\S+",
        f"hybrid_EoS_H2O {PERPLEX_HYBRID_EOS_H2O}",
        text,
    )
    text = re.sub(
        r"(?m)^hybrid_EoS_CO2\s+\S+",
        f"hybrid_EoS_CO2 {PERPLEX_HYBRID_EOS_CO2}",
        text,
    )
    path.write_text(text, encoding="utf-8")


def resolve_perplex_database_file(template_dir: Path) -> Path:
    candidates = [
        template_dir / PERPLEX_DATABASE_FILE,
        PERPLEX_DATAFILES_DIR / PERPLEX_DATABASE_FILE,
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    raise FileNotFoundError(
        f"Required Perple_X thermodynamic database not found: {PERPLEX_DATABASE_FILE}"
    )


def run_meemum(
    meemum_exe: Path, work_dir: Path, project_root: str, p_bar: float, t_k: float
) -> Path:
    # With warn_interactive=F, Perple_X auto-accepts the fluid warning and then
    # prompts only for whether to enter bulk compositions interactively. We keep
    # using the bulk composition embedded in the .dat file by answering "n".
    input_text = f"{project_root}\nn\n{t_k} {p_bar}\n0 0\n"
    input_file = work_dir / f"{project_root}_input.txt"
    input_file.write_text(input_text, encoding="utf-8")

    with input_file.open("r", encoding="utf-8") as stdin_fh:
        res = subprocess.run(
            [str(meemum_exe)],
            cwd=work_dir,
            stdin=stdin_fh,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            timeout=30,
        )
    stdout_text = res.stdout.decode(errors="ignore")
    if res.returncode != 0:
        raise RuntimeError(
            f"meemum failed for {project_root} (rc={res.returncode}):\n{stdout_text}"
        )

    prn_path = work_dir / f"{project_root}.prn"
    if not prn_path.exists():
        raise FileNotFoundError(f"Expected meemum output not found: {prn_path}")
    return prn_path


def parse_perplex_prn(
    prn_path: Path, mineral: str
) -> tuple[dict[str, float], bool, bool]:
    cfg = MINERAL_CONFIGS[mineral]
    text = prn_path.read_text(encoding="utf-8", errors="ignore")

    if not text.strip():
        raise RuntimeError(f"Perple_X output is empty: {prn_path}")

    has_fluid = (
        re.search(
            rf"^\s*{re.escape(PERPLEX_FLUID_PHASE)}(?:\s|$)",
            text,
            re.MULTILINE,
        )
        is not None
    )
    has_mineral = False

    phase_block = re.search(
        r"Phase Compositions \(molar  proportions\):\s*(.*?)\n\s*Phase speciation",
        text,
        re.DOTALL,
    )
    if phase_block:
        lines = [line for line in phase_block.group(1).splitlines() if line.strip()]
        if len(lines) >= 2:
            component = cfg["mineral_component"]
            components = list(cfg["bulk_components"])
            if component in components:
                comp_idx = components.index(component)
                for line in lines[1:]:
                    tokens = line.split()
                    if len(tokens) < 1 + len(components):
                        continue
                    phase_name = tokens[0]
                    if phase_name in {"COH-Fluid", PERPLEX_FLUID_PHASE}:
                        continue
                    try:
                        comp_tokens = tokens[-len(components) :]
                        comp_val = parse_fortran_float(comp_tokens[comp_idx])
                    except ValueError:
                        continue
                    if comp_val > 0.0:
                        has_mineral = True
                        break

    # Perple_X may print multiple immiscible fluid branches, each with its own
    # "Solute endmember properties" block. We aggregate all valid blocks using
    # fluid-phase mol% weights from "Phase Compositions" to obtain a chemically
    # sensible system-scale dissolved amount.
    solute_blocks = list(
        re.finditer(
            r"Solute endmember properties:\s*(.*?)\n\s*(?:Normalized solvent endmember properties:|Solvent endmember properties:)",
            text,
            re.DOTALL,
        )
    )
    if not solute_blocks:
        raise RuntimeError(f"No solute section found in Perple_X output: {prn_path}")

    fluid_weights = []
    if phase_block:
        fluid_mol_percents = []
        for line in phase_block.group(1).splitlines():
            tokens = line.split()
            if len(tokens) < 4:
                continue
            if tokens[0] != PERPLEX_FLUID_PHASE:
                continue
            try:
                fluid_mol_percents.append(parse_fortran_float(tokens[3]))
            except ValueError:
                continue
        fluid_weights = fluid_mol_percents

    species_by_perplex = {
        perplex_name: reaktoro_name
        for reaktoro_name, perplex_name in cfg["species_map"].items()
    }
    # Perple_X may print certain charged species names without explicit '+' suffix.
    species_aliases = {
        "Ca(HCO3)": "Ca(HCO3)+",
    }
    block_values = []
    for match in solute_blocks:
        values_i = {species: 0.0 for species in cfg["species_map"]}
        for line in match.group(1).splitlines():
            m = re.match(r"^\s*(\S+)\s+[-+]?\d+\s+([0-9.Ee+\-]+)", line)
            if not m:
                continue
            name = m.group(1)
            molality = parse_fortran_float(m.group(2))
            reaktoro_name = species_by_perplex.get(name)
            if reaktoro_name is None:
                alias = species_aliases.get(name)
                if alias is not None:
                    reaktoro_name = species_by_perplex.get(alias)
            if reaktoro_name:
                values_i[reaktoro_name] = molality
        block_values.append(values_i)

    nblocks = len(block_values)
    if fluid_weights and len(fluid_weights) >= nblocks:
        weights = fluid_weights[:nblocks]
    else:
        weights = [1.0] * nblocks

    wsum = sum(w for w in weights if w > 0.0)
    if wsum <= 0.0:
        raise RuntimeError(
            f"Could not determine positive fluid weights in Perple_X output: {prn_path}"
        )

    values = {species: 0.0 for species in cfg["species_map"]}
    for species in values:
        accum = 0.0
        for idx, vals in enumerate(block_values):
            accum += weights[idx] * vals.get(species, 0.0)
        values[species] = accum / wsum

    total = 0.0
    for species, coeff in cfg["element_coeffs"].items():
        total += coeff * values.get(species, 0.0)
    values["total"] = total
    return values, has_mineral, has_fluid


def run_perplex_case(
    case: Case, meemum_exe: Path, keep_work_dir: bool
) -> tuple[dict[str, float], bool, bool]:
    if not TEMPLATE_DAT.exists():
        raise FileNotFoundError(f"Template .dat not found: {TEMPLATE_DAT}")

    parent = tempfile.mkdtemp(prefix=f"perplex_{case.name}_", dir=str(SCRIPT_DIR))
    work_dir = Path(parent)
    try:
        template_dir = TEMPLATE_DAT.parent
        db_src = resolve_perplex_database_file(template_dir)
        for fname in [
            "solution_model.dat",
            "gfsm_fluid_probe_perplex_option.dat",
        ]:
            src = template_dir / fname
            if not src.exists():
                raise FileNotFoundError(f"Required Perple_X file not found: {src}")
            shutil.copy(src, work_dir / fname)
        shutil.copy(db_src, work_dir / db_src.name)

        project_root = case.name
        dat_path = work_dir / f"{project_root}.dat"
        option_path = work_dir / f"{project_root}_perplex_option.dat"
        shutil.copy(work_dir / "gfsm_fluid_probe_perplex_option.dat", option_path)
        configure_option_file(option_path)
        rewrite_dat_for_case(TEMPLATE_DAT, dat_path, option_path.name, case)

        prn_path = run_meemum(
            meemum_exe, work_dir, project_root, case.p_kbar * 1000.0, case.t_c + 273.15
        )
        return parse_perplex_prn(prn_path, case.mineral)
    finally:
        if keep_work_dir:
            print(f"Kept Perple_X work dir: {work_dir}")
        else:
            shutil.rmtree(work_dir, ignore_errors=True)


def compare_cases(
    meemum_exe: Path,
    out_csv: Path,
    fail_threshold_log10: float | None,
    keep_work_dir: bool,
    dh_model: str,
    enable_hydrous_species_correction: bool,
) -> int:
    dew_db = DEWDatabase("dew2024-aqueous")
    supcrt_db = SupcrtDatabase("supcrtbl")
    systems = {
        mineral: build_reaktoro_system(
            mineral,
            dew_db,
            supcrt_db,
            dh_model=dh_model,
            enable_hydrous_species_correction=enable_hydrous_species_correction,
        )
        for mineral in {case.mineral for case in CASE_MATRIX}
    }

    print(
        "Parity configuration: "
        f"fluid phase={PERPLEX_FLUID_PHASE}; "
        f"PerplexDEW DH={dh_model}; "
        f"hydrous_correction={'on' if enable_hydrous_species_correction else 'off'}; "
        "water_activity_correction=on"
    )

    fieldnames = [
        "case",
        "mineral",
        "T_C",
        "P_kbar",
        "xco2",
        "stable_in_perplex",
        "fluid_in_perplex",
        "error",
        "metric",
        "perplexdew_m",
        "perplex_m",
        "delta_log10",
    ]
    rows = []
    violations = []

    for case in CASE_MATRIX:
        print(
            f"[{case.mineral}] {case.name}: T={case.t_c:.1f} C  P={case.p_kbar:.1f} kbar  xCO2={case.xco2:.3f}"
        )
        rt_values = solve_reaktoro_case(systems[case.mineral], case)
        case_error = ""
        try:
            px_values, has_mineral, has_fluid = run_perplex_case(
                case, meemum_exe, keep_work_dir
            )
        except Exception as exc:
            case_error = str(exc)
            has_mineral, has_fluid = False, False
            px_values = {
                k: float("nan") for k in MINERAL_CONFIGS[case.mineral]["species_map"]
            }
            px_values["total"] = float("nan")
            print(f"  WARNING: Perple_X case failed: {case_error}")

        metrics = ["total"] + list(MINERAL_CONFIGS[case.mineral]["species_map"].keys())

        for metric in metrics:
            rt_val = rt_values.get(metric, 0.0)
            px_val = px_values.get(metric, float("nan"))
            delta = (
                safe_log10_delta(rt_val, px_val)
                if math.isfinite(px_val)
                else float("nan")
            )
            rows.append(
                {
                    "case": case.name,
                    "mineral": case.mineral,
                    "T_C": case.t_c,
                    "P_kbar": case.p_kbar,
                    "xco2": case.xco2,
                    "stable_in_perplex": has_mineral,
                    "fluid_in_perplex": has_fluid,
                    "error": case_error,
                    "metric": metric,
                    "perplexdew_m": rt_val,
                    "perplex_m": px_val,
                    "delta_log10": delta,
                }
            )
            if (
                fail_threshold_log10 is not None
                and math.isfinite(delta)
                and abs(delta) > fail_threshold_log10
            ):
                violations.append((case.name, metric, delta))

        total_delta = next(
            row["delta_log10"]
            for row in rows[::-1]
            if row["case"] == case.name and row["metric"] == "total"
        )
        if math.isfinite(px_values.get("total", float("nan"))):
            print(
                f"  total dissolved: PerplexDEW={rt_values['total']:.6e}  Perple_X={px_values['total']:.6e}  delta_log10={total_delta:+.3f}"
            )
        else:
            print(
                f"  total dissolved: PerplexDEW={rt_values['total']:.6e}  Perple_X=NaN  delta_log10=NaN"
            )

    with out_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nWrote {out_csv}")
    if violations:
        print("Threshold violations:")
        for case_name, metric, delta in violations:
            print(f"  {case_name}: {metric} delta_log10={delta:+.3f}")
        return 1
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Direct Perple_X vs PerplexDEW mixed-fluid parity checks"
    )
    parser.add_argument("--meemum-exe", default=str(DEFAULT_MEEMUM))
    parser.add_argument("--out-csv", default=str(OUTPUT_CSV))
    parser.add_argument(
        "--dh-model",
        choices=["Davies", "ExtendedDH"],
        default="Davies",
        help="Debye-Huckel model for PerplexDEW aqueous activity.",
    )
    parser.add_argument(
        "--disable-hydrous-species-correction",
        action="store_true",
        help="Disable the PerplexDEW hydrous-species ln(a_H2O) correction.",
    )
    parser.add_argument(
        "--fail-threshold-log10",
        type=float,
        default=None,
        help="Optional failure threshold on absolute log10 difference.",
    )
    parser.add_argument(
        "--keep-work-dir",
        action="store_true",
        help="Keep temporary Perple_X work directories for inspection.",
    )
    args = parser.parse_args()

    meemum_exe = Path(args.meemum_exe)
    if not meemum_exe.exists():
        raise FileNotFoundError(f"meemum executable not found: {meemum_exe}")

    return compare_cases(
        meemum_exe=meemum_exe,
        out_csv=Path(args.out_csv),
        fail_threshold_log10=args.fail_threshold_log10,
        keep_work_dir=args.keep_work_dir,
        dh_model=args.dh_model,
        enable_hydrous_species_correction=not args.disable_hydrous_species_correction,
    )


if __name__ == "__main__":
    raise SystemExit(main())

