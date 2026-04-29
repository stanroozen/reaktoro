#!/usr/bin/env python3
"""
Generate GFSM reference CSV files using Perple_X meemum output.

This script uses a known-working GFSM probe setup (project .dat + option file)
and regenerates reference CSVs for pure-EOS coverage and representative mixtures.

Outputs are intentionally volume/speciation-focused, because meemum .prn output
for this workflow does not expose per-species fugacity columns in the same format
as the legacy fluids.exe tables.
"""

from __future__ import annotations

import argparse
import csv
import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path


SPECIES_ORDER = [
    "H2O",
    "CO2",
    "CH4",
    "H2",
    "CO",
    "H2S",
    "O2",
    "SO2",
    "N2",
    "NH3",
    "HF",
    "C2H6",
    "HCl",
]


@dataclass
class Case:
    name: str
    h2o_eos: int
    co2_eos: int
    ch4_eos: int
    # Bulk elements for the probe file component list (O2, H2, C)
    n_o2: float
    n_h2: float
    n_c: float
    p_bar: float = 1000.0
    t_k: float = 523.15


# All EOS variants available per GFSM-hybrid species (Perple_X iopt 25/26/27).
# Keys are short names used in case/file naming; values are the integer EOS codes.
H2O_EOS: dict[str, int] = {
    "mrk": 0,  # MRK (default)
    "hsmrk": 1,  # Halbach & Chatterjee 1982
    "cork": 2,  # Holland & Powell 1991
    "pseos": 4,  # Pitzer & Sterner 1994
    "haar": 5,  # Haar et al. 1982
    "zd05": 6,  # Zhang & Duan 2005
    "zd09": 7,  # Zhang & Duan 2009
}
CO2_EOS: dict[str, int] = {
    "mrk": 0,  # MRK (default)
    "hsmrk": 1,  # Halbach & Chatterjee 1982
    "cork": 2,  # Holland & Powell 1991
    # "brmrk": 3 — Bottinga & Richet 1981, EXCLUDED: produces IEEE_INVALID_FLAG
    #              (NaN in BRVOL routine) for pure CO2 at all tested T/P conditions.
    "pseos": 4,  # Pitzer & Sterner 1994
    "zd09": 7,  # Zhang & Duan 2009
}
CH4_EOS: dict[str, int] = {
    "mrk": 0,  # MRK (default)
    "hsmrk": 1,  # Kerrick & Jacobs 1981
    "zd09": 7,  # Zhang & Duan 2009
}

# ---------------------------------------------------------------------------
# Bulk compositions (O2-H2-C component amounts) for each mixture type.
# Derived from stoichiometry: H2O ← 0.5 O2 + H2,  CO2 ← O2 + C,  CH4 ← 2H2 + C
# Pure H2O  : n_o2=0.50, n_h2=1.00, n_c=0.00
# Pure CO2  : n_o2=1.00, n_h2=0.00, n_c=1.00
# Pure CH4  : n_o2=0.00, n_h2=2.00, n_c=1.00
# 50% H2O + 50% CO2 : n_o2=0.75, n_h2=0.50, n_c=0.50
# 50% H2O + 50% CH4 : n_o2=0.25, n_h2=1.50, n_c=0.50
# 50% CO2 + 50% CH4 : n_o2=0.50, n_h2=1.00, n_c=1.00
# 33% each           : n_o2=0.50, n_h2=1.00, n_c=2/3
# ---------------------------------------------------------------------------


def _pure_h2o_cases() -> list[Case]:
    return [Case(f"h2o_{n}", c, 0, 0, 0.5, 1.0, 0.0) for n, c in H2O_EOS.items()]


def _pure_co2_cases() -> list[Case]:
    return [Case(f"co2_{n}", 0, c, 0, 1.0, 0.0, 1.0) for n, c in CO2_EOS.items()]


def _pure_ch4_cases() -> list[Case]:
    return [Case(f"ch4_{n}", 0, 0, c, 0.0, 2.0, 1.0) for n, c in CH4_EOS.items()]


def _binary_h2o_co2_cases() -> list[Case]:
    """Vary H2O EOS (CO2=MRK) then vary CO2 EOS (H2O=MRK, skip MRK duplicate)."""
    cases = [
        Case(f"h2o_co2_{h2o_n}_mrk", h2o_c, 0, 0, 0.75, 0.5, 0.5)
        for h2o_n, h2o_c in H2O_EOS.items()
    ]
    cases += [
        Case(f"h2o_co2_mrk_{co2_n}", 0, co2_c, 0, 0.75, 0.5, 0.5)
        for co2_n, co2_c in CO2_EOS.items()
        if co2_n != "mrk"  # mrk+mrk already generated above
    ]
    return cases


def _binary_h2o_ch4_cases() -> list[Case]:
    """Vary H2O EOS (CH4=MRK) then vary CH4 EOS (H2O=MRK, skip MRK duplicate)."""
    cases = [
        Case(f"h2o_ch4_{h2o_n}_mrk", h2o_c, 0, 0, 0.25, 1.5, 0.5)
        for h2o_n, h2o_c in H2O_EOS.items()
    ]
    cases += [
        Case(f"h2o_ch4_mrk_{ch4_n}", 0, 0, ch4_c, 0.25, 1.5, 0.5)
        for ch4_n, ch4_c in CH4_EOS.items()
        if ch4_n != "mrk"
    ]
    return cases


def _binary_co2_ch4_cases() -> list[Case]:
    """Vary CO2 EOS (CH4=MRK) then vary CH4 EOS (CO2=MRK, skip MRK duplicate)."""
    cases = [
        Case(f"co2_ch4_{co2_n}_mrk", 0, co2_c, 0, 0.5, 1.0, 1.0)
        for co2_n, co2_c in CO2_EOS.items()
    ]
    cases += [
        Case(f"co2_ch4_mrk_{ch4_n}", 0, 0, ch4_c, 0.5, 1.0, 1.0)
        for ch4_n, ch4_c in CH4_EOS.items()
        if ch4_n != "mrk"
    ]
    return cases


def _ternary_cases() -> list[Case]:
    """Equal-thirds H2O+CO2+CH4 with matched EOS codes across all three species."""
    # 1/3 each → n_o2=0.5, n_h2=1.0, n_c=2/3
    _n_c = 2.0 / 3.0
    return [
        Case(f"h2o_co2_ch4_{n}", h2o_c, co2_c, ch4_c, 0.5, 1.0, _n_c)
        for n, (h2o_c, co2_c, ch4_c) in {
            "mrk": (0, 0, 0),
            "hsmrk": (1, 1, 1),
            "zd09": (7, 7, 7),
        }.items()
    ]


def _o2_and_redox_stress_cases() -> list[Case]:
    """OHC-only stress states to exercise O2/H2/CO-rich compositions in GFSM output.

    Note: this generator uses an O2-H2-C bulk-component template; therefore these
    cases can exercise O2/H2/CO behavior directly, but not N/S/F/Cl-bearing species
    without a different template and expanded component basis.
    """
    return [
        # Fallback defaults used if no matrix file is provided.
        Case("redox_o2_excess_mrk", 0, 0, 0, 2.0, 0.4, 0.2, p_bar=1000.0, t_k=523.15),
        Case("redox_o2_excess_hot", 0, 0, 0, 2.0, 0.5, 0.2, p_bar=1000.0, t_k=973.15),
        Case("redox_h2_excess_hot", 0, 0, 0, 0.2, 3.0, 0.2, p_bar=1000.0, t_k=973.15),
        Case("redox_co_bias_hot", 0, 0, 0, 0.6, 0.4, 1.0, p_bar=1000.0, t_k=973.15),
        Case(
            "redox_o2_excess_hybrid", 1, 4, 0, 1.8, 0.6, 0.4, p_bar=1500.0, t_k=873.15
        ),
        Case("redox_mixed_hybrid", 7, 4, 1, 0.9, 1.2, 0.8, p_bar=1500.0, t_k=873.15),
    ]


def load_case_matrix(path: Path) -> list[Case]:
    """Load additional/meemum-specific GFSM cases from CSV.

    Required columns:
      case,h2o_eos,co2_eos,ch4_eos,n_o2,n_h2,n_c,p_bar,t_k
    """
    if not path.exists():
        return []

    loaded: list[Case] = []
    with path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row.get("case", "").strip()
            if not name:
                continue
            loaded.append(
                Case(
                    name=name,
                    h2o_eos=int(row["h2o_eos"]),
                    co2_eos=int(row["co2_eos"]),
                    ch4_eos=int(row["ch4_eos"]),
                    n_o2=float(row["n_o2"]),
                    n_h2=float(row["n_h2"]),
                    n_c=float(row["n_c"]),
                    p_bar=float(row.get("p_bar", 1000.0)),
                    t_k=float(row.get("t_k", 523.15)),
                )
            )

    return loaded


CASES: list[Case] = (
    _pure_h2o_cases()  # 7  pure H2O  (mrk + 6 non-MRK EOS variants)
    + _pure_co2_cases()  # 5  pure CO2  (mrk + 4 non-MRK, brmrk excluded)
    + _pure_ch4_cases()  # 3  pure CH4  (mrk + hsmrk + zd09)
    + _binary_h2o_co2_cases()  # 11 binary H2O+CO2
    + _binary_h2o_ch4_cases()  # 9  binary H2O+CH4
    + _binary_co2_ch4_cases()  # 7  binary CO2+CH4
    + _ternary_cases()  # 3  ternary   H2O+CO2+CH4
)  # 45 base cases


def write_option_file(path: Path, h2o_eos: int, co2_eos: int, ch4_eos: int) -> None:
    text = "\n".join(
        [
            "GFSM T",
            "species_output T",
            f"hybrid_EoS_H2O {h2o_eos}",
            f"hybrid_EoS_CO2 {co2_eos}",
            f"hybrid_EoS_CH4 {ch4_eos}",
            # warn_interactive F: auto-accept missing-endmember and unpaired
            # species warnings without requiring a Y/N response at runtime.
            # The warning text is still written to stdout but execution continues.
            "warn_interactive F",
            "pause_on_error F",
            "",
        ]
    )
    path.write_text(text, encoding="utf-8")


def rewrite_dat(
    template: Path,
    target: Path,
    option_filename: str,
    n_o2: float,
    n_h2: float,
    n_c: float,
    title: str,
) -> None:
    lines = template.read_text(encoding="utf-8").splitlines()

    # Replace title line and option file line, plus component amounts.
    out = []
    in_component_list = False
    for line in lines:
        if "GFSM fluid probe" in line:
            out.append(title)
            continue
        if line.strip().endswith("| Perple_X option file"):
            out.append(f"{option_filename}     | Perple_X option file")
            continue

        if line.strip() == "begin thermodynamic component list":
            in_component_list = True
            out.append(line)
            continue
        if line.strip() == "end thermodynamic component list":
            in_component_list = False
            out.append(line)
            continue

        if in_component_list and line.lstrip().startswith("O2"):
            out.append(
                f"O2    1  {n_o2:.6f}      0.00000      0.00000     molar amount"
            )
            continue
        if in_component_list and line.lstrip().startswith("H2"):
            out.append(
                f"H2    1  {n_h2:.6f}      0.00000      0.00000     molar amount"
            )
            continue
        if in_component_list and re.match(r"^\s*C\s", line):
            out.append(f"C     1  {n_c:.6f}      0.00000      0.00000     molar amount")
            continue

        out.append(line)

    target.write_text("\n".join(out) + "\n", encoding="utf-8")


def run_meemum(
    meemum_exe: Path, work_dir: Path, project_root: str, p_bar: float, t_k: float
) -> None:
    # Write input to a file and redirect from it — required on Windows because
    # Fortran programs (Intel runtime) use console APIs that don't work with
    # anonymous pipes (subprocess input= parameter).
    input_text = f"{project_root}\nn\n{p_bar} {t_k}\n0 0\n"
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
    # Log output for debugging — will be empty on clean runs.
    if stdout_text.strip():
        print(f"[meemum {project_root}] {stdout_text[:500]}")


def parse_prn(prn_path: Path) -> tuple[dict[str, float], float, float, float]:
    text = prn_path.read_text(encoding="utf-8", errors="ignore")

    # Extract stable P and T from report header block.
    p_match = re.search(r"P\(bar\)\s*=\s*([0-9.Ee+\-]+)", text)
    t_match = re.search(r"T\(K\)\s*=\s*([0-9.Ee+\-]+)", text)
    if not p_match or not t_match:
        raise RuntimeError(f"Could not parse P/T from {prn_path}")
    p_bar = float(p_match.group(1))
    t_k = float(t_match.group(1))

    # Phase speciation line for COH-Fluid.
    spec_line = None
    for line in text.splitlines():
        if line.strip().startswith("COH-Fluid") and "CO2:" in line and "H2O:" in line:
            spec_line = line
            break
    if spec_line is None:
        raise RuntimeError(f"Could not find COH-Fluid speciation in {prn_path}")

    species_values = {name: 0.0 for name in SPECIES_ORDER}
    for pair in spec_line.split(":"):
        # Recover token by scanning with regex on full line instead.
        pass
    for m in re.finditer(r"([A-Za-z0-9]+):\s*([0-9.Ee+\-]+)", spec_line):
        nm = m.group(1)
        val = float(m.group(2))
        if nm in species_values:
            species_values[nm] = val

    # Molar volume V(J/bar) in COH-Fluid row in Molar Properties block.
    vol = None
    in_molar = False
    for line in text.splitlines():
        if line.strip().startswith("Molar Properties and Density"):
            in_molar = True
            continue
        if in_molar and line.strip().startswith("COH-Fluid"):
            tokens = line.split()
            # Expected columns include: N(g), G(J), S(J/K), V(J/bar), ...
            # V should be token index 4 when split by whitespace.
            if len(tokens) >= 5:
                # V column is in J/bar; convert to cm³/mol (1 J/bar = 10 cm³)
                vol = float(tokens[4]) * 10.0
            break
    if vol is None:
        raise RuntimeError(f"Could not parse COH-Fluid molar volume in {prn_path}")

    return species_values, vol, p_bar, t_k


def write_csv(
    path: Path, species_values: dict[str, float], vol: float, p_bar: float, t_k: float
) -> None:
    headers = ["P(bar)", "T(K)"] + [f"y({s})" for s in SPECIES_ORDER] + ["vol[cm3/mol]"]
    row = [p_bar, t_k] + [species_values[s] for s in SPECIES_ORDER] + [vol]

    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(headers)
        w.writerow(row)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate GFSM reference CSVs via Perple_X meemum"
    )
    parser.add_argument(
        "--meemum-exe", default=r"C:\Program Files (x86)\Perplex\meemum.exe"
    )
    parser.add_argument(
        "--template-dat",
        default=r"C:\Users\stanroozen\Documents\Projects\Perplex\Perple_X\test\weigang\gfsm_fluid_probe.dat",
    )
    parser.add_argument(
        "--out-dir",
        default=str(Path(__file__).resolve().parent / "test" / "gfsm"),
    )
    parser.add_argument(
        "--case-matrix",
        default=str(Path(__file__).resolve().parent / "test" / "gfsm_case_matrix.csv"),
        help="Optional CSV with additional GFSM meemum cases (e.g., O2/redox stress states)",
    )
    args = parser.parse_args()

    meemum_exe = Path(args.meemum_exe)
    template_dat = Path(args.template_dat)
    out_dir = Path(args.out_dir)
    case_matrix = Path(args.case_matrix)

    if not meemum_exe.exists():
        raise FileNotFoundError(f"meemum.exe not found: {meemum_exe}")
    if not template_dat.exists():
        raise FileNotFoundError(f"Template .dat not found: {template_dat}")

    out_dir.mkdir(parents=True, exist_ok=True)

    work_dir = out_dir / "_tmp"
    if work_dir.exists():
        shutil.rmtree(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)

    # Copy the thermodynamic data files that meemum requires in its working dir.
    template_dir = template_dat.parent
    for fname in ["DEW19HP633ver_elements.dat", "solution_model.dat"]:
        src = template_dir / fname
        if not src.exists():
            raise FileNotFoundError(f"Required Perple_X data file not found: {src}")
        shutil.copy(src, work_dir / fname)

    # Load base coverage plus explicit matrix-driven extras (user-extendable).
    matrix_cases = load_case_matrix(case_matrix)
    if not matrix_cases:
        matrix_cases = _o2_and_redox_stress_cases()

    # Preserve insertion order, matrix cases can override by name.
    all_cases_map: dict[str, Case] = {c.name: c for c in CASES}
    for c in matrix_cases:
        all_cases_map[c.name] = c

    for case in all_cases_map.values():
        case_root = case.name
        dat_path = work_dir / f"{case_root}.dat"
        opt_path = work_dir / f"{case_root}_perplex_option.dat"

        write_option_file(opt_path, case.h2o_eos, case.co2_eos, case.ch4_eos)
        rewrite_dat(
            template_dat,
            dat_path,
            opt_path.name,
            case.n_o2,
            case.n_h2,
            case.n_c,
            title=f"GFSM {case.name}",
        )

        run_meemum(meemum_exe, work_dir, case_root, case.p_bar, case.t_k)

        prn_path = work_dir / f"{case_root}.prn"
        if not prn_path.exists():
            raise FileNotFoundError(f"Expected meemum output not found: {prn_path}")

        species_values, vol, p_bar, t_k = parse_prn(prn_path)
        out_csv = out_dir / f"{case.name}.csv"
        write_csv(out_csv, species_values, vol, p_bar, t_k)
        print(f"Wrote {out_csv}")


if __name__ == "__main__":
    main()
