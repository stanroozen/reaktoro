#!/usr/bin/env python3
"""
Generate Perple_X reference CSVs from an executable test matrix.

Reads test/test_matrix.csv and runs Perple_X fluids.exe for each row.
Outputs one CSV per case under test/ (e.g., test/h2o_co2_mrk.csv).

Usage:
    python generate_reference_matrix.py --fluids-exe "C:\\Program Files (x86)\\Perplex\\fluids.exe"
"""

from __future__ import annotations

import argparse
import csv
import os
import subprocess
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple


@dataclass
class MatrixRow:
    case: str
    ifug: int
    P_bar: float
    T_K: float
    XCO2: float


def read_matrix(path: Path) -> List[MatrixRow]:
    rows: List[MatrixRow] = []
    with path.open() as f:
        reader = csv.DictReader(f)
        for r in reader:
            rows.append(
                MatrixRow(
                    case=r["case"].strip(),
                    ifug=int(r["ifug"]),
                    P_bar=float(r["P_bar"]),
                    T_K=float(r["T_K"]),
                    XCO2=float(r["XCO2"]),
                )
            )
    return rows


def build_input_text(
    ifug: int, P_bar: float, T_K: float, XCO2: float, out_root: str
) -> str:
    # Tabulate properties; 2 independent variables (P, T); single point
    lines = []
    lines.append(str(ifug))  # EoS choice
    lines.append("y")  # tabulate output
    lines.append("2")  # two independent vars
    lines.append("1")  # P
    lines.append("2")  # T
    lines.append(f"{P_bar:.6f} {P_bar:.6f} 1.0")
    lines.append(f"{T_K:.6f} {T_K:.6f} 50.0")
    lines.append(f"{XCO2:.6f}")  # sectioning variable X(CO2)
    lines.append(out_root)
    lines.append(out_root)
    lines.append("n")  # log output? no
    lines.append("n")  # more calculations? no
    return "\n".join(lines) + "\n"


def parse_tab_file(path: Path) -> Tuple[List[str], List[List[float]]]:
    lines = path.read_text().splitlines()
    if len(lines) < 6:
        raise ValueError(f"Unexpected .tab format in {path}")

    idx = 0
    idx += 1  # version
    idx += 1  # title
    jpot = int(lines[idx].strip())
    idx += 1

    # Skip independent variable metadata blocks (4 lines each)
    for _ in range(jpot):
        idx += 4

    kount = int(lines[idx].strip())
    idx += 1

    # Read column tags
    tags: List[str] = []
    while len(tags) < kount and idx < len(lines):
        tokens = lines[idx].split()
        tags.extend(tokens)
        idx += 1

    if len(tags) != kount:
        raise ValueError(f"Failed to read {kount} column tags in {path}")

    # Read data values
    values: List[float] = []
    for line in lines[idx:]:
        tokens = line.split()
        for token in tokens:
            values.append(float(token))

    rows: List[List[float]] = []
    for i in range(0, len(values), kount):
        chunk = values[i : i + kount]
        if len(chunk) != kount:
            break
        rows.append(chunk)

    return tags, rows


def run_case(
    fluids_exe: Path, work_dir: Path, row: MatrixRow
) -> Tuple[List[str], List[float]]:
    out_root = f"{row.case}_XCO2_{str(row.XCO2).replace('.', 'p')}"
    input_text = build_input_text(row.ifug, row.P_bar, row.T_K, row.XCO2, out_root)

    result = subprocess.run(
        [str(fluids_exe)],
        input=input_text.encode("utf-8"),
        cwd=work_dir,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    if result.returncode != 0:
        raise RuntimeError(
            f"fluids.exe failed for {row.case} XCO2={row.XCO2}\n"
            f"stdout:\n{result.stdout.decode('utf-8', errors='ignore')}\n"
            f"stderr:\n{result.stderr.decode('utf-8', errors='ignore')}"
        )

    tab_path = work_dir / f"{out_root}.tab"
    if not tab_path.exists():
        raise FileNotFoundError(f"Expected output not found: {tab_path}")

    tags, rows = parse_tab_file(tab_path)
    if not rows:
        raise RuntimeError(f"No data rows in {tab_path}")

    # Single row for single point
    return tags, rows[0]


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate Perple_X reference CSVs from test matrix"
    )
    parser.add_argument(
        "--fluids-exe",
        default=r"C:\Program Files (x86)\Perplex\fluids.exe",
        help="Path to Perple_X fluids executable",
    )
    parser.add_argument(
        "--matrix",
        default=str((Path(__file__).parent / "test/test_matrix.csv").resolve()),
        help="Path to test matrix CSV",
    )
    parser.add_argument(
        "--out-dir",
        default=str((Path(__file__).parent / "test").resolve()),
        help="Output directory for reference CSVs",
    )

    args = parser.parse_args()

    fluids_exe = Path(args.fluids_exe)
    if not fluids_exe.exists():
        raise FileNotFoundError(f"fluids.exe not found: {fluids_exe}")

    matrix_path = Path(args.matrix)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = read_matrix(matrix_path)
    work_dir = out_dir

    grouped: Dict[str, List[MatrixRow]] = defaultdict(list)
    for r in rows:
        grouped[r.case].append(r)

    for case, case_rows in grouped.items():
        combined_tags: List[str] = []
        combined_rows: List[List[float]] = []

        for r in case_rows:
            tags, data = run_case(fluids_exe, work_dir, r)
            if not combined_tags:
                combined_tags = tags
            combined_rows.append(data)

        out_csv = out_dir / f"{case}.csv"
        with out_csv.open("w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(combined_tags)
            writer.writerows(combined_rows)

        print(f"Wrote {out_csv}")


if __name__ == "__main__":
    main()
