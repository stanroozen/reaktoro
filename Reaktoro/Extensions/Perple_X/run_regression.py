#!/usr/bin/env python3
"""
Unified regression runner for Perple_X parity testing.

Workflow:
1. Regenerate GFSM reference CSVs from Perple_X meemum executable.
2. Run Reaktoro regression executable (test_regression.exe) against those references.

Usage examples:
    python run_regression.py --meemum-exe "C:\\Program Files (x86)\\Perplex\\meemum.exe"
    python run_regression.py --meemum-exe "...\\meemum.exe" --skip-generate
    python run_regression.py --meemum-exe "...\\meemum.exe" --build-cmd "cmake --build build --config Release --target test_regression"
"""

from __future__ import annotations

import argparse
import shlex
import subprocess
import sys
from pathlib import Path


REQUIRED_GFSM_MATRIX_CASES = {
    # Pure H2O — all EOS variants (iopt(25))
    "h2o_mrk",
    "h2o_hsmrk",
    "h2o_cork",
    "h2o_pseos",
    "h2o_haar",
    "h2o_zd05",
    "h2o_zd09",
    # Pure CO2 — all non-broken EOS variants (iopt(26); brmrk=3 excluded)
    "co2_mrk",
    "co2_hsmrk",
    "co2_cork",
    "co2_pseos",
    "co2_zd09",
    # Pure CH4 — all EOS variants (iopt(27))
    "ch4_mrk",
    "ch4_hsmrk",
    "ch4_zd09",
}

REQUIRED_GFSM_REFERENCE_CSVS = {
    # Pure-species anchors (same as REQUIRED_GFSM_MATRIX_CASES)
    *REQUIRED_GFSM_MATRIX_CASES,
    # Binary H2O+CO2 — vary H2O EOS (CO2=MRK)
    "h2o_co2_mrk_mrk",
    "h2o_co2_hsmrk_mrk",
    "h2o_co2_cork_mrk",
    "h2o_co2_pseos_mrk",
    "h2o_co2_haar_mrk",
    "h2o_co2_zd05_mrk",
    "h2o_co2_zd09_mrk",
    # Binary H2O+CO2 — vary CO2 EOS (H2O=MRK, skip mrk+mrk duplicate)
    "h2o_co2_mrk_hsmrk",
    "h2o_co2_mrk_cork",
    "h2o_co2_mrk_pseos",
    "h2o_co2_mrk_zd09",
    # Binary H2O+CH4 — vary H2O EOS (CH4=MRK)
    "h2o_ch4_mrk_mrk",
    "h2o_ch4_hsmrk_mrk",
    "h2o_ch4_cork_mrk",
    "h2o_ch4_pseos_mrk",
    "h2o_ch4_haar_mrk",
    "h2o_ch4_zd05_mrk",
    "h2o_ch4_zd09_mrk",
    # Binary H2O+CH4 — vary CH4 EOS (H2O=MRK, skip mrk duplicate)
    "h2o_ch4_mrk_hsmrk",
    "h2o_ch4_mrk_zd09",
    # Binary CO2+CH4 — vary CO2 EOS (CH4=MRK)
    "co2_ch4_mrk_mrk",
    "co2_ch4_hsmrk_mrk",
    "co2_ch4_cork_mrk",
    "co2_ch4_pseos_mrk",
    "co2_ch4_zd09_mrk",
    # Binary CO2+CH4 — vary CH4 EOS (CO2=MRK, skip mrk duplicate)
    "co2_ch4_mrk_hsmrk",
    "co2_ch4_mrk_zd09",
    # Ternary H2O+CO2+CH4 — equal thirds, matched EOS across all three species
    "h2o_co2_ch4_mrk",
    "h2o_co2_ch4_hsmrk",
    "h2o_co2_ch4_zd09",
    # O2-bearing and redox-stress OHC states
    "redox_o2_excess_mrk",
    "redox_o2_excess_hot",
    "redox_h2_excess_hot",
    "redox_co_bias_hot",
    "redox_o2_excess_hybrid",
    "redox_mixed_hybrid",
}


def emit_console_text(text: str) -> None:
    """Write text to stdout without failing on Windows codepage encoding issues."""
    if not text:
        return

    try:
        sys.stdout.write(text)
    except UnicodeEncodeError:
        encoding = sys.stdout.encoding or "utf-8"
        payload = text.encode(encoding, errors="backslashreplace")
        if hasattr(sys.stdout, "buffer"):
            sys.stdout.buffer.write(payload)
        else:
            # Fallback for environments without a writable buffer.
            sys.stdout.write(payload.decode(encoding, errors="ignore"))
    sys.stdout.flush()

def run_cmd(cmd: list[str], cwd: Path, label: str) -> None:
    print(f"\n[{label}] {' '.join(cmd)}")
    result = subprocess.run(
        cmd,
        cwd=str(cwd),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        errors="replace",
    )
    emit_console_text(result.stdout)
    if result.returncode != 0:
        raise RuntimeError(f"Step '{label}' failed with exit code {result.returncode}")


def resolve_test_exe(base_dir: Path, explicit: str | None) -> Path:
    if explicit:
        p = Path(explicit)
        if not p.is_absolute():
            p = base_dir / p
        return p

    candidates = [
        base_dir / "test_regression.exe",
        base_dir / "test_regression",
    ]
    for c in candidates:
        if c.exists():
            return c

    return candidates[0]


def main() -> int:
    base_dir = Path(__file__).resolve().parent

    parser = argparse.ArgumentParser(
        description="Run Perple_X -> Reaktoro parity regression"
    )
    parser.add_argument(
        "--meemum-exe",
        required=True,
        help="Path to Perple_X meemum executable",
    )
    parser.add_argument(
        "--template-dat",
        default=r"C:\Users\stanroozen\Documents\Projects\Perplex\Perple_X\test\weigang\gfsm_fluid_probe.dat",
        help="Path to GFSM meemum template .dat file",
    )
    parser.add_argument(
        "--case-matrix",
        default=str(base_dir / "test" / "gfsm_case_matrix.csv"),
        help="Optional CSV with additional GFSM meemum cases",
    )
    parser.add_argument(
        "--out-dir",
        default=str(base_dir / "test" / "gfsm"),
        help="Output directory for generated GFSM reference CSVs",
    )
    parser.add_argument(
        "--require-full-gfsm-pure",
        action="store_true",
        help="Require all GFSM-callable pure EoS case names and generated CSVs to be present.",
    )
    parser.add_argument(
        "--test-exe",
        default=None,
        help="Path to compiled test_regression executable (default: auto-detect in extension folder)",
    )
    parser.add_argument(
        "--build-cmd",
        default=None,
        help="Optional build command to compile test_regression before running",
    )
    parser.add_argument(
        "--skip-generate",
        action="store_true",
        help="Skip reference regeneration step",
    )
    parser.add_argument(
        "--skip-run",
        action="store_true",
        help="Skip running regression executable",
    )

    args = parser.parse_args()

    meemum_exe = Path(args.meemum_exe)
    if not meemum_exe.exists():
        raise FileNotFoundError(f"meemum executable not found: {meemum_exe}")

    template_dat = Path(args.template_dat)
    if not template_dat.exists() and not args.skip_generate:
        raise FileNotFoundError(f"template dat file not found: {template_dat}")

    case_matrix = Path(args.case_matrix)
    if not case_matrix.exists() and not args.skip_generate:
        print(f"[WARN] case matrix not found, generator will use built-in redox defaults: {case_matrix}")

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    try:
        if not args.skip_generate:
            run_cmd(
                [
                    sys.executable,
                    str(base_dir / "generate_gfsm_meemum_references.py"),
                    "--meemum-exe",
                    str(meemum_exe),
                    "--template-dat",
                    str(template_dat),
                    "--case-matrix",
                    str(case_matrix),
                    "--out-dir",
                    str(out_dir),
                ],
                cwd=base_dir,
                label="Generate GFSM references via Perple_X meemum",
            )

        if args.build_cmd:
            run_cmd(
                shlex.split(args.build_cmd),
                cwd=base_dir,
                label="Build regression executable",
            )

        if not args.skip_run:
            if args.require_full_gfsm_pure:
                gfsm_dir = out_dir
                missing_csv = []
                for case in sorted(REQUIRED_GFSM_REFERENCE_CSVS):
                    p = gfsm_dir / f"{case}.csv"
                    if not p.exists() or p.stat().st_size == 0:
                        missing_csv.append(str(p))
                if missing_csv:
                    raise FileNotFoundError(
                        "Missing/empty required GFSM reference CSVs: "
                        + "; ".join(missing_csv)
                    )

            test_exe = resolve_test_exe(base_dir, args.test_exe)
            if not test_exe.exists():
                raise FileNotFoundError(
                    "test_regression executable not found. "
                    "Provide --build-cmd to compile it or --test-exe to point to an existing binary."
                )

            run_cmd(
                [str(test_exe)], cwd=base_dir, label="Run Reaktoro regression tests"
            )

        print("\n[DONE] Regression workflow completed successfully.")
        return 0

    except Exception as exc:
        print(f"\n[FAILED] {exc}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
