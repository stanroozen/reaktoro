#!/usr/bin/env python3
"""Structured regression suite runner for DEW/PerplexDEW benchmark scripts.

This file centralizes benchmark regression execution in a single entrypoint.
It is intended to be run after a successful local build in build/Reaktoro/Release.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Iterable


BENCHMARK_DIR = Path(__file__).resolve().parent
SUITE_RESULTS_DIR = BENCHMARK_DIR / "regression_results"


@dataclass(frozen=True)
class RegressionCase:
    name: str
    relative_script: str
    tags: tuple[str, ...]
    timeout_seconds: int = 1800

    @property
    def script_path(self) -> Path:
        return BENCHMARK_DIR / self.relative_script


CASES: tuple[RegressionCase, ...] = (
    # Core DEW smoke tests
    RegressionCase(
        "test_build_system",
        "regression/smoke/test_build_system.py",
        ("smoke", "core"),
        300,
    ),
    RegressionCase(
        "test_dew_water_models",
        "regression/smoke/test_dew_water_models.py",
        ("smoke", "water"),
        300,
    ),
    RegressionCase(
        "test_h2o_aq_dew",
        "regression/smoke/test_h2o_aq_dew.py",
        ("smoke", "water"),
        300,
    ),
    RegressionCase(
        "test_h2o_dummy", "regression/smoke/test_h2o_dummy.py", ("smoke", "water"), 300
    ),
    RegressionCase(
        "test_perplex_conditions_nosilence",
        "regression/smoke/test_perplex_conditions_nosilence.py",
        ("smoke", "perplexdew"),
        300,
    ),
    # Quartz backend runs
    RegressionCase(
        "run_quartz_dew",
        "regression/quartz/run_quartz_DEW.py",
        ("regression", "quartz", "dew"),
        900,
    ),
    RegressionCase(
        "run_quartz_perplexdew",
        "regression/quartz/run_quartz_PerplexDEW.py",
        ("regression", "quartz", "perplexdew"),
        1200,
    ),
    RegressionCase(
        "run_quartz_perplexdew_extendeddh",
        "regression/quartz/run_quartz_PerplexDEW_ExtendedDH.py",
        ("regression", "quartz", "perplexdew"),
        1200,
    ),
    RegressionCase(
        "run_quartz_perplexdew_davies_h2oco2",
        "regression/quartz/run_quartz_PerplexDEW_Davies_H2OCO2.py",
        ("regression", "quartz", "perplexdew", "h2oco2"),
        1200,
    ),
    RegressionCase(
        "run_quartz_perplexdew_extendeddh_h2oco2",
        "regression/quartz/run_quartz_PerplexDEW_ExtendedDH_H2OCO2.py",
        ("regression", "quartz", "perplexdew", "h2oco2"),
        1200,
    ),
    # Calcite backend runs
    RegressionCase(
        "run_calcite_dew",
        "regression/calcite/run_calcite_DEW.py",
        ("regression", "calcite", "dew"),
        1200,
    ),
    RegressionCase(
        "run_calcite_dew_co2aq",
        "regression/calcite/run_calcite_DEW_CO2aq.py",
        ("regression", "calcite", "dew", "h2oco2"),
        1200,
    ),
    RegressionCase(
        "run_calcite_perplexdew_davies",
        "regression/calcite/run_calcite_PerplexDEW_Davies.py",
        ("regression", "calcite", "perplexdew"),
        1200,
    ),
    RegressionCase(
        "run_calcite_perplexdew_extendeddh",
        "regression/calcite/run_calcite_PerplexDEW_ExtendedDH.py",
        ("regression", "calcite", "perplexdew"),
        1200,
    ),
    RegressionCase(
        "run_calcite_perplexdew_davies_h2oco2",
        "regression/calcite/run_calcite_PerplexDEW_Davies_H2OCO2.py",
        ("regression", "calcite", "perplexdew", "h2oco2"),
        1200,
    ),
    RegressionCase(
        "run_calcite_perplexdew_extendeddh_h2oco2",
        "regression/calcite/run_calcite_PerplexDEW_ExtendedDH_H2OCO2.py",
        ("regression", "calcite", "perplexdew", "h2oco2"),
        1200,
    ),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run DEW/PerplexDEW regression suite.")
    parser.add_argument(
        "--list", action="store_true", help="List available cases and exit."
    )
    parser.add_argument(
        "--tags",
        nargs="*",
        default=[],
        help="Only run cases that include all provided tags (e.g. smoke perplexdew).",
    )
    parser.add_argument(
        "--python-exe",
        default=sys.executable,
        help="Python executable used to run each regression case.",
    )
    parser.add_argument(
        "--continue-on-fail",
        action="store_true",
        help="Continue executing remaining cases after a failure.",
    )
    parser.add_argument(
        "--max-failures",
        type=int,
        default=0,
        help="Stop after this many failures (0 means no limit).",
    )
    parser.add_argument(
        "--max-cases",
        type=int,
        default=0,
        help="Only execute the first N selected cases (0 means all).",
    )
    parser.add_argument(
        "--reaktoro-lib-dir",
        default="",
        help="Path to the Reaktoro build lib dir (e.g. build/Reaktoro/Release). "
        "Added to PATH and PYTHONPATH so scripts can find reaktoro4py.pyd and Reaktoro.dll.",
    )
    return parser.parse_args()


def select_cases(required_tags: Iterable[str]) -> list[RegressionCase]:
    tags = tuple(t.strip().lower() for t in required_tags if t.strip())
    selected: list[RegressionCase] = []
    for case in CASES:
        case_tags = tuple(t.lower() for t in case.tags)
        if all(tag in case_tags for tag in tags):
            selected.append(case)
    return selected


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if text is None:
        text = ""
    path.write_text(text, encoding="utf-8")


def main() -> int:
    args = parse_args()
    selected = select_cases(args.tags)

    if args.max_cases and args.max_cases > 0:
        selected = selected[: args.max_cases]

    if args.list:
        for case in selected:
            print(
                f"{case.name:45} tags={','.join(case.tags)} script={case.relative_script}"
            )
        return 0

    if not selected:
        print("No regression cases selected. Use --list to inspect available cases.")
        return 2

    run_stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_dir = SUITE_RESULTS_DIR / run_stamp
    run_dir.mkdir(parents=True, exist_ok=True)

    failures = 0
    results: list[dict] = []

    print(f"Running {len(selected)} regression cases...")
    print(f"Results directory: {run_dir}")

    for index, case in enumerate(selected, start=1):
        script_path = case.script_path
        if not script_path.is_file():
            failures += 1
            result = {
                "name": case.name,
                "status": "missing",
                "exit_code": None,
                "duration_seconds": 0.0,
                "script": str(script_path),
                "tags": list(case.tags),
                "stdout_log": None,
                "stderr_log": None,
            }
            results.append(result)
            print(f"[{index}/{len(selected)}] MISSING {case.name}: {script_path}")
            if not args.continue_on_fail:
                break
            if args.max_failures and failures >= args.max_failures:
                break
            continue

        print(f"[{index}/{len(selected)}] RUN {case.name}")
        t0 = time.time()

        run_env = os.environ.copy()
        run_env["PYTHONUTF8"] = "1"
        if args.reaktoro_lib_dir:
            lib_dir = str(Path(args.reaktoro_lib_dir).resolve())
            run_env["PATH"] = lib_dir + os.pathsep + run_env.get("PATH", "")
            run_env["PYTHONPATH"] = lib_dir + os.pathsep + run_env.get("PYTHONPATH", "")

        proc = subprocess.run(
            [args.python_exe, script_path.name],
            cwd=script_path.parent,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=case.timeout_seconds,
            env=run_env,
        )

        duration = time.time() - t0
        status = "passed" if proc.returncode == 0 else "failed"

        stdout_log = run_dir / f"{case.name}.stdout.log"
        stderr_log = run_dir / f"{case.name}.stderr.log"
        write_text(stdout_log, proc.stdout)
        write_text(stderr_log, proc.stderr)

        if status != "passed":
            failures += 1

        result = {
            "name": case.name,
            "status": status,
            "exit_code": proc.returncode,
            "duration_seconds": round(duration, 3),
            "script": str(script_path),
            "tags": list(case.tags),
            "stdout_log": str(stdout_log),
            "stderr_log": str(stderr_log),
        }
        results.append(result)

        print(f"    -> {status.upper()} (exit={proc.returncode}, {duration:.1f}s)")

        if status != "passed" and not args.continue_on_fail:
            print("Stopping on first failure (use --continue-on-fail to override).")
            break
        if args.max_failures and failures >= args.max_failures:
            print(f"Stopping after reaching --max-failures={args.max_failures}.")
            break

    summary = {
        "timestamp": run_stamp,
        "python_executable": args.python_exe,
        "selected_tags": list(args.tags),
        "max_cases": args.max_cases,
        "total": len(results),
        "passed": sum(1 for r in results if r["status"] == "passed"),
        "failed": sum(1 for r in results if r["status"] == "failed"),
        "missing": sum(1 for r in results if r["status"] == "missing"),
        "results": results,
    }

    summary_path = run_dir / "summary.json"
    write_text(summary_path, json.dumps(summary, indent=2))

    print("\nRegression summary")
    print("------------------")
    print(f"Total  : {summary['total']}")
    print(f"Passed : {summary['passed']}")
    print(f"Failed : {summary['failed']}")
    print(f"Missing: {summary['missing']}")
    print(f"JSON   : {summary_path}")

    return 1 if summary["failed"] or summary["missing"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
