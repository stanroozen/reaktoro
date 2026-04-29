"""Runs the two H2O-CO2 calcite benchmark scripts and writes a combined log."""

import subprocess, sys, os, pathlib

py = sys.executable
base = pathlib.Path(__file__).parent
log = pathlib.Path(r"C:\Temp\calcite_h2oco2_runs.log")
log.parent.mkdir(parents=True, exist_ok=True)

scripts = [
    ("PerplexDEW_Davies_H2OCO2", base / "run_calcite_PerplexDEW_Davies_H2OCO2.py"),
    (
        "PerplexDEW_ExtendedDH_H2OCO2",
        base / "run_calcite_PerplexDEW_ExtendedDH_H2OCO2.py",
    ),
]

with log.open("w", encoding="utf-8") as flog:
    for tag, script in scripts:
        header = f"\n{'=' * 60}\n[{tag}]\n{'=' * 60}\n"
        print(header, end="")
        flog.write(header)
        flog.flush()

        proc = subprocess.run(
            [py, "-X", "utf8", str(script)],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            encoding="utf-8",
            errors="replace",
        )
        print(proc.stdout, end="")
        flog.write(proc.stdout)
        flog.write(f"\n[exit code: {proc.returncode}]\n")
        flog.flush()

print("ALL DONE. Log at", log)
