#!/usr/bin/env python3
"""
Analyze G0 (Gibbs energy) differences for all phases
in the H2O + CO2(aq) = H+ + HCO3- reaction
Compares Excel truth vs Reaktoro DEW model predictions
"""

import subprocess
import re
from pathlib import Path
import pandas as pd


def run_test_and_extract_species_errors():
    """Run test and extract per-species G0 differences"""

    exe = Path("build-msvc/Reaktoro/Release/reaktoro-cpptests.exe")

    print("Running test suite to capture per-species G0 values...")
    print("=" * 80)

    result = subprocess.run(
        [str(exe), "[dew][reaction]"],
        capture_output=True,
        text=True,
        timeout=180,
        cwd=".",
    )

    output = result.stdout + result.stderr

    # Parse DEBUG per-species lines
    species_pattern = r"DEBUG per-species at T=([\d.]+) C, P=([\d.]+) kb\n(.*?)\nTruth"

    errors = []

    for match in re.finditer(species_pattern, output, re.DOTALL):
        T_C = float(match.group(1))
        P_kb = float(match.group(2))
        block = match.group(3)

        # Extract individual species data
        h2o_match = re.search(r"H2O: G0=([-\d.]+)", block)
        co2_match = re.search(r"CO2: G0=([-\d.]+)", block)
        h_match = re.search(r"H\+: G0=([-\d.]+)", block)
        hco3_match = re.search(r"HCO3-: G0=([-\d.]+)", block)

        # Extract truth values (after "Truth" line)
        truth_match = re.search(
            r"Truth H2O G0=([-\d.]+) J/mol.*?CO2 G0=([-\d.]+).*?H\+ G0=([-\d.]+).*?HCO3 G0=([-\d.]+)",
            output[match.start() : match.end() + 200],
            re.DOTALL,
        )

        if h2o_match and truth_match:
            g0_h2o_model = float(h2o_match.group(1))
            g0_h2o_truth = float(truth_match.group(1))

            g0_co2_model = float(co2_match.group(1))
            g0_co2_truth = float(truth_match.group(2))

            g0_h_model = float(h_match.group(1))
            g0_h_truth = float(truth_match.group(3))

            g0_hco3_model = float(hco3_match.group(1))
            g0_hco3_truth = float(truth_match.group(4))

            errors.append(
                {
                    "T_C": T_C,
                    "P_kb": P_kb,
                    "H2O_model": g0_h2o_model,
                    "H2O_truth": g0_h2o_truth,
                    "H2O_err": g0_h2o_model - g0_h2o_truth,
                    "CO2_model": g0_co2_model,
                    "CO2_truth": g0_co2_truth,
                    "CO2_err": g0_co2_model - g0_co2_truth,
                    "H_model": g0_h_model,
                    "H_truth": g0_h_truth,
                    "H_err": g0_h_model - g0_h_truth,
                    "HCO3_model": g0_hco3_model,
                    "HCO3_truth": g0_hco3_truth,
                    "HCO3_err": g0_hco3_model - g0_hco3_truth,
                }
            )

    return errors, output


def analyze_species_errors(errors):
    """Analyze per-species error distribution"""

    if not errors:
        print("[!] No per-species data found in test output")
        print("    The test suite may not be printing DEBUG data")
        print(
            "    Run: .\build-msvc\Reaktoro\Release\reaktoro-cpptests.exe '[dew][reaction]' 2>&1 | head -100"
        )
        return

    df = pd.DataFrame(errors)

    print("\n" + "=" * 80)
    print("PER-SPECIES G0 ERROR ANALYSIS")
    print("Reaction: H2O + CO2(aq) = H+ + HCO3-")
    print("=" * 80)

    print(f"\nTotal test points: {len(df)}\n")

    species = ["H2O", "CO2", "H", "HCO3"]

    for species_name, col_prefix in [
        ("H2O", "H2O"),
        ("CO2", "CO2"),
        ("H+", "H"),
        ("HCO3-", "HCO3"),
    ]:
        err_col = f"{col_prefix}_err"

        print(f"\n{species_name}:")
        print("-" * 80)

        abs_errs = df[err_col].abs()

        print(f"  Error statistics (J/mol):")
        print(f"    Min:  {abs_errs.min():>10.2f}")
        print(f"    Max:  {abs_errs.max():>10.2f}")
        print(f"    Mean: {abs_errs.mean():>10.2f}")
        print(f"    Std:  {abs_errs.std():>10.2f}")

        # Find worst conditions
        worst_idx = abs_errs.idxmax()
        worst = df.loc[worst_idx]

        print(f"\n  Worst case:")
        print(f"    Condition: T={worst['T_C']:.0f}°C, P={worst['P_kb']:.0f} kb")
        print(f"    Model: {worst[f'{col_prefix}_model']:>12.2f} J/mol")
        print(f"    Truth: {worst[f'{col_prefix}_truth']:>12.2f} J/mol")
        print(
            f"    Error: {worst[err_col]:>12.2f} J/mol ({100 * worst[err_col] / worst[f'{col_prefix}_truth']:.3f}%)"
        )

        # Top 3 errors
        top_3_idx = abs_errs.nlargest(3).index

        print(f"\n  Top 3 error conditions:")
        for rank, idx in enumerate(top_3_idx, 1):
            row = df.loc[idx]
            print(
                f"    {rank}. T={row['T_C']:.0f}°C, P={row['P_kb']:.0f} kb: "
                f"error={row[err_col]:>8.2f} J/mol"
            )

    print("\n" + "=" * 80)
    print("SUMMARY BY PHASE")
    print("=" * 80)

    phase_errors = {
        "H2O": df["H2O_err"].abs().mean(),
        "CO2": df["CO2_err"].abs().mean(),
        "H+": df["H_err"].abs().mean(),
        "HCO3-": df["HCO3_err"].abs().mean(),
    }

    print("\nAverage absolute error by phase:\n")
    for phase in sorted(
        phase_errors.keys(), key=lambda x: phase_errors[x], reverse=True
    ):
        print(f"  {phase:<10}: {phase_errors[phase]:>8.2f} J/mol")

    # Which phase contributes most to reaction error?
    print("\n" + "=" * 80)
    print("IMPACT ON REACTION GIBBS ENERGY")
    print("=" * 80)
    print("\nΔGr = G(H+) + G(HCO3-) - G(H2O) - G(CO2)")
    print("\nPhase contributions to reaction error:\n")

    for phase, sign in [("H+", "+"), ("HCO3-", "+"), ("H2O", "-"), ("CO2", "-")]:
        col_prefix = "H" if phase == "H+" else phase.replace("-", "")
        err_col = f"{col_prefix}_err"
        avg_contrib = (df[err_col].mean()) * (1 if sign == "+" else -1)
        print(f"  {phase:<10} ({sign}): {avg_contrib:>8.2f} J/mol average contribution")

    total_avg = (
        df["H_err"].mean()
        + df["HCO3_err"].mean()
        - df["H2O_err"].mean()
        - df["CO2_err"].mean()
    )
    print(
        f"\n  Total reaction error: {total_avg:>8.2f} J/mol (should match ~15.36 from test)"
    )


def main():
    print("\n" + "=" * 80)
    print("DETAILED PHASE-BY-PHASE ERROR ANALYSIS")
    print("Reaktoro DEW vs Excel Truth")
    print("=" * 80 + "\n")

    errors, output = run_test_and_extract_species_errors()

    if errors:
        analyze_species_errors(errors)
    else:
        print("\n[!] Could not extract per-species G0 data")
        print("    Looking for test output summary instead...\n")

        # Try to find the test output summary
        lines = output.split("\n")
        for i, line in enumerate(lines):
            if (
                "Tested 180" in line
                or "passed" in line.lower()
                or "error:" in line.lower()
            ):
                print(line)


if __name__ == "__main__":
    main()
