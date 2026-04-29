#!/usr/bin/env python3
"""
Run Perple_X fluids.exe to generate reference data and auto-fill regression tests.

This script:
1. Runs fluids.exe for specific test cases
2. Parses the output
3. Generates reference_data.json
4. Auto-updates test_regression.cpp with reference values
5. Creates a C++ header with all reference data

Usage:
    python run_perplex_and_fill_tests.py --fluids-exe fluids.exe
"""

import subprocess
import json
import argparse
import sys
import re
from pathlib import Path
from typing import List, Dict, Tuple, Optional


def run_fluids_exe(
    fluids_path: str, ifug: int, T: float, P: float, composition: List[float]
) -> Dict:
    """
    Run Perple_X fluids.exe for a single point calculation.

    Args:
        fluids_path: Path to fluids.exe
        ifug: EoS option (0=MRK, 1=HSMRK, 2=Hybrid, etc.)
        T: Temperature (K)
        P: Pressure (bar)
        composition: Mole fractions [H2O, CO2, ...]

    Returns:
        Dict with parsed results
    """

    # Build input for fluids.exe
    # Format for single point calculation:
    # ifug
    # n (no tabulation)
    # y (log output format)
    # P T
    # y1 y2 y3 ...
    # n (no more calculations)

    comp_str = " ".join(f"{c:.10f}" for c in composition)

    input_text = f"""{ifug}
n
y
{P:.6f} {T:.6f}
{comp_str}
n
"""

    try:
        result = subprocess.run(
            [fluids_path], input=input_text, capture_output=True, text=True, timeout=30
        )

        if result.returncode != 0:
            print(f"fluids.exe failed with return code {result.returncode}")
            print("STDOUT:", result.stdout)
            print("STDERR:", result.stderr)
            return None

        # Parse output
        output = result.stdout

        parsed = {
            "T": T,
            "P": P,
            "composition": composition,
            "ifug": ifug,
        }

        # Parse fugacities: look for "ln(fO2)" pattern
        for line in output.split("\n"):
            # Fugacity coefficients: g(H2O), g(CO2), etc.
            if "g(H2O)" in line:
                match = re.search(
                    r"g\(H2O\)\s*=\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)", line
                )
                if match:
                    parsed["ln_f_H2O"] = float(match.group(1))

            if "g(CO2)" in line:
                match = re.search(
                    r"g\(CO2\)\s*=\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)", line
                )
                if match:
                    parsed["ln_f_CO2"] = float(match.group(1))

            # Volume: vol = X.XX cm3/mol
            if "vol" in line and "cm3" in line:
                match = re.search(r"vol\s*=\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)", line)
                if match:
                    parsed["vol"] = float(match.group(1))

            # Mole fractions: y(H2O) = X.XX
            if "y(H2O)" in line:
                match = re.search(
                    r"y\(H2O\)\s*=\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)", line
                )
                if match:
                    parsed["y_H2O"] = float(match.group(1))

            if "y(CO2)" in line:
                match = re.search(
                    r"y\(CO2\)\s*=\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)", line
                )
                if match:
                    parsed["y_CO2"] = float(match.group(1))

        return parsed

    except subprocess.TimeoutExpired:
        print(f"fluids.exe timed out for T={T}, P={P}")
        return None
    except Exception as e:
        print(f"Error running fluids.exe: {e}")
        return None


def generate_test_data(fluids_path: str) -> List[Dict]:
    """
    Generate comprehensive test data by running fluids.exe.

    Returns:
        List of test case results
    """

    test_cases = []

    # Test 1: Pure H2O with HSMRK (ifug=1)
    print("Generating Test 1: Pure H2O HSMRK...")
    for T, P in [(500.0, 1000.0), (600.0, 2000.0), (700.0, 5000.0)]:
        result = run_fluids_exe(fluids_path, ifug=1, T=T, P=P, composition=[1.0])
        if result:
            result["test_name"] = "pure_h2o_hsmrk"
            test_cases.append(result)

    # Test 2: Pure CO2 with CORK (ifug=5)
    print("Generating Test 2: Pure CO2 CORK...")
    for T, P in [(500.0, 1000.0), (600.0, 2000.0)]:
        result = run_fluids_exe(fluids_path, ifug=5, T=T, P=P, composition=[0.0, 1.0])
        if result:
            result["test_name"] = "pure_co2_cork"
            test_cases.append(result)

    # Test 3: H2O-CO2 binary with Hybrid (ifug=2)
    print("Generating Test 3: H2O-CO2 Binary...")
    for xco2 in [0.1, 0.25, 0.5, 0.75, 0.9]:
        result = run_fluids_exe(
            fluids_path, ifug=2, T=500.0, P=1000.0, composition=[1.0 - xco2, xco2]
        )
        if result:
            result["test_name"] = f"h2o_co2_binary_xco2_{int(xco2 * 100)}"
            test_cases.append(result)

    # Test 4: P-T Grid (H2O-CO2 70-30)
    print("Generating Test 4: P-T Grid...")
    pressures = [1.0, 100.0, 500.0, 1000.0, 2000.0, 5000.0]
    temperatures = [298.15, 373.15, 473.15, 573.15, 673.15]

    for P in pressures:
        for T in temperatures:
            result = run_fluids_exe(
                fluids_path, ifug=2, T=T, P=P, composition=[0.7, 0.3]
            )
            if result:
                result["test_name"] = "pt_grid"
                test_cases.append(result)

    return test_cases


def save_reference_json(test_cases: List[Dict], output_path: str):
    """Save test cases to JSON file."""
    with open(output_path, "w") as f:
        json.dump(test_cases, f, indent=2)
    print(f"\nSaved {len(test_cases)} reference points to {output_path}")


def generate_cpp_header(test_cases: List[Dict], output_path: str):
    """
    Generate a C++ header file with all reference data.
    """

    with open(output_path, "w") as f:
        f.write("// Auto-generated reference data from Perple_X\n")
        f.write("// DO NOT EDIT MANUALLY\n\n")
        f.write("#pragma once\n\n")
        f.write("#include <vector>\n")
        f.write("#include <string>\n")
        f.write("#include <map>\n\n")
        f.write("namespace PerpleXReferenceData {\n\n")

        # Group by test name
        tests = {}
        for case in test_cases:
            name = case.get("test_name", "unknown")
            if name not in tests:
                tests[name] = []
            tests[name].append(case)

        # Generate struct for each test
        for test_name, cases in tests.items():
            struct_name = "".join(word.capitalize() for word in test_name.split("_"))

            f.write(f"// {test_name}\n")
            f.write(f"struct {struct_name} {{\n")
            f.write(f"    double T;\n")
            f.write(f"    double P;\n")
            f.write(f"    std::vector<double> composition;\n")
            f.write(f"    double ln_f_H2O;\n")
            f.write(f"    double ln_f_CO2;\n")
            f.write(f"    double vol;\n")
            f.write(f"}};\n\n")

            f.write(f"const std::vector<{struct_name}> {test_name}_data = {{\n")
            for case in cases:
                f.write(f"    {{ // T={case['T']}K P={case['P']}bar\n")
                f.write(f"        {case['T']},  // T\n")
                f.write(f"        {case['P']},  // P\n")
                comp_str = ", ".join(str(c) for c in case["composition"])
                f.write(f"        {{{comp_str}}},  // composition\n")
                f.write(f"        {case.get('ln_f_H2O', 0.0)},  // ln_f_H2O\n")
                f.write(f"        {case.get('ln_f_CO2', 0.0)},  // ln_f_CO2\n")
                f.write(f"        {case.get('vol', 0.0)}   // vol\n")
                f.write(f"    }},\n")
            f.write(f"}};\n\n")

        f.write("}  // namespace PerpleXReferenceData\n")

    print(f"Generated C++ header: {output_path}")


def update_test_regression_cpp(test_cases: List[Dict], cpp_path: str):
    """
    Update test_regression.cpp with actual reference values.
    """

    # Read current file
    with open(cpp_path, "r") as f:
        content = f.read()

    # Find pure H2O HSMRK test and update
    h2o_data = [c for c in test_cases if c.get("test_name") == "pure_h2o_hsmrk"]
    if h2o_data:
        case = h2o_data[0]  # First data point

        # Replace placeholder values
        pattern = r"ref\.ln_f = \{0\.0\};  // Placeholder"
        replacement = f"ref.ln_f = {{{case.get('ln_f_H2O', 0.0)}}};  // From Perple_X"
        content = re.sub(pattern, replacement, content)

        pattern = r"ref\.v = \{18\.0\};    // Placeholder"
        replacement = f"ref.v = {{{case.get('vol', 18.0)}}};    // From Perple_X"
        content = re.sub(pattern, replacement, content)

        pattern = r"ref\.vol = 18\.0;    // Placeholder"
        replacement = f"ref.vol = {case.get('vol', 18.0)};    // From Perple_X"
        content = re.sub(pattern, replacement, content)

    # Update CO2 test
    co2_data = [c for c in test_cases if c.get("test_name") == "pure_co2_cork"]
    if co2_data:
        case = co2_data[0]

        pattern = r"// TODO: Reference data from Perple_X\s+ref\.ln_f = \{0\.0\};"
        replacement = f"// Reference data from Perple_X\n    ref.ln_f = {{{case.get('ln_f_CO2', 0.0)}}};"
        content = re.sub(pattern, replacement, content)

    # Write updated file
    with open(cpp_path, "w") as f:
        f.write(content)

    print(f"Updated {cpp_path} with reference values")


def main():
    parser = argparse.ArgumentParser(
        description="Generate Perple_X reference data and fill regression tests"
    )
    parser.add_argument(
        "--fluids-exe", default="fluids.exe", help="Path to Perple_X fluids executable"
    )
    parser.add_argument(
        "--output-json", default="reference_data.json", help="Output JSON file"
    )
    parser.add_argument(
        "--output-header",
        default="perplex_reference_data.hpp",
        help="Output C++ header file",
    )
    parser.add_argument(
        "--test-cpp",
        default="test_regression.cpp",
        help="Path to test_regression.cpp to update",
    )
    parser.add_argument(
        "--quick", action="store_true", help="Run quick test (fewer data points)"
    )

    args = parser.parse_args()

    # Check if fluids.exe exists
    try:
        result = subprocess.run(
            [args.fluids_exe], input="n\n", capture_output=True, text=True, timeout=5
        )
    except FileNotFoundError:
        print(f"Error: fluids executable not found: {args.fluids_exe}")
        print("\nPlease specify the path to fluids.exe:")
        print("  --fluids-exe /path/to/fluids.exe")
        print("\nOn Windows, it might be:")
        print("  --fluids-exe 'C:\\Program Files (x86)\\Perplex\\fluids.exe'")
        sys.exit(1)
    except Exception as e:
        print(f"Warning: Could not verify fluids.exe: {e}")

    print("=" * 70)
    print("Perple_X Reference Data Generation")
    print("=" * 70)
    print(f"Fluids executable: {args.fluids_exe}")
    print(f"Output JSON: {args.output_json}")
    print(f"Output header: {args.output_header}")
    print()

    # Generate test data
    test_cases = generate_test_data(args.fluids_exe)

    if not test_cases:
        print("\nError: No test cases generated!")
        sys.exit(1)

    # Save to JSON
    save_reference_json(test_cases, args.output_json)

    # Generate C++ header
    generate_cpp_header(test_cases, args.output_header)

    # Update test file if it exists
    if Path(args.test_cpp).exists():
        update_test_regression_cpp(test_cases, args.test_cpp)
    else:
        print(f"\nNote: {args.test_cpp} not found, skipping update")

    print("\n" + "=" * 70)
    print("✓ Reference data generation complete!")
    print("=" * 70)
    print("\nNext steps:")
    print(f"1. Review {args.output_json}")
    print(f"2. Include {args.output_header} in test_regression.cpp")
    print(
        "3. Compile: g++ -std=c++17 test_regression.cpp PerpleX*.cpp -o test_regression"
    )
    print("4. Run: ./test_regression")
    print()


if __name__ == "__main__":
    main()
