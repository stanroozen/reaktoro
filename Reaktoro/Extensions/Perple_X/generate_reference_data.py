#!/usr/bin/env python3
"""
Generate Reference Data from Perple_X for Regression Testing

This script runs Perple_X COHSRK/fluids executable to generate reference
data for validating the Reaktoro implementation.

Usage:
    python generate_reference_data.py --perplex-path /path/to/Perple_X
"""

import subprocess
import json
import os
import sys
from dataclasses import dataclass, asdict
from typing import List, Dict, Any
import argparse


@dataclass
class FluidState:
    """Fluid state computed by Perple_X"""

    T: float  # Temperature (K)
    P: float  # Pressure (bar)
    composition: List[float]  # Mole fractions
    species: List[int]  # Species indices (1=H2O, 2=CO2, etc.)

    # Outputs
    ln_fugacity: List[float]  # ln(fugacity coefficients)
    volumes: List[float]  # Partial molar volumes (cm3/mol)
    total_volume: float  # Total molar volume (cm3/mol)
    epsilon: float = 0.0  # Dielectric constant
    g_function: float = 0.0  # Shock g-function
    adh_factor: float = 0.0  # Debye-Hückel factor


def run_perplex_fluids(
    perplex_path: str,
    T: float,
    P: float,
    species: List[int],
    composition: List[float],
    hybrid_options: Dict[int, int],
) -> FluidState:
    """
    Run Perple_X fluids executable and parse output

    Args:
        perplex_path: Path to Perple_X installation
        T: Temperature (K)
        P: Pressure (bar)
        species: Species indices (1=H2O, 2=CO2, etc.)
        composition: Mole fractions
        hybrid_options: Dict mapping species index to EoS option (ifug)

    Returns:
        FluidState with computed properties
    """

    fluids_exe = os.path.join(perplex_path, "fluids")
    if not os.path.exists(fluids_exe):
        raise FileNotFoundError(f"Perple_X fluids not found at {fluids_exe}")

    # Prepare input for fluids
    # Format:
    # ifug(1) ifug(2) ... (EoS options)
    # T P
    # y(1) y(2) ... (mole fractions)

    ifug_line = " ".join(str(hybrid_options.get(s, 0)) for s in species)
    comp_line = " ".join(str(c) for c in composition)

    input_data = f"{ifug_line}\n{T} {P}\n{comp_line}\n"

    try:
        # Run fluids executable
        result = subprocess.run(
            [fluids_exe], input=input_data, capture_output=True, text=True, timeout=10
        )

        # Parse output
        output_lines = result.stdout.strip().split("\n")

        # TODO: Parse Perple_X output format
        # This depends on the exact output format of fluids
        # For now, placeholder parsing

        state = FluidState(
            T=T,
            P=P,
            composition=composition,
            species=species,
            ln_fugacity=[0.0] * len(species),  # TODO: Parse
            volumes=[0.0] * len(species),  # TODO: Parse
            total_volume=0.0,  # TODO: Parse
        )

        return state

    except subprocess.TimeoutExpired:
        raise RuntimeError(f"Perple_X fluids timed out for T={T} P={P}")
    except Exception as e:
        raise RuntimeError(f"Error running Perple_X: {e}")


def generate_test_cases() -> List[Dict[str, Any]]:
    """
    Generate comprehensive test cases for regression testing

    Returns:
        List of test case specifications
    """

    test_cases = []

    # Test 1: Pure H2O with HSMRK
    test_cases.append(
        {
            "name": "pure_h2o_hsmrk",
            "species": [1],
            "composition": [1.0],
            "hybrid_options": {1: 1},  # ifug=1 for HSMRK
            "conditions": [
                (298.15, 1.0),
                (373.15, 1.0),
                (500.0, 1000.0),
                (600.0, 2000.0),
                (700.0, 5000.0),
            ],
        }
    )

    # Test 2: Pure CO2 with CORK
    test_cases.append(
        {
            "name": "pure_co2_cork",
            "species": [2],
            "composition": [1.0],
            "hybrid_options": {2: 2},  # ifug=2 for CORK
            "conditions": [
                (298.15, 1.0),
                (373.15, 100.0),
                (500.0, 1000.0),
                (600.0, 2000.0),
            ],
        }
    )

    # Test 3: H2O-CO2 binary mixtures
    for xco2 in [0.1, 0.25, 0.5, 0.75, 0.9]:
        test_cases.append(
            {
                "name": f"h2o_co2_xco2_{int(xco2 * 100)}",
                "species": [1, 2],
                "composition": [1.0 - xco2, xco2],
                "hybrid_options": {1: 1, 2: 2},  # HSMRK + CORK
                "conditions": [(500.0, 1000.0), (600.0, 2000.0)],
            }
        )

    # Test 4: P-T grid for H2O-CO2 (70-30)
    pressures = [1.0, 100.0, 500.0, 1000.0, 2000.0, 5000.0]
    temperatures = [298.15, 373.15, 473.15, 573.15, 673.15, 773.15]

    test_cases.append(
        {
            "name": "h2o_co2_pt_grid",
            "species": [1, 2],
            "composition": [0.7, 0.3],
            "hybrid_options": {1: 1, 2: 2},
            "conditions": [(T, P) for T in temperatures for P in pressures],
        }
    )

    # Test 5: Different EoS options for water
    for ifug, name in [(1, "hsmrk"), (3, "brmrk"), (4, "haar")]:
        test_cases.append(
            {
                "name": f"h2o_{name}",
                "species": [1],
                "composition": [1.0],
                "hybrid_options": {1: ifug},
                "conditions": [(500.0, 1000.0), (600.0, 2000.0)],
            }
        )

    # Test 6: Different EoS options for CO2
    for ifug, name in [(2, "cork"), (5, "pseos"), (6, "zhang_duan")]:
        test_cases.append(
            {
                "name": f"co2_{name}",
                "species": [2],
                "composition": [1.0],
                "hybrid_options": {2: ifug},
                "conditions": [(500.0, 1000.0), (600.0, 2000.0)],
            }
        )

    # Test 7: Ternary H2O-CO2-CH4
    test_cases.append(
        {
            "name": "h2o_co2_ch4_ternary",
            "species": [1, 2, 4],
            "composition": [0.6, 0.3, 0.1],
            "hybrid_options": {1: 1, 2: 2, 4: 0},  # HSMRK + CORK + MRK
            "conditions": [(500.0, 1000.0), (600.0, 2000.0)],
        }
    )

    return test_cases


def main():
    parser = argparse.ArgumentParser(
        description="Generate Perple_X reference data for regression testing"
    )
    parser.add_argument(
        "--perplex-path", required=True, help="Path to Perple_X installation directory"
    )
    parser.add_argument(
        "--output",
        default="reference_data.json",
        help="Output JSON file for reference data",
    )
    parser.add_argument(
        "--test-cases", nargs="+", help="Specific test cases to run (default: all)"
    )

    args = parser.parse_args()

    # Check Perple_X installation
    if not os.path.isdir(args.perplex_path):
        print(f"Error: Perple_X path not found: {args.perplex_path}")
        sys.exit(1)

    fluids_exe = os.path.join(args.perplex_path, "fluids")
    if not os.path.exists(fluids_exe):
        print(f"Error: fluids executable not found at {fluids_exe}")
        sys.exit(1)

    print("=" * 70)
    print("Perple_X Reference Data Generation")
    print("=" * 70)
    print(f"Perple_X path: {args.perplex_path}")
    print(f"Output file: {args.output}")
    print()

    # Generate test cases
    all_test_cases = generate_test_cases()

    if args.test_cases:
        test_cases = [tc for tc in all_test_cases if tc["name"] in args.test_cases]
        print(f"Running {len(test_cases)} selected test cases")
    else:
        test_cases = all_test_cases
        print(f"Running all {len(test_cases)} test cases")

    print()

    # Run tests and collect data
    reference_data = []

    for i, test_case in enumerate(test_cases, 1):
        print(f"[{i}/{len(test_cases)}] Running: {test_case['name']}")

        for T, P in test_case["conditions"]:
            try:
                state = run_perplex_fluids(
                    args.perplex_path,
                    T,
                    P,
                    test_case["species"],
                    test_case["composition"],
                    test_case["hybrid_options"],
                )

                reference_data.append(
                    {"test_case": test_case["name"], "state": asdict(state)}
                )

                print(f"  ✓ T={T}K P={P}bar")

            except Exception as e:
                print(f"  ✗ T={T}K P={P}bar - Error: {e}")

    # Save to JSON
    with open(args.output, "w") as f:
        json.dump(reference_data, f, indent=2)

    print()
    print(f"Generated {len(reference_data)} reference data points")
    print(f"Saved to: {args.output}")
    print()
    print("Next steps:")
    print("1. Review reference_data.json")
    print("2. Update test_regression.cpp with reference values")
    print("3. Compile and run: ./test_regression")


if __name__ == "__main__":
    main()
