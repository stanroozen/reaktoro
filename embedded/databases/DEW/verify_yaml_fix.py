#!/usr/bin/env python3
"""Quick verification that the YAML fix worked"""

import yaml

# Read the YAML file
with open("dew2024-aqueous.yaml", "r") as f:
    data = yaml.safe_load(f)

# Find CO2,aq
for name, species in data["Species"].items():
    if name == "CO2_aq":
        print("CO2,aq HKF parameters:")
        print(f"  a1 = {species['HKF']['a1']:.6e} J/(mol·Pa) = m³/mol")
        print(f"  a2 = {species['HKF']['a2']:.6e} J/mol")
        print(f"  a3 = {species['HKF']['a3']:.6e} J·K/(mol·Pa)")
        print(f"  a4 = {species['HKF']['a4']:.6e} J·K/mol")

        # Calculate V0 at 25°C, 1 bar using HKF equation
        a1 = species["HKF"]["a1"]
        a2 = species["HKF"]["a2"]
        a3 = species["HKF"]["a3"]
        a4 = species["HKF"]["a4"]

        T = 298.15  # K
        P = 1.0e5  # Pa
        psi = 2.6e7  # Pa
        theta = 228.0  # K

        # For neutral species, Born terms ≈ 0
        V0 = a1 + a2 / (psi + P) + (a3 + a4 / (psi + P)) / (T - theta)

        print(f"\nCalculated V0 at 25°C, 1 bar:")
        print(f"  V0 = {V0:.6e} m³/mol")
        print(f"  V0 = {V0 * 1e6:.2f} cm³/mol")
        print(f"\nExpected from Excel: V° = 30.0 cm³/mol")

        error_pct = abs(V0 * 1e6 - 30.0) / 30.0 * 100
        print(f"Error: {error_pct:.2f}%")

        if error_pct < 5:
            print("\n✓ SUCCESS: Volume calculation is correct!")
        else:
            print("\n✗ FAILED: Volume calculation is still wrong!")

        break
