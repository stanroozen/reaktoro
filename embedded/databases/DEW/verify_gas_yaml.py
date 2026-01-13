#!/usr/bin/env python3
"""
Verify gas YAML matches Excel after scaling fix.
"""

import pandas as pd
import yaml

CAL2J = 4.184

# Read Excel
raw = pd.read_excel("Latest_DEW2024.xlsm", sheet_name="Gas Table", header=None)
header_names = raw.iloc[2]
df = raw.iloc[4:].copy()
df.columns = header_names
df = df[df["Chemical"].notna()].copy()

co2_excel = df[df["Chemical"].str.contains("CO2", case=False, na=False)].iloc[0]

# Excel values → actual parameters → SI units
excel_a = float(co2_excel["a"]) * CAL2J
excel_b = float(co2_excel["b x 103"]) * 1e3 * CAL2J  # multiply by 1000
excel_c = float(co2_excel["c x 10-5"]) * 1e-5 * CAL2J  # multiply by 1e-5

print("=" * 80)
print("EXCEL GAS TABLE (CO2):")
print("=" * 80)
print(f"  Raw values in Excel:")
print(f"    a       = {float(co2_excel['a'])}")
print(
    f"    b x 103 = {float(co2_excel['b x 103'])} → actual b = {float(co2_excel['b x 103'])} × 10³ = {float(co2_excel['b x 103']) * 1e3}"
)
print(
    f"    c x 10-5 = {float(co2_excel['c x 10-5'])} → actual c = {float(co2_excel['c x 10-5'])} × 10⁻⁵ = {float(co2_excel['c x 10-5']) * 1e-5}"
)
print(f"\n  In SI units:")
print(f"    a = {excel_a:.6e} J/(mol·K)")
print(f"    b = {excel_b:.6e} J/(mol·K²)")
print(f"    c = {excel_c:.6e} J/(mol·K³)")

# Read YAML
with open("dew2024-gas.yaml") as f:
    yaml_data = yaml.safe_load(f)

# Find CO2(g)
co2_yaml = None
for key, species in yaml_data["Species"].items():
    if "CO2" in key:
        co2_yaml = species
        break

if co2_yaml:
    nasa = co2_yaml["StandardThermoModel"]["Nasa"]
    poly = nasa["polynomials"][0]

    # NASA polynomial: Cp/R = a1 + a2*T + a3*T²
    # Our Cp = a_J + b_JK*T + c_JK2*T²
    # So: a_J = a1*R, b_JK = a2*R, c_JK2 = a3*R

    R = 8.314462618
    yaml_a = poly["a1"] * R
    yaml_b = poly["a2"] * R
    yaml_c = poly["a3"] * R

    print("\n" + "=" * 80)
    print("YAML GAS DATABASE (CO2):")
    print("=" * 80)
    print(f"  NASA coefficients:")
    print(f"    a1 = {poly['a1']:.6e}")
    print(f"    a2 = {poly['a2']:.6e}")
    print(f"    a3 = {poly['a3']:.6e}")
    print(f"\n  Reconstructed Cp coefficients:")
    print(f"    a = {yaml_a:.6e} J/(mol·K)")
    print(f"    b = {yaml_b:.6e} J/(mol·K²)")
    print(f"    c = {yaml_c:.6e} J/(mol·K³)")

    print("\n" + "=" * 80)
    print("VERIFICATION:")
    print("=" * 80)
    print(
        f"  a: Excel={excel_a:.6e}, YAML={yaml_a:.6e}, Match={'✓' if abs(excel_a - yaml_a) < 1e-6 else '✗'}"
    )
    print(
        f"  b: Excel={excel_b:.6e}, YAML={yaml_b:.6e}, Match={'✓' if abs(excel_b - yaml_b) < 1e-3 else '✗'}"
    )
    print(
        f"  c: Excel={excel_c:.6e}, YAML={yaml_c:.6e}, Match={'✓' if abs(excel_c - yaml_c) < 1e-9 else '✗'}"
    )
