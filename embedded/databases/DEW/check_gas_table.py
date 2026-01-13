#!/usr/bin/env python3
"""
Check gas table Excel headers to verify scaling is correct.
"""

import pandas as pd

# Read Excel file
raw = pd.read_excel("Latest_DEW2024.xlsm", sheet_name="Gas Table", header=None)

header_names = raw.iloc[2]
df = raw.iloc[4:].copy()
df.columns = header_names

print("Gas Table Columns:")
for i, col in enumerate(df.columns):
    print(f"  {i}: '{col}'")

print("\nFirst gas species:")
df = df[df["Chemical"].notna()].copy()
first_gas = df.iloc[0]
print(f"\nChemical: {first_gas['Chemical']}")

if "b x 103" in first_gas:
    print(f"b x 103 value: {first_gas['b x 103']}")
if "c x 10-5" in first_gas:
    print(f"c x 10-5 value: {first_gas['c x 10-5']}")

# Check CO2 if it exists
co2_gas = df[df["Chemical"].str.contains("CO2", case=False, na=False)]
if not co2_gas.empty:
    print(f"\nCO2 gas:")
    row = co2_gas.iloc[0]
    print(f"  a = {row['a']}")
    print(f"  b x 103 = {row['b x 103']}")
    print(f"  c x 10-5 = {row['c x 10-5']}")

    CAL2J = 4.184

    # Current interpretation (multiplying by scale)
    b_current = float(row["b x 103"]) * 1e-3 * CAL2J
    c_current = float(row["c x 10-5"]) * 1e-5 * CAL2J

    print(f"\n  Current code gives:")
    print(f"    b (coefficient of T):  {b_current:.6e} J/(mol·K)")
    print(f"    c (coefficient of T²): {c_current:.6e} J/(mol·K²)")
