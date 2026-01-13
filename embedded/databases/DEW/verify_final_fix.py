#!/usr/bin/env python3
"""
Final verification that YAML matches Excel after the scaling fix.
"""

import pandas as pd
import yaml

CAL2J = 4.184
BAR2PA = 1.0e5

# Read Excel
raw = pd.read_excel(
    "Latest_DEW2024.xlsm", sheet_name="Aqueous Species Table", header=None
)
header_names = raw.iloc[2]
df = raw.iloc[4:].copy()
df.columns = header_names
df = df[df["Chemical"].notna()].copy()

co2_excel = df[df["Chemical"].str.contains("CO2", case=False, na=False)].iloc[0]

# Excel values → actual parameters → SI units
excel_a1 = float(co2_excel["a1 x 10"]) * 10.0 * CAL2J / BAR2PA
excel_a2 = float(co2_excel["a2 x 10-2"]) * 1.0e-2 * CAL2J
excel_omega = float(co2_excel["ω x 10-5"]) * 1.0e-5 * CAL2J

# Read YAML
with open("dew2024-aqueous.yaml") as f:
    lines = f.readlines()

# Find CO2 section
in_co2 = False
yaml_a1 = None
yaml_a2 = None
yaml_omega = None

for i, line in enumerate(lines):
    if "Name: CO2,aq" in line:
        in_co2 = True
    elif in_co2:
        if "a1:" in line:
            yaml_a1 = float(line.split(":")[1].strip())
        elif "a2:" in line:
            yaml_a2 = float(line.split(":")[1].strip())
        elif "wref:" in line:
            yaml_omega = float(line.split(":")[1].strip())
        elif "Name:" in line and "CO2" not in line:
            break

print("=" * 80)
print("FINAL VERIFICATION: Excel → YAML")
print("=" * 80)

print(f"\na₁ [J/(mol·Pa)]:")
print(f"  Excel:  {excel_a1:.10e}")
print(f"  YAML:   {yaml_a1:.10e}")
print(f"  Match:  {'✓ YES' if abs(excel_a1 - yaml_a1) < 1e-10 else '✗ NO'}")

print(f"\na₂ [J/mol]:")
print(f"  Excel:  {excel_a2:.10e}")
print(f"  YAML:   {yaml_a2:.10e}")
print(f"  Match:  {'✓ YES' if abs(excel_a2 - yaml_a2) < 1e-10 else '✗ NO'}")

print(f"\nω [J/mol]:")
print(f"  Excel:  {excel_omega:.10e}")
print(f"  YAML:   {yaml_omega:.10e}")
print(f"  Match:  {'✓ YES' if abs(excel_omega - yaml_omega) < 1e-10 else '✗ NO'}")

print("\n" + "=" * 80)
if (
    abs(excel_a1 - yaml_a1) < 1e-10
    and abs(excel_a2 - yaml_a2) < 1e-10
    and abs(excel_omega - yaml_omega) < 1e-10
):
    print("✓✓✓ ALL PARAMETERS MATCH PERFECTLY! ✓✓✓")
else:
    print("✗✗✗ MISMATCH DETECTED! ✗✗✗")
print("=" * 80)
