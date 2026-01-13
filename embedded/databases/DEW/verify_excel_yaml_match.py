#!/usr/bin/env python3
"""Verify YAML matches Excel after unit conversions"""

import yaml
import pandas as pd

# Constants
CAL2J = 4.184
BAR2PA = 1.0e5

print("=" * 70)
print("VERIFICATION: Excel → YAML Unit Conversion for CO2(0)")
print("=" * 70)

# Load Excel data
excel_file = "Latest_DEW2024.xlsm"
raw = pd.read_excel(excel_file, sheet_name="Aqueous Species Table", header=None)
header_names = raw.iloc[2]
df = raw.iloc[4:].copy()
df.columns = header_names
df = df[df["Chemical"].notna()].copy()
df = df.iloc[:, 1:].reset_index(drop=True)

# Find CO2
co2_row = df[df["Chemical"] == "CO2,aq"].iloc[0]

print("\n1. EXCEL VALUES (as displayed in cells):")
print("-" * 70)
print(f"  a₁ × 10        = {co2_row['a1 x 10']:.4f}")
print(f"  a₂ × 10⁻²      = {co2_row['a2 x 10-2']:.4f}")
print(f"  a₃             = {co2_row['a3']:.4f}")
print(f"  a₄ × 10⁻⁴      = {co2_row['a4 x 10-4']:.4f}")
print(f"  c₁             = {co2_row['c1']:.4f}")
print(f"  c₂ × 10⁻⁴      = {co2_row['c2 x 10-4']:.4f}")
print(f"  ω × 10⁻⁵       = {co2_row['ω x 10-5']:.4f}")

print("\n2. ACTUAL PHYSICAL VALUES (after removing Excel scaling):")
print("-" * 70)
a1_cal_per_bar = float(co2_row["a1 x 10"]) / 10.0
a2_cal = float(co2_row["a2 x 10-2"]) * 1.0e-2
a3_cal_per_bar = float(co2_row["a3"])
a4_cal = float(co2_row["a4 x 10-4"]) * 1.0e-4
c1_cal = float(co2_row["c1"])
c2_cal = float(co2_row["c2 x 10-4"]) * 1.0e-4
wref_cal = float(co2_row["ω x 10-5"]) * 1.0e-5

print(f"  a₁ = {a1_cal_per_bar:.6f} cal/(mol·bar)")
print(f"  a₂ = {a2_cal:.6f} cal/mol")
print(f"  a₃ = {a3_cal_per_bar:.6f} cal/(K·mol·bar)")
print(f"  a₄ = {a4_cal:.6f} cal·K/mol")
print(f"  c₁ = {c1_cal:.6f} cal/(mol·K)")
print(f"  c₂ = {c2_cal:.6f} cal·K/mol")
print(f"  ω  = {wref_cal:.6e} cal/mol")

print("\n3. PYTHON SCRIPT CONVERSION (to SI units):")
print("-" * 70)
a1_si = a1_cal_per_bar * CAL2J / BAR2PA
a2_si = a2_cal * CAL2J
a3_si = a3_cal_per_bar * CAL2J / BAR2PA
a4_si = a4_cal * CAL2J
c1_si = c1_cal * CAL2J
c2_si = c2_cal * CAL2J
wref_si = wref_cal * CAL2J

print(f"  a₁ = {a1_si:.6e} J/(mol·Pa) = m³/mol")
print(f"  a₂ = {a2_si:.6e} J/mol")
print(f"  a₃ = {a3_si:.6e} J·K/(mol·Pa)")
print(f"  a₄ = {a4_si:.6e} J·K/mol")
print(f"  c₁ = {c1_si:.6e} J/(mol·K)")
print(f"  c₂ = {c2_si:.6e} J·K/mol")
print(f"  ω  = {wref_si:.6e} J/mol")

print("\n4. YAML DATABASE VALUES:")
print("-" * 70)
with open("dew2024-aqueous.yaml", "r") as f:
    data = yaml.safe_load(f)

co2_yaml = data["Species"]["CO2(0)"]["StandardThermoModel"]["HKF"]

print(f"  a₁ = {co2_yaml['a1']:.6e} J/(mol·Pa)")
print(f"  a₂ = {co2_yaml['a2']:.6e} J/mol")
print(f"  a₃ = {co2_yaml['a3']:.6e} J·K/(mol·Pa)")
print(f"  a₄ = {co2_yaml['a4']:.6e} J·K/mol")
print(f"  c₁ = {co2_yaml['c1']:.6e} J/(mol·K)")
print(f"  c₂ = {co2_yaml['c2']:.6e} J·K/mol")
print(f"  ω  = {co2_yaml['wref']:.6e} J/mol")

print("\n5. VERIFICATION (Python conversion == YAML):")
print("-" * 70)
tol = 1e-10
checks = [
    ("a₁", a1_si, co2_yaml["a1"]),
    ("a₂", a2_si, co2_yaml["a2"]),
    ("a₃", a3_si, co2_yaml["a3"]),
    ("a₄", a4_si, co2_yaml["a4"]),
    ("c₁", c1_si, co2_yaml["c1"]),
    ("c₂", c2_si, co2_yaml["c2"]),
    ("ω", wref_si, co2_yaml["wref"]),
]

all_match = True
for name, expected, actual in checks:
    match = abs(expected - actual) < tol * max(abs(expected), abs(actual), 1e-15)
    status = "✓ MATCH" if match else "✗ MISMATCH"
    print(f"  {name}: {status}")
    if not match:
        print(f"      Expected: {expected:.10e}")
        print(f"      Actual:   {actual:.10e}")
        print(f"      Diff:     {expected - actual:.10e}")
        all_match = False

print("\n" + "=" * 70)
if all_match:
    print("✓ SUCCESS: YAML perfectly matches Excel after unit conversions!")
else:
    print("✗ FAILURE: Some values don't match - check conversion formulas!")
print("=" * 70)

print("\n6. THERMODYNAMIC USAGE (what C++ code needs):")
print("-" * 70)
print(f"  ω (YAML):     {co2_yaml['wref']:.6e} J/mol")
print(f"  ω (×1e10):    {co2_yaml['wref'] * 1e10:.6e} J/mol")
print(f"  ω (expected): ~-334,700 J/mol for DEW calculations")
print("\nNOTE: The ×1e10 factor in DEWDatabase.cpp is needed because")
print("      the DEW/HKF formalism uses omega in a scaled form.")
print("=" * 70)
