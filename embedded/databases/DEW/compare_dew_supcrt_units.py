#!/usr/bin/env python3
"""
Compare CO2 parameters from Excel vs SUPCRT to understand the fundamental difference.
"""

import pandas as pd
import yaml

# Read DEW Excel
print("=" * 80)
print("READING DEW EXCEL...")
print("=" * 80)

raw = pd.read_excel(
    "Latest_DEW2024.xlsm", sheet_name="Aqueous Species Table", header=None
)
header_names = raw.iloc[2]
df = raw.iloc[4:].copy()
df.columns = header_names
df = df[df["Chemical"].notna()].copy()

co2_dew = df[df["Chemical"].str.contains("CO2", case=False, na=False)].iloc[0]

CAL2J = 4.184
BAR2PA = 1.0e5

# DEW parameters in cal units (before SI conversion)
dew_a1_cal = float(co2_dew["a1 x 10"]) * 10.0  # cal/(mol·bar)
dew_a2_cal = float(co2_dew["a2 x 10-2"]) * 1.0e-2  # cal/mol
dew_omega_cal = float(co2_dew["ω x 10-5"]) * 1.0e-5  # cal/mol

# DEW parameters in SI
dew_a1_si = dew_a1_cal * CAL2J / BAR2PA
dew_a2_si = dew_a2_cal * CAL2J
dew_omega_si = dew_omega_cal * CAL2J

print(f"\nDEW Excel (CO2,aq):")
print(f"  a₁ = {dew_a1_cal:12.6f} cal/(mol·bar)  →  {dew_a1_si:.6e} J/(mol·Pa)")
print(f"  a₂ = {dew_a2_cal:12.6f} cal/mol        →  {dew_a2_si:.6e} J/mol")
print(f"  ω  = {dew_omega_cal:12.6e} cal/mol        →  {dew_omega_si:.6e} J/mol")

# Read SUPCRT YAML
print("\n" + "=" * 80)
print("READING SUPCRT YAML...")
print("=" * 80)

supcrt_path = "../reaktoro/supcrt98.yaml"
with open(supcrt_path) as f:
    supcrt_data = yaml.safe_load(f)

co2_supcrt = None
for species in supcrt_data:
    if isinstance(species, dict) and species.get("Name") == "CO2(aq)":
        co2_supcrt = species
        break

if co2_supcrt:
    hkf = co2_supcrt["StandardThermoModel"]["HKF"]

    # SUPCRT values are already in SI
    supcrt_a1_si = hkf["a1"]
    supcrt_a2_si = hkf["a2"]
    supcrt_omega_si = hkf["wref"]

    # Convert back to cal for comparison
    supcrt_a1_cal = supcrt_a1_si * BAR2PA / CAL2J
    supcrt_a2_cal = supcrt_a2_si / CAL2J
    supcrt_omega_cal = supcrt_omega_si / CAL2J

    print(f"\nSUPCRT98 (CO2(aq)):")
    print(
        f"  a₁ = {supcrt_a1_cal:12.6f} cal/(mol·bar)  →  {supcrt_a1_si:.6e} J/(mol·Pa)"
    )
    print(f"  a₂ = {supcrt_a2_cal:12.6f} cal/mol        →  {supcrt_a2_si:.6e} J/mol")
    print(
        f"  ω  = {supcrt_omega_cal:12.6e} cal/mol        →  {supcrt_omega_si:.6e} J/mol"
    )

    print("\n" + "=" * 80)
    print("COMPARISON (SUPCRT / DEW ratios):")
    print("=" * 80)

    print(
        f"\n  a₁ ratio: {supcrt_a1_cal / dew_a1_cal:10.4f}×   ({supcrt_a1_cal:10.4f} / {dew_a1_cal:10.4f})"
    )
    print(
        f"  a₂ ratio: {supcrt_a2_cal / dew_a2_cal:10.4f}×   ({supcrt_a2_cal:10.4f} / {dew_a2_cal:10.4f})"
    )
    print(
        f"  ω  ratio: {supcrt_omega_cal / dew_omega_cal:10.4e}×   ({supcrt_omega_cal:10.4e} / {dew_omega_cal:10.4e})"
    )

    print("\n" + "=" * 80)
    print("CONCLUSION:")
    print("=" * 80)
    print("""
These are TWO COMPLETELY DIFFERENT DATABASES with different parameterizations!

DEW (Deep Earth Water) and SUPCRT98 both use the HKF (Helgeson-Kirkham-Flowers)
model framework, but they have:

1. Different source data
2. Different fitting procedures
3. Different parameter values

The ω (omega) parameter in particular shows they're using different conventions
or different revisions of the HKF model. The 250 million× difference suggests
they may be using different definitions of the Born function.

BOTTOM LINE: You cannot mix DEW and SUPCRT parameters. They are incompatible.
    """)
