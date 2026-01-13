#!/usr/bin/env python3
"""
Check what the correct magnitude of wref should be by looking at expected omega values.
We have CO2 with wref = -3.347e-05 J/mol in YAML (after Python conversion).
The Excel shows -0.8 in the "ω × 10⁻⁵" column.

Question: Is the intended ω value:
  A) -0.8 × 10⁻⁵ = -8e-6 cal/mol = -3.347e-5 J/mol (current YAML)
  B) -0.8 cal/mol = -3.347 J/mol (×1e5 larger)
  C) Something else entirely?
"""

import pandas as pd

# Read Excel
raw = pd.read_excel(
    "Latest_DEW2024.xlsm", sheet_name="Aqueous Species Table", header=None
)
header_names = raw.iloc[2]
df = raw.iloc[4:].copy()
df.columns = header_names
df = df[df["Chemical"].notna()].copy()

print("=" * 80)
print("CHECKING OMEGA MAGNITUDE")
print("=" * 80)

# Get several species with different omega values
species_sample = df.head(10)

print("\nSample species from Excel:")
print(
    f"{'Species':<20} {'ω × 10⁻⁵ (cell)':<15} {'Actual ω (cal/mol)':<20} {'ω (J/mol)':<15}"
)
print("-" * 80)

CAL2J = 4.184

for _, row in species_sample.iterrows():
    name = row["Chemical"]
    omega_cell = float(row["ω x 10-5"])

    # Current interpretation: cell × 10⁻⁵ × CAL2J
    omega_cal_current = omega_cell * 1e-5
    omega_J_current = omega_cal_current * CAL2J

    print(
        f"{name:<20} {omega_cell:<15.6f} {omega_cal_current:<20.6e} {omega_J_current:<15.6e}"
    )

print("\n" + "=" * 80)
print("COMPARISON WITH TYPICAL HKF VALUES")
print("=" * 80)
print("""
From SUPCRT98 for CO2(aq): ω = -8,368 J/mol
From DEW 2024 YAML:        ω = -3.347e-05 J/mol
From DEW with ×1e10:       ω = -334,720 J/mol

The ×1e10 scaling makes DEW ω ~40× LARGER than SUPCRT.
This suggests DEW may use a different Born function formulation.

Let's check the DEW constant η = 166,027 Å·cal/mol to understand the scaling.
""")

# Check a charged species
charged = df[df["Z"].notna() & (df["Z"] != 0)].head(5)

print("\nCharged species omega values:")
print(
    f"{'Species':<15} {'Z':<5} {'ω × 10⁻⁵':<12} {'ω (J/mol, current)':<20} {'ω with ×1e10':<15}"
)
print("-" * 85)

for _, row in charged.iterrows():
    name = row["Chemical"]
    Z = float(row["Z"])
    omega_cell = float(row["ω x 10-5"])
    omega_J = omega_cell * 1e-5 * CAL2J
    omega_scaled = omega_J * 1e10

    print(
        f"{name:<15} {Z:<5.0f} {omega_cell:<12.6f} {omega_J:<20.6e} {omega_scaled:<15.2f}"
    )

print("\n" + "=" * 80)
print("The ×1e10 factor seems to be REQUIRED for DEW Born calculations.")
print("Without it, omega values would be ~1e-5 J/mol, which is physically too small.")
print("=" * 80)
