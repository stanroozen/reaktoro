#!/usr/bin/env python3
"""Verify the a1 parameter fix is correct"""

# Expected values from Excel (DEW_2024)
excel_a1_times_10 = 6.962976  # This is "a1 x 10" column value
excel_Vo = 30.0  # cm³/mol - reference volume at 25°C, 1 bar

# Unit conversions
CAL_TO_J = 4.184
BAR_TO_PA = 1.0e5

# WRONG conversion (what we had before):
a1_wrong = excel_a1_times_10 * 10.0 * CAL_TO_J / BAR_TO_PA
print("❌ WRONG conversion (multiply by 10):")
print(f"   a1 = {a1_wrong:.6e} m³/mol = {a1_wrong * 1e6:.2f} cm³/mol")
print(f"   This gives V° ≈ {a1_wrong * 1e6:.0f} cm³/mol at 25°C, 1 bar")
print(f"   Error: {abs(a1_wrong * 1e6 - excel_Vo) / excel_Vo * 100:.1f}%\n")

# CORRECT conversion (what we have now):
a1_correct = excel_a1_times_10 / 10.0 * CAL_TO_J / BAR_TO_PA
print("✅ CORRECT conversion (divide by 10):")
print(f"   a1 = {a1_correct:.6e} m³/mol = {a1_correct * 1e6:.2f} cm³/mol")

# HKF equation at 25°C, 1 bar for neutral species
T = 298.15  # K
P = 1.0e5  # Pa
psi = 2.6e7  # Pa
theta = 228.0  # K

# From YAML file (DEW 2024 CO2):
a2 = 0.09816951818142544  # J/mol
a3 = 0.00012011192468233996  # J·K/(mol·Pa)
a4 = -0.0012033168788162015  # J·K/mol

# Calculate V0 (for neutral species, Born terms = 0)
V0 = a1_correct + a2 / (psi + P) + (a3 + a4 / (psi + P)) / (T - theta)

print(f"\nHKF equation result at 25°C, 1 bar:")
print(f"   V° = {V0:.6e} m³/mol = {V0 * 1e6:.2f} cm³/mol")
print(f"   Expected from Excel: {excel_Vo} cm³/mol")
print(f"   Error: {abs(V0 * 1e6 - excel_Vo) / excel_Vo * 100:.2f}%")

if abs(V0 * 1e6 - excel_Vo) / excel_Vo * 100 < 5:
    print("\n🎉 SUCCESS! The fix is correct!")
else:
    print("\n❌ Still has errors")
