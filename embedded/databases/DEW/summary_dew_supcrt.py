#!/usr/bin/env python3
"""
Summary of DEW vs SUPCRT comparison after fixing the a1 scaling.
"""

print("=" * 80)
print("PARAMETER COMPARISON: DEW vs SUPCRT for CO2(aq)")
print("=" * 80)

# DEW values (from Excel, properly scaled)
dew_a1_cal = 69.63  # cal/(mol·bar)
dew_a2_cal = 0.02346  # cal/mol
dew_omega_cal = -8.0e-6  # cal/mol

# SUPCRT values (from previous analysis, converted back to cal)
supcrt_a1_cal = 62.5  # cal/(mol·bar)
supcrt_a2_cal = 747.0  # cal/mol
supcrt_omega_cal = -2000.0  # cal/mol

print(f"\nIN CALORIES:")
print(f"  Parameter    DEW Value        SUPCRT Value     Ratio")
print(f"  ---------    ---------        ------------     -----")
print(
    f"  a₁           {dew_a1_cal:10.2f}       {supcrt_a1_cal:10.2f}         {supcrt_a1_cal / dew_a1_cal:6.2f}×"
)
print(
    f"  a₂           {dew_a2_cal:10.5f}      {supcrt_a2_cal:10.2f}         {supcrt_a2_cal / dew_a2_cal:6.0f}×"
)
print(
    f"  ω            {dew_omega_cal:10.2e}     {supcrt_omega_cal:10.2e}     {supcrt_omega_cal / dew_omega_cal:6.0e}×"
)

print("\n" + "=" * 80)
print("CONCLUSION:")
print("=" * 80)
print("""
After fixing the Excel header interpretation:
- a₁ is now similar between DEW and SUPCRT (0.9× ratio) ✓
- BUT a₂ differs by 31,000×
- AND ω differs by 250,000,000×

These are DIFFERENT DATABASES using DIFFERENT PARAMETERIZATIONS!

The huge ω difference (250 million×) indicates they're using different
conventions for the Born function. This is NOT a unit conversion issue.

DEW and SUPCRT are incompatible - you cannot mix their parameters.
""")
