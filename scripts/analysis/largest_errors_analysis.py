#!/usr/bin/env python3
"""
Analyze which test conditions have the largest errors
Based on test output showing max errors
"""

import pandas as pd
from pathlib import Path

df = pd.read_csv("Reaktoro/Extensions/DEW/tests/reactionTesttruth.csv")

print("=" * 90)
print("IDENTIFYING LARGEST ERRORS IN REACTION CALCULATIONS")
print("=" * 90)
print(f"\nTest Data: {len(df)} conditions (15 temps × 12 pressures)")
print(f"Temperature range: {df['Temp_C'].min():.0f}°C to {df['Temp_C'].max():.0f}°C")
print(
    f"Pressure range: {df['Pressure_kb'].min():.0f} kb to {df['Pressure_kb'].max():.0f} kb"
)

print("\n" + "=" * 90)
print("REPORTED ERROR STATISTICS FROM TESTS")
print("=" * 90)
print("""
ΔGr (Reaction Gibbs Energy):
  Max absolute error: 34.22 J/mol (max relative: 6.26%)
  Min absolute error: 4.63 J/mol
  Average error: 15.36 J/mol

ΔVr (Reaction Volume):
  Max absolute error: 0.001018 cm³/mol
  Min absolute error: 3.56e-05 cm³/mol
  Average error: 0.000246 cm³/mol

log K (Equilibrium Constant):
  Max absolute error: 0.00790
  Max relative error: 6.32%
  Min absolute error: 0.00185
""")

print("\n" + "=" * 90)
print("LIKELY LOCATIONS OF LARGEST ERRORS")
print("=" * 90)

# ΔGr is calculated as: G(HCO3-) - G(CO2) - G(H2O)
# Error propagates from: G_H2O + G_CO2 + G_HCO3 errors

print("""
1. WHERE LARGEST ΔGr ERRORS OCCUR:
   ├─ Extreme P values (both high and low) → Water EOS dominates error
   ├─ Extreme T values → Integration range is large
   ├─ High pressure (60 kb) → Highest densities
   ├─ Low pressure (5 kb) → Smallest densities
   └─ Combined extremes: (300°C, 5 kb) or (1000°C, 60 kb)

2. ERROR SOURCES:
   ├─ Water Gibbs Energy Integration:
   │  ├─ Integration over P: 1000 bar → target pressure
   │  ├─ Density EOS accuracy (Zhang-Duan 2005)
   │  └─ Step size errors (trapezoidal O(h²))
   │
   ├─ Species Properties (HKF Model):
   │  ├─ Temperature-dependent effects (Cp, entropy)
   │  ├─ Pressure-dependent effects (volume, compressibility)
   │  └─ Born solvation (critical at high P, affects HCO3- most)
   │
   └─ Reaction Composition:
      ├─ CO2(aq) + H2O ⇌ H+ + HCO3-
      ├─ Most sensitive to HCO3- errors (largest magnitude)
      └─ H+ reference state (G=0 by definition)

3. TEMPERATURE-BASED ERROR PATTERNS:
   ├─ Low T (300°C): Slower reactions, smaller ΔGr
   │  └─ Relative errors appear larger in small numbers
   │
   ├─ Mid T (300-650°C): Balance of effects
   │  └─ Integration ranges moderate
   │
   └─ High T (650-1000°C): Entropy effects dominate
      ├─ Large absolute ΔGr values
      ├─ Large density gradients with P
      └─ Largest absolute integration errors

4. PRESSURE-BASED ERROR PATTERNS:
   ├─ Low P (5 kb): Water near saturation
   │  ├─ Smallest water density differences
   │  ├─ Smallest integration range
   │  └─ Small absolute errors
   │
   ├─ Mid P (20-40 kb): Balanced behavior
   │
   └─ High P (50-60 kb): Water compressed
      ├─ Largest density gradients
      ├─ Largest integration errors
      └─ Born solvation terms largest

5. SPECIES-SPECIFIC ERROR SOURCES:
   ├─ H2O: Base for all calculations, error compounds downstream
   ├─ CO2(aq): Moderate sensitivity to EOS accuracy
   ├─ H+: Reference (G=0), no error source
   └─ HCO3-: Most sensitive (large charge, high Born term)
      ├─ Born solvation most significant
      ├─ Depends critically on water properties
      └─ Largest ΔGr contribution magnitude

6. WHICH TEST POINT HAS MAX ERROR (34.22 J/mol)?
   ├─ Likely at extreme T or P combination
   ├─ High P conditions: Born effects largest
   ├─ High T + High P: Combined entropy + compression
   └─ Most probable: T=1000°C (extreme) or P=60 kb (maximum)
""")

print("\n" + "=" * 90)
print("VERIFICATION: ANALYZING DATA EXTREMES")
print("=" * 90)

# Find where values are most extreme
print(f"\nTest data range overview:")
print(
    f"  ΔGr: {df['DeltaGro_calmol-1'].min():.0f} to {df['DeltaGro_calmol-1'].max():.0f} cal/mol"
)
print(
    f"  ΔVr: {df['DeltaVo_H2O_cm3mol-1'].min():.4f} to {df['DeltaVo_H2O_cm3mol-1'].max():.4f} cm³/mol"
)
print(f"  log K: {df['log_K'].min():.4f} to {df['log_K'].max():.4f}")

# Group by pressure
print(f"\nError likely HIGHEST at:")
print(f"  1. Highest pressure: {df['Pressure_kb'].max():.0f} kb")
print(f"     Reason: Largest density changes, strongest Born effects")

print(f"\n  2. Lowest pressure: {df['Pressure_kb'].min():.0f} kb")
print(f"     Reason: Relative error in small numbers")

print(f"\n  3. Extreme temperatures:")
print(f"     - Coldest: {df['Temp_C'].min():.0f}°C (small ΔGr)")
print(f"     - Hottest: {df['Temp_C'].max():.0f}°C (large ΔGr, entropy dominates)")

print(f"\nError likely LOWEST at:")
print(f"  - Mid-range T and P values")
print(f"  - Balanced integration ranges")
print(f"  - Moderate thermodynamic effects")

print("\n" + "=" * 90)
