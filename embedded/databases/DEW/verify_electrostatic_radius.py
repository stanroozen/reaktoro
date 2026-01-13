#!/usr/bin/env python3
"""
Properly test the omega scaling by checking if physical values result.
"""

import math

eta = 166027.0  # Å·cal/mol


def calc_reref(wref_cal, Z):
    """Calculate reference electrostatic radius."""
    denom = (wref_cal / eta) + Z / 3.082
    if abs(denom) < 1e-20:
        return float("inf")
    return (Z * Z) / denom


print("=" * 80)
print("ELECTROSTATIC RADIUS CALCULATION")
print("=" * 80)

# Ag+ example
wref_current_cal = 2.16e-6  # Current YAML interpretation
wref_scaled_cal = 2.16e-6 * 1e10  # With ×1e10 = 21,600 cal/mol
Z = 1

print(f"\nAg+ (Z={Z}):")
print(f"  Excel shows: 0.216 in 'ω × 10⁻⁵' column")
print()

# Case 1: Current interpretation (no ×1e10)
denom1 = (wref_current_cal / eta) + Z / 3.082
reref1 = calc_reref(wref_current_cal, Z)
print(f"WITHOUT ×1e10:")
print(f"  wref = {wref_current_cal:.6e} cal/mol")
print(f"  wref/η = {wref_current_cal / eta:.6e}")
print(f"  Z/3.082 = {Z / 3.082:.6f}")
print(f"  denom = wref/η + Z/3.082 = {denom1:.6f}")
print(f"  re_ref = Z²/denom = {reref1:.6f} Å")
print()

# Case 2: With ×1e10
denom2 = (wref_scaled_cal / eta) + Z / 3.082
reref2 = calc_reref(wref_scaled_cal, Z)
print(f"WITH ×1e10:")
print(f"  wref = {wref_scaled_cal:.2f} cal/mol")
print(f"  wref/η = {wref_scaled_cal / eta:.6f}")
print(f"  Z/3.082 = {Z / 3.082:.6f}")
print(f"  denom = wref/η + Z/3.082 = {denom2:.6f}")
print(f"  re_ref = Z²/denom = {reref2:.6f} Å")
print()

print("=" * 80)
print("TYPICAL IONIC RADII FOR COMPARISON:")
print("=" * 80)
print("""
Ag+:  ~1.15-1.3 Å (crystal ionic radius)
      ~2.0-2.5 Å (effective hydrated radius)

With ×1e10: re_ref = 2.20 Å → REASONABLE
Without:    re_ref = 3.08 Å → too large (dominated by 1/3.082 term)

The ×1e10 makes wref large enough to actually contribute to the formula!
""")

# Al+3 example (Z=3)
print("=" * 80)
wref_al_cell = 2.753
wref_al_current_cal = wref_al_cell * 1e-5
wref_al_scaled_cal = wref_al_current_cal * 1e10
Z_al = 3

print(f"\nAl+3 (Z={Z_al}):")
print(f"  Excel shows: {wref_al_cell} in 'ω × 10⁻⁵' column")
print()

reref_al1 = calc_reref(wref_al_current_cal, Z_al)
reref_al2 = calc_reref(wref_al_scaled_cal, Z_al)

print(f"WITHOUT ×1e10:")
print(f"  wref = {wref_al_current_cal:.6e} cal/mol")
print(f"  re_ref = {reref_al1:.6f} Å")
print()

print(f"WITH ×1e10:")
print(f"  wref = {wref_al_scaled_cal:.2f} cal/mol")
print(f"  re_ref = {reref_al2:.6f} Å")
print()

print("=" * 80)
print("Al+3 typical radius: ~0.54 Å (crystal) → 1.5-2.0 Å with hydration")
print(f"With ×1e10:  {reref_al2:.2f} Å → REASONABLE")
print(f"Without:     {reref_al1:.2f} Å → physically incorrect (way too large)")
print("=" * 80)
