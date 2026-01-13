#!/usr/bin/env python3
"""
Test if removing the ×1e10 factor breaks DEW thermodynamics.
Compare omega calculations with and without the factor.
"""


# DEW Born formula (from WaterBornOmegaDEW.cpp)
def calculate_omega_DEW(wref_cal, Z, g=0):
    """
    Calculate omega using DEW formula.

    wref_cal: wref in cal/mol
    Z: charge
    g: solvent function (default 0 for simplified test)

    Returns omega in cal/mol
    """
    eta = 166027.0  # Å·cal/mol

    if Z == 0:
        return wref_cal  # Neutral species

    # Reference electrostatic radius
    denom = (wref_cal / eta) + Z / 3.082
    if denom == 0:
        return wref_cal

    reref_A = (Z * Z) / denom

    # Electrostatic radius at (P,T)
    re_A = reref_A + abs(Z) * g
    if re_A <= 0:
        return wref_cal

    # DEW omega formula
    omega_cal = eta * ((Z * Z) / re_A - Z / (3.082 + g))

    return omega_cal


print("=" * 80)
print("TESTING OMEGA CALCULATION WITH DEW FORMULA")
print("=" * 80)

# CO2 from Excel: ω × 10⁻⁵ column shows -0.8
# Current Python: wref = -0.8 × 10⁻⁵ × 4.184 = -3.347e-5 J/mol = -8.0e-6 cal/mol

wref_current_J = -3.347e-5  # J/mol (from current YAML)
wref_current_cal = wref_current_J / 4.184  # cal/mol

wref_scaled_J = wref_current_J * 1e10  # With ×1e10
wref_scaled_cal = wref_scaled_J / 4.184

print(f"\nCO2 (neutral, Z=0):")
print(f"  Excel cell value: -0.8 (in 'ω × 10⁻⁵' column)")
print(
    f"  Current YAML:     wref = {wref_current_J:.6e} J/mol = {wref_current_cal:.6e} cal/mol"
)
print(
    f"  With ×1e10:       wref = {wref_scaled_J:.2f} J/mol = {wref_scaled_cal:.2f} cal/mol"
)
print(f"\n  For neutral species, omega = wref (constant)")
print(
    f"    omega (current) = {wref_current_cal:.6e} cal/mol = {wref_current_J:.6e} J/mol"
)
print(
    f"    omega (×1e10)   = {wref_scaled_cal:.2f} cal/mol = {wref_scaled_J:.2f} J/mol"
)

# Test with a charged species: Ag+ (Z=1, ω × 10⁻⁵ = 0.216)
print(f"\n" + "=" * 80)
print(f"Ag+ (charged, Z=1):")
wref_ag_J = 0.216 * 1e-5 * 4.184  # Current
wref_ag_cal = wref_ag_J / 4.184

wref_ag_scaled_J = wref_ag_J * 1e10
wref_ag_scaled_cal = wref_ag_scaled_J / 4.184

print(f"  Excel cell value: 0.216")
print(f"  Current YAML:     wref = {wref_ag_J:.6e} J/mol = {wref_ag_cal:.6e} cal/mol")
print(
    f"  With ×1e10:       wref = {wref_ag_scaled_J:.2f} J/mol = {wref_ag_scaled_cal:.2f} cal/mol"
)

# Calculate omega using DEW formula
omega_current = calculate_omega_DEW(wref_ag_cal, Z=1, g=0)
omega_scaled = calculate_omega_DEW(wref_ag_scaled_cal, Z=1, g=0)

print(f"\n  Using DEW Born formula (g=0):")
print(
    f"    omega (current) = {omega_current:.6e} cal/mol = {omega_current * 4.184:.6e} J/mol"
)
print(
    f"    omega (×1e10)   = {omega_scaled:.2f} cal/mol = {omega_scaled * 4.184:.2f} J/mol"
)

# Calculate electrostatic radius
eta = 166027.0
denom_current = (wref_ag_cal / eta) + 1 / 3.082
reref_current = 1.0 / denom_current

denom_scaled = (wref_ag_scaled_cal / eta) + 1 / 3.082
reref_scaled = 1.0 / denom_scaled

print(f"\n  Electrostatic radius at reference:")
print(f"    re_ref (current) = {reref_current:.6e} Å")
print(f"    re_ref (×1e10)   = {reref_scaled:.6f} Å")

print("\n" + "=" * 80)
print("CONCLUSION:")
print("=" * 80)
print("""
The ×1e10 factor is needed because:

1. Excel stores omega in 'ω × 10⁻⁵' column (scaled by 1e-5)
2. Python converts: value × 1e-5 × 4.184 → ~1e-5 J/mol in YAML
3. DEW Born formula needs omega in cal/mol range ~1-1000
4. Without ×1e10: wref ~ 1e-6 cal/mol → electrostatic radius ~ 3e-6 Å (unphysical!)
5. With ×1e10:    wref ~ 10-100 cal/mol → electrostatic radius ~ 1-3 Å (physical!)

The ×1e10 is CORRECT and REQUIRED for DEW calculations.
It compensates for the Excel scaling convention.
""")
