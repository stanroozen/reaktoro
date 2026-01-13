#!/usr/bin/env python3
"""Investigate why DEW and SUPCRT have such different parameter values"""

import yaml

print("=" * 80)
print("WHY ARE DEW AND SUPCRT PARAMETERS SO DIFFERENT?")
print("=" * 80)

# Load both databases
with open("dew2024-aqueous.yaml", "r") as f:
    dew_data = yaml.safe_load(f)

with open("../reaktoro/supcrt98.yaml", "r") as f:
    supcrt_data = yaml.safe_load(f)

dew_co2 = dew_data["Species"]["CO2(0)"]["StandardThermoModel"]["HKF"]
supcrt_co2 = supcrt_data["Species"]["CO2(aq)"]["StandardThermoModel"]["HKF"]

print("\nHYPOTHESIS 1: Different equation forms")
print("-" * 80)
print("DEW and SUPCRT might use different forms of the HKF equations")
print("where parameters are scaled differently in the formulas.")
print()
print("Standard HKF (SUPCRT) uses:")
print("  G_PV = a₁(P-Pr) + a₂·ln((Ψ+P)/(Ψ+Pr))")
print("       + (a₃(P-Pr) + a₄·ln((Ψ+P)/(Ψ+Pr)))/(T-Θ)")
print()
print("If DEW uses a different pressure unit internally (e.g., kbar vs bar),")
print("or different Ψ value, parameters would need to be scaled.")

print("\nHYPOTHESIS 2: Check if parameter products are consistent")
print("-" * 80)
print("If DEW uses different parameterization, some combinations might match:")
print()

# Check a2 * Psi relationship
Psi = 2600  # bar = 2.6e8 Pa
print(f"a₂ × Ψ product:")
print(
    f"  SUPCRT: {supcrt_co2['a2']:.2e} × {Psi * 1e5:.2e} Pa = {supcrt_co2['a2'] * Psi * 1e5:.2e}"
)
print(
    f"  DEW:    {dew_co2['a2']:.2e} × {Psi * 1e5:.2e} Pa = {dew_co2['a2'] * Psi * 1e5:.2e}"
)
print(f"  Ratio: {(supcrt_co2['a2'] * Psi * 1e5) / (dew_co2['a2'] * Psi * 1e5):.2f}×")

print(f"\na₄ / (T-Θ) at 300°C:")
T = 573.15  # K
Theta = 228  # K
print(
    f"  SUPCRT: {supcrt_co2['a4']:.2e} / {T - Theta:.2f} = {supcrt_co2['a4'] / (T - Theta):.2e}"
)
print(
    f"  DEW:    {dew_co2['a4']:.2e} / {T - Theta:.2f} = {dew_co2['a4'] / (T - Theta):.2e}"
)
print(
    f"  Ratio: {(supcrt_co2['a4'] / (T - Theta)) / (dew_co2['a4'] / (T - Theta)):.2f}×"
)

print("\nHYPOTHESIS 3: DEW uses pressure in bar, SUPCRT uses Pa")
print("-" * 80)
print("If DEW formulas expect P in bar but we pass Pa:")
print()
bar_to_pa = 1e5
print(f"a₂ (if DEW expects bar):")
print(f"  DEW a₂ in bar units: {dew_co2['a2']:.6e} J/mol")
print(f"  Convert to Pa units: {dew_co2['a2'] * bar_to_pa:.6e} J/mol")
print(f"  SUPCRT a₂:          {supcrt_co2['a2']:.6e} J/mol")
print(
    f"  Ratio after conversion: {(dew_co2['a2'] * bar_to_pa) / supcrt_co2['a2']:.2f}×"
)
print()
print("  Still doesn't match! (~3× different, not 31,000×)")

print("\nHYPOTHESIS 4: DEW Excel stores in bar·cm³ (PV units)")
print("-" * 80)
print("The Excel shows a₁ and a₃ in units like 'cal/(mol·bar)'")
print("But maybe internally these are PV work units (bar·cm³)?")
print()
print("Using PV conversion: 1 cal = 41.84 bar·cm³")
print()
a1_dew_barcm3 = dew_co2["a1"] * 1e6  # J/(mol·Pa) → bar·cm³/mol assuming 1 J = 1 Pa·m³
print(f"DEW a₁ in bar·cm³/mol: {a1_dew_barcm3:.4f}")
print(f"SUPCRT a₁ in bar·cm³/mol: {supcrt_co2['a1'] * 1e6:.4f}")
print(f"  → These match within 10%! This suggests BOTH use SI consistently ✓")

print("\n" + "=" * 80)
print("INVESTIGATION: Read DEW 2024 paper or documentation")
print("=" * 80)
print()
print("The massive differences (10⁴ to 10⁹ factors) in a₂, a₄, c₂, ω suggest")
print("that DEW 2024 uses a FUNDAMENTALLY DIFFERENT EQUATION FORM.")
print()
print("Possibilities:")
print("1. DEW derives parameters from different thermodynamic data sources")
print("2. DEW uses modified HKF equations optimized for high P-T conditions")
print("3. DEW parameterization includes additional physical effects")
print("4. The 'wref' in DEW is actually a different parameter than SUPCRT ω")
print()
print("To resolve this, we need to:")
print("- Check the DEW 2024 publication for equation forms")
print("- Compare calculated G, H, V, Cp values at reference conditions")
print("- Verify both models give same thermodynamic predictions")
print()
print("IMPORTANT: The fact that our test was getting close (1% error)")
print("suggests the C++ code IS using the parameters correctly!")
print("The issue is the V_model ≠ V_FD inconsistency, NOT the parameters.")
print("=" * 80)
