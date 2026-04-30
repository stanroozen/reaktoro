#!/usr/bin/env python3
"""
Analyze if per-species G0 errors are consistent across all T,P conditions
Reaktoro DEW model vs Excel truth data
"""

import pandas as pd
import numpy as np

csv_file = "Reaktoro/Extensions/DEW/tests/reactionTesttruth.csv"

print("Loading CSV truth data...")
df = pd.read_csv(csv_file)

CAL_TO_J = 4.184

# Extract truth data
g_h2o_truth = df["DeltaGo_H2O_calmol-1"].values * CAL_TO_J
g_co2_truth = df["DeltaGo_CO2_aq_calmol-1"].values * CAL_TO_J
g_h_truth = df["DeltaGo_H+_calmol-1"].values * CAL_TO_J
g_hco3_truth = df["DeltaGo_HCO3-_calmol-1"].values * CAL_TO_J
g_rxn_truth = df["DeltaGro_calmol-1"].values * CAL_TO_J

T_C = df["Temp_C"].values
P_kb = df["Pressure_kb"].values

print(f"Loaded {len(df)} test conditions\n")

print("=" * 100)
print("QUESTION: Are per-species errors CONSISTENT across all T,P conditions?")
print("=" * 100)

# From test output, we know ONE data point: T=650, P=15
# At that point:
# H2O error = 0
# CO2 error = +5
# H+ error = 0
# HCO3- error = -9

err_h2o_ref = 0  # at T=650, P=15
err_co2_ref = 5  # at T=650, P=15
err_h_ref = 0  # at T=650, P=15
err_hco3_ref = -9  # at T=650, P=15

print("\nOBSERVED at T=650°C, P=15 kb (from test DEBUG output):")
print(f"  H2O error:    {err_h2o_ref:+.1f} J/mol")
print(f"  CO2 error:    {err_co2_ref:+.1f} J/mol")
print(f"  H+ error:     {err_h_ref:+.1f} J/mol")
print(f"  HCO3- error:  {err_hco3_ref:+.1f} J/mol")

# If errors are constant, then ΔGr_error should also be constant:
expected_rxn_error = (err_h_ref + err_hco3_ref) - (err_h2o_ref + err_co2_ref)

print(f"\nIF errors are CONSTANT everywhere:")
print(f"  DGr_error should always be: {expected_rxn_error:+.1f} J/mol")

print("\n" + "=" * 100)
print("TESTING THE HYPOTHESIS")
print("=" * 100)

print("\nActual test results for all 180 conditions:")
print("  DGr absolute error: min=4.63, max=34.22, avg=15.36 J/mol")
print("  Variation: 4.63 to 34.22 J/mol = 7.4x range")

print("\n>>> ERRORS ARE NOT CONSTANT! <<<")
print("    Varies from 4.63 to 34.22 J/mol across different T,P")
print("    Expected constant ~= -14 J/mol if individual species errors didn't vary")

print("\n" + "=" * 100)
print("WHAT THIS MEANS")
print("=" * 100)

print("""
CONCLUSION: Per-species errors MUST VARY with Temperature and Pressure

Evidence:
  - At T=650C, P=15 kb: DGr_error ~= -14 J/mol
  - Full test range: DGr_error varies 4.63 to 34.22 J/mol
  - If H2O, CO2, HCO3-, H+ errors were truly constant:
    -> DGr_error = (-9 + 0) - (5 + 0) = -14 J/mol (would be SAME everywhere)
  - But observed: varies 7.4x
  - Therefore: At least one (or more) of the per-species errors varies with T,P

LIKELY SCENARIOS:

1. CO2 error varies with T and P
   - Model accuracy depends on water density/properties
   - Worse at extreme T or P

2. HCO3- error varies with T and P
   - Ion behavior is temperature/pressure dependent
   - Model assumptions may break down at extremes

3. H2O error stays perfect (0 J/mol everywhere)
   - Most likely given Zhang-Duan is well-validated
   - Water integration already correct

4. H+ error stays 0 (reference species)
   - Fixed by convention, not model-dependent
""")

# Analyze by temperature
temps = sorted(set(T_C))
pressures = sorted(set(P_kb))

print("\n" + "=" * 100)
print("ANALYSIS BY TEMPERATURE (T dependency)")
print("=" * 100)

print(f"\nDGr truth values by temperature:")
print(f"{'T (C)':<10} {'Count':<8} {'DGr min':<15} {'DGr max':<15} {'Range':<15}")
print("-" * 60)

for t in temps:
    mask = T_C == t
    gr_vals = g_rxn_truth[mask]
    print(
        f"{t:<10.0f} {np.sum(mask):<8} {gr_vals.min():<15.0f} {gr_vals.max():<15.0f} {gr_vals.max() - gr_vals.min():<15.0f}"
    )

# Analyze by pressure
print("\n" + "=" * 100)
print("ANALYSIS BY PRESSURE (P dependency)")
print("=" * 100)

print(f"\nDGr truth values by pressure:")
print(f"{'P (kb)':<10} {'Count':<8} {'DGr min':<15} {'DGr max':<15} {'Range':<15}")
print("-" * 60)

for p in pressures:
    mask = P_kb == p
    gr_vals = g_rxn_truth[mask]
    print(
        f"{p:<10.0f} {np.sum(mask):<8} {gr_vals.min():<15.0f} {gr_vals.max():<15.0f} {gr_vals.max() - gr_vals.min():<15.0f}"
    )

print("\n" + "=" * 100)
print("KEY INSIGHT")
print("=" * 100)

print("""
The Excel TRUTH data itself shows large variations in ΔGr with T,P.
This is EXPECTED - thermodynamic properties change dramatically with conditions.

What MATTERS for error consistency:
  - Are model DEVIATIONS from truth consistent?
  - Or do individual species models degrade at certain T,P?

OBSERVED FACT:
  - One test point (T=650, P=15): error ~= 14.46 J/mol
  - Across full range: errors 4.63 to 34.22 J/mol
  - This 7.4x variation in error magnitude indicates:

    EITHER the models have TEMPERATURE/PRESSURE dependent biases
    OR the models become LESS ACCURATE at certain extremes
    OR there's ERROR COMPENSATION (sometimes errors add, sometimes subtract)

NEXT DIAGNOSTIC STEP:
  - Need per-species errors for multiple (T,P) points, not just one
  - Would require enabling DEBUG output for representative subset:
    * Low T (300°C) vs High T (1000°C)
    * Low P (5 kb) vs High P (60 kb)
  - Then compare if CO2 and HCO3- errors vary systematically
  - Check if errors are RANDOM (noise) or SYSTEMATIC (bias)
""")

print("\n" + "=" * 100)
print("CURRENT DATA LIMITATION")
print("=" * 100)

print("""
Available: 1 debug point (T=650°C, P=15 kb)
  - H2O: 0 J/mol error
  - CO2: +5 J/mol error
  - HCO3-: -9 J/mol error

Cannot determine from single point:
  - Temperature sensitivity of errors
  - Pressure sensitivity of errors
  - Whether errors are systematic bias or random scatter
  - Which T,P conditions have worst accuracy

To answer your question fully:
  Would need to modify test to capture DEBUG output for:
  - Low temperature point (e.g., T=300°C, P=5 kb)
  - Medium temperature point (e.g., T=650°C, P=15 kb) [HAVE THIS]
  - High temperature point (e.g., T=1000°C, P=60 kb)
  - Low pressure point (e.g., T=500°C, P=5 kb)
  - High pressure point (e.g., T=500°C, P=60 kb)

Then compare: Are errors consistent across these points?
""")
