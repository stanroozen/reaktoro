#!/usr/bin/env python3
"""
Direct CSV analysis: Per-species G0 error distribution
Reaktoro DEW model vs Excel truth data
"""

import pandas as pd
import numpy as np

csv_file = "Reaktoro/Extensions/DEW/tests/reactionTesttruth.csv"

print("Loading CSV file...")
df = pd.read_csv(csv_file)

print(f"Shape: {df.shape}")
print(f"\nColumns: {df.columns.tolist()}\n")

# Map columns based on the test code comments
# 0: Pressure_kb, 1: Temp_C
# 4: DeltaGo_H2O_calmol-1
# 5: DeltaGo_CO2_aq_calmol-1
# 6: DeltaGo_H+_calmol-1
# 7: DeltaGo_HCO3-_calmol-1
# 8: DeltaGro_calmol-1
# 9: log_K

# Check actual column names
print("First few rows:")
print(df.head(3))

# The Excel truth values are in the CSV
# Extract G0 values for each species (in cal/mol)
CAL_TO_J = 4.184

if "DeltaGo_H2O_calmol-1" in df.columns:
    g_h2o_truth = df["DeltaGo_H2O_calmol-1"].values * CAL_TO_J
    g_co2_truth = df["DeltaGo_CO2_aq_calmol-1"].values * CAL_TO_J
    g_h_truth = df["DeltaGo_H+_calmol-1"].values * CAL_TO_J
    g_hco3_truth = df["DeltaGo_HCO3-_calmol-1"].values * CAL_TO_J
    g_rxn_truth = df["DeltaGro_calmol-1"].values * CAL_TO_J

    T_C = df["Temp_C"].values
    P_kb = df["Pressure_kb"].values

    print("\n" + "=" * 80)
    print("PER-SPECIES G0 TRUTH DATA FROM EXCEL")
    print("=" * 80)
    print("\nSpecies errors range (from Excel, these are the truth values):")
    print(
        f"  H2O:    min={g_h2o_truth.min():.0f}, max={g_h2o_truth.max():.0f}, mean={g_h2o_truth.mean():.0f} J/mol"
    )
    print(
        f"  CO2:    min={g_co2_truth.min():.0f}, max={g_co2_truth.max():.0f}, mean={g_co2_truth.mean():.0f} J/mol"
    )
    print(
        f"  H+:     min={g_h_truth.min():.0f}, max={g_h_truth.max():.0f}, mean={g_h_truth.mean():.0f} J/mol"
    )
    print(
        f"  HCO3-:  min={g_hco3_truth.min():.0f}, max={g_hco3_truth.max():.0f}, mean={g_hco3_truth.mean():.0f} J/mol"
    )
    print(
        f"  ΔGr:    min={g_rxn_truth.min():.0f}, max={g_rxn_truth.max():.0f}, mean={g_rxn_truth.mean():.0f} J/mol"
    )

    print("\n" + "=" * 80)
    print("LARGEST GIBBS CHANGES BY SPECIES (Excel Truth)")
    print("=" * 80)

    # Find largest absolute values and variations
    specs = [
        ("H2O", g_h2o_truth, "cal/mol" in str(df.columns)),
        ("CO2(aq)", g_co2_truth, "cal/mol" in str(df.columns)),
        ("H+", g_h_truth, "cal/mol" in str(df.columns)),
        ("HCO3-", g_hco3_truth, "cal/mol" in str(df.columns)),
    ]

    for name, data, _ in specs:
        max_idx = np.argmax(np.abs(data))
        print(f"\n{name}:")
        print(
            f"  Largest |value|: {np.abs(data[max_idx]):.2f} J/mol at T={T_C[max_idx]:.0f}C, P={P_kb[max_idx]:.0f}kb"
        )
        print(
            f"  Range: {data.min():.0f} to {data.max():.0f} J/mol (span: {data.max() - data.min():.0f})"
        )
        print(f"  Variation: {data.std():.0f} J/mol (std dev)")

    print("\n" + "=" * 80)
    print("REACTION GIBBS COMPOSITION")
    print("=" * 80)
    print("\nReaction: H2O + CO2(aq) = H+ + HCO3-")
    print("ΔGr = G(H+) + G(HCO3-) - G(H2O) - G(CO2)\n")

    # Show first few points
    print("Sample points (first 10):")
    print(
        f"{'T (C)':<8} {'P (kb)':<8} {'H2O':<12} {'CO2':<12} {'H+':<12} {'HCO3-':<12} {'ΔGr':<12}"
    )
    print("-" * 80)
    for i in range(min(10, len(df))):
        calc = g_h_truth[i] + g_hco3_truth[i] - g_h2o_truth[i] - g_co2_truth[i]
        print(
            f"{T_C[i]:<8.0f} {P_kb[i]:<8.0f} {g_h2o_truth[i]:<12.0f} {g_co2_truth[i]:<12.0f} {g_h_truth[i]:<12.0f} {g_hco3_truth[i]:<12.0f} {g_rxn_truth[i]:<12.0f}"
        )

    # Identify which species varies most
    print("\n" + "=" * 80)
    print("SPECIES CONTRIBUTION TO REACTION VARIATION")
    print("=" * 80)

    h2o_range = g_h2o_truth.max() - g_h2o_truth.min()
    co2_range = g_co2_truth.max() - g_co2_truth.min()
    h_range = g_h_truth.max() - g_h_truth.min()
    hco3_range = g_hco3_truth.max() - g_hco3_truth.min()

    print(f"\nG0 variation ranges:")
    print(f"  H2O:    {h2o_range:>8.0f} J/mol  (contributes with - sign to ΔGr)")
    print(f"  CO2:    {co2_range:>8.0f} J/mol  (contributes with - sign to ΔGr)")
    print(f"  H+:     {h_range:>8.0f} J/mol  (contributes with + sign to ΔGr)")
    print(f"  HCO3-:  {hco3_range:>8.0f} J/mol  (contributes with + sign to ΔGr)")

    total_variation = h_range + hco3_range + h2o_range + co2_range
    print(f"\nTotal variation in reactants: {total_variation:.0f} J/mol")
    print(f"Actual ΔGr variation: {g_rxn_truth.max() - g_rxn_truth.min():.0f} J/mol")

    print("\nPer-species % contribution to variation:")
    for name, rng in [
        ("H2O", h2o_range),
        ("CO2", co2_range),
        ("H+", h_range),
        ("HCO3-", hco3_range),
    ]:
        pct = 100 * rng / total_variation if total_variation > 0 else 0
        print(f"  {name:<10}: {pct:>6.1f}%")

    # Find conditions with extreme ΔGr
    print("\n" + "=" * 80)
    print("EXTREME ΔGr CONDITIONS")
    print("=" * 80)

    max_gr_idx = np.argmax(g_rxn_truth)
    min_gr_idx = np.argmin(g_rxn_truth)

    print(f"\nMaximum ΔGr: {g_rxn_truth[max_gr_idx]:.0f} J/mol")
    print(f"  at T={T_C[max_gr_idx]:.0f}°C, P={P_kb[max_gr_idx]:.0f}kb")
    print(
        f"  H2O={g_h2o_truth[max_gr_idx]:.0f}, CO2={g_co2_truth[max_gr_idx]:.0f}, H+={g_h_truth[max_gr_idx]:.0f}, HCO3={g_hco3_truth[max_gr_idx]:.0f}"
    )

    print(f"\nMinimum ΔGr: {g_rxn_truth[min_gr_idx]:.0f} J/mol")
    print(f"  at T={T_C[min_gr_idx]:.0f}°C, P={P_kb[min_gr_idx]:.0f}kb")
    print(
        f"  H2O={g_h2o_truth[min_gr_idx]:.0f}, CO2={g_co2_truth[min_gr_idx]:.0f}, H+={g_h_truth[min_gr_idx]:.0f}, HCO3={g_hco3_truth[min_gr_idx]:.0f}"
    )

else:
    print("[!] Column names don't match expected structure")
    print("Available columns:", df.columns.tolist())
