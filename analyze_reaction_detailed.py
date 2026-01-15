import pandas as pd
import numpy as np

# Load comprehensive diffs
df = pd.read_csv("reaction_column_diffs.csv")

print("=" * 80)
print("COMPREHENSIVE ERROR ANALYSIS - ALL 180 TEST CONDITIONS")
print("=" * 80)

# Per-species G0 errors
print("\n=== GIBBS ENERGY ERRORS BY SPECIES (cal/mol) ===")
species_g = {
    "H2O": df["dG0_H2O_cal"],
    "CO2": df["dG0_CO2_cal"],
    "H+": df["dG0_H+_cal"],
    "HCO3-": df["dG0_HCO3-_cal"],
}

for sp, errors in species_g.items():
    print(f"\n{sp}:")
    print(f"  Min error:  {errors.min():>12.6f} cal/mol")
    print(f"  Max error:  {errors.max():>12.6f} cal/mol")
    print(f"  Mean error: {errors.mean():>12.6f} cal/mol")
    print(f"  Std dev:    {errors.std():>12.6f} cal/mol")

# Per-species V0 errors
print("\n=== VOLUME ERRORS BY SPECIES (cm³/mol) ===")
species_v = {
    "H2O": df["dV0_H2O_cm3"],
    "CO2": df["dV0_CO2_cm3"],
    "H+": df["dV0_H+_cm3"],
    "HCO3-": df["dV0_HCO3-_cm3"],
}

for sp, errors in species_v.items():
    print(f"\n{sp}:")
    print(f"  Min error:  {errors.min():>14.8f} cm³/mol")
    print(f"  Max error:  {errors.max():>14.8f} cm³/mol")
    print(f"  Mean error: {errors.mean():>14.8f} cm³/mol")
    print(f"  Std dev:    {errors.std():>14.8f} cm³/mol")

# Error distribution by T and P
print("\n=== ERROR DISTRIBUTION BY TEMPERATURE ===")
g_by_t = df.groupby("T_C")["dDeltaGro_cal"].agg(["mean", "std", "min", "max"])
print(g_by_t.to_string())

print("\n=== ERROR DISTRIBUTION BY PRESSURE ===")
g_by_p = df.groupby("P_kb")["dDeltaGro_cal"].agg(["mean", "std", "min", "max"])
print(g_by_p.to_string())

print("\n=== TOP 10 WORST GIBBS ERRORS (by |ΔΔGr|) ===")
df["abs_dG"] = df["dDeltaGro_cal"].abs()
worst = df.nlargest(10, "abs_dG")[
    ["T_C", "P_kb", "dDeltaGro_cal", "dG0_H2O_cal", "dG0_CO2_cal", "dG0_HCO3-_cal"]
]
print(worst.to_string(index=False))

print("\n=== SPECIES CONTRIBUTION TO TOTAL ERROR ===")
df["H2O_contrib"] = df["dG0_H2O_cal"].abs()
df["CO2_contrib"] = df["dG0_CO2_cal"].abs()
df["HCO3_contrib"] = df["dG0_HCO3-_cal"].abs()
print(f"Average absolute contribution to ΔΔGr:")
print(f"  H2O:   {df['H2O_contrib'].mean():>8.4f} cal/mol")
print(f"  CO2:   {df['CO2_contrib'].mean():>8.4f} cal/mol")
print(f"  H+:    {df['dG0_H+_cal'].abs().mean():>8.4f} cal/mol")
print(f"  HCO3-: {df['HCO3_contrib'].mean():>8.4f} cal/mol")

# Correlation between species errors
print("\n=== CORRELATION BETWEEN SPECIES ERRORS ===")
corr_cols = ["dG0_H2O_cal", "dG0_CO2_cal", "dG0_HCO3-_cal"]
corr = df[corr_cols].corr()
print(corr.to_string())

# Pressure-Temperature heatmap data
print("\n=== WORST ERROR LOCATIONS (T, P) ===")
worst_loc = df.nlargest(20, "abs_dG")[["T_C", "P_kb", "dDeltaGro_cal"]]
print(worst_loc.to_string(index=False))
