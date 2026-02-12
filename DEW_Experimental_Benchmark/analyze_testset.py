"""Analyze distinct T-P conditions in quartz solubility dataset"""

import pandas as pd

df = pd.read_csv("quartz_DEW_testset.csv")

print("=" * 80)
print("QUARTZ SOLUBILITY TEST SET ANALYSIS")
print("=" * 80)
print(f"\nTotal data points: {len(df)}")
print(f"Distinct temperatures: {df['T_C'].nunique()}")
print(f"Distinct pressures (kbar): {df['P_kbar'].nunique()}")
print(f"\nTemperature range: {df['T_C'].min():.1f} - {df['T_C'].max():.1f} °C")
print(f"Pressure range: {df['P_kbar'].min():.3f} - {df['P_kbar'].max():.3f} kbar")
print(f"\nUnique T-P combinations: {df[['T_C', 'P_kbar']].drop_duplicates().shape[0]}")

print("\n" + "=" * 80)
print("DISTINCT TEMPERATURES (°C):")
print("=" * 80)
temps = sorted(df["T_C"].unique())
for i in range(0, len(temps), 10):
    print(", ".join([f"{t:.1f}" for t in temps[i : i + 10]]))

print("\n" + "=" * 80)
print("DISTINCT PRESSURES (kbar):")
print("=" * 80)
pressures = sorted(df["P_kbar"].unique())
for i in range(0, len(pressures), 10):
    print(", ".join([f"{p:.3f}" for p in pressures[i : i + 10]]))

# Analyze experiments
print("\n" + "=" * 80)
print("EXPERIMENTS:")
print("=" * 80)
if "experiment_id" in df.columns:
    experiments = df.groupby("experiment_id").agg(
        {"T_C": ["min", "max", "count"], "P_kbar": ["min", "max"]}
    )
    print(experiments)
