"""Identify which data points lie on or near water saturation pressure curve"""

import pandas as pd
import numpy as np


# Simple Psat approximation (Antoine equation for water, valid ~0-100°C)
# For better accuracy above 100°C, we'd need IAPWS formulation
def psat_bar_simple(T_C):
    """Calculate saturation pressure in bar (simple Antoine equation)"""
    if T_C < 0 or T_C > 374:  # Beyond critical point
        return np.nan
    T_K = T_C + 273.15
    # Antoine parameters for water (bar units)
    A, B, C = 5.40221, 1838.675, -31.737
    log10_P = A - B / (T_K + C)
    return 10**log10_P


def psat_kbar_simple(T_C):
    """Calculate saturation pressure in kbar"""
    P_bar = psat_bar_simple(T_C)
    return P_bar / 1000.0 if not np.isnan(P_bar) else np.nan


# Read data
df = pd.read_csv("quartz_DEW_testset.csv")

# Calculate Psat for each temperature
df["Psat_kbar"] = df["T_C"].apply(psat_kbar_simple)

# Check for NaN pressures
nan_pressure = df[df["P_kbar"].isna()]

# Calculate relative difference from Psat
df["P_Psat_ratio"] = df["P_kbar"] / df["Psat_kbar"]
df["P_minus_Psat"] = df["P_kbar"] - df["Psat_kbar"]

# Identify points near Psat (within ±5% or ±0.01 kbar, whichever is larger)
tolerance_rel = 0.05
tolerance_abs = 0.01  # kbar

df["near_Psat"] = (np.abs(df["P_Psat_ratio"] - 1.0) < tolerance_rel) | (
    np.abs(df["P_minus_Psat"]) < tolerance_abs
)

near_psat = df[df["near_Psat"] == True].copy()

print("=" * 80)
print("SATURATION PRESSURE ANALYSIS")
print("=" * 80)

print(f"\nTotal data points: {len(df)}")
print(f"Points with NaN pressure: {len(nan_pressure)}")
print(f"Points near Psat (±5% or ±0.01 kbar): {len(near_psat)}")

if len(nan_pressure) > 0:
    print("\n" + "=" * 80)
    print("DATA POINTS WITH NaN PRESSURE (likely Psat conditions):")
    print("=" * 80)
    print(
        nan_pressure[
            ["reference", "T_C", "P_kbar", "molality_m", "experiment_type"]
        ].to_string(index=False)
    )

if len(near_psat) > 0:
    print("\n" + "=" * 80)
    print("DATA POINTS NEAR SATURATION PRESSURE:")
    print("=" * 80)
    near_psat_sorted = near_psat.sort_values("T_C")
    print(
        near_psat_sorted[
            [
                "reference",
                "T_C",
                "P_kbar",
                "Psat_kbar",
                "P_minus_Psat",
                "molality_m",
            ]
        ].to_string(index=False)
    )

# Low pressure points that might be Psat
low_p = df[df["P_kbar"] < 0.05].copy()
low_p = low_p.sort_values("T_C")

print("\n" + "=" * 80)
print(f"LOW PRESSURE POINTS (P < 0.05 kbar = 50 bar):")
print("=" * 80)
print(f"Total low-pressure points: {len(low_p)}")
if len(low_p) > 0:
    print("\nThese are likely at or near saturation pressure:")
    print(
        low_p[
            [
                "reference",
                "T_C",
                "P_kbar",
                "Psat_kbar",
                "P_Psat_ratio",
                "molality_m",
            ]
        ].to_string(index=False)
    )

# Check critical point region (374.15°C)
critical = df[(df["T_C"] > 370) & (df["T_C"] < 380)]
if len(critical) > 0:
    print("\n" + "=" * 80)
    print("POINTS NEAR CRITICAL POINT (370-380°C):")
    print("=" * 80)
    print(
        critical[["reference", "T_C", "P_kbar", "Psat_kbar", "molality_m"]]
        .sort_values("T_C")
        .to_string(index=False)
    )

print("\n" + "=" * 80)
print("PRESSURE DISTRIBUTION SUMMARY:")
print("=" * 80)
print(
    f"Sub-critical pressure (P < Psat): {len(df[df['P_kbar'] < df['Psat_kbar']])} points"
)
print(f"At saturation (±5%): {len(near_psat)} points")
print(
    f"Super-critical pressure (P > Psat): {len(df[df['P_kbar'] > df['Psat_kbar']])} points"
)
print(f"Unknown (NaN pressure): {len(nan_pressure)} points")
