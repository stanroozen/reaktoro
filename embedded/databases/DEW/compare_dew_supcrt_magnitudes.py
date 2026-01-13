#!/usr/bin/env python3
"""Compare parameter magnitudes between DEW and HKF/SUPCRT databases"""

import yaml

print("=" * 80)
print("PARAMETER MAGNITUDE COMPARISON: DEW vs SUPCRT for CO2(aq)")
print("=" * 80)

# Load DEW database
with open("dew2024-aqueous.yaml", "r") as f:
    dew_data = yaml.safe_load(f)

# Load SUPCRT database
supcrt_path = "../reaktoro/supcrt98.yaml"
with open(supcrt_path, "r") as f:
    supcrt_data = yaml.safe_load(f)

# Find CO2 in both databases
dew_co2 = dew_data["Species"]["CO2(0)"]["StandardThermoModel"]["HKF"]

# Find CO2 in SUPCRT (might be named differently)
supcrt_co2 = None
for name, species in supcrt_data["Species"].items():
    if "CO2" in name and "StandardThermoModel" in species:
        if "HKF" in species["StandardThermoModel"]:
            supcrt_co2 = species["StandardThermoModel"]["HKF"]
            supcrt_name = name
            break

if not supcrt_co2:
    print("ERROR: Could not find CO2 in SUPCRT database")
    exit(1)

print(f"\nDEW database: CO2(0)")
print(f"SUPCRT database: {supcrt_name}")
print("\n" + "-" * 80)

params = ["a1", "a2", "a3", "a4", "c1", "c2", "wref"]
param_names = {
    "a1": "a₁",
    "a2": "a₂",
    "a3": "a₃",
    "a4": "a₄",
    "c1": "c₁",
    "c2": "c₂",
    "wref": "ω",
}

print(
    f"\n{'Parameter':<12} {'DEW Value':<20} {'SUPCRT Value':<20} {'Ratio (SUPCRT/DEW)':<20}"
)
print("-" * 80)

for param in params:
    dew_val = dew_co2.get(param, 0.0)
    supcrt_val = supcrt_co2.get(param, 0.0)

    if dew_val != 0 and supcrt_val != 0:
        ratio = supcrt_val / dew_val
    elif dew_val == 0 and supcrt_val == 0:
        ratio = 1.0
    else:
        ratio = float("inf") if dew_val == 0 else 0.0

    name = param_names.get(param, param)
    print(f"{name:<12} {dew_val:<20.6e} {supcrt_val:<20.6e} {ratio:<20.2f}")

print("\n" + "=" * 80)
print("OBSERVATIONS:")
print("=" * 80)

# Calculate key ratios
if dew_co2["a1"] != 0:
    a1_ratio = supcrt_co2["a1"] / dew_co2["a1"]
    print(f"\na₁ ratio: {a1_ratio:.1f}×")
    print(f"  DEW:    {dew_co2['a1']:.6e} J/(mol·Pa)")
    print(f"  SUPCRT: {supcrt_co2['a1']:.6e} J/(mol·Pa)")

if dew_co2["wref"] != 0:
    w_ratio = supcrt_co2["wref"] / dew_co2["wref"]
    print(f"\nω ratio: {w_ratio:.1f}×")
    print(f"  DEW:    {dew_co2['wref']:.6e} J/mol")
    print(f"  SUPCRT: {supcrt_co2['wref']:.6e} J/mol")
    print(f"\n  DEW with ×1e10:    {dew_co2['wref'] * 1e10:.1f} J/mol")
    print(f"  SUPCRT (no scale): {supcrt_co2['wref']:.1f} J/mol")
    print(f"  → DEW needs ×1e10 to match SUPCRT magnitude!")

print("\n" + "=" * 80)
