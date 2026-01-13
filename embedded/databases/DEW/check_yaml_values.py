#!/usr/bin/env python3
import yaml

with open("dew2024-aqueous.yaml") as f:
    data = yaml.safe_load(f)

# Find CO2(0) species
for species in data:
    if species.get("Name") == "CO2(0)":
        print("CO2(0) in YAML after fix:")
        print(f"  a1   = {species['HKF']['a1']:.6e} J/(mol·Pa)")
        print(f"  a2   = {species['HKF']['a2']:.6e} J/mol")
        print(f"  wref = {species['HKF']['wref']:.6e} J/mol")
        break
else:
    print("CO2(0) not found in YAML")
