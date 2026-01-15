#!/usr/bin/env python3
"""
NEXT INVESTIGATION: Does C++ YAML parsing lose precision?

After verifying Python extraction is EXACT, need to check if C++ loses precision
when reading the YAML file back into Reaktoro.
"""

import yaml
import json

print("=" * 100)
print("YAML PRECISION ANALYSIS: Python vs C++")
print("=" * 100)

# Load the YAML file
with open("dew2024-aqueous.yaml", "r") as f:
    data = yaml.safe_load(f)

# Check specific values
print("\n" + "=" * 100)
print("CHECKING YAML VALUES FOR FULL PRECISION")
print("=" * 100)

species_to_check = ["CO2(0)", "HCO3(-)"]

for species_key in species_to_check:
    if species_key not in data["Species"]:
        print(f"\n⚠ {species_key} not found")
        continue

    species = data["Species"][species_key]
    hkf = species["StandardThermoModel"]["HKF"]

    print(f"\n{species_key}:")
    print("-" * 100)

    params = {
        "Gf": hkf["Gf"],
        "Hf": hkf["Hf"],
        "Sr": hkf["Sr"],
        "a1": hkf["a1"],
        "a2": hkf["a2"],
        "a3": hkf["a3"],
        "a4": hkf["a4"],
        "c1": hkf["c1"],
        "c2": hkf["c2"],
        "wref": hkf["wref"],
    }

    for param, value in params.items():
        if isinstance(value, float):
            # Show the exact representation
            print(f"  {param:6s}: {value:20.15g}  (repr: {repr(value)})")
        else:
            print(f"  {param:6s}: {value}")

print("\n" + "=" * 100)
print("RAW YAML TEXT: Check for truncation")
print("=" * 100)

with open("dew2024-aqueous.yaml", "r") as f:
    lines = f.readlines()

# Find CO2 section
in_co2 = False
in_hco3 = False
print("\nCO2(0) section from YAML file:")
for i, line in enumerate(lines):
    if "CO2(0):" in line:
        in_co2 = True
        start = i
    elif in_co2 and ("wref:" in line):
        # Print around wref
        for j in range(max(0, i - 2), min(len(lines), i + 5)):
            print(f"  {lines[j]}", end="")
        in_co2 = False
        break

print("\nHCO3(-) section from YAML file:")
for i, line in enumerate(lines):
    if "HCO3(-):" in line:
        in_hco3 = True
    elif in_hco3 and ("wref:" in line):
        for j in range(max(0, i - 1), min(len(lines), i + 3)):
            print(f"  {lines[j]}", end="")
        in_hco3 = False
        break

print("\n" + "=" * 100)
print("WHAT C++ NEEDS TO DO")
print("=" * 100)

print("""
When C++ reads the YAML file (using a YAML parser like nlohmann/yaml or similar):

1. YAML TEXT EXAMPLE:
   wref: -334720.0
   wref: 532748.7200000001

2. C++ PARSING OPTIONS:

   Option A: Parse as string, convert to double
   - Result: Full precision (matches Python)
   - Precision: ~15-17 significant digits

   Option B: Parse as YAML number type
   - Result: Full precision (depends on parser)
   - Precision: ~15-17 significant digits

   Option C: Fixed parsing (e.g., only 6 decimals)
   - Result: LOSS of precision ✗
   - Example: 532748.7200000001 → 532748.720000

3. REAKTORO'S YAML PARSER:
   File: cpp/Reaktoro/Core/Database.cpp or similar
   Need to verify:
   - Uses full double precision when reading floats
   - Doesn't truncate or round values
   - Preserves all digits from YAML file

4. TESTING HYPOTHESIS:
   If C++ reads YAML correctly:
   - CO2 wref should be: -334720.0 exactly
   - HCO3- wref should be: 532748.7200000001 exactly

   If we see different values in debug output:
   - C++ parser might be truncating
   - Or rounding to N decimal places
""")

print("\n" + "=" * 100)
print("RECOMMENDATION FOR NEXT STEP")
print("=" * 100)

print("""
To verify C++ is reading YAML values exactly:

1. Add debug output to C++ code that reads HKF parameters
   Location: Reaktoro/Core/Database.cpp or similar
   Code: Print each parameter value read from YAML to 15+ decimal places

   Example:
   std::cout << std::setprecision(15) << "wref = " << wref << std::endl;

2. Compare C++ output with YAML values:
   YAML:  wref: 532748.7200000001
   C++:   wref = 532748.7200000001  ← Should match exactly

3. If they differ, identify:
   - Which YAML parser is being used?
   - Does it have precision limits?
   - Are there any custom conversion functions?

4. Files to check:
   - cpp/Reaktoro/Core/Database.cpp
   - cpp/Reaktoro/Core/StandardThermoModel*.cpp
   - cpp/Reaktoro/Core/SpeciesData.hpp
   - Any YAML parser integration files
""")

print("\n" + "=" * 100)
