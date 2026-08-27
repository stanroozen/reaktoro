#!/usr/bin/env python3
"""
Analyze gas species integration in DEW databases.
"""

import json
from pathlib import Path

db_dir = Path(__file__).parent.parent / "embedded" / "databases" / "perplex"

print("=" * 100)
print("GAS SPECIES INTEGRATION IN DEW DATABASES")
print("=" * 100)
print()

dbs_to_check = [
    "DEW24HP633ver_elements-reaktoro.json",
    "DEW17HP622_Zn_2025-reaktoro.json",
    "DEW19HP622ver_elements-reaktoro.json",
]

for db_name in dbs_to_check:
    db_file = db_dir / db_name
    if not db_file.exists():
        continue

    with open(db_file) as f:
        data = json.load(f)

    print(f"Database: {db_name}")
    print("-" * 100)

    gases = {}
    for sp_name, sp in data["Species"].items():
        if sp.get("AggregateState") == "Gas":
            gases[sp_name] = sp

    print(f"Total gas species: {len(gases)}\n")

    if not gases:
        print("  ⚠️  No gas species found\n\n")
        continue

    for gas_name, gas_data in sorted(gases.items()):
        print(f"Species: {gas_name}")
        print(f"  Formula: {gas_data.get('Formula')}")
        print(f"  Charge: {gas_data.get('Charge', 0)}")

        # Check thermo model
        if "StandardThermoModel" in gas_data:
            models = list(gas_data["StandardThermoModel"].keys())
            print(f"  Models: {models}")

            for model in models:
                model_data = gas_data["StandardThermoModel"][model]
                if model == "PerplexGFSM":
                    print(f"    PerplexGFSM:")
                    print(f"      - speciesIndex: {model_data.get('speciesIndex')}")
                    print(f"      - G0: {model_data.get('G0')} J/mol")
                    print(f"      - H0: {model_data.get('H0')} J/mol")
                    print(f"      - V0: {model_data.get('V0')} m³/mol")
                    print(f"      - Tmax: {model_data.get('Tmax')} K")
        else:
            print(f"  ⚠️  No StandardThermoModel")

        # Check metadata
        if "Metadata" in gas_data:
            meta = gas_data["Metadata"]
            print(f"  EoS from PerpleX: {meta.get('PerpleX_EoS')}")
            print(f"  Conversion mode: {meta.get('GfConversionMode')}")

        print()

    print()

print("=" * 100)
print("ANALYSIS OF GAS INTEGRATION")
print("=" * 100)
print()

# Now check Reaktoro's support for PerplexGFSM
reaktoro_dir = Path(__file__).parent.parent.parent / "Reaktoro"

gfsm_support_files = [
    "Extensions/Perple_X/StandardThermoModelPerplexGFSM.cpp",
    "Extensions/Perple_X/StandardThermoModelPerplexGFSM.hpp",
    "Models/StandardThermoModels.hpp",
]

print("PerplexGFSM Model Support in Reaktoro:")
print("-" * 100)

for fname in gfsm_support_files:
    fpath = reaktoro_dir / fname
    if fpath.exists():
        print(f"✅ {fname}")
    else:
        print(f"❌ {fname} (not found)")

print()
print("Summary:")
print("  - PerplexGFSM is the gas phase model used in DEW databases")
print(
    "  - Gases are marked with EoS codes (e.g., 101, 102, 103, etc. for different fluids)"
)
print("  - Each gas has its own speciesIndex (unique identifier)")
print("  - G0, H0, V0 thermodynamic properties are provided")
print()
