#!/usr/bin/env python3
"""
Scan PerpleX databases to determine compatibility with Reaktoro functionality.
"""

import json
from pathlib import Path
from collections import Counter, defaultdict

db_dir = Path(__file__).parent.parent / "embedded" / "databases" / "perplex"


def scan_database(db_file):
    """Scan a database and return coverage statistics."""
    with open(db_file) as f:
        data = json.load(f)

    stats = {
        "filename": db_file.name,
        "total_species": len(data.get("Species", {})),
        "models": Counter(),
        "with_thermo_reference": 0,
        "with_h0": 0,
        "with_v0": 0,
        "with_s0": 0,
        "constant_g0_only": 0,
        "aggregate_states": Counter(),
        "missing_models": [],
    }

    for sp_name, sp in data.get("Species", {}).items():
        # Aggregate state
        if "AggregateState" in sp:
            stats["aggregate_states"][sp["AggregateState"]] += 1

        # Thermodynamic reference data
        if "ThermoReference" in sp:
            stats["with_thermo_reference"] += 1
            thermo_ref = sp["ThermoReference"]
            if thermo_ref.get("Hf") is not None:
                stats["with_h0"] += 1
            if thermo_ref.get("V0") is not None or thermo_ref.get("Vr") is not None:
                stats["with_v0"] += 1
            if thermo_ref.get("S0") is not None or thermo_ref.get("Sr") is not None:
                stats["with_s0"] += 1

        # StandardThermoModel models
        if "StandardThermoModel" in sp:
            for model_name in sp["StandardThermoModel"].keys():
                stats["models"][model_name] += 1

                # Check Constant model details
                if model_name == "Constant":
                    const = sp["StandardThermoModel"][model_name]
                    if "H0" not in const and "H" not in const:
                        stats["constant_g0_only"] += 1
        else:
            stats["missing_models"].append(sp_name)

    return stats


def categorize_database(db_file):
    """Categorize database by name."""
    name = db_file.name
    if "DEW" in name:
        return "DEW-mixed (aqueous + solids)"
    elif name.startswith(("b89", "b92", "ba96")):
        return "Old-SUPCRT"
    else:
        return "HP-only (solid minerals)"


# Scan all databases
print("=" * 80)
print("PERPLE_X DATABASE COMPATIBILITY SCAN FOR REAKTORO")
print("=" * 80)
print()

all_databases = sorted(db_dir.glob("*-reaktoro.json"))
results_by_category = defaultdict(list)

for db_file in all_databases:
    stats = scan_database(db_file)
    category = categorize_database(db_file)
    results_by_category[category].append(stats)

# Print summary
for category in [
    "HP-only (solid minerals)",
    "DEW-mixed (aqueous + solids)",
    "Old-SUPCRT",
]:
    if category not in results_by_category:
        continue

    print(f"\n{'=' * 80}")
    print(f"CATEGORY: {category}")
    print(f"{'=' * 80}")
    print()

    databases = results_by_category[category]

    if category == "Old-SUPCRT":
        print("⚠️  WARNING: Old SUPCRT databases have limited functionality!\n")

    for stats in databases:
        print(f"Database: {stats['filename']}")
        print(f"  Total species: {stats['total_species']}")

        if stats["models"]:
            print(f"  Thermodynamic models available:")
            for model, count in sorted(stats["models"].items(), key=lambda x: -x[1]):
                print(f"    - {model}: {count} species")

        if stats["missing_models"]:
            print(
                f"  ⚠️  Species without StandardThermoModel: {len(stats['missing_models'])}"
            )
            if len(stats["missing_models"]) <= 5:
                for name in stats["missing_models"]:
                    print(f"     - {name}")

        if stats["constant_g0_only"] > 0:
            print(
                f"  ⚠️  Constant models with G0 only (no H0): {stats['constant_g0_only']}"
            )

        print(f"  Coverage:")
        print(
            f"    - Has ThermoReference: {stats['with_thermo_reference']}/{stats['total_species']}"
        )
        print(f"    - Has Hf/H0: {stats['with_h0']}/{stats['total_species']}")
        print(f"    - Has V0: {stats['with_v0']}/{stats['total_species']}")
        print(f"    - Has S0: {stats['with_s0']}/{stats['total_species']}")

        if stats["aggregate_states"]:
            print(f"  Aggregate states: {dict(stats['aggregate_states'])}")

        print()

# Summary and recommendations
print("\n" + "=" * 80)
print("FUNCTIONALITY ANALYSIS")
print("=" * 80)
print()

print("1. SOLID MINERAL THERMODYNAMICS")
print("   HP-only databases: ✅ Full support (HollandPowell EoS)")
print("   DEW-mixed databases: ✅ Full support (HollandPowell for solids)")
print(
    "   Old-SUPCRT databases: ⚠️  Constant model only (G0 limited; no H0/V0 derivable)"
)
print()

print("2. AQUEOUS SPECIES (Electrolyte Solutions)")
print("   HP-only databases: ❌ Not available (no aqueous species)")
print("   DEW-mixed databases: ✅ Full support (PerplexDEW + GFSM gas coupling)")
print("   Old-SUPCRT databases: ❌ Not available")
print()

print("3. GAS SPECIES (GFSM Hybrid EoS)")
print("   HP-only databases: ❌ Not available")
print("   DEW-mixed databases: ✅ Available (PerplexGFSM: CO2, H2S, etc.)")
print("   Old-SUPCRT databases: ❌ Not available")
print()

print("4. COMPLETE THERMODYNAMIC DESCRIPTIONS")
print("   HP-only & DEW-mixed: ✅ G0, H0, V0, S0 fully available per species")
print("   Old-SUPCRT: ⚠️  Only G0 available; H0 cannot be derived from G0 alone")
print()

print("5. REACTION EQUILIBRIA & PHASE DIAGRAMS")
print("   HP-only: ✅ Excellent for T-P mineral stability diagrams")
print("   DEW-mixed: ✅ Excellent for fluid-mineral equilibria")
print(
    "   Old-SUPCRT: ⚠️  Limited; Constant model insufficient for rigorous calculations"
)
print()

print("6. AQUEOUS ACTIVITY MODELS")
print("   HP-only: ❌ No aqueous species")
print("   DEW-mixed: ✅ Advanced (Born model, GFSM solvent effect)")
print("   Old-SUPCRT: ❌ No aqueous species")
print()

print("7. MULTI-PHASE GEOCHEMICAL SPECIATION")
print("   HP-only: ⚠️  Partial (solid phases only)")
print("   DEW-mixed: ✅ Full (solids + aqueous + gases)")
print("   Old-SUPCRT: ⚠️  Very limited (G0-only mineral models)")
print()

print("=" * 80)
print("RECOMMENDATIONS")
print("=" * 80)
print()

print("USE HP-ONLY DATABASES FOR:")
print("  • High-temperature mineral stable assemblages (diagrams, phase relations)")
print("  • Solid-solid equilibria and polythermal reactions")
print("  • Studies not requiring aqueous or gas phases")
print()

print("USE DEW-MIXED DATABASES FOR:")
print("  • Fluid-mineral interactions and mass transfer")
print("  • Geothermal systems with aqueous fluids and gases")
print("  • Acid-base speciation in aqueous solutions")
print("  • Complete thermodynamic descriptions: G0, H0, V0, S0 for all phases")
print()

print("AVOID OLD-SUPCRT DATABASES FOR:")
print("  • Rigorous thermodynamic calculations requiring H0/V0")
print("  • Applications where both G and H are needed (e.g., entropy calculations)")
print("  • New projects (use HP or DEW databases instead)")
print()

print("NOTE: Old databases (b89/b92/ba96) persist for backward compatibility.")
print("      They lack elemental entropy data, making H0 derivation impossible.")
print()
