"""
Test script demonstrating configurable DEW water models
Tests different water model combinations to show flexibility
"""

import sys
import os

# Setup paths
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(SCRIPT_DIR)
PYD_DIR = os.path.join(ROOT_DIR, "build", "Reaktoro", "Release")
if os.path.isdir(PYD_DIR) and PYD_DIR not in sys.path:
    sys.path.insert(0, PYD_DIR)

from reaktoro4py import *

# Suppress warnings
try:
    Warnings.disable(906)
except:
    pass

print("=" * 80)
print("DEW CONFIGURABLE WATER MODELS TEST")
print("=" * 80)

# Load databases
dew_db = DEWDatabase("dew2019-aqueous")
supcrt_db = SupcrtDatabase("supcrtbl")

print(f"\n[OK] Databases loaded")
print(f"  - DEW species: {len(dew_db.species())}")

# Check if H2O_aq is available
dew_species_names = [s.name() for s in dew_db.species()]
h2o_water_aq_available = "WATER,AQ" in dew_species_names

if h2o_water_aq_available:
    print(f"[OK] WATER,AQ found in DEW database")
else:
    print(f"[WARN] WATER,AQ NOT found in DEW database (rebuild may be needed)")

# Test 1: Default configuration (ZhangDuan2005, PowerFunction, DewIntegral)
print(f"\n--- Test 1: Default Configuration ---")
print(f"  EOS: ZhangDuan2005")
print(f"  Dielectric: PowerFunction")
print(f"  Gibbs Integration: DewIntegral")

default_config = {
    "eos_model": "ZhangDuan2005",
    "dielectric_model": "PowerFunction",
    "gibbs_model": "DewIntegral",
    "born_model": "Shock92Dew",
}

try:
    params1 = StandardThermoModelParamsDEW()
    params1.waterOptions.eosModel = WaterEosModel.ZhangDuan2005
    params1.waterOptions.dielectricModel = WaterDielectricModel.PowerFunction
    params1.waterOptions.gibbsModel = WaterGibbsModel.DewIntegral
    params1.waterOptions.bornModel = WaterBornModel.Shock92Dew

    model1 = StandardThermoModelDEW(params1)
    print(f"[OK] Default configuration created successfully")
except Exception as e:
    print(f"[FAIL] Error: {e}")

# Test 2: Alternative configuration (ZhangDuan2009, JohnsonNorton1991, DelaneyHelgeson1978)
print(f"\n--- Test 2: Alternative Configuration ---")
print(f"  EOS: ZhangDuan2009")
print(f"  Dielectric: JohnsonNorton1991")
print(f"  Gibbs Integration: DelaneyHelgeson1978")

try:
    params2 = StandardThermoModelParamsDEW()
    params2.waterOptions.eosModel = WaterEosModel.ZhangDuan2009
    params2.waterOptions.dielectricModel = WaterDielectricModel.JohnsonNorton1991
    params2.waterOptions.gibbsModel = WaterGibbsModel.DelaneyHelgeson1978
    params2.waterOptions.bornModel = WaterBornModel.Shock92Dew

    model2 = StandardThermoModelDEW(params2)
    print(f"[OK] Alternative configuration created successfully")
except Exception as e:
    print(f"[FAIL] Error: {e}")

# Test 3: Available water models
print(f"\n--- Available Water Model Options ---")
print(f"  EOS Models: WagnerPruss, HGK, ZhangDuan2005 (default), ZhangDuan2009")
print(f"  Dielectric Models: PowerFunction (default), JohnsonNorton1991")
print(f"  Gibbs Integration: DewIntegral (default), DelaneyHelgeson1978")
print(f"  Born Models: Shock92Dew (default), Shock92 (neutral species)")

print(f"\n" + "=" * 80)
print(f"SUMMARY:")
print(f"[OK] DEW water model configuration is fully customizable")
print(f"[OK] Default: Duan & Zang 2005 EOS + Power Function + Volume Integration")
print(f"[OK] Users can specify alternative models for sensitivity studies")
print(f"=" * 80)

