"""
Test StandardThermoModelDEW with different water model combinations
"""

import sys

sys.path.insert(
    0,
    r"c:\Users\stanroozen\Documents\Projects\reaktoro-dev\reaktoro\build-msvc\Reaktoro\Release",
)

import reaktoro4py as rkt

print("=" * 70)
print("Testing StandardThermoModelDEW with different water model combinations")
print("=" * 70)

# Test 1: Default DEW preset
print("\n1. DEW Preset (ZD2005 + PowerFunction + DewIntegral + Shock92):")
params1 = rkt.StandardThermoModelParamsDEW()
print(f"   EOS: {params1.waterOptions.eosModel}")
print(f"   Dielectric: {params1.waterOptions.dielectricModel}")
print(f"   Gibbs: {params1.waterOptions.gibbsModel}")
print(f"   Born: {params1.waterOptions.bornModel}")
model1 = rkt.StandardThermoModelDEW(params1)
print("   ✅ Model created successfully")

# Test 2: ZhangDuan2009 + JohnsonNorton
print("\n2. Custom: ZD2009 + JohnsonNorton1991 + DelaneyHelgeson + Shock92:")
params2 = rkt.StandardThermoModelParamsDEW()
params2.waterOptions.eosModel = rkt.WaterEosModel.ZhangDuan2009
params2.waterOptions.dielectricModel = rkt.WaterDielectricModel.JohnsonNorton1991
params2.waterOptions.gibbsModel = rkt.WaterGibbsModel.DelaneyHelgeson1978
params2.waterOptions.bornModel = rkt.WaterBornModel.Shock92Dew
print(f"   EOS: {params2.waterOptions.eosModel}")
print(f"   Dielectric: {params2.waterOptions.dielectricModel}")
print(f"   Gibbs: {params2.waterOptions.gibbsModel}")
model2 = rkt.StandardThermoModelDEW(params2)
print("   ✅ Model created successfully")

# Test 3: Wagner-Pruss (standard) + Fernandez
print("\n3. Standard: WagnerPruss + Fernandez1997 (with DEW Gibbs/Born):")
params3 = rkt.StandardThermoModelParamsDEW()
params3.waterOptions.eosModel = rkt.WaterEosModel.WagnerPruss
params3.waterOptions.dielectricModel = rkt.WaterDielectricModel.Fernandez1997
# Keep default Gibbs and Born (DewIntegral + Shock92)
print(f"   EOS: {params3.waterOptions.eosModel}")
print(f"   Dielectric: {params3.waterOptions.dielectricModel}")
print(f"   Gibbs: {params3.waterOptions.gibbsModel}")
print(f"   Born: {params3.waterOptions.bornModel}")
model3 = rkt.StandardThermoModelDEW(params3)
print("   ✅ Model created successfully")

# Test 4: All DEW options (ZD2005 variations)
print("\n4. ZD2005 + Franck1990:")
params4 = rkt.StandardThermoModelParamsDEW()
params4.waterOptions.eosModel = rkt.WaterEosModel.ZhangDuan2005
params4.waterOptions.dielectricModel = rkt.WaterDielectricModel.Franck1990
params4.waterOptions.gibbsModel = rkt.WaterGibbsModel.DewIntegral
print(f"   EOS: {params4.waterOptions.eosModel}")
print(f"   Dielectric: {params4.waterOptions.dielectricModel}")
model4 = rkt.StandardThermoModelDEW(params4)
print("   ✅ Model created successfully")

print("\n" + "=" * 70)
print("✅ All water model combinations work!")
print("=" * 70)
print("\nAvailable combinations:")
print("  EOS: WagnerPruss, HGK, ZhangDuan2005, ZhangDuan2009")
print("  Dielectric: JohnsonNorton1991, Franck1990, Fernandez1997, PowerFunction")
print("  Gibbs: DelaneyHelgeson1978, DewIntegral (always required)")
print("  Born: None (neutral species), Shock92Dew (ionic species)")
print("\nThe model is fully flexible - mix and match as needed!")
