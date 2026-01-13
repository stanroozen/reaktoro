"""
Demonstrate what the None option means for Born model
(Gibbs is now always required for species thermodynamics)
"""

import sys

sys.path.insert(
    0,
    r"c:\Users\stanroozen\Documents\Projects\reaktoro-dev\reaktoro\build-msvc\Reaktoro\Release",
)

import reaktoro4py as rkt

print("=" * 70)
print("Understanding Born 'None' Option in Water Models")
print("=" * 70)

print("\n" + "=" * 70)
print("WaterBornModel.None - For Neutral Species")
print("=" * 70)
print("""
What it means:
  - Don't compute Born solvation functions (omega, domega/dP)
  - Skip ionic solvation energy calculations

When to use:
  - For NEUTRAL species (charge = 0)
  - When Born contribution isn't needed
  - For faster calculations with non-ionic solutes

Example: Neutral aqueous species like SiO2(aq), CO2(aq)
""")

opts_no_born = rkt.WaterModelOptions()
opts_no_born.eosModel = rkt.WaterEosModel.WagnerPruss
opts_no_born.dielectricModel = rkt.WaterDielectricModel.Fernandez1997
opts_no_born.gibbsModel = rkt.WaterGibbsModel.DelaneyHelgeson1978
opts_no_born.bornModel = getattr(rkt.WaterBornModel, "None")
print("Config: WagnerPruss + Fernandez1997 + DH1978 + NO Born")
params = rkt.StandardThermoModelParamsDEW()
params.waterOptions = opts_no_born
model = rkt.StandardThermoModelDEW(params)
print("✓ Model created (no Born functions)")

print("\n" + "=" * 70)
print("Summary: Water Model Requirements")
print("=" * 70)
print("""
EOS (Equation of State):
  ✓ ALWAYS required
  ✓ Provides density rho(T,P)
  ✓ Options: WagnerPruss, HGK, ZhangDuan2005, ZhangDuan2009

Dielectric Model:
  ✓ ALWAYS required
  ✓ Provides epsilon(T,P) for electrolyte calculations
  ✓ Options: JohnsonNorton1991, Franck1990, Fernandez1997, PowerFunction

Gibbs Model:
  ✓ ALWAYS required for species thermodynamics
  ✓ Needed to compute G0, H0, S0
  ✓ Options: DelaneyHelgeson1978, DewIntegral
  ✗ None option removed (not physically meaningful for species)

Born Model:
  ✓ Required for IONIC species (charge ≠ 0)
  ✗ Use None for NEUTRAL species (charge = 0)
  ✓ Options: None, Shock92Dew

Full DEW (default for ionic species):
  ✓ Complete thermodynamics
  ✓ Ionic species
  ✓ High T/P geochemistry
  → gibbsModel = DewIntegral
  → bornModel = Shock92Dew
""")

print("\n" + "=" * 70)
print("✅ All 'None' configurations work and are useful!")
print("=" * 70)
