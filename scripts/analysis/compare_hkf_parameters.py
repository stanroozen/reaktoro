#!/usr/bin/env python3
"""
Compare HKF parameters between DEW database and Excel truth to identify discrepancy sources
"""

import yaml
import pandas as pd

# Load DEW database
with open("embedded/databases/DEW/dew2024-aqueous.yaml", "r") as f:
    dew_db = yaml.safe_load(f)

# Get species data
co2_data = dew_db["Species"]["CO2(0)"]["StandardThermoModel"]["HKF"]
hco3_data = dew_db["Species"]["HCO3(-)"]["StandardThermoModel"]["HKF"]
h_data = None  # H+ is reference, usually in database with Gf=0
h2o_data = None  # H2O uses different model

print("=" * 100)
print("HKF PARAMETERS IN DEW DATABASE")
print("=" * 100)

print("\nCO2(aq) - DEW 2024 database:")
for key, val in co2_data.items():
    print(f"  {key:<10}: {val}")

print("\nHCO3- - DEW 2024 database:")
for key, val in hco3_data.items():
    print(f"  {key:<10}: {val}")

# Convert from cal/mol to J/mol for comparison
CAL_TO_J = 4.184

print("\n" + "=" * 100)
print("UNIT CONVERSIONS (cal/mol → J/mol)")
print("=" * 100)

print("\nCO2(aq):")
print(f"  Gf: {co2_data['Gf']} cal/mol = {co2_data['Gf'] * CAL_TO_J:.2f} J/mol")
print(f"  Sr: {co2_data['Sr']} cal/mol-K = {co2_data['Sr'] * CAL_TO_J:.6f} J/mol-K")
print(f"  c1: {co2_data['c1']} cal/mol-K = {co2_data['c1'] * CAL_TO_J:.6f} J/mol-K")
print(f"  c2: {co2_data['c2']} cal-K/mol = {co2_data['c2'] * CAL_TO_J:.2f} J-K/mol")
print(f"  wref: {co2_data['wref']} cal/mol = {co2_data['wref'] * CAL_TO_J:.2f} J/mol")

print("\nHCO3-:")
print(f"  Gf: {hco3_data['Gf']} cal/mol = {hco3_data['Gf'] * CAL_TO_J:.2f} J/mol")
print(f"  Sr: {hco3_data['Sr']} cal/mol-K = {hco3_data['Sr'] * CAL_TO_J:.6f} J/mol-K")
print(f"  c1: {hco3_data['c1']} cal/mol-K = {hco3_data['c1'] * CAL_TO_J:.6f} J/mol-K")
print(f"  c2: {hco3_data['c2']} cal-K/mol = {hco3_data['c2'] * CAL_TO_J:.2f} J-K/mol")
print(f"  wref: {hco3_data['wref']} cal/mol = {hco3_data['wref'] * CAL_TO_J:.2f} J/mol")

print("\n" + "=" * 100)
print("POSSIBLE DISCREPANCY SOURCES")
print("=" * 100)

print(
    """
Given that:
- CO2 error: +5 J/mol (model higher than Excel)
- HCO3- error: -9 J/mol (model lower than Excel)
- H2O error: 0 J/mol (perfect match)

Potential sources of discrepancies:

1. **Reference State Parameters (Gf, Hf, Sr at 298.15K, 1 bar)**
   - Excel might use different values than DEW 2024 database
   - Source databases: SUPCRT92, SUPCRT07, DEW 2019, DEW 2024
   - Small differences in Gf propagate directly to G0

2. **HKF Volume Parameters (a1, a2, a3, a4)**
   - These are pressure-dependent terms
   - Affect both V0 and pressure corrections to G0
   - CO2: a1={co2_a1:.6e}, a2={co2_a2:.2f}, a3={co2_a3:.6e}, a4={co2_a4:.2f}
   - HCO3-: a1={hco3_a1:.6e}, a2={hco3_a2:.2f}, a3={hco3_a3:.6e}, a4={hco3_a4:.2f}

3. **HKF Heat Capacity Parameters (c1, c2)**
   - Affect temperature-dependent corrections
   - CO2: c1={co2_c1:.2f}, c2={co2_c2:.2f} cal/mol-K
   - HCO3-: c1={hco3_c1:.2f}, c2={hco3_c2:.2f} cal/mol-K

4. **Born Function wref (omega at reference state)**
   - Critical for solvation energy
   - Depends on effective Born radius
   - CO2: wref={co2_wref:.2f} cal/mol (neutral species)
   - HCO3-: wref={hco3_wref:.2f} cal/mol (charged species)
   - For ions: ω scales with charge²/radius
   - Small changes in ionic radius → large changes in wref

5. **Born Function Model Differences**
   - Code uses: waterBornOmegaDEW() in WaterBornOmegaDEW.cpp
   - Excel might use: Different Born coefficient calculation
   - Different effective ionic radius for HCO3-
   - Different solvent dielectric treatment

6. **Dielectric Constant Calculation**
   - Code uses: Johnson-Norton (1991) via WaterElectroPropsJohnsonNorton.cpp
   - Excel might use: Different dielectric model (Shock-Helgeson? Franck?)
   - Small ε differences → different Born Z,Y,Q → different solvation

7. **Water Density Values**
   - Both use Zhang-Duan 2005 EOS (confirmed from test code)
   - However: iteration tolerance matters
   - Code uses: 0.001 bar density tolerance
   - Excel might use: different tolerance or iteration method

8. **Temperature/Pressure Interpolation**
   - Excel VBA might have rounding/truncation differences
   - Logarithm precision: log((ψ+P)/(ψ+Pr))
   - Power calculations: (T-θ)⁻¹ terms

MOST LIKELY CULPRITS (ranked):

A. **Different HKF parameter SOURCE** (most likely)
   - Excel might use SUPCRT92/SUPCRT07 parameters
   - Reaktoro uses DEW 2024 parameters
   - Even small Gf differences (few cal/mol) = observable errors

B. **Born omega calculation differences**
   - HCO3- is charged → large Born contribution
   - Different ionic radius or Z,Y calculation method
   - This explains why HCO3- error (9 J/mol) > CO2 error (5 J/mol)

C. **Dielectric constant model**
   - Different ε(T,P) model used in Excel
   - Affects Born functions for both species
   - Neutral CO2 still affected via Born term

D. **Numerical precision/rounding**
   - Excel VBA vs C++ double precision
   - Unlikely to cause systematic 5-9 J/mol errors
""".format(
        co2_a1=co2_data["a1"],
        co2_a2=co2_data["a2"],
        co2_a3=co2_data["a3"],
        co2_a4=co2_data["a4"],
        hco3_a1=hco3_data["a1"],
        hco3_a2=hco3_data["a2"],
        hco3_a3=hco3_data["a3"],
        hco3_a4=hco3_data["a4"],
        co2_c1=co2_data["c1"],
        co2_c2=co2_data["c2"],
        hco3_c1=hco3_data["c1"],
        hco3_c2=hco3_data["c2"],
        co2_wref=co2_data["wref"],
        hco3_wref=hco3_data["wref"],
    )
)

print("\n" + "=" * 100)
print("RECOMMENDED INVESTIGATION STEPS")
print("=" * 100)

print("""
1. VERIFY EXCEL PARAMETER SOURCE:
   - Check which HKF database Excel uses (SUPCRT92? SUPCRT07? DEW 2019?)
   - Compare Gf, Hf, Sr, a1-a4, c1-c2, wref values with DEW 2024
   - Even 1-2 cal/mol difference in Gf explains observed errors

2. CHECK BORN OMEGA CALCULATION:
   - Compare waterBornOmegaDEW() output with Excel
   - Test at T=650°C, P=15 kb for both CO2 and HCO3-
   - Verify ionic radius used for HCO3- (typical: 2.06 Å)
   - Check Born Z,Y,Q calculation methods

3. COMPARE DIELECTRIC CONSTANTS:
   - Extract ε(T,P) from Excel at test conditions
   - Compare with Johnson-Norton (1991) from Reaktoro
   - Check if Excel uses Helgeson-Kirkham or other model

4. SENSITIVITY ANALYSIS:
   - Vary Gf by ±10 cal/mol, see impact on G0
   - Vary wref by ±100 cal/mol, see impact
   - Vary a1-a4 by 1%, see impact at high pressure
   - Identify which parameter changes reproduce observed errors

5. CROSS-CHECK WITH PUBLISHED DATA:
   - Compare CO2(aq) G0 at test conditions with literature
   - Compare HCO3- G0 with Shock & Helgeson (1988)
   - Determine which implementation (Excel or Reaktoro) matches literature
""")
