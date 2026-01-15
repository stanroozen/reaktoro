#!/usr/bin/env python3
"""
Investigation: Physical Constants, Density, Epsilon, and Convention Differences

These could explain 5-9 J/mol discrepancies even with identical parameters:
1. Standard constants (R, CAL2J)
2. Water density (Zhang-Duan EOS)
3. Dielectric constant epsilon (Johnson-Norton model)
4. Thermodynamic conventions (apparent vs standard Gibbs)
"""

print("=" * 100)
print("INVESTIGATION: CONSTANTS, DENSITY, EPSILON, AND CONVENTIONS")
print("=" * 100)

print("""
The 5-9 J/mol species errors could come from:

A) Different PHYSICAL CONSTANTS
B) Different WATER DENSITY calculation
C) Different DIELECTRIC CONSTANT calculation
D) Different THERMODYNAMIC CONVENTIONS

Let's investigate each possibility...
""")

print("\n" + "=" * 100)
print("A. PHYSICAL CONSTANTS")
print("=" * 100)

# Constants from build_dew_reaktoro_db.py
CAL2J_python = 4.184
R_python = 8.31446261815324

print(f"""
Python extraction code (build_dew_reaktoro_db.py):
  CAL2J = {CAL2J_python}         [J/cal] - conversion factor
  R = {R_python} [J/(mol·K)] - gas constant

Common values in literature:
  CAL2J (thermochemical calorie):
    - 4.184 J/cal (exact, by definition since 1956) ← STANDARD
    - 4.1868 J/cal (old 15°C calorie)
    - 4.18580 J/cal (mean calorie)

  R (gas constant):
    - 8.314462618 J/(mol·K) (CODATA 2018) ← CURRENT STANDARD
    - 8.31451 J/(mol·K) (older value)
    - 1.98720 cal/(mol·K) = 8.314462618 J/(mol·K) / 4.184

POTENTIAL ISSUE:
If Excel uses CAL2J = 4.1868 instead of 4.184:
  - Error = (4.1868 - 4.184) / 4.184 = 0.065% relative
  - For Gf ~ 500,000 J/mol: Error ~ 325 J/mol (HUGE!)
  - For wref ~ 300,000 J/mol: Error ~ 195 J/mol

But this is TOO LARGE to explain 5-9 J/mol errors.
So constants are likely the same.

VERDICT: Constants probably NOT the issue (would cause >100 J/mol errors)
""")

print("\n" + "=" * 100)
print("B. WATER DENSITY (Zhang-Duan 2005 EOS)")
print("=" * 100)

print("""
Water density affects:
1. Volume integration for H2O Gibbs energy
2. Born omega calculation for aqueous species (via solvent function g)

Reaktoro uses: Zhang-Duan (2005) EOS
Excel might use:
  - Zhang-Duan (2005) ✓ SAME
  - IAPWS-95 (different formulation)
  - Wagner & Pruss (2002)
  - Older Haar-Gallagher-Kell (1984)

TEST: H2O Gibbs error = 0 J/mol at T=650°C, P=15 kb
  → This proves density calculation is CORRECT and MATCHES Excel

Evidence:
  - H2O G0 error = 0 (perfect match)
  - Volume integration = ∫(M/ρ) dP
  - If ρ were wrong, H2O G0 would be wrong
  - Since H2O G0 is perfect, ρ is correct

VERDICT: Density is CORRECT (proven by H2O error = 0)
""")

print("\n" + "=" * 100)
print("C. DIELECTRIC CONSTANT (Johnson-Norton 1991)")
print("=" * 100)

print("""
Dielectric constant ε(T,P) affects:
1. Born function Z = -1/ε (enters omega calculation)
2. Solvent function g(T,P,ρ,ε) for charged species

Reaktoro uses: Johnson-Norton (1991) correlation
Excel might use:
  - Johnson-Norton (1991) ✓ SAME (likely, it's the DEW standard)
  - Franck (1990) model
  - Different polynomial fit
  - Experimental data interpolation

Formula comparison:
Johnson-Norton (1991):
  ε(T,P) = A + B/T + C/T² + ... (pressure-dependent coefficients)

Test point: T=650°C = 923 K, P=15 kb = 1500 MPa
  Expected ε ≈ 3-5 (very low at high T)

Born function: Z = -1/ε
  If ε = 4: Z = -0.25
  If ε = 5: Z = -0.20
  Difference: ΔZ = 0.05 (20% relative)

Impact on omega for HCO3- (Z=-1, wref=532748 J/mol):
  Born term = -ω*(Z+1) where ω depends on Z

  For ions, omega calculation:
  ω(T,P) = η*(Z²/re - Z/(3.082+g))

  Where g = solvent function that depends on ε

  If g changes by 1% → ω changes by ~0.5%
  For wref=532748 J/mol: Δω ≈ 2600 J/mol

  This affects Born solvation term:
  ΔG_Born = Δω * (Z+1) = 2600 * 0 = 0 for HCO3- (charge=-1, Z+1=0)

  Wait, that's not quite right. Let me reconsider...

Actually, the Born term in HKF is:
  G0_Born = -(w + wref*Yr*(T-Tr))*(Z+1)

Where:
  w = omega(T,P) calculated from Born model
  Z = Born function = charge effect, NOT dielectric function

The dielectric ε enters through:
  - Solvent function g(T,P,ρ,ε)
  - Which affects omega calculation
  - Which affects G0

For neutral species (CO2):
  - Small Born contribution
  - Less sensitive to ε

For charged species (HCO3-):
  - Large Born contribution
  - Very sensitive to ε through omega

POTENTIAL ISSUE:
If Excel uses different ε model:
  → Different g(T,P)
  → Different ω(T,P)
  → Different G0 by 1-10 J/mol for ions ⚠

But H2O error = 0 suggests ε is correct too (H2O also uses ε indirectly).

VERDICT: Epsilon likely CORRECT but could differ slightly for Born calculations
        (would need to compare Excel's dielectric model directly)
""")

print("\n" + "=" * 100)
print("D. THERMODYNAMIC CONVENTIONS")
print("=" * 100)

print("""
CRITICAL QUESTION: Are we comparing the same thermodynamic quantities?

Standard Gibbs Free Energy (G0):
  - Referenced to elements in standard states at T,P
  - G0 = Gf + ∫[Tr→T] Cp dT + ∫[Pr→P] V dP + other terms

Apparent Gibbs Free Energy (G_apparent):
  - Referenced to ions of infinite dilution
  - Used in equilibrium calculations
  - Accounts for activity coefficients implicitly

Formation Energy (Gf):
  - At reference state (25°C, 1 bar)
  - Gf(elements) = 0 by definition

POTENTIAL ISSUE:
If Excel reports G_apparent and Reaktoro calculates G0:
  - Could differ by activity coefficient terms
  - Typically 1-10 kJ/mol = 1000-10000 J/mol (TOO LARGE)

If Excel reports G0_bar (partial molar) vs G0_infinity (infinite dilution):
  - Could differ by a few J/mol at high P,T
  - This is more reasonable ⚠

HKF MODEL CONVENTION:
The HKF model calculates G0 at infinite dilution:
  G0_∞(T,P) = Gf + ΔG(T,Pr) + ΔG(Pr,P) + Born terms

Where:
  - Gf is formation energy at Tr=25°C, Pr=1 bar
  - All properties are for infinite dilution
  - No activity coefficients involved

EXCEL MIGHT USE:
A) Same convention (G0_∞) ✓ Most likely
B) Apparent properties with γ=1 assumption ⚠ Possible
C) Different reference pressure (Pref = 1 bar vs 1000 bar?) ⚠ Possible

TEST CASE:
For our reaction: H2O + CO2(aq) ⇌ H+ + HCO3-

Excel log_K vs Reaktoro log_K:
  - If conventions differ, log_K should differ systematically
  - If log_K matches but species G0 differ → likely pressure reference difference

VERDICT: Convention differences possible but unlikely to be the main issue
         (would need to check Excel's exact calculation method)
""")

print("\n" + "=" * 100)
print("E. REFERENCE STATE DIFFERENCES")
print("=" * 100)

print("""
CRITICAL: What reference pressure does each code use?

HKF MODEL STANDARD:
  - Standard state: Pr = 1 bar (0.1 MPa)
  - But water properties often referenced to Psat(T)

DEW MODEL MODIFICATION:
  - Reference pressure: Pr = 1000 bar (100 MPa)?
  - Or pressure-dependent reference?
  - This affects the integration limits

Integration formula in Reaktoro:
  G(T,P) = G(T,Pref) + ∫[Pref→P] V(T,P') dP'

From WaterGibbsModel.cpp:
  - Integration from 1000 bar (1 kb) to target P
  - Uses GAtOneKb_cal_per_mol() as reference

If Excel uses different reference:
  - Pref = 1 bar instead of 1000 bar
  - Integration limits differ
  - Result differs by ∫[1→1000 bar] V dP
  - For water: V ~ 18 cm³/mol, ΔP = 999 bar
  - ΔG ~ 18e-6 m³/mol * 999e5 Pa = 1800 J/mol

But we have H2O error = 0, so reference pressure must match.

VERDICT: Reference states appear consistent (H2O error = 0 proves it)
""")

print("\n" + "=" * 100)
print("F. SUMMARY OF LIKELY CAUSES")
print("=" * 100)

print("""
Ranking by likelihood to explain 5-9 J/mol errors:

1. ✓ DIFFERENT HKF PARAMETERS (HIGH - most likely)
   - Excel uses DEW 2019, Reaktoro uses DEW 2024
   - Or Excel uses SUPCRT database
   - Gf, a1-a4, c1-c2, wref could differ by 0.01-0.1%
   - Would cause 1-10 J/mol errors ✓ MATCHES OBSERVATIONS

2. ⚠ DIFFERENT DIELECTRIC MODEL (MEDIUM - possible)
   - If Excel uses different ε(T,P) correlation
   - Affects Born omega for ions
   - Could cause 1-5 J/mol errors for charged species
   - Need to verify Excel's dielectric model

3. ⚠ OMEGA CALCULATION DIFFERENCES (MEDIUM - possible)
   - If Excel uses different Born model formulation
   - Shock et al. vs Helgeson vs DEW
   - Could cause 1-10 J/mol errors
   - Need to verify Excel's omega formula

4. ✗ PHYSICAL CONSTANTS (LOW - unlikely)
   - CAL2J and R mismatch would cause >100 J/mol errors
   - Observed errors are only 5-9 J/mol
   - Almost certainly constants match

5. ✗ DENSITY CALCULATION (NONE - proven correct)
   - H2O error = 0 proves density is correct
   - Both codes use same Zhang-Duan EOS
   - This is NOT the issue ✓

6. ✗ REFERENCE STATE (NONE - proven correct)
   - H2O error = 0 proves reference state matches
   - Integration limits are the same
   - This is NOT the issue ✓

7. ✗ THERMODYNAMIC CONVENTIONS (LOW - unlikely)
   - Both use G0 at infinite dilution (standard HKF)
   - Convention differences would cause larger errors
   - Almost certainly same convention
""")

print("\n" + "=" * 100)
print("G. RECOMMENDATION FOR NEXT INVESTIGATION")
print("=" * 100)

print("""
To identify the exact source of 5-9 J/mol errors:

1. CHECK DIELECTRIC CONSTANT:
   File: Reaktoro/Models/WaterDielectricConstant*.cpp
   - Which model is used? Johnson-Norton 1991?
   - Calculate ε at T=650°C, P=15 kb
   - Compare with Excel's ε value
   - If different: This could explain ion errors (HCO3-)

2. CHECK OMEGA FORMULA:
   File: Reaktoro/Models/WaterBornOmegaDEW.cpp
   - Formula: ω = η*(Z²/re - Z/(3.082+g))
   - Is this the same as Excel's formula?
   - Check η = 166027 cal·Å/mol constant
   - Check solvent function g calculation

3. CHECK HKF PARAMETER SOURCE:
   - What database does Excel use?
   - SUPCRT92, SUPCRT07, DEW 2019, or DEW 2024?
   - Compare Gf, a1-a4, c1-c2, wref values one-by-one
   - Look for 0.01-0.1% differences

4. VERIFY CONSTANTS IN C++:
   File: Reaktoro/Common/Constants.hpp (or similar)
   - R = 8.314462618... J/(mol·K)?
   - CAL2J = 4.184 exactly?
   - Other physical constants

5. CHECK BORN FUNCTION CALCULATION:
   - Is Z calculated correctly from ε?
   - Formula: Z = (some function of ε, ρ, T, P)
   - Compare Reaktoro vs Excel Born function values

MOST LIKELY: Parameter source differences (DEW version mismatch)
NEXT LIKELY: Dielectric or Born model formulation differences
LEAST LIKELY: Constants, density, reference state (proven correct by H2O=0)
""")

print("\n" + "=" * 100)
