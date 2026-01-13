"""
Verify HKF calculations against Excel formulas from the Calculations sheet.

This script checks if our C++ implementation matches the Excel VBA logic
for aqueous species Gibbs energy and volume calculations.
"""

import math

# Constants (from Excel Calculations sheet)
Tr = 298.15  # K (C8)
Pr = 1.0  # bar (C9)
eps_r = 78.47  # dimensionless (C10)
R = 1.9858775  # cal/(mol·K) (C11)
Psi = 2600.0  # bar (C12)
Theta = 228.0  # K (C13)
Y = -5.79865e-5  # K^-1 (C14)
Zr = -1.278055636e-2  # dimensionless

# Test conditions (from reactionTesttruth.csv, row 1)
T_C = 300.0  # °C
P_bar = 5000.0  # bar
T_K = T_C + 273.15

# Water properties at 300°C, 5000 bar (from test output or calculations)
rho_H2O = 1.0  # g/cm³ (approximate, need actual value)
eps = 10.0  # dimensionless (approximate, need actual value)
Q_born = 0.0  # Born Q coefficient (need actual value)

print("=" * 80)
print("HKF Calculation Verification: HCO3- at 300°C, 5000 bar")
print("=" * 80)

# HCO3- parameters from dew2024-aqueous.yaml
params_HCO3 = {
    "Gf": -586939.888,  # cal/mol
    "Hf": -689933.2320000001,  # cal/mol
    "Sr": 98.44952,  # cal/(mol·K)
    "a1": 3.20076e-05,  # (dimensionless)
    "a2": 0.0384928,  # (dimensionless)
    "a3": 2.5104e-05,  # (dimensionless)
    "a4": -0.0011798880000000002,  # (dimensionless)
    "c1": 46.024,  # cal/(mol·K)
    "c2": -0.00158992,  # (dimensionless)
    "wref": 5.327487200000001e-05,  # (dimensionless)
    "charge": -1.0,
}

# Note: Excel stores a1, a2, a3, a4 in different units
# Need to check if these are already converted or in original HKF units
# Original HKF units: a1 (cal/bar), a2 (cal), a3 (cal·K/bar), a4 (cal·K)

print("\nHCO3- Parameters:")
for key, value in params_HCO3.items():
    print(f"  {key:8s} = {value}")

# Calculate temperature-dependent terms
T_diff = T_K - Tr
Tth = T_K - Theta
Trth = Tr - Theta

print(f"\nTemperature terms:")
print(f"  T = {T_K} K")
print(f"  T - Tr = {T_diff} K")
print(f"  T - Theta = {Tth} K")

# Calculate pressure-dependent terms
P_diff = P_bar - Pr
psiP = Psi + P_bar
psiPr = Psi + Pr

print(f"\nPressure terms:")
print(f"  P = {P_bar} bar")
print(f"  P - Pr = {P_diff} bar")
print(f"  Psi + P = {psiP} bar")
print(f"  Psi + Pr = {psiPr} bar")

# Excel HKF formula for volume (from your description):
# V0 = a1*(P-Pr) + a2*ln((Psi+P)/(Psi+Pr))
#      + (1/(T-Theta)) * [a3*(P-Pr) + a4*ln((Psi+P)/(Psi+Pr))]
#      + Born terms

# Check units: Excel may have different conventions
# Let's calculate with assumption that a1, a2, a3, a4 are in proper units

print("\n" + "=" * 80)
print("VOLUME CALCULATION")
print("=" * 80)

# Non-Born volume terms
V_term1 = params_HCO3["a1"] * P_diff
V_term2 = params_HCO3["a2"] * math.log(psiP / psiPr)
V_term3 = (1.0 / Tth) * (
    params_HCO3["a3"] * P_diff + params_HCO3["a4"] * math.log(psiP / psiPr)
)

print(f"\nNon-Born volume components:")
print(f"  a1*(P-Pr) = {V_term1:.6e}")
print(f"  a2*ln((Psi+P)/(Psi+Pr)) = {V_term2:.6e}")
print(f"  (1/(T-Theta))*[a3*(P-Pr) + a4*ln(...)] = {V_term3:.6e}")
print(f"  Sum (non-Born) = {V_term1 + V_term2 + V_term3:.6e}")

# Note: Born terms require omega, Q, and their derivatives
# These depend on water properties (eps, rho) which we need from actual calculations

print("\n" + "=" * 80)
print("GIBBS ENERGY CALCULATION")
print("=" * 80)

# Excel HKF formula for Gibbs (from your description):
# G0 = Gf - Sr*(T-Tr) - c1*(T*ln(T/Tr) - T + Tr)
#      + a1*(P-Pr) + a2*ln((Psi+P)/(Psi+Pr))
#      - c2*[(1/(T-Theta) - 1/(Tr-Theta))*(Theta-T)/Theta
#            - T/(Theta^2)*ln(Tr/T * (T-Theta)/(Tr-Theta))]
#      + (1/(T-Theta))*[a3*(P-Pr) + a4*ln((Psi+P)/(Psi+Pr))]
#      + Born terms

G_term1 = params_HCO3["Gf"]
G_term2 = -params_HCO3["Sr"] * T_diff
G_term3 = -params_HCO3["c1"] * (T_K * math.log(T_K / Tr) - T_K + Tr)
G_term4 = params_HCO3["a1"] * P_diff
G_term5 = params_HCO3["a2"] * math.log(psiP / psiPr)

# c2 term (complex heat capacity correction)
c2_bracket = (1.0 / Tth - 1.0 / Trth) * (Theta - T_K) / Theta - T_K / (
    Theta * Theta
) * math.log((Tr / T_K) * (Tth / Trth))
G_term6 = -params_HCO3["c2"] * c2_bracket

G_term7 = (1.0 / Tth) * (
    params_HCO3["a3"] * P_diff + params_HCO3["a4"] * math.log(psiP / psiPr)
)

print(f"\nNon-Born Gibbs components (cal/mol):")
print(f"  Gf = {G_term1:.2f}")
print(f"  -Sr*(T-Tr) = {G_term2:.2f}")
print(f"  -c1*(T*ln(T/Tr) - T + Tr) = {G_term3:.2f}")
print(f"  a1*(P-Pr) = {G_term4:.2f}")
print(f"  a2*ln((Psi+P)/(Psi+Pr)) = {G_term5:.2f}")
print(f"  -c2*[...complex...] = {G_term6:.2f}")
print(f"  (1/(T-Theta))*[a3*(P-Pr) + a4*ln(...)] = {G_term7:.2f}")

G_non_born = G_term1 + G_term2 + G_term3 + G_term4 + G_term5 + G_term6 + G_term7
print(f"\n  Sum (non-Born) = {G_non_born:.2f} cal/mol")

# Convert to J/mol for comparison
G_non_born_J = G_non_born * 4.184
print(f"  Sum (non-Born) = {G_non_born_J:.2f} J/mol")

print("\n" + "=" * 80)
print("EXPECTED vs ACTUAL (from test output)")
print("=" * 80)
print("\nHCO3- at 300°C, 5000 bar:")
print(f"  Truth G0:  -143644 cal/mol (-601056 J/mol)")
print(f"  Model G0:  -144015 cal/mol (-602608 J/mol)")
print(f"  Difference: -371 cal/mol (-1552 J/mol)")
print(f"\n  Non-Born contribution calculated: {G_non_born:.2f} cal/mol")
print(f"  Born term needed: {-143644 - G_non_born:.2f} cal/mol (to match truth)")
print(f"  Born term in model: {-144015 - G_non_born:.2f} cal/mol")
print(f"  Born term difference: {-371:.2f} cal/mol")

print("\n" + "=" * 80)
print("CONCLUSION")
print("=" * 80)
print("""
The error appears to be in the Born omega calculations, which depend on:
1. Water dielectric constant ε(T, P, ρ)
2. Born coefficient Q = (1/ε²)(∂ε/∂P)_T
3. Born omega ω(P, T) via the Shock et al. (1992) g-function

To identify the exact source, we need to:
1. Verify water properties (ρ, ε, Q) match Excel at 300°C, 5000 bar
2. Check the g-function implementation
3. Verify omega calculation and its derivatives
""")
