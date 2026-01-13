import yaml

with open("dew2024-aqueous.yaml", "r") as f:
    data = yaml.safe_load(f)

co2 = data["Species"]["CO2(0)"]["StandardThermoModel"]["HKF"]

print("CO2(0) HKF parameters from YAML:")
print(f"  a1 = {co2['a1']:.6e} J/(mol·Pa) = m³/mol")
print(f"  a2 = {co2['a2']:.6e} J/mol")
print(f"  a3 = {co2['a3']:.6e} J·K/(mol·Pa)")
print(f"  a4 = {co2['a4']:.6e} J·K/mol")
print()

# Calculate V0 at 300°C, 5000 bar
T = 573.15  # K
P = 5e8  # Pa (5000 bar)
psi = 2.6e8  # Pa (2600 bar)
theta = 228.0  # K

a1, a2, a3, a4 = co2["a1"], co2["a2"], co2["a3"], co2["a4"]

# HKF V formula (without Born terms for neutral species)
V_pv = a1 + a2 / (psi + P) + (a3 + a4 / (psi + P)) / (T - theta)

print(f"Manual V calculation at T={T}K, P={P / 1e8:.0f} kbar:")
print(f"  psi + P = {(psi + P) / 1e8:.2f}e8 Pa")
print(f"  T - theta = {T - theta:.2f} K")
print()
print(f"  Term a1                = {a1 * 1e6:.4f} cm³/mol")
print(f"  Term a2/(psi+P)        = {(a2 / (psi + P)) * 1e6:.4f} cm³/mol")
print(
    f"  Term (a3+a4/(psi+P))/Tth = {((a3 + a4 / (psi + P)) / (T - theta)) * 1e6:.4f} cm³/mol"
)
print(f"  Total V_PV             = {V_pv * 1e6:.4f} cm³/mol")
print()
print(f"Expected from Excel: V° ≈ 33.40 cm³/mol")
print(f"Expected from test: V_FD ≈ 33.55 cm³/mol")
print(f"Actual model output: V_model ≈ 32.57 cm³/mol")
