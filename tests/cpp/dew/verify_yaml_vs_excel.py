#!/usr/bin/env python3
"""Verify base thermodynamic properties in YAML match Excel image"""

import yaml

# Load YAML database
with open(
    r"C:\Users\stanroozen\Documents\Projects\reaktoro-dev\reaktoro\embedded\databases\DEW\dew2024-aqueous.yaml",
    "r",
) as f:
    data = yaml.safe_load(f)

# Expected values from Excel image (with units as shown)
excel_data = {
    "CO2(0)": {  # YAML name is CO2(0), not CO2_aq
        "ΔGf°": -92200,  # cal/mol
        "ΔHf°": None,  # not shown
        "S°": 35.800,  # cal/(mol·K)
        "V°": 30.0,  # cm³/mol
        "Cp": 16.9,  # cal/(mol·K)
        "a1×10": 6.9630,  # cal/(mol·bar) ALREADY ×10!
        "a2×10⁻²": 2.3463,  # cal/mol ALREADY ×10⁻²
        "a3": 2.8707,  # cal·K/(mol·bar)
        "a4×10⁻⁴": -2.8760,  # cal·K/mol ALREADY ×10⁻⁴
        "c1": 8.7008,  # cal/(K·mol)
        "c2×10⁻⁴": 0.4079,  # cal·K/mol ALREADY ×10⁻⁴
        "ω×10⁻⁵": -0.8000,  # cal/mol ALREADY ×10⁻⁵
        "Z": 0,
    },
    "H(+)": {  # YAML name is H(+)
        "ΔGf°": 0,
        "ΔHf°": 0,
        "S°": 0.000,
        "V°": 0.0,
        "Cp": 0.0,
        "a1×10": 0.0000,
        "a2×10⁻²": 0.0000,
        "a3": 0.0000,
        "a4×10⁻⁴": 0.0000,
        "c1": 0.0000,
        "c2×10⁻⁴": 0.0000,
        "ω×10⁻⁵": 0.0000,
        "Z": 1,
    },
    "HCO3(-)": {
        "ΔGf°": -140282,
        "ΔHf°": -164898,
        "S°": 23.530,
        "V°": 24.2,
        "Cp": -8.5,
        "a1×10": 7.6500,
        "a2×10⁻²": 0.9200,
        "a3": 0.6000,
        "a4×10⁻⁴": -2.8200,
        "c1": 11.0000,
        "c2×10⁻⁴": -3.8000,
        "ω×10⁻⁵": 1.2733,
        "Z": -1,
    },
}

# Unit conversion constants
CAL_TO_J = 4.184
BAR_TO_PA = 1.0e5

print("=" * 100)
print("VERIFICATION: YAML Database vs Excel Image")
print("=" * 100)

for species_name, excel_vals in excel_data.items():
    print(f"\n{'=' * 100}")
    print(f"Species: {species_name}")
    print(f"{'=' * 100}")

    # Find species in YAML
    yaml_species = data["Species"].get(species_name)
    if not yaml_species:
        print(f"ERROR: Species {species_name} not found in YAML!")
        continue

    hkf = yaml_species["StandardThermoModel"]["HKF"]

    # Check Gf (cal/mol -> J/mol)
    yaml_Gf_cal = hkf["Gf"] / CAL_TO_J
    excel_Gf = excel_vals["ΔGf°"]
    print(f"\nΔGf°:")
    print(f"  Excel:  {excel_Gf:12.4f} cal/mol")
    print(f"  YAML:   {yaml_Gf_cal:12.4f} cal/mol")
    print(f"  Match:  {'✓' if abs(yaml_Gf_cal - excel_Gf) < 0.1 else '✗ MISMATCH'}")

    # Check Hf if available
    if (
        excel_vals["ΔHf°"] is not None
        and not str(hkf.get("Hf", "nan")).lower() == "nan"
    ):
        yaml_Hf_cal = hkf["Hf"] / CAL_TO_J
        excel_Hf = excel_vals["ΔHf°"]
        print(f"\nΔHf°:")
        print(f"  Excel:  {excel_Hf:12.4f} cal/mol")
        print(f"  YAML:   {yaml_Hf_cal:12.4f} cal/mol")
        print(f"  Match:  {'✓' if abs(yaml_Hf_cal - excel_Hf) < 0.1 else '✗ MISMATCH'}")

    # Check S°
    yaml_S_cal = hkf["Sr"] / CAL_TO_J
    excel_S = excel_vals["S°"]
    print(f"\nS°:")
    print(f"  Excel:  {excel_S:12.4f} cal/(mol·K)")
    print(f"  YAML:   {yaml_S_cal:12.4f} cal/(mol·K)")
    print(f"  Match:  {'✓' if abs(yaml_S_cal - excel_S) < 0.1 else '✗ MISMATCH'}")

    # Check c1
    yaml_c1_cal = hkf["c1"] / CAL_TO_J
    excel_c1 = excel_vals["c1"]
    print(f"\nc1:")
    print(f"  Excel:  {excel_c1:12.4f} cal/(K·mol)")
    print(f"  YAML:   {yaml_c1_cal:12.4f} cal/(K·mol)")
    print(f"  Match:  {'✓' if abs(yaml_c1_cal - excel_c1) < 0.1 else '✗ MISMATCH'}")

    # Check c2×10⁻⁴ (Excel value is ALREADY ×10⁻⁴)
    yaml_c2_cal = hkf["c2"] / CAL_TO_J
    excel_c2_scaled = excel_vals["c2×10⁻⁴"]
    excel_c2_actual = excel_c2_scaled * 1e-4  # Convert from ×10⁻⁴ to actual value
    print(f"\nc2×10⁻⁴:")
    print(f"  Excel:  {excel_c2_scaled:12.4f} (×10⁻⁴ cal·K/mol)")
    print(f"  Excel actual: {excel_c2_actual:12.4e} cal·K/mol")
    print(f"  YAML:   {yaml_c2_cal:12.4e} cal·K/mol")
    print(f"  YAML scaled: {yaml_c2_cal * 1e4:12.4f} (×10⁻⁴)")
    if abs(excel_c2_actual) > 1e-10:
        print(
            f"  Match:  {'✓' if abs(yaml_c2_cal - excel_c2_actual) / abs(excel_c2_actual) < 0.01 else '✗ MISMATCH'}"
        )
    else:
        print(f"  Match:  {'✓' if abs(yaml_c2_cal) < 1e-10 else '✗ MISMATCH'}")

    # Check ω×10⁻⁵ (Excel value is ALREADY ×10⁻⁵)
    yaml_w_cal = hkf["wref"] / CAL_TO_J
    excel_w_scaled = excel_vals["ω×10⁻⁵"]
    excel_w_actual = excel_w_scaled * 1e-5
    print(f"\nω×10⁻⁵:")
    print(f"  Excel:  {excel_w_scaled:12.4f} (×10⁻⁵ cal/mol)")
    print(f"  Excel actual: {excel_w_actual:12.4e} cal/mol")
    print(f"  YAML:   {yaml_w_cal:12.4e} cal/mol")
    print(f"  YAML scaled: {yaml_w_cal * 1e5:12.4f} (×10⁻⁵)")
    if abs(excel_w_actual) > 1e-10:
        print(
            f"  Match:  {'✓' if abs(yaml_w_cal - excel_w_actual) / abs(excel_w_actual) < 0.01 else '✗ MISMATCH'}"
        )
    else:
        print(f"  Match:  {'✓' if abs(yaml_w_cal) < 1e-10 else '✗ MISMATCH'}")

    # Check charge
    yaml_Z = hkf["charge"]
    excel_Z = excel_vals["Z"]
    print(f"\nCharge (Z):")
    print(f"  Excel:  {excel_Z}")
    print(f"  YAML:   {yaml_Z}")
    print(f"  Match:  {'✓' if yaml_Z == excel_Z else '✗ MISMATCH'}")

    # Check HKF parameters
    print(f"\nHKF Volume Parameters:")

    # a1×10 (Excel is ALREADY ×10, need to DIVIDE by 10)
    excel_a1_scaled = excel_vals["a1×10"]
    excel_a1_actual_cal_bar = excel_a1_scaled / 10.0  # Divide by 10!
    excel_a1_actual_SI = excel_a1_actual_cal_bar * CAL_TO_J / BAR_TO_PA
    yaml_a1 = hkf["a1"]
    print(f"  a1×10 (Excel): {excel_a1_scaled:12.4f}")
    print(f"  a1 actual (Excel): {excel_a1_actual_cal_bar:12.4e} cal/(mol·bar)")
    print(f"  a1 expected SI: {excel_a1_actual_SI:12.4e} J/(mol·Pa) = m³/mol")
    print(f"  a1 YAML:        {yaml_a1:12.4e} J/(mol·Pa) = m³/mol")
    if abs(excel_a1_actual_SI) > 1e-10:
        print(
            f"  Match:  {'✓' if abs((yaml_a1 - excel_a1_actual_SI) / excel_a1_actual_SI) < 0.001 else '✗ MISMATCH'}"
        )
    else:
        print(f"  Match:  {'✓' if abs(yaml_a1) < 1e-10 else '✗ MISMATCH'}")

    # a2×10⁻² (Excel is ALREADY ×10⁻², multiply by 10⁻²)
    excel_a2_scaled = excel_vals["a2×10⁻²"]
    excel_a2_actual_cal = excel_a2_scaled * 1e-2
    excel_a2_actual_SI = excel_a2_actual_cal * CAL_TO_J
    yaml_a2 = hkf["a2"]
    print(f"  a2×10⁻² (Excel): {excel_a2_scaled:12.4f}")
    print(f"  a2 actual (Excel): {excel_a2_actual_cal:12.4e} cal/mol")
    print(f"  a2 expected SI: {excel_a2_actual_SI:12.4e} J/mol")
    print(f"  a2 YAML:          {yaml_a2:12.4e} J/mol")
    if abs(excel_a2_actual_SI) > 1e-10:
        print(
            f"  Match:  {'✓' if abs((yaml_a2 - excel_a2_actual_SI) / excel_a2_actual_SI) < 0.001 else '✗ MISMATCH'}"
        )
    else:
        print(f"  Match:  {'✓' if abs(yaml_a2) < 1e-10 else '✗ MISMATCH'}")

    # a3 (no scaling)
    excel_a3_cal_bar = excel_vals["a3"]
    excel_a3_actual_SI = excel_a3_cal_bar * CAL_TO_J / BAR_TO_PA
    yaml_a3 = hkf["a3"]
    print(f"  a3 (Excel): {excel_a3_cal_bar:12.4f} cal·K/(mol·bar)")
    print(f"  a3 expected SI: {excel_a3_actual_SI:12.4e} J·K/(mol·Pa)")
    print(f"  a3 YAML:          {yaml_a3:12.4e} J·K/(mol·Pa)")
    if abs(excel_a3_actual_SI) > 1e-10:
        print(
            f"  Match:  {'✓' if abs((yaml_a3 - excel_a3_actual_SI) / excel_a3_actual_SI) < 0.001 else '✗ MISMATCH'}"
        )
    else:
        print(f"  Match:  {'✓' if abs(yaml_a3) < 1e-10 else '✗ MISMATCH'}")

    # a4×10⁻⁴ (Excel is ALREADY ×10⁻⁴)
    excel_a4_scaled = excel_vals["a4×10⁻⁴"]
    excel_a4_actual_cal = excel_a4_scaled * 1e-4
    excel_a4_actual_SI = excel_a4_actual_cal * CAL_TO_J
    yaml_a4 = hkf["a4"]
    print(f"  a4×10⁻⁴ (Excel): {excel_a4_scaled:12.4f}")
    print(f"  a4 actual (Excel): {excel_a4_actual_cal:12.4e} cal·K/mol")
    print(f"  a4 expected SI: {excel_a4_actual_SI:12.4e} J·K/mol")
    print(f"  a4 YAML:         {yaml_a4:12.4e} J·K/mol")
    if abs(excel_a4_actual_SI) > 1e-10:
        print(
            f"  Match:  {'✓' if abs((yaml_a4 - excel_a4_actual_SI) / excel_a4_actual_SI) < 0.001 else '✗ MISMATCH'}"
        )
    else:
        print(f"  Match:  {'✓' if abs(yaml_a4) < 1e-10 else '✗ MISMATCH'}")

print(f"\n{'=' * 100}")
print("VERIFICATION COMPLETE")
print(f"{'=' * 100}")
