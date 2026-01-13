import pandas as pd

# Check both Excel files
for excel_file in ["DEW_2019.xlsm", "Latest_DEW2024.xlsm"]:
    print(f"\n{'=' * 70}")
    print(f"Checking: {excel_file}")
    print("=" * 70)

    # Read the Excel file
    df = pd.read_excel(excel_file, sheet_name="Aqueous Species Table", header=None)

    # Get headers from row 2 (index 2)
    headers = df.iloc[2, :].tolist()

    # Find CO2
    for idx in range(len(df)):
        chemical = df.iloc[idx, 1]  # Column 1 is Chemical name
        symbol = df.iloc[idx, 2]  # Column 2 is Symbol
        if pd.notna(chemical) and str(chemical).strip() == "CO2,aq":
            print(f"\nFound CO2,aq at row {idx}:")
            print(f"  Chemical: {chemical}")
            print(f"  Symbol: {symbol}")
            print(f"  ΔGfo: {df.iloc[idx, 3]} cal/mol")
            print(f"  Vo: {df.iloc[idx, 6]} cm³/mol")
            print(f"  a1 x 10: {df.iloc[idx, 8]} cal mol-1 bar-1")
            print(f"  a2 x 10-2: {df.iloc[idx, 9]} cal mol-1")
            print(f"  a3: {df.iloc[idx, 10]} cal K mol-1 bar-1")
            print(f"  a4 x 10-4: {df.iloc[idx, 11]} cal K mol-1")

            # Calculate what should be in YAML
            a1_yaml = df.iloc[idx, 8] * 10.0 * 4.184 / 1e5
            a2_yaml = df.iloc[idx, 9] * 0.01 * 4.184
            a3_yaml = df.iloc[idx, 10] * 4.184 / 1e5
            a4_yaml = df.iloc[idx, 11] * 0.0001 * 4.184

            print(f"\nExpected YAML values:")
            print(f"  a1: {a1_yaml}")
            print(f"  a2: {a2_yaml}")
            print(f"  a3: {a3_yaml}")
            print(f"  a4: {a4_yaml}")
            break

# Now check what's actually in the YAML
print(f"\n{'=' * 70}")
print("Actual YAML values:")
print("=" * 70)
import yaml

for yaml_file in ["dew2019-aqueous.yaml", "dew2024-aqueous.yaml"]:
    with open(yaml_file, "r") as f:
        data = yaml.safe_load(f)

    # Find CO2
    for key, species in data["Species"].items():
        if species["Name"] == "CO2,aq":
            print(f"\n{yaml_file}:")
            print(f"  a1: {species['StandardThermoModel']['HKF']['a1']}")
            print(f"  a2: {species['StandardThermoModel']['HKF']['a2']}")
            print(f"  a3: {species['StandardThermoModel']['HKF']['a3']}")
            print(f"  a4: {species['StandardThermoModel']['HKF']['a4']}")
            break
