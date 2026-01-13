import pandas as pd

# Read the Excel file
df = pd.read_excel("DEW_2019.xlsm", sheet_name="Aqueous Species Table", header=None)

# Get headers from row 2 (index 2)
headers = df.iloc[2, :].tolist()
print("Column headers (row 2):")
for i, h in enumerate(headers[:20]):
    print(f"  Col {i}: {h}")

# Get units from row 3
units = df.iloc[3, :].tolist()
print("\nColumn units (row 3):")
for i, u in enumerate(units[:20]):
    print(f"  Col {i}: {u}")

# Find CO2
print("\n" + "=" * 60)
print("Searching for CO2...")
for idx in range(len(df)):
    chemical = df.iloc[idx, 1]  # Column 1 is Chemical name
    symbol = df.iloc[idx, 2]  # Column 2 is Symbol
    if pd.notna(chemical) and "CO2" in str(chemical):
        print(f"\nFound at row {idx}:")
        print(f"  Chemical: {chemical}")
        print(f"  Symbol: {symbol}")
        print("Values:")
        for i in range(17):
            if pd.notna(df.iloc[idx, i]) and i < len(headers):
                header = headers[i]
                print(f"  {header}: {df.iloc[idx, i]}")
