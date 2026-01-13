#!/usr/bin/env python3
import pandas as pd

excel_file = "Latest_DEW2024.xlsm"
df = pd.read_excel(excel_file, sheet_name="Aqueous Species Table", header=1)

print("Columns in 'Aqueous Species Table':")
for col in df.columns:
    print(f"  - '{col}'")

print("\nFirst few rows:")
print(df.head())
