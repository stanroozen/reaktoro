#!/usr/bin/env python3
import pandas as pd

excel_file = "Latest_DEW2024.xlsm"
xl = pd.ExcelFile(excel_file)
print("Sheet names in", excel_file)
for name in xl.sheet_names:
    print(f"  - {name}")
