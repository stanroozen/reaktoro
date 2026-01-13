#!/usr/bin/env python3
"""
Check the raw values in Excel for CO2 to understand scaling.
"""

import pandas as pd

# Read Excel file - use header=None and manually select header row
excel_file = "Latest_DEW2024.xlsm"
raw = pd.read_excel(excel_file, sheet_name="Aqueous Species Table", header=None)

# row 2 has the header names
header_names = raw.iloc[2]
# data starts at row 4
df = raw.iloc[4:].copy()
df.columns = header_names

# Drop rows without a Chemical name
df = df[df["Chemical"].notna()].copy()

print("Columns available:")
for i, col in enumerate(df.columns):
    print(f"  {i}: '{col}'")
print()

# Find CO2(0) - look in Chemical column
co2_row = df[df["Chemical"].str.contains("CO2", case=False, na=False)]
if co2_row.empty:
    print("CO2 not found!")
    print("\nFirst few species:")
    print(df["Chemical"].head(10).tolist())
else:
    print(f"Found {len(co2_row)} CO2 species:")
    for idx, row_data in co2_row.iterrows():
        print(f"  - {row_data['Chemical']}")

    # Use the first one
    row = co2_row.iloc[0]
    print(f"\nUsing: {row['Chemical']}")
    print()
    row = co2_row.iloc[0]

    print("=" * 70)
    print("RAW EXCEL VALUES FOR CO2(0)")
    print("=" * 70)
    print()

    # Show the raw values as they appear in Excel
    print(f"a1 x 10        = {row['a1 x 10']}")
    print(f"a2 x 10-2      = {row['a2 x 10-2']}")
    print(f"a3             = {row['a3']}")
    print(f"a4 x 10-4      = {row['a4 x 10-4']}")
    print(f"c1             = {row['c1']}")
    print(f"c2 x 10-4      = {row['c2 x 10-4']}")
    print(f"ω x 10-5       = {row['ω x 10-5']}")
    print()

    print("=" * 70)
    print("INTERPRETATION (what the actual parameters should be):")
    print("=" * 70)
    print()

    # Apply the scaling factors
    CAL2J = 4.184
    BAR2PA = 1.0e5

    a1_raw = float(row["a1 x 10"])
    a1_actual = a1_raw * 10.0  # Header says "× 10", so multiply by 10
    a1_SI = a1_actual * CAL2J / BAR2PA

    a2_raw = float(row["a2 x 10-2"])
    a2_actual = a2_raw * 1.0e-2  # Header says "× 10⁻²"
    a2_SI = a2_actual * CAL2J

    omega_raw = float(row["ω x 10-5"])
    omega_actual = omega_raw * 1.0e-5  # Header says "× 10⁻⁵"
    omega_SI = omega_actual * CAL2J

    print(f"a₁:")
    print(f"  Excel shows:    {a1_raw} (in column 'a1 x 10')")
    print(f"  Actual a₁:      {a1_actual} cal/(mol·bar)")
    print(f"  In SI units:    {a1_SI:.6e} J/(mol·Pa)")
    print()

    print(f"a₂:")
    print(f"  Excel shows:    {a2_raw} (in column 'a2 x 10-2')")
    print(f"  Actual a₂:      {a2_actual} cal/mol")
    print(f"  In SI units:    {a2_SI:.6e} J/mol")
    print()

    print(f"ω:")
    print(f"  Excel shows:    {omega_raw} (in column 'ω x 10-5')")
    print(f"  Actual ω:       {omega_actual} cal/mol")
    print(f"  In SI units:    {omega_SI:.6e} J/mol")
    print()
