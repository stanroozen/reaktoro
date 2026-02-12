import pandas as pd

df = pd.read_csv("quartz_DEW_testset.csv")
morey = df[df["reference"] == "Morey_Fournier_Rowe_1962"]
hemley = df[df["reference"] == "Hemley_1980"]

print("=== Morey_Fournier_Rowe_1962 ===")
print(f"Total: {len(morey)}")
print(f"P < 1 kbar: {len(morey[morey['P_kbar'] < 1.0])}")
print(f"P >= 1 kbar: {len(morey[morey['P_kbar'] >= 1.0])}")
print(f"P NaN: {morey['P_kbar'].isna().sum()}")
print(f"Unique pressures: {sorted(morey['P_kbar'].dropna().unique())}")

print("\n=== Hemley_1980 ===")
print(f"Total: {len(hemley)}")
print(f"P < 1 kbar: {len(hemley[hemley['P_kbar'] < 1.0])}")
print(f"P >= 1 kbar: {len(hemley[hemley['P_kbar'] >= 1.0])}")
print(f"P NaN: {hemley['P_kbar'].isna().sum()}")
print(f"Unique pressures: {sorted(hemley['P_kbar'].dropna().unique())}")

# Show NaN rows
if hemley["P_kbar"].isna().sum() > 0:
    print("\nHemley rows with NaN pressure:")
    print(hemley[hemley["P_kbar"].isna()][["T_C", "P_kbar", "experiment_type"]])
