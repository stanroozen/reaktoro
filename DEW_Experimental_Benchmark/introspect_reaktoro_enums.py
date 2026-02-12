import sys, os

ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
PYD_DIR = os.path.join(os.path.dirname(ROOT_DIR), "build-msvc", "Reaktoro", "Release")

print("Python:", sys.version)

mod = None
try:
    import reaktoro as r

    print("Using module: reaktoro")
    mod = r
except Exception as e:
    print("Failed to import reaktoro:", e)
    sys.path.insert(0, PYD_DIR)
    try:
        import reaktoro4py as r

        print("Using module: reaktoro4py")
        mod = r
    except Exception as e2:
        print("Failed to import reaktoro4py:", e2)
        sys.exit(1)

names = [
    "WaterBornModel",
    "WaterEosModel",
    "WaterDielectricModel",
    "WaterGibbsModel",
    "ActivityModelDEW",
    "StandardThermoModelDEW",
]
for name in names:
    obj = getattr(mod, name, None)
    print(f"\n{name}:")
    if obj is None:
        print("  (missing)")
        continue
    attrs = [a for a in dir(obj) if not a.startswith("_")]
    # Cut long outputs
    if len(attrs) > 0:
        print("  members:", ", ".join(attrs))
    else:
        print("  (no public members)")

# Try accessing potential enum values that often exist
candidates = [
    ("WaterBornModel", ["Shock92", "Shock92Dew", "ReaktoroShock92", "ShockBorn92"]),
    ("WaterEosModel", ["ZhangDuan2005", "ZhangDuan", "IAPWS95"]),
    (
        "WaterDielectricModel",
        ["Fernandez2006", "Fernandez", "WaterDielectricFernandez2006"],
    ),
    (
        "WaterGibbsModel",
        ["WagnerPruss2002", "WagnerPruss", "WaterGibbsWagnerPruss2002"],
    ),
]

for enum_name, values in candidates:
    enum = getattr(mod, enum_name, None)
    print(f"\nProbing {enum_name} values:")
    if enum is None:
        print("  (enum missing)")
        continue
    for v in values:
        ok = hasattr(enum, v)
        print(f"  {v}: {'OK' if ok else 'missing'}")
