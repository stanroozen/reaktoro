# Quick Reference: DEW Water Model Configuration

## The Issue Solved
Previously, water thermodynamics in quartz solubility calculations came from SUPCRT92 (via IAPWS-95), which is limited to 1000 MPa (10 kbar). The DEW model is parameterized for 0-60 kbar, so we needed a way to use Duan & Zang 2005 water EOS instead.

## The Solution
✅ **Now available:** Configurable water models with DEW database integration

## Default Behavior (No Changes Needed)
```python
# Your existing code works as-is - uses Duan & Zang 2005 by default
system = build_system(dew_db, supcrt_db)
molality = calculate_quartz_solubility(300, 1000, dew_db, supcrt_db)
```

**Default Configuration:**
- **EOS:** Duan & Zang 2005 (extends to 60 kbar)
- **Dielectric:** Power Function
- **Gibbs Integration:** Volume-based (DewIntegral)
- **Data Source:** H2O_aq from DEW database (not SUPCRT/IAPWS-95)

## Custom Water Models (For Sensitivity Studies)

### Option 1: Change Globally
```python
from quartz_solubility_analysis import DEW_CONFIG, build_system

# Modify the global config
DEW_CONFIG["eos_model"] = "ZhangDuan2009"
DEW_CONFIG["dielectric_model"] = "JohnsonNorton1991"

# All subsequent calls use this configuration
system = build_system(dew_db, supcrt_db)
```

### Option 2: Override Per-Call
```python
custom_water_model = {
    "eos_model": "ZhangDuan2009",
    "dielectric_model": "JohnsonNorton1991",
    "gibbs_model": "DewIntegral",
    "born_model": "Shock92Dew",
}

# One-off calculation with different water model
molality = calculate_quartz_solubility(
    T_C=300,
    P_bar=5000,  # Up to 50 kbar works now!
    dew_db=dew_db,
    supcrt_db=supcrt_db,
    water_config=custom_water_model
)
```

## Available Water Model Options

### EOS (Water Density Model)
Choose based on pressure range and accuracy requirements:
- **ZhangDuan2005** ← **DEFAULT** (Duan & Zang 2005, covers 0-60 kbar)
- **ZhangDuan2009** (Updated Duan & Zang 2009)
- **WagnerPruss** (Limited to ~1 kbar)
- **HGK** (Holton-Goodwin-Kerley, limited range)

### Dielectric Model
- **PowerFunction** ← **DEFAULT** (Best for DEW compatibility)
- **JohnsonNorton1991** (Alternative correlations)

### Gibbs Integration
- **DewIntegral** ← **DEFAULT** (Volume-based, recommended for DEW)
- **DelaneyHelgeson1978** (Legacy, limited P-T range)

### Born Model (Ionic Contributions)
- **Shock92Dew** ← **DEFAULT** (For charged species)
- **Shock92** (For neutral species only)

## What Changed in the Code

### Database
- `H2O_aq` entry added to DEW database (2019 and 2024 versions)
- This is just a placeholder; actual water EOS is selected via configuration

### Script Functions
- `build_system(dew_db, supcrt_db, water_config=None)` - Now configurable
- `calculate_quartz_solubility(..., water_config=None)` - Now configurable

### Species Names
- ✅ Use `H2O_aq` (from DEW) instead of `H2O(aq)` (SUPCRT)
- ✅ Use `SiO2(aq)` (standard naming)

## Common Workflows

### Workflow 1: Standard Solubility Curve (0.01-10 kbar)
```python
# No changes needed - default config works great
system = build_system(dew_db, supcrt_db)
for T_C in range(150, 551, 50):
    m = calculate_quartz_solubility(T_C, 1000, dew_db, supcrt_db)
```

### Workflow 2: High-Pressure Study (Up to 50 kbar)
```python
# Uses default Duan & Zang 2005 which extends to 60 kbar
for P_kbar in [1, 5, 10, 20, 30, 40, 50]:
    m = calculate_quartz_solubility(300, P_kbar*1000, dew_db, supcrt_db)
    # This now works! (Previously failed above 10 kbar)
```

### Workflow 3: Model Sensitivity Study
```python
models = {
    "DZ2005 (default)": {"eos_model": "ZhangDuan2005"},
    "DZ2009": {"eos_model": "ZhangDuan2009"},
    "WP": {"eos_model": "WagnerPruss"},
}

for model_name, config in models.items():
    results[model_name] = []
    for T in temps:
        m = calculate_quartz_solubility(T, 5000, dew_db, supcrt_db, config)
        results[model_name].append(m)
    print(f"{model_name}: {results[model_name]}")
```

## Troubleshooting

### "Could not find H2O_aq"
→ Library needs to be rebuilt to embed the updated database
→ Run: `cmake --build build-msvc --config Release`

### Water models not taking effect
→ Check that you're passing `water_config` parameter
→ Verify keys match exactly: `"eos_model"`, `"dielectric_model"`, etc.

### Still getting IAPWS-95 errors above 10 kbar
→ Make sure you're using the rebuilt library
→ Verify H2O_aq species is in system (check output log)

## Architecture Summary

```
┌─────────────────────────────────────────────────┐
│ User Code                                       │
│ - build_system(dew_db, supcrt_db, config)      │
│ - calculate_quartz_solubility(..., config)     │
└──────────────────┬──────────────────────────────┘
                   │
                   ▼
┌─────────────────────────────────────────────────┐
│ DEW_CONFIG / User Config                        │
│ - eos_model, dielectric_model, gibbs_model     │
│ - Maps strings to enum values                  │
└──────────────────┬──────────────────────────────┘
                   │
                   ▼
┌─────────────────────────────────────────────────┐
│ StandardThermoModelParamsDEW                    │
│ - waterOptions.eosModel = ZhangDuan2005        │
│ - waterOptions.dielectricModel = PowerFunction  │
│ - waterOptions.gibbsModel = DewIntegral        │
└──────────────────┬──────────────────────────────┘
                   │
                   ▼
┌─────────────────────────────────────────────────┐
│ ActivityModelDEW                                │
│ - Species: HKF model (from DEW database)       │
│ - Solvent: Configured water EOS (NOT IAPWS-95)│
│ - Pressure range: 0-60 kbar (with DZ2005)     │
└─────────────────────────────────────────────────┘
```

## Key Takeaway
✨ **You now have full control over water thermodynamics while maintaining DEW species model consistency and extending to high pressures!**
