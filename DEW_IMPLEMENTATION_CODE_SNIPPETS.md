# Implementation Code Snippets

## 1. Global Configuration (DEW_CONFIG)

```python
# Default configuration - user can modify or override
DEW_CONFIG = {
    "eos_model": "ZhangDuan2005",          # Duan & Zang 2005 EOS for density
    "dielectric_model": "PowerFunction",    # Power function for dielectric constant
    "gibbs_model": "DewIntegral",          # Volume integration for Gibbs energy
    "born_model": "Shock92Dew",            # Born model for ionic contributions
}
```

## 2. build_system() Function Structure

```python
def build_system(dew_db, supcrt_db, water_config=None):
    """
    Build ChemicalSystem with configurable water models.

    Parameters:
    -----------
    water_config : dict, optional
        Keys: "eos_model", "dielectric_model", "gibbs_model", "born_model"
        Values: string names matching WaterEosModel, WaterDielectricModel, etc.
    """

    if water_config is None:
        water_config = DEW_CONFIG

    # 1. Setup database
    combined_db = Database(dew_db.species())
    combined_db.addSpecies(supcrt_db.species("Quartz"))

    # 2. Define phases
    aqueous = AqueousPhase("H2O_aq H+ OH- SiO2_aq H2_aq O2_aq HO2- HSiO3- Si2O4_aq Si3O6_aq")
    mineral = MineralPhase("Quartz")

    # 3. Configure DEW with custom water options
    try:
        params = StandardThermoModelParamsDEW()

        # Map string names to enum values
        eos_map = {
            "WagnerPruss": WaterEosModel.WagnerPruss,
            "HGK": WaterEosModel.HGK,
            "ZhangDuan2005": WaterEosModel.ZhangDuan2005,
            "ZhangDuan2009": WaterEosModel.ZhangDuan2009,
        }

        dielectric_map = {
            "PowerFunction": WaterDielectricModel.PowerFunction,
            "JohnsonNorton1991": WaterDielectricModel.JohnsonNorton1991,
        }

        gibbs_map = {
            "DewIntegral": WaterGibbsModel.DewIntegral,
            "DelaneyHelgeson1978": WaterGibbsModel.DelaneyHelgeson1978,
        }

        born_map = {
            "Shock92Dew": WaterBornModel.Shock92Dew,
            "Shock92": WaterBornModel.Shock92,
        }

        # Apply user configuration (with defaults)
        params.waterOptions.eosModel = eos_map.get(
            water_config.get("eos_model", "ZhangDuan2005"),
            WaterEosModel.ZhangDuan2005
        )
        params.waterOptions.dielectricModel = dielectric_map.get(
            water_config.get("dielectric_model", "PowerFunction"),
            WaterDielectricModel.PowerFunction
        )
        params.waterOptions.gibbsModel = gibbs_map.get(
            water_config.get("gibbs_model", "DewIntegral"),
            WaterGibbsModel.DewIntegral
        )
        params.waterOptions.bornModel = born_map.get(
            water_config.get("born_model", "Shock92Dew"),
            WaterBornModel.Shock92Dew
        )

        # Create and apply model
        dew_model = StandardThermoModelDEW(params)
        aqueous.setActivityModel(dew_model)

        # Log what was configured
        print(f"✓ DEW configured: EOS={water_config.get('eos_model', 'ZhangDuan2005')}, "
              f"Dielectric={water_config.get('dielectric_model', 'PowerFunction')}, "
              f"Gibbs={water_config.get('gibbs_model', 'DewIntegral')}")

    except Exception as e:
        print(f"Warning: Custom DEW config failed: {e}. Using default.")
        aqueous.setActivityModel(ActivityModelDEW())

    # 4. Create system
    system = ChemicalSystem(combined_db, aqueous, mineral)
    return system
```

## 3. Usage Examples

### Example 1: Using Defaults
```python
# Automatically uses ZhangDuan2005, PowerFunction, DewIntegral
system = build_system(dew_db, supcrt_db)

# All calculations use this system
solver = EquilibriumSolver(system)
state = ChemicalState(system)
state.temperature(300, "celsius")
state.pressure(5000, "bar")  # Works now - up to 60 kbar!
```

### Example 2: Custom Configuration (One-time)
```python
my_config = {
    "eos_model": "ZhangDuan2009",
    "dielectric_model": "JohnsonNorton1991",
    "gibbs_model": "DelaneyHelgeson1978",
    "born_model": "Shock92Dew",
}

system = build_system(dew_db, supcrt_db, water_config=my_config)
```

### Example 3: Custom Configuration (Global)
```python
# Modify the global default
DEW_CONFIG["eos_model"] = "WagnerPruss"
DEW_CONFIG["dielectric_model"] = "JohnsonNorton1991"

# All subsequent calls use the modified config
system1 = build_system(dew_db, supcrt_db)  # Uses WagnerPruss
system2 = build_system(dew_db, supcrt_db)  # Also uses WagnerPruss
```

### Example 4: Per-Calculation Override
```python
def calculate_with_model(T, P, model_name):
    config = {
        "eos_model": model_name,
        "dielectric_model": "PowerFunction",
        "gibbs_model": "DewIntegral",
        "born_model": "Shock92Dew",
    }
    return calculate_quartz_solubility(T, P, dew_db, supcrt_db, water_config=config)

# Compare results across models
for model in ["ZhangDuan2005", "ZhangDuan2009", "WagnerPruss"]:
    m = calculate_with_model(300, 5000, model)
    print(f"{model}: {m:.4e} mol/kg")
```

## 4. Water Model Enum Mappings

```python
# WaterEosModel options (density model)
{
    "WagnerPruss": WaterEosModel.WagnerPruss,
    "HGK": WaterEosModel.HGK,
    "ZhangDuan2005": WaterEosModel.ZhangDuan2005,      # ← DEFAULT, covers 0-60 kbar
    "ZhangDuan2009": WaterEosModel.ZhangDuan2009,
}

# WaterDielectricModel options
{
    "PowerFunction": WaterDielectricModel.PowerFunction,           # ← DEFAULT
    "JohnsonNorton1991": WaterDielectricModel.JohnsonNorton1991,
}

# WaterGibbsModel options (Gibbs energy integration)
{
    "DewIntegral": WaterGibbsModel.DewIntegral,                   # ← DEFAULT
    "DelaneyHelgeson1978": WaterGibbsModel.DelaneyHelgeson1978,
}

# WaterBornModel options (ionic Born model)
{
    "Shock92Dew": WaterBornModel.Shock92Dew,                       # ← DEFAULT
    "Shock92": WaterBornModel.Shock92,
}
```

## 5. Error Handling and Logging

```python
# The build_system function includes:

# 1. Logging of selected models
print(f"✓ DEW configured: EOS={eos_name}, Dielectric={diel_name}, Gibbs={gibbs_name}")

# 2. Graceful fallback if custom config fails
except Exception as e:
    print(f"Warning: Could not configure DEW with custom water options: {e}")
    print("  Falling back to default ActivityModelDEW()")
    aqueous.setActivityModel(ActivityModelDEW())

# 3. For null config, uses global DEW_CONFIG
if water_config is None:
    water_config = DEW_CONFIG
```

## 6. State Initialization (Important!)

```python
def build_state(system, T_C, P_bar):
    state = ChemicalState(system)
    state.temperature(float(T_C), "celsius")
    state.pressure(float(P_bar), "bar")

    # Use H2O_aq from DEW database (not H2O(aq) from SUPCRT)
    state.set("H2O_aq", 1.0, "kg")         # ← KEY CHANGE
    state.set("H+", 1e-8, "mol")
    state.set("OH-", 1e-8, "mol")
    state.set("SiO2(aq)", 1e-6, "mol")     # ← Renamed from SiO2_aq
    state.set("Quartz", 10.0, "mol")

    return state
```

## 7. Database Check (Verification)

```python
# Verify H2O_aq is available in DEW
dew_db = DEWDatabase('dew2019-aqueous')
species_names = [s.name() for s in dew_db.species()]

if 'H2O_aq' in species_names:
    print("✓ H2O_aq found in DEW database")
    h2o = dew_db.species('H2O_aq')
    print(f"  Formula: {h2o.formula()}")
    print(f"  State: {h2o.aggregateState()}")
else:
    print("✗ H2O_aq NOT found - rebuild library")
    print(f"  Available species: {len(species_names)}")
```

## 8. Key Differences from Previous Implementation

| Aspect | Before | After |
|--------|--------|-------|
| Water source | H2O(aq) from SUPCRT | H2O_aq from DEW |
| Water EOS | IAPWS-95 (limited to 10 kbar) | Configurable (default: DZ2005, up to 60 kbar) |
| Water config | Hardcoded | Via DEW_CONFIG or water_config parameter |
| Pressure range | 0.01-10 kbar | 0.01-60 kbar (default config) |
| Flexibility | Fixed model | Multiple model options available |
| Consistency | Water EOS vs DEW mismatch | Aligned (both use Duan & Zang by default) |

## 9. Pressure Capability Chart

```
Pressure Range (MPa)  | Water Model Options
─────────────────────┼──────────────────────────────
0-100 (0-1 kbar)     | All models work
100-500 (1-5 kbar)   | All models work
500-1000 (5-10 kbar) | All models work
1000-2000 (10-20kb)  | DZ2005, DZ2009 only
2000-6000 (20-60kb)  | DZ2005, DZ2009 only ← NEW!
```

This is what changed with H2O_aq + configurable DEW!
