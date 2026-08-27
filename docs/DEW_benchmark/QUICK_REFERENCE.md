# Quick Reference: Mineral Solubility Analysis

## 🎯 3-Step Setup for Any Mineral

### Step 1: Edit Configuration (Lines 50-68)
```python
MINERAL_CONFIG = {
    "mineral_name": "YourMineral",        # e.g., "Calcite", "Halite"
    "mineral_formula": "Formula",         # e.g., "CaCO3", "NaCl"
    "solute_species": "Species_aq",       # e.g., "Ca+2", "Na+"
    "aqueous_species": "Complex1 Complex2",  # e.g., "CO3-2 HCO3-"
    "csv_file": "your_data.csv",
    "output_prefix": "yourmineral",
    "plot_title": "Your Mineral Solubility",
    "y_label": "Your Mineral Solubility (mol/kg-H₂O)",
}
```

### Step 2: Prepare CSV File
Must have these columns:
- `T_C` - Temperature (°C)
- `P_kbar` - Pressure (kbar)
- `molality_m` - Measured molality
- `reference` - Citation
- `experiment_type` - Experiment type

### Step 3: Run
```bash
python quartz_solubility_analysis_v2_dew24.py
```

## 📋 Common Minerals Cheat Sheet

| Mineral | mineral_name | solute_species | Key aqueous_species |
|---------|-------------|----------------|---------------------|
| Quartz | `"Quartz"` | `"SiO2_aq"` | `"HSiO3- Si2O4_aq Si3O6_aq"` |
| Calcite | `"Calcite"` | `"Ca+2"` | `"CO2_aq CO3-2 HCO3- CaCO3_aq"` |
| Halite | `"Halite"` | `"Na+"` | `"Cl- NaCl_aq"` |
| Corundum | `"Corundum"` | `"Al+3"` | `"AlOH+2 Al(OH)3_aq Al(OH)4-"` |
| Periclase | `"Periclase"` | `"Mg+2"` | `"MgOH+ Mg(OH)2_aq"` |
| Rutile | `"Rutile"` | `"Ti(OH)4_aq"` | `"TiO2_aq Ti(OH)3+"` |

## 🔍 Find Species Names

```python
from reaktoro import *

# Check mineral name
db_supcrt = SupcrtDatabase("supcrtbl")
print(db_supcrt.species("YourMineral"))

# Find aqueous species
db_dew = DEWDatabase("dew2024-aqueous")
for species in db_dew.species():
    if "Ca" in species.name():  # Search for element
        print(species.name())
```

## 📊 Output Files

Automatically generated with your `output_prefix`:

1. `{prefix}_solubility_comparison_low_P_dew24.png`
2. `{prefix}_solubility_comparison_high_P_dew24.png`
3. `{prefix}_solubility_residuals_dew24.png`

## ⚙️ Advanced Settings

### Temperature/Pressure Range
```python
T_MIN, T_MAX = 150, 550    # Change default range
N_POINTS = 100             # Number of calculation points
```

### Water Model
```python
DEW_CONFIG = {
    "eos_model": "ZhangDuan2005",        # or "WagnerPruss", "HGK"
    "dielectric_model": "PowerFunction",  # or "JohnsonNorton1991"
    "gibbs_model": "DewIntegral",        # or "DelaneyHelgeson1978"
    "born_model": "Shock92Dew",
}
```

## 🐛 Troubleshooting

| Problem | Solution |
|---------|----------|
| "Species not found" | Check exact spelling in database |
| Convergence failures | Verify T-P range is appropriate |
| Empty plots | Check CSV file path and format |
| Wrong solubility | Verify solute_species is correct |

## 📚 Documentation

- **Full Guide**: `MINERAL_SOLUBILITY_CONFIG_GUIDE.md`
- **Template**: `mineral_template.py`
- **Summary**: `STANDARDIZATION_SUMMARY.md`

## 💡 Example: Switch from Quartz to Calcite

**Change this:**
```python
"mineral_name": "Quartz",
"mineral_formula": "SiO2",
"solute_species": "SiO2_aq",
"aqueous_species": "H2_aq O2_aq HO2- HSiO3- Si2O4_aq Si3O6_aq",
"csv_file": "quartz_DEW_testset.csv",
"output_prefix": "quartz",
```

**To this:**
```python
"mineral_name": "Calcite",
"mineral_formula": "CaCO3",
"solute_species": "Ca+2",
"aqueous_species": "CO2_aq HCO3- CO3-2 CaCO3_aq CaHCO3+",
"csv_file": "calcite_DEW_testset.csv",
"output_prefix": "calcite",
```

**That's it! 🎉**

No other code changes needed anywhere!
