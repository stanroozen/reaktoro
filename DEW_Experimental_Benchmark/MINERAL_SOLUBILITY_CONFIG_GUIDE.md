# Example Mineral Configurations for Solubility Analysis

This document shows how to adapt `quartz_solubility_analysis_v2_dew24.py` for different minerals.

## How to Use

Simply modify the `MINERAL_CONFIG` dictionary at the top of the script (lines 47-68).

## Example Configurations

### 1. Quartz (SiO2) - Current Default

```python
MINERAL_CONFIG = {
    "mineral_name": "Quartz",
    "mineral_formula": "SiO2",
    "solute_species": "SiO2_aq",
    "aqueous_species": "H2_aq O2_aq HO2- HSiO3- Si2O4_aq Si3O6_aq",
    "csv_file": "quartz_DEW_testset.csv",
    "output_prefix": "quartz",
    "plot_title": "Quartz Solubility",
    "y_label": "Quartz Solubility (mol/kg-H₂O)",
}
```

### 2. Calcite (CaCO3)

```python
MINERAL_CONFIG = {
    "mineral_name": "Calcite",
    "mineral_formula": "CaCO3",
    "solute_species": "Ca+2",              # Primary dissolved ion
    "aqueous_species": "CO2_aq HCO3- CO3-2 CaCO3_aq CaHCO3+ CaOH+",
    "csv_file": "calcite_DEW_testset.csv",
    "output_prefix": "calcite",
    "plot_title": "Calcite Solubility",
    "y_label": "Calcite Solubility (mol/kg-H₂O)",
}
```

### 3. Halite (NaCl)

```python
MINERAL_CONFIG = {
    "mineral_name": "Halite",
    "mineral_formula": "NaCl",
    "solute_species": "Na+",               # Primary dissolved ion
    "aqueous_species": "Cl- NaCl_aq NaOH_aq",
    "csv_file": "halite_DEW_testset.csv",
    "output_prefix": "halite",
    "plot_title": "Halite Solubility",
    "y_label": "Halite Solubility (mol/kg-H₂O)",
}
```

### 4. Corundum (Al2O3)

```python
MINERAL_CONFIG = {
    "mineral_name": "Corundum",
    "mineral_formula": "Al2O3",
    "solute_species": "Al+3",
    "aqueous_species": "AlOH+2 Al(OH)2+ Al(OH)3_aq Al(OH)4- AlO2-",
    "csv_file": "corundum_DEW_testset.csv",
    "output_prefix": "corundum",
    "plot_title": "Corundum Solubility",
    "y_label": "Corundum Solubility (mol/kg-H₂O)",
}
```

### 5. Periclase (MgO)

```python
MINERAL_CONFIG = {
    "mineral_name": "Periclase",
    "mineral_formula": "MgO",
    "solute_species": "Mg+2",
    "aqueous_species": "MgOH+ Mg(OH)2_aq",
    "csv_file": "periclase_DEW_testset.csv",
    "output_prefix": "periclase",
    "plot_title": "Periclase Solubility",
    "y_label": "Periclase Solubility (mol/kg-H₂O)",
}
```

## Configuration Parameters Explained

| Parameter | Description | Example |
|-----------|-------------|---------|
| `mineral_name` | Mineral name as it appears in the SUPCRT database | `"Quartz"` |
| `mineral_formula` | Chemical formula for display purposes | `"SiO2"` |
| `solute_species` | Primary aqueous species to track solubility | `"SiO2_aq"` |
| `aqueous_species` | Additional aqueous species to include (space-separated) | `"HSiO3- Si2O4_aq"` |
| `csv_file` | Experimental data CSV filename | `"quartz_DEW_testset.csv"` |
| `output_prefix` | Prefix for output figure filenames | `"quartz"` |
| `plot_title` | Title for plots | `"Quartz Solubility"` |
| `y_label` | Y-axis label for plots | `"Quartz Solubility (mol/kg-H₂O)"` |

## Required CSV File Format

Your experimental data CSV file must contain these columns:

- `T_C`: Temperature in Celsius
- `P_kbar`: Pressure in kbar
- `molality_m`: Measured molality (mol/kg H₂O)
- `reference`: Citation/reference for the data
- `experiment_type`: Type of experiment (e.g., "hydrothermal", "saturation")

## Steps to Analyze a New Mineral

1. **Prepare experimental data CSV** with the required columns
2. **Check species names** in the DEW and SUPCRT databases:
   ```python
   from reaktoro import *
   db_dew = DEWDatabase("dew2024-aqueous")
   db_supcrt = SupcrtDatabase("supcrtbl")

   # Check if your mineral exists
   print(db_supcrt.species("YourMineralName"))

   # Check aqueous species
   for species in db_dew.species():
       print(species.name())
   ```
3. **Update MINERAL_CONFIG** in the script
4. **Run the analysis**:
   ```bash
   python quartz_solubility_analysis_v2_dew24.py
   ```

## Tips for Choosing Aqueous Species

- **Simple minerals**: Include the primary ion and common complexes
- **Check speciation**: Look at dominant species at your T-P conditions
- **DEW database**: Use species available in `dew2024-aqueous`
- **Start minimal**: Begin with just the primary species, add complexes if needed

## Output Files

The script will generate three figures with your `output_prefix`:

1. `{output_prefix}_solubility_comparison_low_P_dew24.png` - Low pressure (<1 kbar)
2. `{output_prefix}_solubility_comparison_high_P_dew24.png` - High pressure (≥1 kbar)
3. `{output_prefix}_solubility_residuals_dew24.png` - Model-experiment residuals

## Common Issues

### Mineral not found
**Error**: `"Species not found: YourMineralName"`
**Solution**: Check the exact spelling in SUPCRTBL database. Try variations like "Calcite" vs "calcite"

### Species not found
**Error**: `"Species not found: YourSpecies_aq"`
**Solution**: Verify the species exists in DEW2024 database. Some complex species may not be included.

### Convergence failures
**Warning**: Many equilibrium calculations fail
**Solution**:
- Check if T-P range is appropriate for your mineral
- Verify aqueous species are compatible
- Adjust initial state amounts if needed

## Need Help?

- Check available species: Look at database documentation
- Test single point: Run one T-P calculation to debug
- Simplify first: Start with minimal species list, add complexity gradually
