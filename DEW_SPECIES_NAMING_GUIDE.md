# DEW vs SUPCRT Database - Silica Species Comparison

## Executive Summary

**The issue:** Your Python script is trying to use `H4SiO4(aq)` with the DEW database, but this species does NOT exist in DEW. The DEW database uses `SiO2_aq` for neutral dissolved silica.

## Key Findings

### 1. Species Naming Conventions

**DEW Database:**
- Commas in species names are replaced with underscores
- Example: `SiO2,aq` (original) → `SiO2_aq` (in DEW database)
- Example: `NaHSiO3,aq` → `NaHSiO3_aq`

**SUPCRT Database:**
- Uses standard aqueous nomenclature with parentheses
- Example: `H4SiO4(aq)`, `SiO2(aq)`, `HSiO3-`

### 2. Available Silica Species

#### DEW2019-aqueous Database (10 species)

| Species Name | Formula | Type |
|-------------|---------|------|
| **SiO2_aq** | SiO2(0) | **Neutral dissolved silica (primary)** |
| Si2O4_aq | Si2O4(0) | Dimeric silica |
| Si3O6_aq | Si3O6(0) | Trimeric silica |
| HSiO3- | HSiO3(-) | Deprotonated silicic acid |
| NaHSiO3_aq | NaHSiO3(0) | Sodium silicate complex |
| Ca(HSiO3)+ | Ca(HSiO3)+ | Calcium silicate complex |
| Mg(HSiO3)+ | Mg(HSiO3)(+) | Magnesium silicate complex |
| Fe(HSiO3)+ | Fe(HSiO3)+ | Iron silicate complex |
| AlO2(SiO2)- | AlO2(SiO2)- | Aluminum silicate complex |
| Mg(SiO2)(HCO3)+ | MgSiC+ | Magnesium carbonate silicate |

#### SUPCRT (supcrtbl) Database (8 aqueous silica species)

| Species Name | Formula | Type |
|-------------|---------|------|
| **H4SiO4(aq)** | H4SiO4 | **Silicic acid (primary)** |
| SiO2(aq) | SiO2 | Dissolved silica |
| HSiO3- | HSiO3- | Deprotonated silicic acid |
| NaHSiO3(aq) | NaHSiO3 | Sodium silicate complex |
| Ca(HSiO3)+ | Ca(HSiO3)+ | Calcium silicate complex |
| Mg(HSiO3)+ | Mg(HSiO3)+ | Magnesium silicate complex |
| AlH3SiO4+2 | AlH3SiO4+2 | Aluminum silicate complex |
| SiF6-2 | SiF6-2 | Hexafluorosilicate |

### 3. Chemical Equivalence

Note that `H4SiO4(aq)` and `SiO2(aq)` represent the same chemical species:
- **H4SiO4** = Si(OH)4 = silicic acid
- **SiO2·2H2O** = Si(OH)4 = same structure
- The difference is in how the formula is written, not the actual molecular structure

In solution at low temperatures and neutral pH, dissolved silica primarily exists as **monomeric orthosilicic acid Si(OH)4**, which can be written as either:
- H4SiO4 (SUPCRT convention)
- SiO2 + 2H2O (DEW convention)

### 4. DEW Database Unique Features

The DEW database includes polymeric silica species that are important at high temperatures:
- **Si2O4_aq**: Dimeric silica
- **Si3O6_aq**: Trimeric silica

These species become increasingly important at higher temperatures where silica polymerization occurs.

## Solution for Your Code

### Current Code (INCORRECT for DEW):
```python
aqueous = AqueousPhase("H+ OH- H4SiO4(aq)")  # ❌ H4SiO4(aq) not in DEW!
```

### Corrected Code for DEW Database:
```python
# For simple quartz solubility (low T, neutral pH)
aqueous = AqueousPhase("H+ OH- SiO2_aq")

# For high-temperature modeling (include polymerization)
aqueous = AqueousPhase("H+ OH- SiO2_aq HSiO3- Si2O4_aq Si3O6_aq")
```

### If Using SUPCRT Database:
```python
aqueous = AqueousPhase("H+ OH- H4SiO4(aq)")  # ✓ This works with SUPCRT
```

## Database Files Location

**DEW Databases:**
- `embedded/databases/DEW/dew2019-aqueous.yaml` (229 species)
- `embedded/databases/DEW/dew2024-aqueous.yaml` (updated version)
- `embedded/databases/DEW/dew2019-gas.yaml`
- `embedded/databases/DEW/dew2024-gas.yaml`

**SUPCRT Databases:**
- `embedded/databases/reaktoro/supcrtbl.yaml` (1108 species)
- `embedded/databases/reaktoro/supcrtbl-organics.yaml`

## Database Selection Guidance

### Use DEW Database when:
- Working at high temperatures (>100°C)
- High pressures (deep earth, geothermal)
- Need DEW-specific activity model
- Studying polymerization of dissolved species
- Temperature range: 25-1000°C
- Pressure range: 1-60 kbar

### Use SUPCRT Database when:
- Working at standard conditions or moderate T/P
- Need broad species coverage (1108 vs 229 species)
- Using standard HKF activity model
- Following established geochemical modeling conventions

## References

1. **DEW Database Sources:**
   - Sverjensky et al. (2014) - Fit to Raman speciation and quartz solubility data
   - Huang & Sverjensky (2019) - Updated parameters for various species
   - Based on Zhang-Duan water EOS and DEW electrostatic models

2. **Code Implementation:**
   - `Reaktoro/Extensions/DEW/DEWDatabase.hpp` - DEW database class
   - `Reaktoro/Extensions/DEW/DEWDatabase.cpp` - Implementation with YAML parsing
   - `embedded/databases/DEW/build_dew_reaktoro_db.py` - Database builder script

3. **Species Name Conversion:**
   - Line 115 in `DEWDatabase.cpp`: Comma replacement logic
   ```cpp
   String modified_name = name;
   auto comma_in_name = modified_name.find(",aq");
   if(comma_in_name != String::npos)
       modified_name = modified_name.substr(0, comma_in_name) + "_aq";
   ```

## Testing Commands

To verify available species in your installation:

```python
from reaktoro import *

# Check DEW database
dew_db = DEWDatabase('dew2019-aqueous')
print("DEW species:", [s.name() for s in dew_db.species() if 'Si' in str(s.formula())])

# Check SUPCRT database
supcrt_db = SupcrtDatabase('supcrtbl')
aqueous_species = [s for s in supcrt_db.species() if s.aggregateState() == AggregateState.Aqueous]
print("SUPCRT silica:", [s.name() for s in aqueous_species if 'Si' in str(s.formula())])
```

## Recommended Action

**For your quartz solubility analysis script:**

1. **Option A: Use DEW database correctly**
   ```python
   db = DEWDatabase("dew2019-aqueous")
   aqueous = AqueousPhase("H+ OH- SiO2_aq HSiO3-")  # Fixed species names
   ```

2. **Option B: Switch to SUPCRT**
   ```python
   db = SupcrtDatabase("supcrtbl")
   aqueous = AqueousPhase("H+ OH- H4SiO4(aq) HSiO3-")  # Keep original names
   ```

3. **Option C: Hybrid approach** (if you need both databases)
   ```python
   # Use SUPCRT for basic aqueous chemistry
   supcrt_db = SupcrtDatabase("supcrtbl")

   # Extend with DEW for high T/P water properties
   # (This requires careful integration - see DEW_INTEGRATION_STRATEGY.md)
   ```
