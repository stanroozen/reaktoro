"""
Template: Copy this file and modify MINERAL_CONFIG for your mineral analysis

Example: Analyzing Calcite Solubility with DEW2024
"""

# =============================================================================
# CHANGE THIS SECTION FOR YOUR MINERAL
# =============================================================================

MINERAL_CONFIG = {
    # Mineral identification
    "mineral_name": "Calcite",  # Name as it appears in SUPCRTBL database
    "mineral_formula": "CaCO3",  # Chemical formula (for display)
    "solute_species": "Ca+2",  # Primary aqueous species to track
    # Additional aqueous species (besides WATER,AQ, H+, OH-)
    # Include complexes and related species
    "aqueous_species": "CO2_aq HCO3- CO3-2 CaCO3_aq CaHCO3+ CaOH+",
    # File paths
    "csv_file": "calcite_DEW_testset.csv",  # Your experimental data CSV
    "output_prefix": "calcite",  # Prefix for output figure files
    # Plot settings
    "plot_title": "Calcite Solubility",
    "y_label": "Calcite Solubility (mol/kg-H₂O)",
}

# =============================================================================
# NO CHANGES NEEDED BELOW - Just copy the rest from the original script
# =============================================================================

# The rest of the code is completely generic and works for any mineral!
# Just update MINERAL_CONFIG above and run.

"""
Quick Start:
1. Prepare your experimental data CSV with columns:
   - T_C (temperature in Celsius)
   - P_kbar (pressure in kbar)
   - molality_m (measured molality in mol/kg H2O)
   - reference (citation)
   - experiment_type (e.g., "hydrothermal")

2. Check species availability:
   from reaktoro import *
   db_dew = DEWDatabase("dew2024-aqueous")
   db_supcrt = SupcrtDatabase("supcrtbl")

   # Verify mineral exists
   print(db_supcrt.species("Calcite"))

   # Check available aqueous species
   for s in db_dew.species():
       if "Ca" in s.name() or "CO" in s.name():
           print(s.name())

3. Update MINERAL_CONFIG above

4. Copy the remaining code from quartz_solubility_analysis_v2_dew24.py
   (everything after line 68)

5. Run: python your_mineral_analysis.py
"""

# =============================================================================
# Example Configurations for Common Minerals
# =============================================================================

# Halite (NaCl)
HALITE_CONFIG = {
    "mineral_name": "Halite",
    "mineral_formula": "NaCl",
    "solute_species": "Na+",
    "aqueous_species": "Cl- NaCl_aq NaOH_aq",
    "csv_file": "halite_DEW_testset.csv",
    "output_prefix": "halite",
    "plot_title": "Halite Solubility",
    "y_label": "Halite Solubility (mol/kg-H₂O)",
}

# Corundum (Al2O3)
CORUNDUM_CONFIG = {
    "mineral_name": "Corundum",
    "mineral_formula": "Al2O3",
    "solute_species": "Al+3",
    "aqueous_species": "AlOH+2 Al(OH)2+ Al(OH)3_aq Al(OH)4- AlO2-",
    "csv_file": "corundum_DEW_testset.csv",
    "output_prefix": "corundum",
    "plot_title": "Corundum Solubility",
    "y_label": "Corundum Solubility (mol/kg-H₂O)",
}

# Periclase (MgO)
PERICLASE_CONFIG = {
    "mineral_name": "Periclase",
    "mineral_formula": "MgO",
    "solute_species": "Mg+2",
    "aqueous_species": "MgOH+ Mg(OH)2_aq",
    "csv_file": "periclase_DEW_testset.csv",
    "output_prefix": "periclase",
    "plot_title": "Periclase Solubility",
    "y_label": "Periclase Solubility (mol/kg-H₂O)",
}

# Rutile (TiO2)
RUTILE_CONFIG = {
    "mineral_name": "Rutile",
    "mineral_formula": "TiO2",
    "solute_species": "Ti(OH)4_aq",
    "aqueous_species": "TiO2_aq Ti(OH)3+ Ti(OH)3-",
    "csv_file": "rutile_DEW_testset.csv",
    "output_prefix": "rutile",
    "plot_title": "Rutile Solubility",
    "y_label": "Rutile Solubility (mol/kg-H₂O)",
}
