"""Trace build_system step by step"""

import sys
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(SCRIPT_DIR)
PYD_DIR = os.path.join(ROOT_DIR, "build-msvc", "Reaktoro", "Release")
if os.path.isdir(PYD_DIR) and PYD_DIR not in sys.path:
    sys.path.insert(0, PYD_DIR)

print("[1] Importing reaktoro4py...")
from reaktoro4py import *

print("[2] Loading databases...")
dew_db = DEWDatabase("dew2019-aqueous")
supcrt_db = SupcrtDatabase("supcrtbl")
print(f"    DEW species: {len(dew_db.species())}")

print("[3] Getting Quartz...")
quartz = supcrt_db.species("Quartz")
print(f"    Quartz: {quartz.name()}")

print("[4] Creating combined database...")
combined_db = Database(dew_db.species())
combined_db.addSpecies(quartz)
print(f"    Combined DB created")

print("[5] Creating AqueousPhase...")
aqueous = AqueousPhase(
    "WATER,AQ H+ OH- SiO2_aq H2_aq O2_aq HO2- HSiO3- Si2O4_aq Si3O6_aq"
)
print(f"    AqueousPhase created")

print("[6] Setting ActivityModelDEW...")
try:
    aqueous.setActivityModel(ActivityModelDEW())
    print(f"    ActivityModelDEW set")
except Exception as e:
    print(f"    [WARN] ActivityModelDEW failed: {e}")

print("[7] Creating MineralPhase...")
mineral = MineralPhase("Quartz")
print(f"    MineralPhase created")

print("[8] Creating ChemicalSystem...")
system = ChemicalSystem(combined_db, aqueous, mineral)
print(f"    ChemicalSystem created with {len(system.species())} species")

print("[OK] All steps completed successfully!")
