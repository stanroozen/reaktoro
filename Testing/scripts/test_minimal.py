import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "build" / "Reaktoro" / "Release"))
print("step 1: sys.path set", flush=True)
import autodiff

print("step 2: autodiff imported", flush=True)
from reaktoro4py import *

print("step 3: reaktoro4py imported", flush=True)

db = DEWDatabase("dew2024-aqueous")
print("step 4: DEW db loaded", flush=True)

hp_db = Database.fromFile(
    str(
        REPO_ROOT / "embedded" / "databases" / "hollandpowell" / "tc-ds62-reaktoro.json"
    )
)
print("step 5: HP db loaded", flush=True)

combined_db = Database(db.species())
combined_db.addSpecies(hp_db.species("br"))
print("step 6: combined_db built", flush=True)

aqueous = AqueousPhase("H2O(aq) H+(aq) OH-(aq) Mg+2(aq) MgOH+(aq)")
aqueous.setActivityModel(ActivityModelDEW())
mineral = MineralPhase("br")
system = ChemicalSystem(combined_db, aqueous, mineral)
print("step 7: system built", flush=True)

state = ChemicalState(system)
print("step 8: state created", flush=True)

# Test plain float
try:
    state.set("H2O(aq)", 55.5, "mol")
    print("step 9: state.set float OK", flush=True)
except TypeError as e:
    print(f"step 9: state.set float TypeError (expected): {e}", flush=True)
    state.set("H2O(aq)", autodiff.real(55.5), "mol")
    print("step 9b: state.set autodiff.real OK", flush=True)
except Exception as e:
    print(f"step 9: UNEXPECTED error {type(e).__name__}: {e}", flush=True)
