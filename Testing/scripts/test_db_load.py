import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "build" / "Reaktoro" / "Release"))
from reaktoro4py import *

print("Imports OK", flush=True)

# Try loading supcrtbl from file
db_path = str(REPO_ROOT / "embedded" / "databases" / "reaktoro" / "supcrtbl.yaml")
try:
    db = Database.fromFile(db_path)
    print("Loaded supcrtbl.yaml from file OK", flush=True)
    names = [s.name() for s in db.species()]
    brucite_matches = [n for n in names if "rucit" in n.lower()]
    print("Brucite matches:", brucite_matches, flush=True)
except Exception as e:
    print(f"ERROR loading yaml: {type(e).__name__}: {e}", flush=True)

# Also try supcrt07 (which we know doesn't crash SupcrtDatabase)
try:
    db2 = SupcrtDatabase("supcrt07")
    print("SupcrtDatabase supcrt07 OK", flush=True)
    names2 = [s.name() for s in db2.species()]
    brucite2 = [n for n in names2 if "rucit" in n.lower()]
    print("supcrt07 Brucite matches:", brucite2, flush=True)
except Exception as e:
    print(f"supcrt07 ERROR: {type(e).__name__}: {e}", flush=True)
