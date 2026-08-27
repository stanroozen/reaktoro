import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "build" / "Reaktoro" / "Release"))
from reaktoro4py import *

print("Imports OK", flush=True)
try:
    db = SupcrtDatabase("supcrtbl")
    print("SupcrtDatabase supcrtbl loaded OK", flush=True)
    names = [s.name() for s in db.species()]
    brucite_matches = [n for n in names if "rucit" in n.lower()]
    print("Brucite matches:", brucite_matches, flush=True)
except Exception as e:
    print(f"ERROR: {type(e).__name__}: {e}", flush=True)
