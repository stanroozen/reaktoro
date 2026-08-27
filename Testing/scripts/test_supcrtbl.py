import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "build" / "Reaktoro" / "Release"))
import autodiff
from reaktoro4py import *

print("Imports OK", flush=True)
db = SupcrtDatabase("supcrtbl")
print("SupcrtDatabase supcrtbl loaded OK", flush=True)
names = [s.name() for s in db.species()]
brucite_matches = [n for n in names if "rucite" in n or "rucit" in n]
print("Brucite matches:", brucite_matches, flush=True)
print("'Brucite' in db:", "Brucite" in names, flush=True)
