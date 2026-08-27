import os
import sys
from pathlib import Path

p = (
    Path(__file__).resolve().parents[1]
    / "build"
    / "python"
    / "package"
    / "build"
    / "lib"
    / "reaktoro"
)
sys.path.insert(0, str(p))
os.add_dll_directory(str(p))

import reaktoro4py as r

root = Path(__file__).resolve().parents[1]
f = root / "embedded" / "databases" / "DEW" / "dew2024-aqueous.yaml"
print("file", f)
db = r.Database.fromFile(str(f))
print("ok", len(db.species()))
print("has H2O", any(s.name() == "H2O(aq)" for s in db.species()))
