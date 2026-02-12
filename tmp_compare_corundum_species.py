import os
import re
import sys
import importlib.util
import importlib

ROOT_DIR = r"C:\Users\stanroozen\Documents\Projects\reaktoro-dev\reaktoro"
PYD_DIR = os.path.join(ROOT_DIR, "build-msvc", "Reaktoro", "Release")
if os.path.isdir(PYD_DIR):
    if PYD_DIR not in sys.path:
        sys.path.insert(0, PYD_DIR)
    try:
        os.add_dll_directory(PYD_DIR)
    except Exception:
        pass

spec = importlib.util.find_spec("reaktoro")
if spec is None:
    rkt = importlib.import_module("reaktoro4py")
else:
    rkt = importlib.import_module("reaktoro")

dew_db = rkt.DEWDatabase("dew2024-aqueous")
allowed = {"Al", "O", "H"}
pattern = re.compile(r"[A-Z][a-z]?")

aloh_species = set()
for sp in dew_db.species():
    formula = str(sp.formula())
    elems = set(pattern.findall(formula))
    if elems and elems.issubset(allowed):
        aloh_species.add(sp.name())

included = {
    "WATER,AQ",
    "H+",
    "OH-",
    "Al+3",
    "Al(OH)+2",
    "AlO2-",
    "HAlO2_aq",
    "H2_aq",
    "O2_aq",
}

missing = sorted(aloh_species - included)
extra = sorted(included - aloh_species)
print("MISSING:")
print("\n".join(missing) if missing else "<none>")
print("\nEXTRA:")
print("\n".join(extra) if extra else "<none>")
