import os
import sys
import re
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

names = []
for sp in dew_db.species():
    formula = str(sp.formula())
    elems = set(pattern.findall(formula))
    if elems and elems.issubset(allowed):
        names.append(sp.name())

names = sorted(set(names))
print("\n".join(names))
