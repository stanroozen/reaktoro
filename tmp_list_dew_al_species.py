import os
import re
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PYD_DIR = os.path.join(SCRIPT_DIR, "build-msvc", "Reaktoro", "Release")
if os.path.isdir(PYD_DIR) and PYD_DIR not in sys.path:
    sys.path.insert(0, PYD_DIR)

try:
    from reaktoro import *  # noqa: F401,F403
except ModuleNotFoundError:
    from reaktoro4py import *  # noqa: F401,F403


dew_db = DEWDatabase("dew2024-aqueous")
pattern = re.compile(r"[A-Z][a-z]?")

al_species = []
for sp in dew_db.species():
    formula = str(sp.formula())
    elems = set(pattern.findall(formula))
    if "Al" in elems:
        al_species.append(sp.name())

al_species = sorted(set(al_species))
print("Count:", len(al_species))
print("\n".join(al_species))
