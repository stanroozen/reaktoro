"""Check SiO2_aq definition and reaction stoichiometry in DEW database."""

import os, sys

_pfx = sys.prefix
_env_paths = [
    _pfx,
    os.path.join(_pfx, "Library", "mingw-w64", "bin"),
    os.path.join(_pfx, "Library", "bin"),
    os.path.join(_pfx, "Scripts"),
]
_sys_root = os.environ.get("SystemRoot", r"C:\Windows")
os.environ["PATH"] = ";".join(
    p for p in _env_paths + [os.path.join(_sys_root, "System32")] if os.path.isdir(p)
)

import numpy as np
import autodiff  # noqa

_pyd_dir = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "..",
    "..",
    "build",
    "Reaktoro",
    "Release",
)
sys.path.insert(0, os.path.normpath(_pyd_dir))
import reaktoro4py as r4

r4.Warnings.disable(906)

dew_db = r4.DEWDatabase("dew2024-aqueous")
supcrt_db = r4.SupcrtDatabase("supcrtbl")

for db_name, db, sp_name in [
    ("DEW", dew_db, "SiO2_aq"),
    ("DEW", dew_db, "HSiO3-"),
    ("DEW", dew_db, "Si2O4_aq"),
    ("DEW", dew_db, "Si3O6_aq"),
    ("SUPCRTBL", supcrt_db, "SiO2(aq)"),
    ("SUPCRTBL", supcrt_db, "HSiO3-"),
]:
    try:
        sp = db.species(sp_name)
        print(f"\n{db_name}: {sp_name}")
        print(f"  formula: {sp.formula()}")
        elem = sp.elements()
        for i in range(elem.size()):
            print(f"  element: {elem[i].symbol()} = {elem.coefficient(i)}")
        print(f"  charge: {sp.charge()}")
    except Exception as e:
        print(f"\n{db_name}: {sp_name} -> {e}")

# Check if SiO2_aq in DEW includes water (O count)
print("\n=== Element analysis ===")
dew_sp = dew_db.species("SiO2_aq")
formula_obj = dew_sp.formula()
print(f"SiO2_aq formula: {formula_obj}")
print("Elements:")
elems = dew_sp.elements()
for i in range(elems.size()):
    e = elems[i]
    print(f"  {e.symbol()}: {elems.coefficient(i)}")

