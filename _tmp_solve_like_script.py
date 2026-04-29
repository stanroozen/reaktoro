import os
import sys
import importlib

if os.name == "nt":
    env_prefix = sys.prefix
    env_paths = [
        env_prefix,
        os.path.join(env_prefix, "Library", "mingw-w64", "bin"),
        os.path.join(env_prefix, "Library", "usr", "bin"),
        os.path.join(env_prefix, "Library", "bin"),
        os.path.join(env_prefix, "Scripts"),
        os.path.join(env_prefix, "bin"),
    ]
    system_root = os.environ.get("SystemRoot", r"C:\Windows")
    system_paths = [
        os.path.join(system_root, "System32"),
        system_root,
        os.path.join(system_root, "System32", "Wbem"),
    ]
    os.environ["PATH"] = ";".join([p for p in env_paths + system_paths if os.path.isdir(p)])

try:
    import autodiff
except ModuleNotFoundError:
    class _AutoDiffShim:
        @staticmethod
        def real(value):
            return value

    autodiff = _AutoDiffShim()

try:
    from reaktoro import *
except ModuleNotFoundError:
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    ROOT_DIR = SCRIPT_DIR
    pyd_candidates = [
        os.path.join(ROOT_DIR, "build-dew", "Reaktoro", "Release"),
        os.path.join(ROOT_DIR, "build-msvc", "Reaktoro", "Release"),
    ]
    for pyd_dir in pyd_candidates:
        if not os.path.isdir(pyd_dir):
            continue
        sys.path.insert(0, pyd_dir)
        local_mod = importlib.import_module("reaktoro4py")
        globals().update({name: getattr(local_mod, name) for name in dir(local_mod) if not name.startswith("_")})
        break


def to_real(value):
    try:
        return autodiff.real(value)
    except Exception:
        return value

from pathlib import Path
root = Path.cwd()

m = Database.fromFile(str(root / "embedded" / "databases" / "hollandpowell" / "tc-ds62-reaktoro.json"))
dew = DEWDatabase("dew2024-aqueous")
comb = Database(dew.species())
comb.addSpecies(m.species("q"))

aq = AqueousPhase("WATER,AQ H+ OH- SiO2_aq H2_aq O2_aq HO2- HSiO3- Si2O4_aq Si3O6_aq")
aq.setActivityModel(ActivityModelPerplexDEW(ActivityDHModel.Davies))
system = ChemicalSystem(comb, aq, MineralPhase("q"))

specs = EquilibriumSpecs(system)
specs.temperature()
specs.pressure()
solver = EquilibriumSolver(specs)
conds = EquilibriumConditions(specs)
state = ChemicalState(system)
state.set("WATER,AQ", to_real(1.0), "kg")
state.set("H+", to_real(1e-8), "mol")
state.set("OH-", to_real(1e-8), "mol")
state.set("SiO2_aq", to_real(1e-6), "mol")
state.set("q", to_real(10.0), "mol")
conds.temperature(300.0, "celsius")
conds.pressure(1000.0, "bar")
res = solver.solve(state, conds)
print("ok", res.succeeded())
if res.succeeded():
    print(float(AqueousProps(state).elementMolality("Si")))
