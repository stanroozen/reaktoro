import json
import os
import sys
import importlib.util

# Local Reaktoro loading
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
rel = os.path.join(REPO, "build", "Reaktoro", "Release")
if rel not in sys.path:
    sys.path.insert(0, rel)

import reaktoro4py as rkt

DB_PATH = os.path.join(
    REPO, "embedded", "databases", "perplex", "DEW17HP622_Zn_2025-reaktoro.json"
)

db_json = json.load(open(DB_PATH, encoding="utf-8"))
species_data = db_json["Species"]

# Start from tutorial aqueous set, then add Ca/Mg/Fe species needed for
# hematite/calcite(-dolomite proxy) buffering.
TUTORIAL_PATH = os.path.join(
    REPO,
    "DEW_Experimental_Benchmark",
    "Tutorial",
    "willemite_solubility_tutorial_dew17hp622_zn.py",
)
spec = importlib.util.spec_from_file_location("wmod", TUTORIAL_PATH)
wmod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(wmod)

extra_aqueous = [
    "Ca+2",
    "Ca(OH)+",
    "CaCO3,aq",
    "Ca(HCO3)",
    "CaCl+",
    "CaCl2,aq",
    "CaSO4,aq",
    "Mg+2",
    "MgOH+",
    "MgCO3,aq",
    "Mg(HCO3)",
    "MgCl+",
    "MgSO4,aq",
    "Fe+2",
    "Fe+3",
    "Fe(OH)+",
    "FeO,aq",
    "HFeO2-",
    "FeCl+",
    "FeCl2,aq",
    "FeCl2+",
]

aqueous_species = list(dict.fromkeys(list(wmod.AQUEOUS_SPECIES) + extra_aqueous))

# Observed minerals mapped to DEW names.
HEMATITE = "hem"
CALCITE = "cc"
QUARTZ = "q"
ZINCITE = "Znc"
WILLEMITE = "Wlm"
# Dolomite is not present in this DEW variant; use magnesite as Mg-carbonate proxy.
MAGNESITE = "mag"

required = [HEMATITE, CALCITE, QUARTZ, ZINCITE, WILLEMITE]
proxy_set = [HEMATITE, CALCITE, MAGNESITE, QUARTZ, ZINCITE, WILLEMITE]

mineral_names = sorted(set(required + [MAGNESITE, "arag", "trd", "crst", "coe", "stv"]))

database = rkt.Database.fromFile(DB_PATH)
aq = rkt.AqueousPhase(" ".join(aqueous_species))
aq.setActivityModel(rkt.ActivityModelPerplexDEW())
mins = rkt.MineralPhases(rkt.StringList(mineral_names))
gas = rkt.GaseousPhase("O2 CO2")
gas.setActivityModel(rkt.ActivityModelIdealGas())
system = rkt.ChemicalSystem(database, aq, mins, gas)

specs = rkt.EquilibriumSpecs(system)
specs.temperature()
specs.pressure()
specs.pH()
specs.fugacity("O2")
specs.fugacity("CO2")

solver = rkt.EquilibriumSolver(specs)
opts = rkt.EquilibriumOptions()
if hasattr(opts, "hessian") and hasattr(rkt, "GibbsHessian"):
    opts.hessian = rkt.GibbsHessian.Exact
solver.setOptions(opts)
conditions = rkt.EquilibriumConditions(specs)

# Coarse search domain; adjust as needed.
T_values = [150.0, 200.0, 250.0, 300.0, 350.0, 400.0]
pH_values = [4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
logfO2_values = [-30.0, -25.0, -20.0, -15.0]
logfCO2_values = [-4.0, -3.0, -2.0, -1.0, 0.0]

THRESH = 1.0e-8

name_to_idx = {system.species(i).name(): i for i in range(system.species().size())}


def seed_state():
    st = rkt.ChemicalState(system)
    amounts = [1.0e-20] * system.species().size()
    # fluid seeds
    for nm, val in [
        ("H2O", 55.5),
        ("H+", 1e-7),
        ("OH-", 1e-7),
        ("O2", 1e-20),
        ("CO2", 1e-20),
    ]:
        i = name_to_idx.get(nm)
        if i is not None:
            amounts[i] = val
    # seed minerals to allow appearance/disappearance
    for nm in mineral_names:
        i = name_to_idx.get(nm)
        if i is not None:
            amounts[i] = 1.0
    st.setSpeciesAmounts(amounts)
    return st


hits_required = []
hits_proxy = []
total = 0
ok = 0

for T in T_values:
    for pH in pH_values:
        for lfO2 in logfO2_values:
            for lfCO2 in logfCO2_values:
                total += 1
                st = seed_state()
                conditions.temperature(T, "celsius")
                conditions.pressure(2000.0, "bar")
                conditions.pH(pH)
                conditions.fugacity("O2", 10.0**lfO2, "bar")
                conditions.fugacity("CO2", 10.0**lfCO2, "bar")
                res = solver.solve(st, conditions)
                if not res.succeeded():
                    continue
                ok += 1

                def present(name):
                    try:
                        return float(st.speciesAmount(name)) > THRESH
                    except Exception:
                        return False

                if all(present(n) for n in required):
                    hits_required.append((T, pH, lfO2, lfCO2))
                if all(present(n) for n in proxy_set):
                    hits_proxy.append((T, pH, lfO2, lfCO2))

print(f"TOTAL_POINTS={total} SOLVED={ok}")
print(f"HITS_REQUIRED(hem+cc+q+Znc+Wlm)={len(hits_required)}")
for h in hits_required[:40]:
    print("REQ", h)
print(f"HITS_PROXY(+mag for dolomite proxy)={len(hits_proxy)}")
for h in hits_proxy[:40]:
    print("PRX", h)
