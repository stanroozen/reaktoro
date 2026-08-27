import importlib.util
import os
import sys

import numpy as np


REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REL = os.path.join(REPO, "build", "Reaktoro", "Release")
if REL not in sys.path:
    sys.path.insert(0, REL)

TUTORIAL = os.path.join(
    REPO,
    "DEW_Experimental_Benchmark",
    "Tutorial",
    "willemite_solubility_tutorial_dew17hp622_zn.py",
)
spec = importlib.util.spec_from_file_location("w", TUTORIAL)
w = importlib.util.module_from_spec(spec)
spec.loader.exec_module(w)


# Zn-focused assemblage with key gangue candidates.
w.USE_COMPETING_ZN_MINERALS = True
w.INCLUDE_COMPATIBLE_GANGUE_MINERALS = False
w.COMPETING_ZN_MINERALS = list(
    dict.fromkeys(w.COMPETING_ZN_MINERALS + ["hem", "cc", "q"])
)

# Add dissolved CO2 for fluid-rock C transfer under mixed-fluid conditions.
if "CO2,aq" not in w.AQUEOUS_SPECIES:
    w.AQUEOUS_SPECIES = list(w.AQUEOUS_SPECIES) + ["CO2,aq"]

# Seed gangue minerals so they can stabilize if favored.
for m in ["hem", "cc", "q"]:
    w.INITIAL_SPECIES_AMOUNTS_MOL[m] = 2.0

w.validate_user_inputs()

# Build a mixed database: DEW Zn system + H2O(g)/CO2(g) gas species from SUPCRT.
dew_db = w.Database.fromFile(w.PERPLEX_DATABASE_FILE)
supcrt_db = w.Database.fromFile(
    os.path.join(REPO, "embedded", "databases", "reaktoro", "supcrt07.yaml")
)
dew_db.addSpecies(supcrt_db.species("H2O(g)"))
dew_db.addSpecies(supcrt_db.species("CO2(g)"))

# PerplexDEW aqueous model with default conflict handling (warn-only).
perplexdew_params = w.ActivityModelParamsPerplexDEW()
perplexdew_params.errorOnConflictingStandardState = False
aq = w.AqueousPhase(" ".join(w.AQUEOUS_SPECIES))
aq.setActivityModel(w.ActivityModelPerplexDEW(perplexdew_params))

mins = w.make_mineral_phases()

gasm = w.GaseousPhase("H2O(g) CO2(g)")
gfsm_params = w.ActivityModelParamsPerplexGFSM()
gasm.setActivityModel(w.ActivityModelPerplexGFSM(gfsm_params))

system = w.ChemicalSystem(dew_db, aq, mins, gasm)

specs = w.EquilibriumSpecs(system)
specs.temperature()
specs.pressure()
specs.pH()

solver = w.make_equilibrium_solver(system, specs)
opts = w.EquilibriumOptions()
if hasattr(opts, "hessian") and hasattr(w, "GibbsHessian"):
    opts.hessian = w.GibbsHessian.Exact
if hasattr(opts, "warmstart"):
    opts.warmstart = True
solver.setOptions(opts)
conditions = w.EquilibriumConditions(specs)

Tvals = [250.0, 300.0, 350.0]
pHvals = [5.0, 7.0, 9.0]
XCO2_vals = [0.0, 0.1, 0.3, 0.5, 0.7, 0.9]
Pbar = 2000.0
TH = 1.0e-8
N_GAS = 1.0


def present(state, name):
    try:
        return float(state.speciesAmount(name)) > TH
    except Exception:
        return False


hits = {x: [] for x in XCO2_vals}
solved = {x: 0 for x in XCO2_vals}

for xco2 in XCO2_vals:
    for T in Tvals:
        for pH in pHvals:
            state = w.make_base_state(system)

            # Initialize mixed-fluid composition.
            nco2 = xco2 * N_GAS
            nh2o = (1.0 - xco2) * N_GAS
            if nh2o > 0.0:
                state.set("H2O(g)", nh2o, "mol")
            if nco2 > 0.0:
                state.set("CO2(g)", nco2, "mol")

            # Seed dissolved CO2 to help speciation pathway.
            state.set("CO2,aq", max(1.0e-12, 1.0e-6 * xco2), "mol")

            conditions.temperature(float(T), "celsius")
            conditions.pressure(float(Pbar), "bar")
            conditions.pH(float(pH))

            res = solver.solve(state, conditions)
            if not res.succeeded():
                continue

            solved[xco2] += 1
            if present(state, "Wlm"):
                hits[xco2].append((T, pH, float(state.speciesAmount("Wlm"))))


total_per_x = len(Tvals) * len(pHvals)
print("MODEL=PerplexDEW(aq)+PerplexGFSM(H2O-CO2 gas)")
print(f"GRID= T:{Tvals} pH:{pHvals} Pbar:{Pbar}")
for xco2 in XCO2_vals:
    print(f"\nXCO2={xco2:.1f}")
    print(f"  solved={solved[xco2]}/{total_per_x}")
    print(f"  wlm_present={len(hits[xco2])}")
    for T, pH, wlm in hits[xco2][:10]:
        print(f"    T={T:.1f} C, pH={pH:.1f}, Wlm={wlm:.3e} mol")
