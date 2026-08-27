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

try:
    w.Warnings.disable(906)
except Exception:
    pass


# Configure a tighter geologic test system.
w.USE_COMPETING_ZN_MINERALS = True
w.INCLUDE_COMPATIBLE_GANGUE_MINERALS = False
w.COMPETING_ZN_MINERALS = list(
    dict.fromkeys(w.COMPETING_ZN_MINERALS + ["hem", "cc", "q", "mag"])
)

# Add major cation aqueous species to permit Fe/Ca/Mg bulk sweeps.
extra_aq = [
    "CO2,aq",
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
w.AQUEOUS_SPECIES = list(dict.fromkeys(list(w.AQUEOUS_SPECIES) + extra_aq))

# Seed gangue solids.
for m in ["hem", "cc", "q", "mag"]:
    w.INITIAL_SPECIES_AMOUNTS_MOL[m] = 2.0

w.validate_user_inputs()

# Build mixed database: DEW17 Zn + gas species from SUPCRT.
dew_db = w.Database.fromFile(w.PERPLEX_DATABASE_FILE)
supcrt_db = w.Database.fromFile(
    os.path.join(REPO, "embedded", "databases", "reaktoro", "supcrt07.yaml")
)
dew_db.addSpecies(supcrt_db.species("H2O(g)"))
dew_db.addSpecies(supcrt_db.species("CO2(g)"))

# Aqueous and gas models.
perplexdew_params = w.ActivityModelParamsPerplexDEW()
perplexdew_params.errorOnConflictingStandardState = False
aq = w.AqueousPhase(" ".join(w.AQUEOUS_SPECIES))
aq.setActivityModel(w.ActivityModelPerplexDEW(perplexdew_params))

mins = w.make_mineral_phases()

gas = w.GaseousPhase("H2O(g) CO2(g)")
gfsm_params = w.ActivityModelParamsPerplexGFSM()
gas.setActivityModel(w.ActivityModelPerplexGFSM(gfsm_params))

system = w.ChemicalSystem(dew_db, aq, mins, gas)

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

# Tighter geologic window requested: lower pH + different pressure + Fe/Ca/Mg bulk sweep.
Tvals = [240.0, 320.0]
pHvals = [3.0, 4.5, 6.0]
Pbars = [1000.0, 3000.0]
XCO2_vals = [0.1, 0.5]

# Fe/Ca/Mg bulk seed levels (mol added to initial state).
bulk_levels = [0.0, 1.0e-3]

TH = 1.0e-8
N_GAS = 1.0


def present(state, name):
    try:
        return float(state.speciesAmount(name)) > TH
    except Exception:
        return False


def add_bulk_cations(state, fe_seed, ca_seed, mg_seed):
    if fe_seed > 0.0:
        state.set("Fe+2", fe_seed, "mol")
    if ca_seed > 0.0:
        state.set("Ca+2", ca_seed, "mol")
    if mg_seed > 0.0:
        state.set("Mg+2", mg_seed, "mol")


cases_tested = 0
cases_solved = 0
wlm_hits = []
strict_hits = []

for fe_seed in bulk_levels:
    for ca_seed in bulk_levels:
        for mg_seed in bulk_levels:
            for Pbar in Pbars:
                for xco2 in XCO2_vals:
                    for T in Tvals:
                        for pH in pHvals:
                            cases_tested += 1
                            state = w.make_base_state(system)

                            # Mixed fluid setup.
                            nco2 = xco2 * N_GAS
                            nh2o = (1.0 - xco2) * N_GAS
                            if nh2o > 0.0:
                                state.set("H2O(g)", nh2o, "mol")
                            if nco2 > 0.0:
                                state.set("CO2(g)", nco2, "mol")
                            state.set("CO2,aq", max(1.0e-12, 1.0e-6 * xco2), "mol")

                            add_bulk_cations(state, fe_seed, ca_seed, mg_seed)

                            conditions.temperature(float(T), "celsius")
                            conditions.pressure(float(Pbar), "bar")
                            conditions.pH(float(pH))

                            res = solver.solve(state, conditions)
                            if not res.succeeded():
                                continue

                            cases_solved += 1

                            wlm_amt = float(state.speciesAmount("Wlm"))
                            hem_amt = (
                                float(state.speciesAmount("hem"))
                                if present(state, "hem")
                                else 0.0
                            )
                            cc_amt = (
                                float(state.speciesAmount("cc"))
                                if present(state, "cc")
                                else 0.0
                            )
                            q_amt = (
                                float(state.speciesAmount("q"))
                                if present(state, "q")
                                else 0.0
                            )

                            if wlm_amt > TH:
                                wlm_hits.append(
                                    {
                                        "T": T,
                                        "pH": pH,
                                        "Pbar": Pbar,
                                        "XCO2": xco2,
                                        "Fe": fe_seed,
                                        "Ca": ca_seed,
                                        "Mg": mg_seed,
                                        "Wlm": wlm_amt,
                                        "hem": hem_amt,
                                        "cc": cc_amt,
                                        "q": q_amt,
                                    }
                                )

                            # Strict assemblage hit requested earlier.
                            if wlm_amt > TH and hem_amt > TH and cc_amt > TH:
                                strict_hits.append(
                                    {
                                        "T": T,
                                        "pH": pH,
                                        "Pbar": Pbar,
                                        "XCO2": xco2,
                                        "Fe": fe_seed,
                                        "Ca": ca_seed,
                                        "Mg": mg_seed,
                                        "Wlm": wlm_amt,
                                        "hem": hem_amt,
                                        "cc": cc_amt,
                                        "q": q_amt,
                                    }
                                )

print("MODEL=PerplexDEW(aq)+PerplexGFSM(H2O-CO2 gas), tighter geologic window")
print(f"CASES_TESTED={cases_tested}")
print(f"CASES_SOLVED={cases_solved}")
print(f"WLM_PRESENCE_HITS={len(wlm_hits)}")
print(f"STRICT_WLM+HEM+CC_HITS={len(strict_hits)}")

if wlm_hits:
    print("\nTOP_WLM_HITS (first 25):")
    for h in wlm_hits[:25]:
        print(
            "T={T:.1f}C pH={pH:.1f} P={Pbar:.0f}bar XCO2={XCO2:.1f} "
            "Fe={Fe:.1e} Ca={Ca:.1e} Mg={Mg:.1e} Wlm={Wlm:.3e} hem={hem:.3e} cc={cc:.3e} q={q:.3e}".format(
                **h
            )
        )

if strict_hits:
    print("\nSTRICT_HITS (first 25):")
    for h in strict_hits[:25]:
        print(
            "T={T:.1f}C pH={pH:.1f} P={Pbar:.0f}bar XCO2={XCO2:.1f} "
            "Fe={Fe:.1e} Ca={Ca:.1e} Mg={Mg:.1e} Wlm={Wlm:.3e} hem={hem:.3e} cc={cc:.3e} q={q:.3e}".format(
                **h
            )
        )
