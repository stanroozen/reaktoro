import importlib.util
import os
import sys


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


w.USE_COMPETING_ZN_MINERALS = True
w.INCLUDE_COMPATIBLE_GANGUE_MINERALS = False
w.COMPETING_ZN_MINERALS = list(
    dict.fromkeys(w.COMPETING_ZN_MINERALS + ["hem", "cc", "q", "mag"])
)

extra_aq = [
    "CO2,aq",
    "Ca+2",
    "Ca(OH)+",
    "CaCO3,aq",
    "Ca(HCO3)",
    "Mg+2",
    "MgOH+",
    "MgCO3,aq",
    "Mg(HCO3)",
    "Fe+2",
    "Fe+3",
    "Fe(OH)+",
    "FeO,aq",
    "HFeO2-",
]
w.AQUEOUS_SPECIES = list(dict.fromkeys(list(w.AQUEOUS_SPECIES) + extra_aq))

for m in ["hem", "cc", "q", "mag"]:
    w.INITIAL_SPECIES_AMOUNTS_MOL[m] = 2.0

w.validate_user_inputs()

# Mixed database with H2O(g), CO2(g).
dew_db = w.Database.fromFile(w.PERPLEX_DATABASE_FILE)
supcrt_db = w.Database.fromFile(
    os.path.join(REPO, "embedded", "databases", "reaktoro", "supcrt07.yaml")
)
dew_db.addSpecies(supcrt_db.species("H2O(g)"))
dew_db.addSpecies(supcrt_db.species("CO2(g)"))

perplexdew_params = w.ActivityModelParamsPerplexDEW()
perplexdew_params.errorOnConflictingStandardState = False
aq = w.AqueousPhase(" ".join(w.AQUEOUS_SPECIES))
aq.setActivityModel(w.ActivityModelPerplexDEW(perplexdew_params))

mins = w.make_mineral_phases()
gas = w.GaseousPhase("H2O(g) CO2(g)")
gas.setActivityModel(w.ActivityModelPerplexGFSM(w.ActivityModelParamsPerplexGFSM()))

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

# Compact geologic window.
Tvals = [240.0, 320.0]
pHvals = [3.0, 5.0]
Pbars = [1000.0, 3000.0]
XCO2_vals = [0.1, 0.5]

# Representative Fe/Ca/Mg bulk seeds: low and enriched.
bulk_cases = [
    (0.0, 0.0, 0.0),
    (1.0e-3, 1.0e-3, 1.0e-3),
    (1.0e-3, 0.0, 0.0),
    (0.0, 1.0e-3, 1.0e-3),
]

TH = 1.0e-8
N_GAS = 1.0


def amt(state, name):
    try:
        return float(state.speciesAmount(name))
    except Exception:
        return 0.0


cases = 0
solved = 0
wlm_hits = []
strict_hits = []

for fe_seed, ca_seed, mg_seed in bulk_cases:
    for pbar in Pbars:
        for xco2 in XCO2_vals:
            for t in Tvals:
                for ph in pHvals:
                    cases += 1

                    st = w.make_base_state(system)
                    st.set("H2O(g)", (1.0 - xco2) * N_GAS, "mol")
                    st.set("CO2(g)", xco2 * N_GAS, "mol")
                    st.set("CO2,aq", max(1.0e-12, 1.0e-6 * xco2), "mol")

                    if fe_seed > 0.0:
                        st.set("Fe+2", fe_seed, "mol")
                    if ca_seed > 0.0:
                        st.set("Ca+2", ca_seed, "mol")
                    if mg_seed > 0.0:
                        st.set("Mg+2", mg_seed, "mol")

                    conditions.temperature(t, "celsius")
                    conditions.pressure(pbar, "bar")
                    conditions.pH(ph)

                    res = solver.solve(st, conditions)
                    if not res.succeeded():
                        continue

                    solved += 1

                    wlm = amt(st, "Wlm")
                    hem = amt(st, "hem")
                    cc = amt(st, "cc")
                    q = amt(st, "q")

                    if wlm > TH:
                        hit = (
                            t,
                            ph,
                            pbar,
                            xco2,
                            fe_seed,
                            ca_seed,
                            mg_seed,
                            wlm,
                            hem,
                            cc,
                            q,
                        )
                        wlm_hits.append(hit)
                    if wlm > TH and hem > TH and cc > TH:
                        strict_hits.append(
                            (
                                t,
                                ph,
                                pbar,
                                xco2,
                                fe_seed,
                                ca_seed,
                                mg_seed,
                                wlm,
                                hem,
                                cc,
                                q,
                            )
                        )

print("MODEL=PerplexDEW+PerplexGFSM H2O-CO2 mixed-fluid")
print(f"CASES_TESTED={cases}")
print(f"CASES_SOLVED={solved}")
print(f"WLM_HITS={len(wlm_hits)}")
print(f"STRICT_WLM_HEM_CC_HITS={len(strict_hits)}")

for h in wlm_hits[:20]:
    print(
        "WLM T={:.1f} pH={:.1f} Pbar={:.0f} XCO2={:.1f} Fe={:.1e} Ca={:.1e} Mg={:.1e} "
        "Wlm={:.3e} hem={:.3e} cc={:.3e} q={:.3e}".format(*h)
    )
