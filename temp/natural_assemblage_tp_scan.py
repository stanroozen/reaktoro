import importlib.util
import os
import sys
import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
rel = os.path.join(REPO, "build", "Reaktoro", "Release")
if rel not in sys.path:
    sys.path.insert(0, rel)

TUTORIAL = os.path.join(
    REPO,
    "DEW_Experimental_Benchmark",
    "Tutorial",
    "willemite_solubility_tutorial_dew17hp622_zn.py",
)
spec = importlib.util.spec_from_file_location("w", TUTORIAL)
w = importlib.util.module_from_spec(spec)
spec.loader.exec_module(w)

# User-observed minerals (dolomite unavailable in this DB; use calcite + magnesite proxy)
obs = ["Wlm", "Znc", "hem", "cc", "q", "mag"]

# Force mineral set to include Zn phases + observed gangue
w.USE_COMPETING_ZN_MINERALS = True
w.INCLUDE_COMPATIBLE_GANGUE_MINERALS = False
w.COMPETING_ZN_MINERALS = list(
    dict.fromkeys(w.COMPETING_ZN_MINERALS + ["hem", "cc", "q", "mag"])
)

# Add aqueous species needed to allow Fe/Ca/Mg exchange with fluid
extra_aq = [
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

# Seed observed gangue in initial state
for m in ["hem", "cc", "q", "mag"]:
    w.INITIAL_SPECIES_AMOUNTS_MOL[m] = 2.0

w.validate_user_inputs()
db = w.Database.fromFile(w.PERPLEX_DATABASE_FILE)
aq = w.AqueousPhase(" ".join(w.AQUEOUS_SPECIES))
aq.setActivityModel(w.AQUEOUS_ACTIVITY_MODEL())
mins = w.make_mineral_phases()
gas = w.GaseousPhase("O2 CO2")
gas.setActivityModel(w.ActivityModelIdealGas())
sysm = w.ChemicalSystem(db, aq, mins, gas)

specs = w.EquilibriumSpecs(sysm)
specs.temperature()
specs.pressure()
specs.pH()
specs.fugacity("O2")
specs.fugacity("CO2")
solver = w.make_equilibrium_solver(sysm, specs)
opts = w.EquilibriumOptions()
if hasattr(opts, "hessian") and hasattr(w, "GibbsHessian"):
    opts.hessian = w.GibbsHessian.Exact
solver.setOptions(opts)
conditions = w.EquilibriumConditions(specs)

Tvals = np.linspace(150.0, 400.0, 7)
pHvals = np.linspace(4.0, 10.0, 9)
logfO2_vals = [-30.0, -25.0, -20.0, -15.0]
logfCO2_vals = [-4.0, -3.0, -2.0, -1.0]
Pbar = 2000.0
TH = 1e-8

hits_strict = []
hits_relaxed = []
hits_relaxed_with_znc = []
solved = 0
for T in Tvals:
    for pH in pHvals:
        for lfO2 in logfO2_vals:
            for lfCO2 in logfCO2_vals:
                st = w.make_base_state(sysm)
                try:
                    st.set("O2", 1.0e-20, "mol")
                    st.set("CO2", 1.0e-20, "mol")
                except Exception:
                    pass
                conditions.temperature(float(T), "celsius")
                conditions.pressure(float(Pbar), "bar")
                conditions.pH(float(pH))
                conditions.fugacity("O2", float(10.0**lfO2), "bar")
                conditions.fugacity("CO2", float(10.0**lfCO2), "bar")
                r = solver.solve(st, conditions)
                if not r.succeeded():
                    continue
                solved += 1

        def present(name):
            try:
                return float(st.speciesAmount(name)) > TH
            except Exception:
                return False

                # Strict observed set with mag as dolomite proxy
                if all(present(m) for m in obs):
                    hits_strict.append((float(T), float(pH), float(lfO2), float(lfCO2)))

                # Relaxed: Wlm + Hem + Q + (Calcite or Magnesite proxy)
                if (
                    present("Wlm")
                    and present("hem")
                    and present("q")
                    and (present("cc") or present("mag"))
                ):
                    hits_relaxed.append(
                        (float(T), float(pH), float(lfO2), float(lfCO2))
                    )
                    if present("Znc"):
                        hits_relaxed_with_znc.append(
                            (float(T), float(pH), float(lfO2), float(lfCO2))
                        )


print(
    f"SOLVED={solved}/{len(Tvals) * len(pHvals) * len(logfO2_vals) * len(logfCO2_vals)}"
)
print(f"HITS_STRICT_WLM+ZNC+HEM+CC+Q+MAG={len(hits_strict)}")
for h in hits_strict:
    print(f"T_C={h[0]:.1f}, pH={h[1]:.2f}, logfO2={h[2]:.1f}, logfCO2={h[3]:.1f}")

print(f"HITS_RELAXED_WLM+HEM+Q+(CC|MAG)={len(hits_relaxed)}")
for h in hits_relaxed:
    print(
        f"RELAXED T_C={h[0]:.1f}, pH={h[1]:.2f}, logfO2={h[2]:.1f}, logfCO2={h[3]:.1f}"
    )

print(f"HITS_RELAXED_WITH_ZNC_OPTIONAL={len(hits_relaxed_with_znc)}")
for h in hits_relaxed_with_znc:
    print(
        f"RELAXED_ZNC T_C={h[0]:.1f}, pH={h[1]:.2f}, logfO2={h[2]:.1f}, logfCO2={h[3]:.1f}"
    )
