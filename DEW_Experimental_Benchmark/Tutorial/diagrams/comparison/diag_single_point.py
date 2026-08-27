"""Single-point equilibrium check to diagnose why Reaktoro shows only Goethite."""

import sys, os, numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
PYD = os.path.join(REPO, "build", "python", "package", "build", "lib", "reaktoro")
sys.path.insert(0, PYD)
os.add_dll_directory(PYD)
if "autodiff" in sys.modules:
    del sys.modules["autodiff"]
import autodiff as ad
import reaktoro4py as rkt

rkt.Warnings.disable(906)

db = rkt.SupcrtDatabase("supcrtbl")
aq = rkt.AqueousPhase(
    rkt.speciate(
        [
            "H2O(aq)",
            "H+",
            "OH-",
            "Fe+2",
            "Fe+3",
            "FeO+",
            "FeO2-",
            "FeOH+",
            "FeOH+2",
            "HFeO2(aq)",
            "HFeO2-",
        ]
    )
)
system = rkt.ChemicalSystem(
    db,
    aq,
    rkt.MineralPhases(rkt.StringList(["Goethite", "Hematite", "Iron", "Magnetite"])),
)

specs = rkt.EquilibriumSpecs(system)
specs.temperature()
specs.pressure()
specs.pH()
specs.Eh()
# Use sweep solver over a tiny diagnostic grid
state = rkt.ChemicalState(system)
solver = rkt.EquilibriumSweepSolver(specs)

pH_test = np.array([0.0, 2.0, 5.0, 5.0, 1.0], dtype=float)
Eh_test = np.array([0.0, 0.0, -0.5, 0.5, -0.5], dtype=float)

aq_sp = ["Fe+2", "Fe+3", "FeOH+", "HFeO2-"]
min_sp = ["Goethite", "Hematite", "Iron", "Magnetite"]
species_all = aq_sp + min_sp

# sweepPHEhGrid expects 2-D grid; use 1-D pH/Eh arrays for 5 points along pH axis
# Workaround: sweep on a 5x1 grid (constant Eh = 0) won't work for variable Eh.
# Instead, run 5 separate single-row sweeps.
for i, (pH, Eh) in enumerate(zip(pH_test, Eh_test)):
    ph_arr = np.array([pH], dtype=float)
    eh_arr = np.array([Eh], dtype=float)
    fresh = rkt.ChemicalState(system)
    fresh.set("H2O(aq)", ad.real(55.5), "mol")
    fresh.set("Fe+2", ad.real(1e-6), "mol")
    grid = solver.sweepPHEhGrid(fresh, ph_arr, eh_arr, "V")
    states = list(grid.states)
    s = states[0]
    props = rkt.ChemicalProps(s)
    print(f"pH={pH:5.1f} Eh={Eh:5.2f}")
    for sp in species_all:
        try:
            amt = float(props.speciesAmount(sp))
            if amt > 1e-20:
                print(f"  {sp:<18s}: {amt:.3e} mol")
        except Exception:
            pass
