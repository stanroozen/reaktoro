"""
Deep-fluid log10(fO2)-pH Fe stability diagram using Reaktoro.

Uses a per-point equilibrium solve (EquilibriumSolver, not the sweep solver)
to avoid warm-start path-dependence which caused the sweep to report only
Hematite across the entire domain.
"""

import os
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
BASE = os.path.abspath(os.path.dirname(__file__))

_release_dir = os.path.join(REPO, "build", "Reaktoro", "Release")
_package_dir = os.path.join(
    REPO, "build", "python", "package", "build", "lib", "reaktoro"
)

PYD_DIR = None
for _cand in [_release_dir, _package_dir]:
    if os.path.isfile(os.path.join(_cand, "reaktoro4py.cp312-win_amd64.pyd")):
        PYD_DIR = _cand
        break

if PYD_DIR is None:
    raise FileNotFoundError(
        "Could not find reaktoro4py.cp312-win_amd64.pyd in build folders."
    )

if PYD_DIR not in sys.path:
    sys.path.insert(0, PYD_DIR)
os.add_dll_directory(PYD_DIR)

import autodiff as ad
import reaktoro4py as rkt

sys.modules.setdefault("reaktoro", rkt)


# ---------------------------------------------------------------------------
# Species and colours – must match CHNOSZ order as closely as possible
# ---------------------------------------------------------------------------
SPECIES = [
    "Fe+2",
    "Fe+3",
    "FeOH+",
    "FeOH+2",
    "HFeO2(aq)",
    "FeO2-",
    "Hematite",
    "Magnetite",
    "Goethite",
    "Iron",
]

# Colour palette that roughly matches CHNOSZ terrain fill
COLOURS = {
    "Iron": "#f0f0f0",  # near-white (very reducing)
    "Magnetite": "#d4a373",  # tan-orange
    "Hematite": "#e9c46a",  # golden
    "Goethite": "#a8c686",  # green
    "Fe+2": "#4caf50",  # green (aqueous)
    "Fe+3": "#2196f3",  # blue
    "FeOH+": "#9c27b0",  # purple
    "FeOH+2": "#ff9800",  # orange
    "HFeO2(aq)": "#00bcd4",  # cyan
    "FeO2-": "#e91e63",  # pink
}


def make_system_and_specs():
    db = rkt.SupcrtDatabase("supcrtbl")
    aq = rkt.AqueousPhase(
        rkt.speciate(
            [
                "H2O(aq)",
                "H+",
                "OH-",
                "Fe+2",
                "Fe+3",
                "FeOH+",
                "FeOH+2",
                "HFeO2(aq)",
                "FeO2-",
            ]
        )
    )
    gas = rkt.GaseousPhase("O2(g)")
    mins = rkt.MineralPhases(
        rkt.StringList(["Hematite", "Magnetite", "Goethite", "Iron"])
    )
    system = rkt.ChemicalSystem(db, aq, gas, mins)

    specs = rkt.EquilibriumSpecs(system)
    specs.temperature()
    specs.pressure()
    specs.fugacity("O2")
    specs.pH()

    return system, specs


def fresh_state(system, T_C, P_bar):
    state = rkt.ChemicalState(system)
    state.temperature(ad.real(T_C), "celsius")
    state.pressure(ad.real(P_bar), "bar")
    state.set("H2O(aq)", ad.real(1.0), "kg")
    state.set("Fe+2", ad.real(1e-3), "mol")
    state.set("O2(g)", ad.real(1e-8), "mol")
    # Small seed amounts so the solver can easily pick the right phase
    for m in ["Hematite", "Magnetite", "Iron"]:
        state.set(m, ad.real(1e-8), "mol")
    return state


def predominant_fe_species(state):
    """Return the name of the species carrying the most moles of Fe."""
    best_name = None
    best_amount = -1.0
    fe_stoich = {
        "Fe+2": 1,
        "Fe+3": 1,
        "FeOH+": 1,
        "FeOH+2": 1,
        "HFeO2(aq)": 1,
        "FeO2-": 1,
        "Hematite": 2,
        "Magnetite": 3,
        "Goethite": 1,
        "Iron": 1,
    }
    for sp, stoich in fe_stoich.items():
        try:
            amt = float(state.speciesAmount(sp)) * stoich
        except Exception:
            continue
        if amt > best_amount:
            best_amount = amt
            best_name = sp
    return best_name


def main():
    print(f"Using reaktoro4py from: {PYD_DIR}")
    try:
        rkt.Warnings.disable(906)
    except Exception:
        pass

    T_C = 350.0
    P_BAR = 2000.0

    logfO2_vals = np.linspace(-48.0, -18.0, 100)
    pH_vals = np.linspace(3.0, 8.5, 80)

    system, specs = make_system_and_specs()
    conds = rkt.EquilibriumConditions(specs)

    nx = len(logfO2_vals)
    ny = len(pH_vals)

    predominance = np.full((nx, ny), -1, dtype=int)
    n_ok = 0
    n_fail = 0

    total = nx * ny
    print(f"Grid: {nx} x {ny} = {total} points")

    for i, lfo2 in enumerate(logfO2_vals):
        for j, pH in enumerate(pH_vals):
            state = fresh_state(system, T_C, P_BAR)
            solver = rkt.EquilibriumSolver(specs)
            conds.temperature(T_C, "celsius")
            conds.pressure(P_BAR, "bar")
            conds.fugacity("O2", 10.0**lfo2, "bar")
            conds.pH(pH)
            try:
                result = solver.solve(state, conds)
                if result.succeeded():
                    name = predominant_fe_species(state)
                    if name in SPECIES:
                        predominance[i, j] = SPECIES.index(name)
                    n_ok += 1
                else:
                    n_fail += 1
            except Exception:
                n_fail += 1

        if (i + 1) % 10 == 0 or i == nx - 1:
            pct = 100.0 * (i + 1) / nx
            print(f"  row {i + 1}/{nx} ({pct:.0f}%)  ok={n_ok}  fail={n_fail}")

    print(f"Done. ok={n_ok}  fail={n_fail}")

    # ------------------------------------------------------------------
    # Build figure using PredominancePlot from diagrams.py
    # ------------------------------------------------------------------
    import importlib.util as _ilu

    _diagrams_file = os.path.join(
        REPO, "python", "package", "reaktoro", "extensions", "diagrams.py"
    )
    _spec_mod = _ilu.spec_from_file_location("reaktoro_diagrams", _diagrams_file)
    _dmod = _ilu.module_from_spec(_spec_mod)
    _spec_mod.loader.exec_module(_dmod)
    PredominancePlot = _dmod.PredominancePlot

    pp = PredominancePlot(
        logfO2_vals,
        pH_vals,
        predominance,
        SPECIES,
        xlabel=r"$\log_{10} f_{\mathrm{O_2}}$",
        ylabel="pH",
        palette="Set3",
    )
    fig, ax = pp.plot(
        figsize=(8.6, 6.0),
        label_min_fraction=0.003,
        boundary_color="black",
        boundary_linewidth=1.0,
    )

    for label, lfo2 in (("HM-ish", -24.0), ("FMQ-ish", -30.0), ("IW-ish", -36.0)):
        if logfO2_vals[0] <= lfo2 <= logfO2_vals[-1]:
            ax.axvline(lfo2, color="gray", linestyle=":", linewidth=0.9)
            ax.text(lfo2 + 0.2, pH_vals[-1] - 0.15, label, fontsize=8, color="gray")

    ax.set_title(
        f"Reaktoro deep-fluid Fe stability: log$_{{10}}$(f$_{{O_2}}$) vs pH\n"
        f"T = {T_C:.0f} C, P = {P_BAR:.0f} bar"
    )
    plt.tight_layout()

    out_png = os.path.join(BASE, "Reaktoro_DeepFluid_LogfO2_pH_Fe.png")
    fig.savefig(out_png, dpi=180)
    plt.close(fig)

    out_txt = os.path.join(BASE, "deepfluid_logfo2_ph_setup.txt")
    with open(out_txt, "w", encoding="utf-8") as f:
        f.write("Deep-fluid potential diagram setup (Reaktoro)\n")
        f.write(f"- Binary path: {PYD_DIR}\n")
        f.write(f"- T = {T_C} C\n")
        f.write(f"- P = {P_BAR} bar\n")
        f.write(f"- log10(fO2) range: [{logfO2_vals[0]}, {logfO2_vals[-1]}] bar\n")
        f.write(f"- pH range: [{pH_vals[0]}, {pH_vals[-1]}]\n")
        f.write(f"- Grid: {nx} x {ny} points\n")
        f.write(f"- Solved: ok={n_ok}, fail={n_fail}\n")
        f.write("- Method: per-point EquilibriumSolver (no sweep warm-start)\n")
        f.write(f"- Output: {out_png}\n")

    print("Wrote:", out_png)
    print("Wrote:", out_txt)


if __name__ == "__main__":
    main()
