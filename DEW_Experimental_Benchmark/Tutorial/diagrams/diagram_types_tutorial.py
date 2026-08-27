"""
Diagram Types Tutorial — Reaktoro extensions.diagrams
======================================================

Demonstrates all seven diagram classes:
    1. SpeciationPlot       — Carbonate speciation vs pH
    2. PredominancePlot     — Iron Pourbaix (Eh–pH)
    3. SolubilityPlot       — Calcite Ca solubility on T–P grid
    4. ActivityDiagram      — Calcite stability in log a(Ca²⁺)–log a(CO₃²⁻)
    5. MosaicPlot           — Fe minerals + Fe aq species overlaid
    6. LogfO2pHDiagram      — Fe oxide stability vs log fO2 and pH
    7. TPDiagram            — Calcite/Aragonite polymorphism

All calculations use SupcrtDatabase("supcrtbl") and updated per-point equilibrium recipes.
Output PNG files are saved in the same folder as this script.

Usage:
    cd <repo-root>/DEW_Experimental_Benchmark/Tutorial/diagrams
    python diagram_types_tutorial.py
"""

import os
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")  # headless backend — works without a display
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# ---------------------------------------------------------------------------
# Locate repository root and add local build / local diagrams package.
# ---------------------------------------------------------------------------
THIS_DIR = os.path.dirname(os.path.abspath(__file__))  # .../diagrams/
REPO_ROOT = os.path.abspath(os.path.join(THIS_DIR, "..", "..", ".."))  # repo root

# The local build keeps Reaktoro.dll and reaktoro4py.pyd in the same directory.
# We try that first, registering the DLL directory so Windows can load it.
_LOCAL_PYD_DIR = os.path.join(
    REPO_ROOT, "build", "python", "package", "build", "lib", "reaktoro"
)
_FALLBACK_PYD_DIR = os.path.join(REPO_ROOT, "build", "Reaktoro", "Release")

_reaktoro_loaded = False
for _candidate in [_FALLBACK_PYD_DIR, _LOCAL_PYD_DIR]:
    if not os.path.isdir(_candidate):
        continue
    if _candidate not in sys.path:
        sys.path.insert(0, _candidate)
    # Register the directory so Windows finds Reaktoro.dll.
    try:
        os.add_dll_directory(_candidate)
    except (AttributeError, OSError):
        pass  # Python < 3.8 or directory not found

    try:
        import autodiff  # shipped alongside the pyd in some build layouts
    except ImportError:
        pass

    try:
        import reaktoro4py as _r4p
        from reaktoro4py import *

        # diagrams.py expects `import reaktoro` to succeed when building
        # ChemicalProps from sweep states. Alias local reaktoro4py for this run.
        sys.modules.setdefault("reaktoro", _r4p)

        print(f"Using local Reaktoro build from: {_candidate}")
        _reaktoro_loaded = True
        break
    except (ModuleNotFoundError, ImportError):
        pass

if not _reaktoro_loaded:
    from reaktoro import *

    print("Using installed 'reaktoro' package.")
    print("WARNING: EquilibriumSweepSolver is only available in the local build.")
    print("         Build Reaktoro from source to run this script.")

try:
    Warnings.disable(906)
except Exception:
    pass

# Load diagrams.py directly by path so the installed reaktoro package doesn't
# shadow the local extension module.
import importlib.util as _ilu

_diagrams_file = os.path.join(
    REPO_ROOT, "python", "package", "reaktoro", "extensions", "diagrams.py"
)
if not os.path.isfile(_diagrams_file):
    raise FileNotFoundError(
        f"diagrams.py not found at: {_diagrams_file}\n"
        "Check that REPO_ROOT points to the repository root."
    )

_spec = _ilu.spec_from_file_location("reaktoro_diagrams", _diagrams_file)
_dmod = _ilu.module_from_spec(_spec)
_spec.loader.exec_module(_dmod)

SpeciationPlot = _dmod.SpeciationPlot
PredominancePlot = _dmod.PredominancePlot
SolubilityPlot = _dmod.SolubilityPlot
ActivityDiagram = _dmod.ActivityDiagram
MosaicPlot = _dmod.MosaicPlot
LogfO2pHDiagram = _dmod.LogfO2pHDiagram
TPDiagram = _dmod.TPDiagram
water_lines = _dmod.water_lines
saturation_curve = _dmod.saturation_curve

print("Diagram classes imported successfully.\n")

# ---------------------------------------------------------------------------
# Shared database
# ---------------------------------------------------------------------------
db = SupcrtDatabase("supcrtbl")


def env_int(name, default):
    try:
        return int(os.environ.get(name, str(default)))
    except Exception:
        return int(default)


def list_species(element_symbol, aggregate_states=None):
    """Return species names in db that contain element_symbol."""
    results = []
    for sp in db.species():
        syms = [pair[0].symbol() for pair in sp.elements()]
        if element_symbol in syms:
            if aggregate_states is None or sp.aggregateState() in aggregate_states:
                results.append(sp.name())
    return sorted(results)


class GridResultProxy:
    def __init__(self, xvalues, yvalues, grid):
        self.xvalues = np.asarray(xvalues, dtype=float)
        self.yvalues = np.asarray(yvalues, dtype=float)
        self._grid = np.asarray(grid, dtype=float)

    def predominantSpeciesGrid(self, species):
        return self._grid

    def elementMolalityGrid(self, element):
        return self._grid


def _fresh_state(system, seed_amounts):
    state = ChemicalState(system)
    for species_name, amount, unit in seed_amounts:
        try:
            state.set(species_name, amount, unit)
        except Exception:
            pass
    return state


def _species_indices(system, species_names):
    species_list = system.species()
    lookup = {}
    for i in range(species_list.size()):
        lookup[species_list[i].name()] = i

    indices = Indices()
    for species_name in species_names:
        if species_name in lookup:
            indices.push_back(lookup[species_name])
    return indices


def _predominant_by_amount(state, species_names):
    props = ChemicalProps(state)
    best_idx = -1
    best_amount = -np.inf
    for idx, species_name in enumerate(species_names):
        try:
            amount = float(props.speciesAmount(species_name))
        except Exception:
            continue
        if np.isfinite(amount) and amount > best_amount:
            best_amount = amount
            best_idx = idx
    return best_idx


def _predominant_fe_field(state, aqueous_species, mineral_species):
    props = ChemicalProps(state)

    best_min_idx = -1
    best_min_amount = 1e-12
    for idx, species_name in enumerate(mineral_species):
        try:
            amount = float(props.speciesAmount(species_name))
        except Exception:
            amount = 0.0
        if amount > best_min_amount:
            best_min_amount = amount
            best_min_idx = len(aqueous_species) + idx

    if best_min_idx >= 0:
        return best_min_idx

    best_aq_idx = -1
    best_aq_lga = -np.inf
    for idx, species_name in enumerate(aqueous_species):
        try:
            lga = float(props.speciesActivityLg(species_name))
        except Exception:
            lga = np.nan
        if np.isfinite(lga) and lga > best_aq_lga:
            best_aq_lga = lga
            best_aq_idx = idx

    return best_aq_idx


# ===========================================================================
# Section 1 — SpeciationPlot: Carbonate speciation vs pH
# ===========================================================================
print("Section 1: SpeciationPlot — carbonate speciation vs pH")

carbonate_aq = ["H2O(aq)", "H+", "OH-", "CO2(aq)", "HCO3-", "CO3-2"]

system_carb = ChemicalSystem(db, AqueousPhase(speciate(carbonate_aq)))

specs_ph = EquilibriumSpecs(system_carb)
specs_ph.temperature()
specs_ph.pressure()
specs_ph.pH()

state_carb = ChemicalState(system_carb)
state_carb.temperature(25.0, "celsius")
state_carb.pressure(1.0, "bar")
state_carb.set("H2O(aq)", 1.0, "kg")
state_carb.set("HCO3-", 0.01, "mol")

solver_ph = EquilibriumSweepSolver(specs_ph)

pH_values = np.linspace(2.0, 12.0, max(30, env_int("DIAG_SEC1_NPH", 140)))
result_carb = solver_ph.sweepPH(
    state_carb,
    pH_values,
)

sp = SpeciationPlot.from_sweep_result(
    result_carb,
    species=["CO2(aq)", "HCO3-", "CO3-2"],
    xlabel="pH",
    xvalues=pH_values,
    palette="tab10",
)

fig1, ax1 = sp.plot(figsize=(7, 4), ylim=(-14, 0))
ax1.set_title("Carbonate speciation at 25 °C, 1 bar (10 mM total C)", fontsize=11)
ax1.set_ylabel(r"$\log_{10}(a_i)$ of aqueous species $i$ (dimensionless)")
ax1.text(
    0.01,
    0.02,
    r"$a_i$ is activity (effective concentration) relative to the standard state",
    transform=ax1.transAxes,
    fontsize=8,
    alpha=0.9,
)
plt.tight_layout()
out1 = os.path.join(THIS_DIR, "Section1_SpeciationPlot_carbonate.png")
fig1.savefig(out1, dpi=150)
plt.close(fig1)
print(f"  Saved: {out1}\n")


# ===========================================================================
# Section 2 — PredominancePlot: Iron Pourbaix (Eh–pH)
# ===========================================================================
print("Section 2: PredominancePlot — Fe Pourbaix diagram")

fe_aq_oh = [
    "H2O(aq)",
    "H+",
    "OH-",
    "e-",
    "Fe+2",
    "Fe+3",
    "FeO(aq)",
    "FeO+",
    "FeO2-",
    "FeOH+",
    "FeOH+2",
    "HFeO2(aq)",
    "HFeO2-",
]
fe_min = ["Hematite", "Magnetite", "Goethite", "Iron"]
fe_min_with_aqueous_only = fe_min + ["Aqueous_only"]
aq_species = fe_aq_oh[4:]

system_fe = ChemicalSystem(
    db,
    AqueousPhase(speciate(fe_aq_oh)),
    GaseousPhase("O2(g)"),
    MineralPhases(StringList(fe_min)),
)

specs_pourbaix = EquilibriumSpecs(system_fe)
specs_pourbaix.temperature()
specs_pourbaix.pressure()
specs_pourbaix.lgActivity("H+")
specs_pourbaix.Eh()

conds_pourbaix = EquilibriumConditions(specs_pourbaix)
conds_pourbaix.set("T", 298.15)
conds_pourbaix.set("P", 1.0e5)
solver_pourbaix = EquilibriumSolver(specs_pourbaix)

seed_amounts = [
    ("H2O(aq)", 1.0, "kg"),
    ("Fe+2", 1e-6, "mol"),
    ("O2(g)", 1e-8, "mol"),
]

pH_grid = np.linspace(0.0, 14.0, max(30, env_int("DIAG_SEC2_NPH", 70)))
Eh_grid = np.linspace(-1.0, 1.2, max(30, env_int("DIAG_SEC2_NEH", 70)))
pred_fe = np.full((len(pH_grid), len(Eh_grid)), -1, dtype=int)
pred_fe_mineral_only = np.full((len(pH_grid), len(Eh_grid)), -1, dtype=int)
pred_fe_aq_only = np.full((len(pH_grid), len(Eh_grid)), -1, dtype=int)
sec2_mineral_stability_tol = 1e-14
for i, pH in enumerate(pH_grid):
    for j, Eh in enumerate(Eh_grid):
        state = _fresh_state(system_fe, seed_amounts)
        conds_pourbaix.set("T", 298.15)
        conds_pourbaix.set("P", 1.0e5)
        conds_pourbaix.set("ln(a[H+])", float(np.log(10.0) * -pH))
        conds_pourbaix.set("Eh", float(Eh))
        try:
            result = solver_pourbaix.solve(state, conds_pourbaix)
            if result.succeeded():
                pred_fe[i, j] = _predominant_fe_field(state, aq_species, fe_min)
                pred_fe_aq_only[i, j] = _predominant_by_amount(state, aq_species)
                mineral_amounts = [state.speciesAmount(s).val() for s in fe_min]
                if max(mineral_amounts) <= sec2_mineral_stability_tol:
                    pred_fe_mineral_only[i, j] = len(fe_min)
                else:
                    pred_fe_mineral_only[i, j] = int(np.argmax(mineral_amounts))
        except Exception:
            pass

pourbaix_species = aq_species + fe_min
pp = PredominancePlot(
    pH_grid,
    Eh_grid,
    pred_fe,
    pourbaix_species,
    xlabel="pH",
    ylabel="Eh",
    palette="Set3",
)

fig2, ax2 = pp.plot(
    figsize=(8, 5),
    label_min_fraction=0.004,
    boundary_color="black",
    boundary_linewidth=1.0,
)
mineral_flag = np.where(
    pred_fe >= 0, (pred_fe >= len(aq_species)).astype(float), np.nan
)
if np.any(np.isfinite(mineral_flag)):
    X2, Y2 = np.meshgrid(pH_grid, Eh_grid)
    try:
        ax2.contour(
            X2,
            Y2,
            mineral_flag.T,
            levels=[0.5],
            colors="0.35",
            linewidths=1.2,
            linestyles="--",
        )
    except Exception:
        pass

pp.add_water_lines(ax2, T_K=298.15, color="black", linestyle="-.", linewidth=1.0)
ax2.axhline(0.0, color="0.55", linestyle=":", linewidth=1.0)
ax2.axvline(7.0, color="0.55", linestyle=":", linewidth=1.0)
ax2.legend(
    handles=[
        Line2D(
            [0],
            [0],
            color="black",
            lw=1.4,
            linestyle="-",
            label="Predominance boundary (solid)",
        ),
        Line2D(
            [0],
            [0],
            color="0.35",
            lw=1.2,
            linestyle="--",
            label="Mineral-stability frontier (dashed)",
        ),
        Line2D(
            [0],
            [0],
            color="black",
            lw=1.0,
            linestyle="-.",
            label="Water-stability lines (dash-dot)",
        ),
        Line2D(
            [0],
            [0],
            color="0.55",
            lw=1.0,
            linestyle=":",
            label="Reference guides: Eh = 0, pH = 7",
        ),
    ],
    loc="lower right",
    fontsize=8,
    framealpha=0.9,
)
ax2.set_title("Fe–O–H Pourbaix diagram at 25 °C, 1 bar", fontsize=11)
plt.tight_layout()
out2 = os.path.join(THIS_DIR, "Section2_PredominancePlot_Fe_Pourbaix.png")
fig2.savefig(out2, dpi=150)
plt.close(fig2)
print(f"  Saved: {out2}\n")

pp2_min_only = PredominancePlot(
    pH_grid,
    Eh_grid,
    pred_fe_mineral_only,
    fe_min_with_aqueous_only,
    xlabel="pH",
    ylabel="Eh",
    palette="YlOrBr",
)
fig2_min, ax2_min = pp2_min_only.plot(
    figsize=(8, 5),
    label_min_fraction=0.004,
    boundary_color="black",
    boundary_linewidth=1.0,
)
pp2_min_only.add_water_lines(
    ax2_min, T_K=298.15, color="black", linestyle="-.", linewidth=1.0
)
ax2_min.axhline(0.0, color="0.55", linestyle=":", linewidth=1.0)
ax2_min.axvline(7.0, color="0.55", linestyle=":", linewidth=1.0)
ax2_min.set_title("Fe–O–H Pourbaix mineral-only at 25 °C, 1 bar", fontsize=11)
plt.tight_layout()
out2_min = os.path.join(
    THIS_DIR, "Section2_PredominancePlot_Fe_Pourbaix_Mineral_only.png"
)
fig2_min.savefig(out2_min, dpi=150)
plt.close(fig2_min)
print(f"  Saved: {out2_min}\n")

pp2_aq_only = PredominancePlot(
    pH_grid,
    Eh_grid,
    pred_fe_aq_only,
    aq_species,
    xlabel="pH",
    ylabel="Eh",
    palette="tab20",
)
fig2_aq, ax2_aq = pp2_aq_only.plot(
    figsize=(8, 5),
    label_min_fraction=0.004,
    boundary_color="midnightblue",
    boundary_linewidth=1.0,
)
pp2_aq_only.add_water_lines(
    ax2_aq, T_K=298.15, color="black", linestyle="-.", linewidth=1.0
)
ax2_aq.axhline(0.0, color="0.55", linestyle=":", linewidth=1.0)
ax2_aq.axvline(7.0, color="0.55", linestyle=":", linewidth=1.0)
ax2_aq.set_title("Fe–O–H Pourbaix aqueous-only at 25 °C, 1 bar", fontsize=11)
plt.tight_layout()
out2_aq = os.path.join(
    THIS_DIR, "Section2_PredominancePlot_Fe_Pourbaix_Aqueous_only.png"
)
fig2_aq.savefig(out2_aq, dpi=150)
plt.close(fig2_aq)
print(f"  Saved: {out2_aq}\n")


# ===========================================================================
# Section 3 — SolubilityPlot: Calcite Ca solubility on T–P grid
# ===========================================================================
print("Section 3: SolubilityPlot — calcite Ca solubility vs T and P")

calcite_aq = [
    "H2O(aq)",
    "H+",
    "OH-",
    "Ca+2",
    "CaOH+",
    "Ca(HCO3)+",
    "CaCO3(aq)",
    "CO2(aq)",
    "HCO3-",
    "CO3-2",
]

system_calc = ChemicalSystem(
    db,
    AqueousPhase(speciate(calcite_aq)),
    MineralPhases(StringList(["Calcite"])),
)

specs_tp = EquilibriumSpecs(system_calc)
specs_tp.temperature()
specs_tp.pressure()

conds_tp = EquilibriumConditions(specs_tp)
solver_tp = EquilibriumSolver(specs_tp)

seed_amounts = [
    ("H2O(aq)", 1.0, "kg"),
    ("Calcite", 10.0, "mol"),
]

T_vals = np.linspace(25.0, 300.0, max(25, env_int("DIAG_SEC3_NT", 55)))
P_vals = np.linspace(1.0, 1000.0, max(25, env_int("DIAG_SEC3_NP", 55)))
ca_grid = np.full((len(T_vals), len(P_vals)), np.nan, dtype=float)
water_indices = _species_indices(system_calc, ["H2O(aq)"])
ca_indices = _species_indices(system_calc, ["Ca+2", "CaOH+", "Ca(HCO3)+", "CaCO3(aq)"])
for i, T_C in enumerate(T_vals):
    for j, P_bar in enumerate(P_vals):
        state = _fresh_state(system_calc, seed_amounts)
        conds_tp.set("T", float(T_C + 273.15))
        conds_tp.set("P", float(P_bar * 1e5))
        try:
            result = solver_tp.solve(state, conds_tp)
            if result.succeeded():
                props = ChemicalProps(state)
                ca_amount = float(props.elementAmountAmongSpecies("Ca", ca_indices))
                water_moles = (
                    float(props.elementAmountAmongSpecies("H", water_indices)) / 2.0
                )
                water_mass_kg = water_moles * 18.01528 / 1000.0
                if water_mass_kg > 0:
                    ca_grid[i, j] = ca_amount / water_mass_kg
        except Exception:
            pass

grid_calc_tp = GridResultProxy(T_vals, P_vals, ca_grid)
sol = SolubilityPlot.from_grid_result(
    grid_calc_tp,
    element="Ca",
    xlabel="T (°C)",
    ylabel="P (bar)",
)

fig3, ax3 = sol.plot(
    figsize=(7, 5),
    levels=25,
    iso_levels=[-3.0, -2.5, -2.0, -1.5, -1.0],
    cmap="plasma_r",
)
ax3.set_title(
    "Calcite solubility: log₁₀ [Ca] (mol/kg), 25–300 °C, 1–1000 bar", fontsize=10
)
plt.tight_layout()
out3 = os.path.join(THIS_DIR, "Section3_SolubilityPlot_calcite_TP.png")
fig3.savefig(out3, dpi=150)
plt.close(fig3)
print(f"  Saved: {out3}\n")


# ===========================================================================
# Section 4 — ActivityDiagram: Calcite in log a(Ca²⁺)–log a(CO₃²⁻) space
# ===========================================================================
print("Section 4: ActivityDiagram — calcite stability (log-activity axes)")

specs_act = EquilibriumSpecs(system_calc)
specs_act.temperature()
specs_act.pressure()
specs_act.lgActivity("Ca+2")
specs_act.lgActivity("CO3-2")

conds_act = EquilibriumConditions(specs_act)
conds_act.set("T", 298.15)
conds_act.set("P", 1.0e5)
solver_act = EquilibriumSolver(specs_act)

# Include all aqueous carbonate/calcium species in this system plus the mineral
# so the activity diagram reflects the full modeled assemblage.
stability_species = [
    "Ca+2",
    "CaOH+",
    "Ca(HCO3)+",
    "CaCO3(aq)",
    "CO2(aq)",
    "HCO3-",
    "CO3-2",
    "Calcite",
]
seed_amounts = [
    ("H2O(aq)", 1.0, "kg"),
    ("Ca+2", 1e-4, "mol"),
    ("CO3-2", 1e-4, "mol"),
]

lga_ca = np.linspace(-6.0, 0.0, max(25, env_int("DIAG_SEC4_NX", 55)))
lga_co3 = np.linspace(-6.0, 0.0, max(25, env_int("DIAG_SEC4_NY", 55)))
pred_act = np.full((len(lga_ca), len(lga_co3)), -1, dtype=int)
for i, lga_x in enumerate(lga_ca):
    for j, lga_y in enumerate(lga_co3):
        state = _fresh_state(system_calc, seed_amounts)
        conds_act.set("T", 298.15)
        conds_act.set("P", 1.0e5)
        conds_act.set("ln(a[Ca+2])", float(np.log(10.0) * lga_x))
        conds_act.set("ln(a[CO3-2])", float(np.log(10.0) * lga_y))
        try:
            result = solver_act.solve(state, conds_act)
            if result.succeeded():
                pred_act[i, j] = _predominant_by_amount(state, stability_species)
        except Exception:
            pass

ad_proxy = GridResultProxy(lga_ca, lga_co3, pred_act)
ad = ActivityDiagram.from_grid_result(
    ad_proxy,
    stability_species,
    speciesX="Ca+2",
    speciesY="CO3-2",
    palette="Pastel1",
)

fig4, ax4 = ad.plot(figsize=(6, 5))
ax4.set_xlabel("log10 a(Ca+2)")
ax4.set_ylabel("log10 a(CO3-2)")
ax4.set_title("Calcite-carbonate activity predominance at 25 °C, 1 bar", fontsize=11)
plt.tight_layout()
out4 = os.path.join(THIS_DIR, "Section4_ActivityDiagram_calcite.png")
fig4.savefig(out4, dpi=150)
plt.close(fig4)
print(f"  Saved: {out4}\n")


# ===========================================================================
# Section 5 — MosaicPlot: Fe minerals + Fe aq species overlaid
# ===========================================================================
print("Section 5: MosaicPlot — Fe minerals + aqueous Fe species overlay")

fe_aq_mosaic = [
    "H2O(aq)",
    "H+",
    "OH-",
    "e-",
    "Fe+2",
    "Fe+3",
    "FeO(aq)",
    "FeO+",
    "FeO2-",
    "FeOH+",
    "FeOH+2",
    "HFeO2(aq)",
    "HFeO2-",
]
mineral_layer_species = ["Hematite", "Magnetite", "Goethite", "Iron"]
mineral_layer_species_with_aqueous_only = mineral_layer_species + ["Aqueous_only"]
aqueous_layer_species = fe_aq_mosaic[4:]

system_mosaic = ChemicalSystem(
    db,
    AqueousPhase(speciate(fe_aq_mosaic)),
    MineralPhases(StringList(mineral_layer_species)),
)

specs_mosaic = EquilibriumSpecs(system_mosaic)
specs_mosaic.temperature()
specs_mosaic.pressure()
specs_mosaic.lgActivity("H+")
specs_mosaic.Eh()

conds_mosaic = EquilibriumConditions(specs_mosaic)
conds_mosaic.set("T", 298.15)
conds_mosaic.set("P", 1.0e5)
solver_mosaic = EquilibriumSolver(specs_mosaic)

seed_amounts = [
    ("H2O(aq)", 1.0, "kg"),
    ("Fe+2", 1e-6, "mol"),
]

pH_grid = np.linspace(0.0, 14.0, max(30, env_int("DIAG_SEC2_NPH", 70)))
Eh_grid = np.linspace(-1.0, 1.2, max(30, env_int("DIAG_SEC2_NEH", 70)))
min_predominance = np.full((len(pH_grid), len(Eh_grid)), -1, dtype=int)
aq_predominance = np.full((len(pH_grid), len(Eh_grid)), -1, dtype=int)
mineral_stability_tol = 1e-14
for i, pH in enumerate(pH_grid):
    for j, Eh in enumerate(Eh_grid):
        state = _fresh_state(system_mosaic, seed_amounts)
        conds_mosaic.set("T", 298.15)
        conds_mosaic.set("P", 1.0e5)
        conds_mosaic.set("ln(a[H+])", float(np.log(10.0) * -pH))
        conds_mosaic.set("Eh", float(Eh))
        try:
            result = solver_mosaic.solve(state, conds_mosaic)
            if result.succeeded():
                mineral_amounts = [
                    state.speciesAmount(s).val() for s in mineral_layer_species
                ]
                if max(mineral_amounts) <= mineral_stability_tol:
                    min_predominance[i, j] = len(mineral_layer_species)
                else:
                    min_predominance[i, j] = int(np.argmax(mineral_amounts))
                aq_predominance[i, j] = _predominant_by_amount(
                    state, aqueous_layer_species
                )
        except Exception:
            pass

layers = [
    {
        "label": "Minerals",
        "species": mineral_layer_species_with_aqueous_only,
        "predominance": min_predominance,
        "palette": "YlOrBr",
        "alpha": 0.92,
        "boundary_color": "black",
        "boundary_linewidth": 1.5,
        "boundary_linestyle": "-",
        "draw_edges": True,
    },
    {
        "label": "Aqueous species",
        "species": aqueous_layer_species,
        "predominance": aq_predominance,
        "palette": "tab20",
        "alpha": 0.60,
        "boundary_color": "midnightblue",
        "boundary_linewidth": 0.9,
        "boundary_linestyle": "--",
    },
]

mp = MosaicPlot(pH_grid, Eh_grid, layers, xlabel="pH", ylabel="Eh")
fig5, ax5 = mp.plot(figsize=(8, 5))
mp.add_water_lines(ax5, T_K=298.15, color="navy", linestyle="--", linewidth=1.2)
ax5.set_title("Fe–O–H mosaic diagram at 25 °C, 1 bar", fontsize=11)
plt.tight_layout()
out5 = os.path.join(THIS_DIR, "Section5_MosaicPlot_Fe.png")
fig5.savefig(out5, dpi=150)
plt.close(fig5)
print(f"  Saved: {out5}\n")

layers_mineral_only = [
    {
        "label": "Minerals",
        "species": mineral_layer_species_with_aqueous_only,
        "predominance": min_predominance,
        "palette": "YlOrBr",
        "alpha": 0.95,
        "boundary_color": "black",
        "boundary_linewidth": 1.3,
        "boundary_linestyle": "-",
        "draw_edges": True,
    }
]

mp_min_only = MosaicPlot(
    pH_grid, Eh_grid, layers_mineral_only, xlabel="pH", ylabel="Eh"
)
fig5_min, ax5_min = mp_min_only.plot(figsize=(8, 5))
mp_min_only.add_water_lines(
    ax5_min, T_K=298.15, color="navy", linestyle="--", linewidth=1.2
)
ax5_min.set_title("Fe–O–H mineral-only mosaic at 25 °C, 1 bar", fontsize=11)
plt.tight_layout()
out5_min = os.path.join(THIS_DIR, "Section5_MosaicPlot_Fe_Mineral_only.png")
fig5_min.savefig(out5_min, dpi=150)
plt.close(fig5_min)
print(f"  Saved: {out5_min}\n")

layers_aqueous_only = [
    {
        "label": "Aqueous species",
        "species": aqueous_layer_species,
        "predominance": aq_predominance,
        "palette": "tab20",
        "alpha": 0.92,
        "boundary_color": "midnightblue",
        "boundary_linewidth": 1.0,
        "boundary_linestyle": "--",
    }
]

mp_aq_only = MosaicPlot(pH_grid, Eh_grid, layers_aqueous_only, xlabel="pH", ylabel="Eh")
fig5_aq, ax5_aq = mp_aq_only.plot(figsize=(8, 5))
mp_aq_only.add_water_lines(
    ax5_aq, T_K=298.15, color="navy", linestyle="--", linewidth=1.2
)
ax5_aq.set_title("Fe–O–H aqueous-only mosaic at 25 °C, 1 bar", fontsize=11)
plt.tight_layout()
out5_aq = os.path.join(THIS_DIR, "Section5_MosaicPlot_Fe_Aqueous_only.png")
fig5_aq.savefig(out5_aq, dpi=150)
plt.close(fig5_aq)
print(f"  Saved: {out5_aq}\n")


# ===========================================================================
# Section 6 — LogfO2pHDiagram: Fe oxide stability vs log fO2 and pH
# ===========================================================================
print("Section 6: LogfO2pHDiagram — Fe stability vs log fO2 and pH")

fe_aq_fo2 = [
    "H2O(aq)",
    "H+",
    "OH-",
    "e-",
    "Fe+2",
    "Fe+3",
    "FeO(aq)",
    "FeO+",
    "FeO2-",
    "FeOH+",
    "FeOH+2",
    "HFeO2(aq)",
    "HFeO2-",
]
fe_min_fo2 = ["Hematite", "Magnetite", "Goethite", "Iron"]
fe_min_fo2_with_aqueous_only = fe_min_fo2 + ["Aqueous_only"]
aq_species_fo2 = fe_aq_fo2[4:]

system_fo2 = ChemicalSystem(
    db,
    AqueousPhase(speciate(fe_aq_fo2)),
    GaseousPhase("O2(g)"),
    MineralPhases(StringList(fe_min_fo2)),
)

specs_fo2 = EquilibriumSpecs(system_fo2)
specs_fo2.temperature()
specs_fo2.pressure()
specs_fo2.fugacity("O2")
specs_fo2.pH()

conds_fo2 = EquilibriumConditions(specs_fo2)
conds_fo2.set("T", 298.15)
conds_fo2.set("P", 1.0e5)
solver_fo2 = EquilibriumSolver(specs_fo2)

seed_amounts = [
    ("H2O(aq)", 1.0, "kg"),
    ("Fe+2", 1e-3, "mol"),
    ("O2(g)", 1e-6, "mol"),
]

logfO2_vals = np.linspace(-80.0, 0.0, max(30, env_int("DIAG_SEC6_NFO2", 70)))
pH_fo2_vals = np.linspace(0.0, 14.0, max(30, env_int("DIAG_SEC6_NPH", 70)))
pred_fo2 = np.full((len(logfO2_vals), len(pH_fo2_vals)), -1, dtype=int)
pred_fo2_mineral_only = np.full((len(logfO2_vals), len(pH_fo2_vals)), -1, dtype=int)
pred_fo2_aq_only = np.full((len(logfO2_vals), len(pH_fo2_vals)), -1, dtype=int)
sec6_mineral_stability_tol = 1e-14
for i, lfo2 in enumerate(logfO2_vals):
    for j, pH in enumerate(pH_fo2_vals):
        state = _fresh_state(system_fo2, seed_amounts)
        conds_fo2.set("T", 298.15)
        conds_fo2.set("P", 1.0e5)
        conds_fo2.fugacity("O2", 10.0 ** float(lfo2), "bar")
        conds_fo2.pH(float(pH))
        try:
            result = solver_fo2.solve(state, conds_fo2)
            if result.succeeded():
                pred_fo2[i, j] = _predominant_fe_field(
                    state, aq_species_fo2, fe_min_fo2
                )
                pred_fo2_aq_only[i, j] = _predominant_by_amount(state, aq_species_fo2)
                mineral_amounts_fo2 = [state.speciesAmount(s).val() for s in fe_min_fo2]
                if max(mineral_amounts_fo2) <= sec6_mineral_stability_tol:
                    pred_fo2_mineral_only[i, j] = len(fe_min_fo2)
                else:
                    pred_fo2_mineral_only[i, j] = int(np.argmax(mineral_amounts_fo2))
        except Exception:
            pass

fo2_proxy = GridResultProxy(logfO2_vals, pH_fo2_vals, pred_fo2)
fd = LogfO2pHDiagram.from_grid_result(
    fo2_proxy, aq_species_fo2 + fe_min_fo2, palette="Set3"
)

fig6, ax6 = fd.plot(figsize=(8, 5))
for label, lfo2 in [("HM", -70.6), ("FMQ", -85.0)]:
    if logfO2_vals[0] <= lfo2 <= logfO2_vals[-1]:
        ax6.axvline(lfo2, color="gray", linestyle=":", linewidth=1.0)
        ax6.text(lfo2 + 0.5, 13.0, label, fontsize=8, color="gray")

ax6.set_title(r"Fe–O–H stability: $\log f_{O_2}$ vs pH at 25 °C, 1 bar", fontsize=11)
plt.tight_layout()
out6 = os.path.join(THIS_DIR, "Section6_LogfO2pHDiagram_Fe.png")
fig6.savefig(out6, dpi=150)
plt.close(fig6)
print(f"  Saved: {out6}\n")

fo2_proxy_min = GridResultProxy(logfO2_vals, pH_fo2_vals, pred_fo2_mineral_only)
fd_min = LogfO2pHDiagram.from_grid_result(
    fo2_proxy_min, fe_min_fo2_with_aqueous_only, palette="YlOrBr"
)
fig6_min, ax6_min = fd_min.plot(figsize=(8, 5))
for label, lfo2 in [("HM", -70.6), ("FMQ", -85.0)]:
    if logfO2_vals[0] <= lfo2 <= logfO2_vals[-1]:
        ax6_min.axvline(lfo2, color="gray", linestyle=":", linewidth=1.0)
        ax6_min.text(lfo2 + 0.5, 13.0, label, fontsize=8, color="gray")
ax6_min.set_title(
    r"Fe–O–H mineral-only: $\log f_{O_2}$ vs pH at 25 °C, 1 bar", fontsize=11
)
plt.tight_layout()
out6_min = os.path.join(THIS_DIR, "Section6_LogfO2pHDiagram_Fe_Mineral_only.png")
fig6_min.savefig(out6_min, dpi=150)
plt.close(fig6_min)
print(f"  Saved: {out6_min}\n")

fo2_proxy_aq = GridResultProxy(logfO2_vals, pH_fo2_vals, pred_fo2_aq_only)
fd_aq = LogfO2pHDiagram.from_grid_result(fo2_proxy_aq, aq_species_fo2, palette="tab20")
fig6_aq, ax6_aq = fd_aq.plot(figsize=(8, 5))
for label, lfo2 in [("HM", -70.6), ("FMQ", -85.0)]:
    if logfO2_vals[0] <= lfo2 <= logfO2_vals[-1]:
        ax6_aq.axvline(lfo2, color="gray", linestyle=":", linewidth=1.0)
        ax6_aq.text(lfo2 + 0.5, 13.0, label, fontsize=8, color="gray")
ax6_aq.set_title(
    r"Fe–O–H aqueous-only: $\log f_{O_2}$ vs pH at 25 °C, 1 bar", fontsize=11
)
plt.tight_layout()
out6_aq = os.path.join(THIS_DIR, "Section6_LogfO2pHDiagram_Fe_Aqueous_only.png")
fig6_aq.savefig(out6_aq, dpi=150)
plt.close(fig6_aq)
print(f"  Saved: {out6_aq}\n")


# ===========================================================================
# Section 7 — TPDiagram: Calcite / Aragonite polymorphism
# ===========================================================================
print("Section 7: TPDiagram — Calcite vs Aragonite polymorphism")

system_poly = ChemicalSystem(
    db,
    AqueousPhase("H2O(aq) H+ OH- Ca+2 HCO3- CO3-2 CO2(aq)"),
    MineralPhases(StringList(["Calcite", "Aragonite"])),
)

specs_poly = EquilibriumSpecs(system_poly)
specs_poly.temperature()
specs_poly.pressure()

conds_poly = EquilibriumConditions(specs_poly)
solver_poly = EquilibriumSolver(specs_poly)

seed_amounts = [
    ("H2O(aq)", 1.0, "kg"),
    ("Calcite", 5.0, "mol"),
    ("Aragonite", 5.0, "mol"),
]

T_poly = np.linspace(25.0, 500.0, max(30, env_int("DIAG_SEC7_NT", 70)))
P_poly = np.linspace(1.0, 5000.0, max(30, env_int("DIAG_SEC7_NP", 70)))
pred_poly = np.full((len(T_poly), len(P_poly)), -1, dtype=int)
for i, T_C in enumerate(T_poly):
    for j, P_bar in enumerate(P_poly):
        state = _fresh_state(system_poly, seed_amounts)
        conds_poly.set("T", float(T_C + 273.15))
        conds_poly.set("P", float(P_bar * 1e5))
        try:
            result = solver_poly.solve(state, conds_poly)
            if result.succeeded():
                pred_poly[i, j] = _predominant_by_amount(
                    state, ["Calcite", "Aragonite"]
                )
        except Exception:
            pass

td_proxy = GridResultProxy(T_poly, P_poly, pred_poly)
td = TPDiagram.from_grid_result(
    td_proxy,
    ["Calcite", "Aragonite"],
    T_unit="C",
    P_unit="bar",
    palette="Set1",
)

fig7, ax7 = td.plot(figsize=(7, 5))
ax7.set_title("CaCO₃ polymorph stability (Calcite vs Aragonite)", fontsize=11)
ax7.set_xlabel("T (°C)")
ax7.set_ylabel("P (bar)")
plt.tight_layout()
out7 = os.path.join(THIS_DIR, "Section7_TPDiagram_CaCO3_polymorphs.png")
fig7.savefig(out7, dpi=150)
plt.close(fig7)
print(f"  Saved: {out7}\n")

print("All diagrams complete.")
