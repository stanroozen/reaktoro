"""
Willemite-focused stability diagram suite for the DEW17HP622_Zn tutorial system.

These diagrams are intended to answer where Willemite is stable more directly
than the H+/HS- activity map. Because the local binding path does not reliably
support direct SiO2,aq activity/chemical-potential constraints from this script,
it uses the initial SiO2,aq inventory as a practical silica-control proxy.

Generated figures:
- pH vs silica-inventory proxy dominant-mineral map + Willemite-only solubility
- Temperature vs pH dominant-mineral map at fixed silica proxy
- pH vs log10(fH2S) dominant-mineral map at fixed silica proxy
"""

import importlib.util
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.lines import Line2D


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(os.path.dirname(SCRIPT_DIR))
TUTORIAL_PATH = os.path.join(
    SCRIPT_DIR,
    "willemite_solubility_tutorial_dew17hp622_zn.py",
)
LOCAL_REAKTORO_RELEASE_DIR = os.path.join(REPO_ROOT, "build", "Reaktoro", "Release")
LOCAL_REAKTORO_PYD_DIR = os.path.join(
    REPO_ROOT,
    "build",
    "python",
    "package",
    "build",
    "lib",
    "reaktoro",
)

OUTPUT_PH_SIO2 = os.path.join(
    SCRIPT_DIR,
    "willemite_stability_ph_sio2_suite.png",
)
OUTPUT_T_PH = os.path.join(
    SCRIPT_DIR,
    "willemite_stability_temperature_ph_fixed_sio2.png",
)
OUTPUT_PH_FH2S = os.path.join(
    SCRIPT_DIR,
    "willemite_stability_ph_fh2s_fixed_sio2.png",
)

TEMPERATURE_C = 300.0
PRESSURE_KBAR = 2.0
FIXED_LOG_SIO2_PROXY_MOL = -3.0

PH_MIN = 0.0
PH_MAX = 14.0
LOG_SIO2_PROXY_MOL_MIN = -8.0
LOG_SIO2_PROXY_MOL_MAX = -1.0
LOG_FH2S_MAX = -6.0
LOG_FH2S_MIN = -20.0
TEMP_MIN_C = 100.0
TEMP_MAX_C = 400.0

GRID_NX = 32
GRID_NY = 32

# Solver option controls exposed by this local reaktoro4py binding.
# These do not change convergence in this system unless ideal activity models
# are forced on, which improves robustness but changes thermodynamic physics.
SOLVER_EPSILON = None  # e.g., 1.0e-13
SOLVER_LOG_BARRIER_FACTOR = None  # e.g., 1.0e-2
SOLVER_USE_IDEAL_ACTIVITY_MODELS = False
SOLVER_HESSIAN_MODE = "Exact"  # One of: Exact, PartiallyExact, Approx, ApproxDiagonal

# pH domain expanded to cover full Earth-relevant aqueous conditions.
# Retain the marker for reporting/comparison purposes only.
pH_CONVERGENCE_LIMIT = PH_MAX
AQUEOUS_ONLY_SOLID_THRESHOLD_MOL = 1.0e-14
WILLEMITE_NAME = "Wlm"


def prepare_state_with_silica_proxy(module, system, log_sio2_proxy_mol, seed_state=None):
    if seed_state is None:
        return make_state_with_silica_proxy(module, system, log_sio2_proxy_mol)

    silica_moles = float(10.0**log_sio2_proxy_mol)
    try:
        seed_state.set(module.SIO2_PROXY_SPECIES, silica_moles, "mol")
    except Exception:
        return make_state_with_silica_proxy(module, system, log_sio2_proxy_mol)
    return seed_state


def configure_solver_options(module, solver):
    try:
        options = module.EquilibriumOptions()
    except Exception:
        return

    if SOLVER_EPSILON is not None and hasattr(options, "epsilon"):
        options.epsilon = float(SOLVER_EPSILON)
    if (
        SOLVER_LOG_BARRIER_FACTOR is not None
        and hasattr(options, "logarithm_barrier_factor")
    ):
        options.logarithm_barrier_factor = float(SOLVER_LOG_BARRIER_FACTOR)
    if hasattr(options, "use_ideal_activity_models"):
        options.use_ideal_activity_models = bool(SOLVER_USE_IDEAL_ACTIVITY_MODELS)

    # Configure Hessian mode when exposed by the local Python binding.
    if hasattr(options, "hessian") and hasattr(module, "GibbsHessian"):
        hessian_mode = str(SOLVER_HESSIAN_MODE).strip()
        if hasattr(module.GibbsHessian, hessian_mode):
            options.hessian = getattr(module.GibbsHessian, hessian_mode)

    try:
        solver.setOptions(options)
    except Exception:
        pass


def load_tutorial_module(path):
    if sys.platform.startswith("win") and os.path.isdir(LOCAL_REAKTORO_RELEASE_DIR):
        if LOCAL_REAKTORO_RELEASE_DIR not in sys.path:
            sys.path.insert(0, LOCAL_REAKTORO_RELEASE_DIR)
    if os.path.isdir(LOCAL_REAKTORO_PYD_DIR) and LOCAL_REAKTORO_PYD_DIR not in sys.path:
        sys.path.append(LOCAL_REAKTORO_PYD_DIR)

    spec = importlib.util.spec_from_file_location("willemite_tutorial", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def plotted_zn_mineral_names(module):
    # Keep calculations free to include gangue minerals, but plots should report
    # only Zn-bearing phase stability fields.
    if getattr(module, "USE_COMPETING_ZN_MINERALS", False):
        return list(dict.fromkeys(getattr(module, "COMPETING_ZN_MINERALS", [])))
    return [module.MINERAL_NAME]


def dominant_mineral_label(module, state, plotted_minerals):
    dominant_name = "AqueousOnly"
    dominant_amount = 0.0

    for mineral_name in plotted_minerals:
        try:
            amount = float(state.speciesAmount(mineral_name))
        except Exception:
            amount = 0.0

        if amount > dominant_amount:
            dominant_amount = amount
            dominant_name = mineral_name

    if dominant_amount < AQUEOUS_ONLY_SOLID_THRESHOLD_MOL:
        return "AqueousOnly"

    return dominant_name


def infer_solvent_species_name(module):
    return module.infer_solvent_species_name(
        module.AQUEOUS_SPECIES,
        module.SOLVENT_SPECIES_NAME,
        module.INITIAL_SPECIES_AMOUNTS_MOL,
    )


def make_state_with_silica_proxy(module, system, log_sio2_proxy_mol):
    state = module.make_base_state(system)
    silica_moles = float(10.0**log_sio2_proxy_mol)
    try:
        state.set(module.SIO2_PROXY_SPECIES, silica_moles, "mol")
    except Exception:
        overrides = dict(module.INITIAL_SPECIES_AMOUNTS_MOL)
        overrides[module.SIO2_PROXY_SPECIES] = silica_moles
        module.apply_species_amount_overrides(state, overrides)
    return state


def make_phase_style(categories):
    cmap = ListedColormap(plt.cm.tab20(np.linspace(0.0, 1.0, len(categories))))
    norm = BoundaryNorm(np.arange(len(categories) + 1) - 0.5, cmap.N)
    return cmap, norm


def add_phase_legend(ax, categories, present_indices, cmap, norm):
    handles = [
        Line2D(
            [0],
            [0],
            marker="s",
            linestyle="",
            markersize=8,
            markerfacecolor=cmap(norm(idx)),
            markeredgecolor="black",
            label=categories[idx],
        )
        for idx in present_indices
    ]
    ax.legend(
        handles=handles,
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        frameon=False,
    )


def add_willemite_boundary(ax, x_values, y_values, label_grid, willemite_index):
    x_mesh, y_mesh = np.meshgrid(x_values, y_values)
    mask = (label_grid == willemite_index).astype(float)
    if np.any(mask > 0.0):
        ax.contour(x_mesh, y_mesh, mask, levels=[0.5], colors="black", linewidths=1.8)
        return True
    return False


def save_ph_sio2_suite(module, categories):
    system = module.build_tutorial_system()
    specs = module.EquilibriumSpecs(system)
    specs.temperature()
    specs.pressure()
    specs.lgActivity("H+")

    solver = module.make_equilibrium_solver(system, specs)
    configure_solver_options(module, solver)
    conditions = module.EquilibriumConditions(specs)
    pressure_bar = PRESSURE_KBAR * 1000.0

    solvent_species_name = infer_solvent_species_name(module)
    plotted_minerals = plotted_zn_mineral_names(module)
    category_to_index = {name: idx for idx, name in enumerate(categories)}
    pH_values = np.linspace(PH_MIN, PH_MAX, GRID_NX)
    log_sio2_proxy_values = np.linspace(
        LOG_SIO2_PROXY_MOL_MIN, LOG_SIO2_PROXY_MOL_MAX, GRID_NY
    )

    label_grid = np.full(
        (len(log_sio2_proxy_values), len(pH_values)),
        category_to_index["NoConvergence"],
    )
    zn_grid = np.full((len(log_sio2_proxy_values), len(pH_values)), np.nan, dtype=float)
    stats = {"converged": 0, "failed": 0, "warm_seed_used": 0}
    column_warm_states = [None] * len(pH_values)

    for j, log_sio2_proxy_mol in enumerate(log_sio2_proxy_values):
        print(f"  pH-Si row {j + 1}/{len(log_sio2_proxy_values)}")
        row_warm_state = None
        x_indices = range(len(pH_values)) if (j % 2 == 0) else range(len(pH_values) - 1, -1, -1)
        for i in x_indices:
            pH_value = pH_values[i]
            seed_state = row_warm_state if row_warm_state is not None else column_warm_states[i]
            if seed_state is not None:
                stats["warm_seed_used"] += 1
            state = prepare_state_with_silica_proxy(
                module, system, log_sio2_proxy_mol, seed_state=seed_state
            )
            conditions.temperature(float(TEMPERATURE_C), "celsius")
            conditions.pressure(float(pressure_bar), "bar")
            conditions.lgActivity("H+", float(-pH_value))
            result = solver.solve(state, conditions)
            if not result.succeeded():
                state = make_state_with_silica_proxy(module, system, log_sio2_proxy_mol)
                result = solver.solve(state, conditions)

            if result.succeeded():
                stats["converged"] += 1
                row_warm_state = state
                column_warm_states[i] = state
                label_grid[j, i] = category_to_index[
                    dominant_mineral_label(module, state, plotted_minerals)
                ]
                zn_grid[j, i] = module.dissolved_element_molality(
                    state, solvent_species_name
                )
            else:
                stats["failed"] += 1
                row_warm_state = None

    cmap, norm = make_phase_style(categories)
    x_mesh, y_mesh = np.meshgrid(pH_values, log_sio2_proxy_values)
    willemite_idx = category_to_index[WILLEMITE_NAME]
    present_indices = sorted(int(idx) for idx in np.unique(label_grid))

    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(14, 6.2))

    axes[0].imshow(
        label_grid,
        origin="lower",
        extent=(PH_MIN, PH_MAX, LOG_SIO2_PROXY_MOL_MIN, LOG_SIO2_PROXY_MOL_MAX),
        aspect="auto",
        cmap=cmap,
        norm=norm,
    )
    add_willemite_boundary(
        axes[0], pH_values, log_sio2_proxy_values, label_grid, willemite_idx
    )
    axes[0].set_xlabel("pH")
    axes[0].set_ylabel("log10(initial SiO2,aq mol)")
    axes[0].set_title("Dominant Zn Mineral")
    axes[0].grid(True, alpha=0.2, linestyle="--")
    add_phase_legend(axes[0], categories, present_indices, cmap, norm)

    valid = np.isfinite(zn_grid) & (zn_grid > 0.0) & (label_grid == willemite_idx)
    if np.any(valid):
        logz = np.full_like(zn_grid, np.nan, dtype=float)
        logz[valid] = np.log10(zn_grid[valid])
        levels = np.linspace(float(np.nanmin(logz)), float(np.nanmax(logz)), 12)
        cf = axes[1].contourf(x_mesh, y_mesh, logz, levels=levels, cmap="viridis")
        axes[1].contour(
            x_mesh,
            y_mesh,
            logz,
            levels=levels[::2],
            colors="white",
            linewidths=0.7,
        )
        add_willemite_boundary(
            axes[1], pH_values, log_sio2_proxy_values, label_grid, willemite_idx
        )
        cbar = fig.colorbar(cf, ax=axes[1])
        cbar.set_label("log10(total dissolved Zn molality)")
    else:
        axes[1].text(
            0.5,
            0.5,
            "Willemite not stable\nin sampled window",
            ha="center",
            va="center",
            transform=axes[1].transAxes,
        )

    axes[1].set_xlabel("pH")
    axes[1].set_ylabel("log10(initial SiO2,aq mol)")
    axes[1].set_title("Willemite-Stable Zn Solubility")
    axes[1].grid(True, alpha=0.2, linestyle="--")

    fig.suptitle(
        "Willemite Stability Search: pH vs Silica Inventory Proxy\n"
        f"T={TEMPERATURE_C:.0f} C, P={PRESSURE_KBAR:.1f} kbar",
        fontsize=12,
    )
    fig.tight_layout()
    fig.savefig(OUTPUT_PH_SIO2, dpi=250)
    plt.close(fig)

    return stats


def save_temperature_ph_map(module, categories):
    system = module.build_tutorial_system()
    specs = module.EquilibriumSpecs(system)
    specs.temperature()
    specs.pressure()
    specs.lgActivity("H+")

    solver = module.make_equilibrium_solver(system, specs)
    configure_solver_options(module, solver)
    conditions = module.EquilibriumConditions(specs)
    pressure_bar = PRESSURE_KBAR * 1000.0
    plotted_minerals = plotted_zn_mineral_names(module)

    category_to_index = {name: idx for idx, name in enumerate(categories)}
    pH_values = np.linspace(PH_MIN, PH_MAX, GRID_NX)
    temperature_values = np.linspace(TEMP_MIN_C, TEMP_MAX_C, GRID_NY)
    label_grid = np.full(
        (len(temperature_values), len(pH_values)),
        category_to_index["NoConvergence"],
    )
    stats = {"converged": 0, "failed": 0, "warm_seed_used": 0}
    column_warm_states = [None] * len(pH_values)

    for j, temperature_c in enumerate(temperature_values):
        print(f"  T-pH row {j + 1}/{len(temperature_values)}")
        row_warm_state = None
        x_indices = range(len(pH_values)) if (j % 2 == 0) else range(len(pH_values) - 1, -1, -1)
        for i in x_indices:
            pH_value = pH_values[i]
            seed_state = row_warm_state if row_warm_state is not None else column_warm_states[i]
            if seed_state is not None:
                stats["warm_seed_used"] += 1
            state = prepare_state_with_silica_proxy(
                module,
                system,
                FIXED_LOG_SIO2_PROXY_MOL,
                seed_state=seed_state,
            )
            conditions.temperature(float(temperature_c), "celsius")
            conditions.pressure(float(pressure_bar), "bar")
            conditions.lgActivity("H+", float(-pH_value))
            result = solver.solve(state, conditions)
            if not result.succeeded():
                state = make_state_with_silica_proxy(
                    module, system, FIXED_LOG_SIO2_PROXY_MOL
                )
                result = solver.solve(state, conditions)
            if result.succeeded():
                stats["converged"] += 1
                row_warm_state = state
                column_warm_states[i] = state
                label_grid[j, i] = category_to_index[
                    dominant_mineral_label(module, state, plotted_minerals)
                ]
            else:
                stats["failed"] += 1
                row_warm_state = None

    cmap, norm = make_phase_style(categories)
    present_indices = sorted(int(idx) for idx in np.unique(label_grid))
    willemite_idx = category_to_index[WILLEMITE_NAME]

    fig, ax = plt.subplots(figsize=(10, 6.2))
    ax.imshow(
        label_grid,
        origin="lower",
        extent=(PH_MIN, PH_MAX, TEMP_MIN_C, TEMP_MAX_C),
        aspect="auto",
        cmap=cmap,
        norm=norm,
    )
    add_willemite_boundary(ax, pH_values, temperature_values, label_grid, willemite_idx)
    ax.set_xlabel("pH")
    ax.set_ylabel("Temperature (C)")
    ax.set_title(
        "Dominant Zn Mineral: Temperature vs pH\n"
        f"P={PRESSURE_KBAR:.1f} kbar, log10(initial SiO2,aq mol)={FIXED_LOG_SIO2_PROXY_MOL:.2f}"
    )
    ax.grid(True, alpha=0.2, linestyle="--")
    add_phase_legend(ax, categories, present_indices, cmap, norm)

    fig.tight_layout()
    fig.savefig(OUTPUT_T_PH, dpi=250)
    plt.close(fig)

    return stats


def save_ph_fh2s_map(module, categories):
    system = module.build_system_with_h2s_gas_phase()
    specs = module.EquilibriumSpecs(system)
    specs.temperature()
    specs.pressure()
    specs.lgActivity("H+")
    specs.fugacity(module.SULFUR_FUGACITY_SPECIES)

    solver = module.make_equilibrium_solver(system, specs)
    configure_solver_options(module, solver)
    conditions = module.EquilibriumConditions(specs)
    pressure_bar = PRESSURE_KBAR * 1000.0
    plotted_minerals = plotted_zn_mineral_names(module)

    category_to_index = {name: idx for idx, name in enumerate(categories)}
    pH_values = np.linspace(PH_MIN, PH_MAX, GRID_NX)
    logf_values = np.linspace(LOG_FH2S_MIN, LOG_FH2S_MAX, GRID_NY)
    label_grid = np.full(
        (len(logf_values), len(pH_values)),
        category_to_index["NoConvergence"],
    )
    stats = {"converged": 0, "failed": 0, "warm_seed_used": 0}
    column_warm_states = [None] * len(pH_values)

    for j, logf in enumerate(logf_values):
        print(f"  pH-fH2S row {j + 1}/{len(logf_values)}")
        row_warm_state = None
        x_indices = range(len(pH_values)) if (j % 2 == 0) else range(len(pH_values) - 1, -1, -1)
        for i in x_indices:
            pH_value = pH_values[i]
            seed_state = row_warm_state if row_warm_state is not None else column_warm_states[i]
            if seed_state is not None:
                stats["warm_seed_used"] += 1
            state = prepare_state_with_silica_proxy(
                module,
                system,
                FIXED_LOG_SIO2_PROXY_MOL,
                seed_state=seed_state,
            )
            try:
                state.set(module.SULFUR_FUGACITY_SPECIES, 1.0e-20, "mol")
            except Exception:
                pass
            conditions.temperature(float(TEMPERATURE_C), "celsius")
            conditions.pressure(float(pressure_bar), "bar")
            conditions.lgActivity("H+", float(-pH_value))
            conditions.fugacity(
                module.SULFUR_FUGACITY_SPECIES, float(10.0**logf), "bar"
            )
            result = solver.solve(state, conditions)
            if not result.succeeded():
                state = make_state_with_silica_proxy(
                    module, system, FIXED_LOG_SIO2_PROXY_MOL
                )
                try:
                    state.set(module.SULFUR_FUGACITY_SPECIES, 1.0e-20, "mol")
                except Exception:
                    pass
                result = solver.solve(state, conditions)
            if result.succeeded():
                stats["converged"] += 1
                row_warm_state = state
                column_warm_states[i] = state
                label_grid[j, i] = category_to_index[
                    dominant_mineral_label(module, state, plotted_minerals)
                ]
            else:
                stats["failed"] += 1
                row_warm_state = None

    cmap, norm = make_phase_style(categories)
    present_indices = sorted(int(idx) for idx in np.unique(label_grid))
    willemite_idx = category_to_index[WILLEMITE_NAME]

    fig, ax = plt.subplots(figsize=(10, 6.2))
    ax.imshow(
        label_grid,
        origin="lower",
        extent=(PH_MIN, PH_MAX, LOG_FH2S_MIN, LOG_FH2S_MAX),
        aspect="auto",
        cmap=cmap,
        norm=norm,
    )
    add_willemite_boundary(ax, pH_values, logf_values, label_grid, willemite_idx)
    ax.set_xlabel("pH")
    ax.set_ylabel("log10(fH2S / bar)")
    ax.set_title(
        "Dominant Zn Mineral: pH vs Sulfur Fugacity\n"
        f"T={TEMPERATURE_C:.0f} C, P={PRESSURE_KBAR:.1f} kbar, log10(initial SiO2,aq mol)={FIXED_LOG_SIO2_PROXY_MOL:.2f}"
    )
    ax.grid(True, alpha=0.2, linestyle="--")
    add_phase_legend(ax, categories, present_indices, cmap, norm)

    fig.tight_layout()
    fig.savefig(OUTPUT_PH_FH2S, dpi=250)
    plt.close(fig)

    return stats


def main():
    module = load_tutorial_module(TUTORIAL_PATH)

    try:
        module.Warnings.disable(906)
    except Exception:
        pass

    module.validate_user_inputs()
    module.USE_COMPETING_ZN_MINERALS = True

    categories = plotted_zn_mineral_names(module) + [
        "AqueousOnly",
        "NoConvergence",
    ]

    print("Computing pH vs silica-proxy map...")
    ph_sio2_stats = save_ph_sio2_suite(module, categories)
    print("Computing temperature vs pH map...")
    t_ph_stats = save_temperature_ph_map(module, categories)
    print("Computing pH vs fH2S map...")
    ph_fh2s_stats = save_ph_fh2s_map(module, categories)

    print("Generated Willemite-focused stability diagrams:")
    print(OUTPUT_PH_SIO2)
    print(OUTPUT_T_PH)
    print(OUTPUT_PH_FH2S)
    print(
        "pH-silica-proxy grid convergence: "
        f"{ph_sio2_stats['converged']}/{ph_sio2_stats['converged'] + ph_sio2_stats['failed']}"
    )
    print(f"pH-silica-proxy warm starts used: {ph_sio2_stats['warm_seed_used']}")
    print(
        "T-pH grid convergence: "
        f"{t_ph_stats['converged']}/{t_ph_stats['converged'] + t_ph_stats['failed']}"
    )
    print(f"T-pH warm starts used: {t_ph_stats['warm_seed_used']}")
    print(
        "pH-fH2S grid convergence: "
        f"{ph_fh2s_stats['converged']}/{ph_fh2s_stats['converged'] + ph_fh2s_stats['failed']}"
    )
    print(f"pH-fH2S warm starts used: {ph_fh2s_stats['warm_seed_used']}")
    print(
        "Solver options: "
        f"epsilon={SOLVER_EPSILON}, "
        f"logarithm_barrier_factor={SOLVER_LOG_BARRIER_FACTOR}, "
        f"use_ideal_activity_models={SOLVER_USE_IDEAL_ACTIVITY_MODELS}"
    )


if __name__ == "__main__":
    main()
