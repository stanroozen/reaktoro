"""
Activity-Activity Phase Diagram for the DEW17HP622_Zn Willemite system.

This script performs a 2D sweep over two imposed aqueous activities:
- log10(a(H+))
- log10(a(HS-))

At each grid point, equilibrium is solved and the dominant Zn-bearing mineral
is classified. The result is a stable-phase map for an open-system style setup
with externally controlled activities.
"""

import importlib.util
import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap
from matplotlib.patches import Rectangle


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
TUTORIAL_PATH = os.path.join(
    SCRIPT_DIR, "willemite_solubility_tutorial_dew17hp622_zn.py"
)

OUTPUT_FIGURE = os.path.join(
    SCRIPT_DIR,
    "willemite_activity_activity_phase_diagram_hplus_hs.png",
)
OUTPUT_SOLUBILITY_CONTOUR_FIGURE = os.path.join(
    SCRIPT_DIR,
    "willemite_activity_activity_solubility_contours_hplus_hs.png",
)
OUTPUT_OVERLAY_FIGURE = os.path.join(
    SCRIPT_DIR,
    "willemite_activity_activity_phase_solubility_overlay_hplus_hs.png",
)
OUTPUT_CROSS_SECTION_FIGURE = os.path.join(
    SCRIPT_DIR,
    "willemite_activity_activity_solubility_sections_hs.png",
)
OUTPUT_REDOX_SULFUR_FIGURE = os.path.join(
    SCRIPT_DIR,
    "willemite_solubility_redox_sulfur_sweeps.png",
)

# Sweep controls
ACTIVITY_SPECIES_X = "H+"
ACTIVITY_SPECIES_Y = "HS-"
LOG_A_X_MIN = -10.0
LOG_A_X_MAX = -3.0
LOG_A_Y_MIN = -12.0
LOG_A_Y_MAX = -2.0

# Adaptive mesh settings
COARSE_NX = 6
COARSE_NY = 6
MAX_REFINEMENT_LEVEL = 1

# Classification settings
AQUEOUS_ONLY_SOLID_THRESHOLD_MOL = 1.0e-14
NONCONV_INFERENCE_X_SHIFTS = (0.02, 0.05, 0.10, 0.20)
NONCONV_INFERENCE_Y_SHIFTS = (0.40,)

# Optional strict-domain clipping to avoid known hard-to-converge low-log a(H+) region.
USE_STRICT_DOMAIN_CLIP = True
STRICT_MIN_LOG_A_X = -6.5

# Reference state
TEMPERATURE_C = 300.0
PRESSURE_KBAR = 2.0

# Regular grid used for solubility contours and section plots.
GRID_NX = 15
GRID_NY = 15
SECTION_LOG_A_HS_VALUES = (-11.0, -8.0, -5.5)


def load_tutorial_module(path):
    spec = importlib.util.spec_from_file_location("willemite_tutorial", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def dominant_mineral_label(module, state, threshold):
    dominant_name = "AqueousOnly"
    dominant_amount = 0.0

    for mineral_name in module.selected_mineral_names():
        try:
            amount = float(state.speciesAmount(mineral_name))
        except Exception:
            amount = 0.0

        if amount > dominant_amount:
            dominant_amount = amount
            dominant_name = mineral_name

    if dominant_amount < threshold:
        return "AqueousOnly"

    return dominant_name


def make_label_evaluator(module, system, solver, conditions, category_to_index):
    cache = {}
    inference_cache = {}
    converged_points = []
    stats = {
        "domain_clipped_points": 0,
        "domain_clipped_inferred_points": 0,
        "domain_clipped_unresolved_points": 0,
        "inferred_from_shift_points": 0,
        "true_noconvergence_points": 0,
    }

    def _nearest_converged_key(loga_x, loga_y):
        if not converged_points:
            return None
        best_key = None
        best_d2 = None
        for px, py in converged_points:
            d2 = (px - loga_x) ** 2 + (py - loga_y) ** 2
            if best_d2 is None or d2 < best_d2:
                best_d2 = d2
                best_key = (px, py)
        return best_key

    def _solve_at(state, loga_x, loga_y):
        conditions.lgActivity(ACTIVITY_SPECIES_X, loga_x)
        conditions.lgActivity(ACTIVITY_SPECIES_Y, loga_y)
        return solver.solve(state, conditions)

    def _infer_label_from_shift(loga_x, loga_y, use_clip_anchor=False):
        if use_clip_anchor:
            cache_key = ("clip", round(float(loga_y), 6))
        else:
            cache_key = (round(float(loga_x), 6), round(float(loga_y), 6))
        if cache_key in inference_cache:
            return inference_cache[cache_key]

        x_candidates = []
        if use_clip_anchor:
            for dx in NONCONV_INFERENCE_X_SHIFTS:
                x_candidates.append(float(STRICT_MIN_LOG_A_X + dx))
        else:
            for dx in NONCONV_INFERENCE_X_SHIFTS:
                x_candidates.append(float(loga_x + dx))

        for shifted_x in x_candidates:
            if shifted_x > LOG_A_X_MAX:
                continue
            shifted_state = module.make_base_state(system)
            shifted_result = _solve_at(shifted_state, shifted_x, loga_y)
            if shifted_result.succeeded():
                inferred = dominant_mineral_label(
                    module,
                    shifted_state,
                    threshold=AQUEOUS_ONLY_SOLID_THRESHOLD_MOL,
                )
                inference_cache[cache_key] = inferred
                return inferred

        for shifted_x in x_candidates:
            if shifted_x > LOG_A_X_MAX:
                continue
            for dy in NONCONV_INFERENCE_Y_SHIFTS:
                for sign in (-1.0, 1.0):
                    shifted_y = float(loga_y + sign * dy)
                    if shifted_y < LOG_A_Y_MIN or shifted_y > LOG_A_Y_MAX:
                        continue
                    shifted_state = module.make_base_state(system)
                    shifted_result = _solve_at(shifted_state, shifted_x, shifted_y)
                    if shifted_result.succeeded():
                        inferred = dominant_mineral_label(
                            module,
                            shifted_state,
                            threshold=AQUEOUS_ONLY_SOLID_THRESHOLD_MOL,
                        )
                        inference_cache[cache_key] = inferred
                        return inferred

        inference_cache[cache_key] = None
        return None

    def evaluate(loga_x, loga_y):
        key = (float(loga_x), float(loga_y))
        if key in cache:
            return cache[key]

        if USE_STRICT_DOMAIN_CLIP and key[0] < STRICT_MIN_LOG_A_X:
            stats["domain_clipped_points"] += 1
            inferred_label = _infer_label_from_shift(
                key[0], key[1], use_clip_anchor=True
            )
            if inferred_label is not None:
                idx = category_to_index.get(
                    inferred_label,
                    category_to_index["AqueousOnly"],
                )
                stats["domain_clipped_inferred_points"] += 1
            else:
                stats["domain_clipped_unresolved_points"] += 1
                stats["true_noconvergence_points"] += 1
                idx = category_to_index["NoConvergence"]
            cache[key] = idx
            return idx

        state = module.make_base_state(system)
        result = _solve_at(state, key[0], key[1])

        if not result.succeeded():
            nearest = _nearest_converged_key(key[0], key[1])
            if nearest is not None:
                warm_state = module.make_base_state(system)
                seed_result = _solve_at(warm_state, nearest[0], nearest[1])
                if seed_result.succeeded():
                    result = _solve_at(warm_state, key[0], key[1])
                    if result.succeeded():
                        state = warm_state

        if result.succeeded():
            converged_points.append(key)
            label = dominant_mineral_label(
                module,
                state,
                threshold=AQUEOUS_ONLY_SOLID_THRESHOLD_MOL,
            )
        else:
            inferred_label = _infer_label_from_shift(
                key[0], key[1], use_clip_anchor=False
            )
            if inferred_label is not None:
                label = inferred_label
                stats["inferred_from_shift_points"] += 1
            else:
                label = "NoConvergence"
                stats["true_noconvergence_points"] += 1

        idx = category_to_index.get(label, category_to_index["AqueousOnly"])
        cache[key] = idx
        return idx

    return evaluate, cache, stats


def refine_cell(x0, x1, y0, y1, depth, max_depth, evaluate, leaf_cells):
    xm = 0.5 * (x0 + x1)
    ym = 0.5 * (y0 + y1)

    sample_points = [
        (x0, y0),
        (x1, y0),
        (x0, y1),
        (x1, y1),
        (xm, ym),
    ]
    labels = [evaluate(x, y) for x, y in sample_points]
    unique_labels = set(labels)

    if depth >= max_depth or len(unique_labels) == 1:
        # Use center label for the final cell classification.
        leaf_cells.append((x0, x1, y0, y1, labels[-1], depth))
        return

    # Subdivide boundary cells into four quadrants.
    refine_cell(x0, xm, y0, ym, depth + 1, max_depth, evaluate, leaf_cells)
    refine_cell(xm, x1, y0, ym, depth + 1, max_depth, evaluate, leaf_cells)
    refine_cell(x0, xm, ym, y1, depth + 1, max_depth, evaluate, leaf_cells)
    refine_cell(xm, x1, ym, y1, depth + 1, max_depth, evaluate, leaf_cells)


def _plot_phase_rectangles(ax, leaf_cells, cmap):
    for x0, x1, y0, y1, label_idx, _ in leaf_cells:
        ax.add_patch(
            Rectangle(
                (x0, y0),
                x1 - x0,
                y1 - y0,
                facecolor=cmap(label_idx),
                edgecolor="none",
            )
        )


def _solve_regular_grid(module, categories, category_to_index):
    system = module.build_tutorial_system()

    specs = module.EquilibriumSpecs(system)
    specs.temperature()
    specs.pressure()
    specs.lgActivity(ACTIVITY_SPECIES_X)
    specs.lgActivity(ACTIVITY_SPECIES_Y)

    solver = module.make_equilibrium_solver(system, specs)
    conditions = module.EquilibriumConditions(specs)

    pressure_bar = PRESSURE_KBAR * 1000.0
    conditions.temperature(float(TEMPERATURE_C), "celsius")
    conditions.pressure(float(pressure_bar), "bar")

    solvent_species_name = module.infer_solvent_species_name(
        module.AQUEOUS_SPECIES,
        module.SOLVENT_SPECIES_NAME,
        module.INITIAL_SPECIES_AMOUNTS_MOL,
    )

    x_vals = np.linspace(LOG_A_X_MIN, LOG_A_X_MAX, GRID_NX)
    y_vals = np.linspace(LOG_A_Y_MIN, LOG_A_Y_MAX, GRID_NY)

    label_grid = np.full(
        (GRID_NY, GRID_NX), category_to_index["NoConvergence"], dtype=int
    )
    zn_grid = np.full((GRID_NY, GRID_NX), np.nan, dtype=float)

    converged_points = []
    point_cache = {}
    inference_cache = {}

    stats = {
        "grid_points": int(GRID_NX * GRID_NY),
        "grid_converged": 0,
        "grid_inferred": 0,
        "grid_true_noconvergence": 0,
        "grid_domain_clipped": 0,
    }

    def solve_at(state, loga_x, loga_y):
        conditions.lgActivity(ACTIVITY_SPECIES_X, float(loga_x))
        conditions.lgActivity(ACTIVITY_SPECIES_Y, float(loga_y))
        return solver.solve(state, conditions)

    def nearest_converged_key(loga_x, loga_y):
        if not converged_points:
            return None
        best_key = None
        best_d2 = None
        for px, py in converged_points:
            d2 = (px - loga_x) ** 2 + (py - loga_y) ** 2
            if best_d2 is None or d2 < best_d2:
                best_d2 = d2
                best_key = (px, py)
        return best_key

    def classify_state(state):
        label = dominant_mineral_label(
            module,
            state,
            threshold=AQUEOUS_ONLY_SOLID_THRESHOLD_MOL,
        )
        idx = category_to_index.get(label, category_to_index["AqueousOnly"])
        zn_molality = module.dissolved_element_molality(state, solvent_species_name)
        return idx, float(zn_molality)

    def infer_from_shift(loga_x, loga_y, use_clip_anchor=False):
        if use_clip_anchor:
            cache_key = ("clip", round(float(loga_y), 6))
        else:
            cache_key = (round(float(loga_x), 6), round(float(loga_y), 6))
        if cache_key in inference_cache:
            return inference_cache[cache_key]

        x_candidates = []
        if use_clip_anchor:
            for dx in NONCONV_INFERENCE_X_SHIFTS:
                x_candidates.append(float(STRICT_MIN_LOG_A_X + dx))
        else:
            for dx in NONCONV_INFERENCE_X_SHIFTS:
                x_candidates.append(float(loga_x + dx))

        for shifted_x in x_candidates:
            if shifted_x > LOG_A_X_MAX:
                continue
            shifted_state = module.make_base_state(system)
            if solve_at(shifted_state, shifted_x, loga_y).succeeded():
                inferred = classify_state(shifted_state)
                inference_cache[cache_key] = inferred
                return inferred

        for shifted_x in x_candidates:
            if shifted_x > LOG_A_X_MAX:
                continue
            for dy in NONCONV_INFERENCE_Y_SHIFTS:
                for sign in (-1.0, 1.0):
                    shifted_y = float(loga_y + sign * dy)
                    if shifted_y < LOG_A_Y_MIN or shifted_y > LOG_A_Y_MAX:
                        continue
                    shifted_state = module.make_base_state(system)
                    if solve_at(shifted_state, shifted_x, shifted_y).succeeded():
                        inferred = classify_state(shifted_state)
                        inference_cache[cache_key] = inferred
                        return inferred

        inference_cache[cache_key] = None
        return None

    for j, y in enumerate(y_vals):
        for i, x in enumerate(x_vals):
            key = (float(x), float(y))
            if key in point_cache:
                label_grid[j, i], zn_grid[j, i] = point_cache[key]
                continue

            if USE_STRICT_DOMAIN_CLIP and x < STRICT_MIN_LOG_A_X:
                stats["grid_domain_clipped"] += 1
                inferred = infer_from_shift(x, y, use_clip_anchor=True)
                if inferred is None:
                    point_cache[key] = (category_to_index["NoConvergence"], np.nan)
                    stats["grid_true_noconvergence"] += 1
                else:
                    point_cache[key] = inferred
                    stats["grid_inferred"] += 1
                label_grid[j, i], zn_grid[j, i] = point_cache[key]
                continue

            state = module.make_base_state(system)
            result = solve_at(state, x, y)

            if not result.succeeded():
                nearest = nearest_converged_key(x, y)
                if nearest is not None:
                    warm_state = module.make_base_state(system)
                    seed_result = solve_at(warm_state, nearest[0], nearest[1])
                    if seed_result.succeeded():
                        retry_result = solve_at(warm_state, x, y)
                        if retry_result.succeeded():
                            state = warm_state
                            result = retry_result

            if result.succeeded():
                converged_points.append(key)
                point_cache[key] = classify_state(state)
                stats["grid_converged"] += 1
            else:
                inferred = infer_from_shift(x, y, use_clip_anchor=False)
                if inferred is None:
                    point_cache[key] = (category_to_index["NoConvergence"], np.nan)
                    stats["grid_true_noconvergence"] += 1
                else:
                    point_cache[key] = inferred
                    stats["grid_inferred"] += 1

            label_grid[j, i], zn_grid[j, i] = point_cache[key]

    return x_vals, y_vals, label_grid, zn_grid, stats


def _save_solubility_contour_plot(x_vals, y_vals, zn_grid):
    valid = np.isfinite(zn_grid) & (zn_grid > 0.0)
    if not np.any(valid):
        return False

    x_mesh, y_mesh = np.meshgrid(x_vals, y_vals)
    logz = np.full_like(zn_grid, np.nan, dtype=float)
    logz[valid] = np.log10(zn_grid[valid])

    min_log = float(np.nanmin(logz))
    max_log = float(np.nanmax(logz))
    levels = np.linspace(min_log, max_log, 14)

    fig, ax = plt.subplots(figsize=(10, 7.5))
    cf = ax.contourf(x_mesh, y_mesh, logz, levels=levels, cmap="viridis")
    ax.contour(
        x_mesh,
        y_mesh,
        logz,
        levels=levels[::2],
        colors="white",
        linewidths=0.6,
        alpha=0.7,
    )

    cbar = fig.colorbar(cf, ax=ax)
    cbar.set_label("log10(total dissolved Zn molality [mol/kg-H2O])")

    ax.set_xlabel(f"log10(a({ACTIVITY_SPECIES_X}))")
    ax.set_ylabel(f"log10(a({ACTIVITY_SPECIES_Y}))")
    ax.set_title(
        "Zn Solubility Contours in Activity-Activity Space\n"
        f"T={TEMPERATURE_C:.0f} C, P={PRESSURE_KBAR:.1f} kbar"
    )
    ax.set_xlim(LOG_A_X_MIN, LOG_A_X_MAX)
    ax.set_ylim(LOG_A_Y_MIN, LOG_A_Y_MAX)
    ax.grid(True, alpha=0.2, linestyle="--")

    fig.tight_layout()
    fig.savefig(OUTPUT_SOLUBILITY_CONTOUR_FIGURE, dpi=250)
    plt.close(fig)
    return True


def _save_overlay_plot(x_vals, y_vals, label_grid, zn_grid, categories, cmap):
    valid = np.isfinite(zn_grid) & (zn_grid > 0.0)
    if not np.any(valid):
        return False

    x_mesh, y_mesh = np.meshgrid(x_vals, y_vals)
    logz = np.full_like(zn_grid, np.nan, dtype=float)
    logz[valid] = np.log10(zn_grid[valid])

    fig, ax = plt.subplots(figsize=(10, 7.5))
    overlay_alpha = 0.55
    phase_norm = plt.Normalize(vmin=-0.5, vmax=len(categories) - 0.5)
    im = ax.imshow(
        label_grid,
        origin="lower",
        interpolation="nearest",
        extent=(LOG_A_X_MIN, LOG_A_X_MAX, LOG_A_Y_MIN, LOG_A_Y_MAX),
        aspect="auto",
        cmap=cmap,
        norm=phase_norm,
        alpha=overlay_alpha,
    )

    levels = np.linspace(float(np.nanmin(logz)), float(np.nanmax(logz)), 11)
    contours = ax.contour(
        x_mesh,
        y_mesh,
        logz,
        levels=levels,
        colors="black",
        linewidths=0.8,
    )
    ax.clabel(contours, inline=True, fontsize=7, fmt="%.2f")

    present_indices = sorted({int(idx) for idx in np.unique(label_grid)})
    handles = [
        plt.Line2D(
            [0],
            [0],
            marker="s",
            linestyle="",
            markersize=8,
            markerfacecolor=cmap(phase_norm(idx)),
            markeredgecolor="black",
            alpha=overlay_alpha,
            label=categories[idx],
        )
        for idx in present_indices
    ]
    ax.legend(
        handles=handles, loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False
    )

    ax.set_xlabel(f"log10(a({ACTIVITY_SPECIES_X}))")
    ax.set_ylabel(f"log10(a({ACTIVITY_SPECIES_Y}))")
    ax.set_title(
        "Zn Mineral Stability + Solubility Contours\n"
        f"T={TEMPERATURE_C:.0f} C, P={PRESSURE_KBAR:.1f} kbar"
    )
    ax.set_xlim(LOG_A_X_MIN, LOG_A_X_MAX)
    ax.set_ylim(LOG_A_Y_MIN, LOG_A_Y_MAX)
    ax.grid(True, alpha=0.2, linestyle="--")

    fig.tight_layout()
    fig.savefig(OUTPUT_OVERLAY_FIGURE, dpi=250)
    plt.close(fig)
    return True


def _save_solubility_sections_plot(x_vals, y_vals, zn_grid):
    fig, ax = plt.subplots(figsize=(9.5, 6.0))

    pH_vals = -x_vals
    plotted = 0
    for target_y in SECTION_LOG_A_HS_VALUES:
        j = int(np.argmin(np.abs(y_vals - target_y)))
        y_used = float(y_vals[j])
        row = zn_grid[j, :]
        valid = np.isfinite(row) & (row > 0.0)
        if not np.any(valid):
            continue
        ax.plot(
            pH_vals[valid],
            row[valid],
            linewidth=2.0,
            label=f"log10(a(HS-))={y_used:.2f}",
        )
        plotted += 1

    if plotted == 0:
        plt.close(fig)
        return False

    ax.set_yscale("log")
    ax.set_xlabel("pH (=-log10(a(H+)))")
    ax.set_ylabel("Total dissolved Zn molality (mol/kg-H2O)")
    ax.set_title(
        "Zn Solubility Sections at Fixed HS- Activity\n"
        f"T={TEMPERATURE_C:.0f} C, P={PRESSURE_KBAR:.1f} kbar"
    )
    ax.grid(True, which="both", alpha=0.25, linestyle="--")
    ax.legend(frameon=False)

    fig.tight_layout()
    fig.savefig(OUTPUT_CROSS_SECTION_FIGURE, dpi=250)
    plt.close(fig)
    return True


def _save_redox_sulfur_sweeps(module):
    system_o2 = module.build_system_with_o2_gas_phase()
    system_h2s = module.build_system_with_h2s_gas_phase()
    solvent_species_name = module.infer_solvent_species_name(
        module.AQUEOUS_SPECIES,
        module.SOLVENT_SPECIES_NAME,
        module.INITIAL_SPECIES_AMOUNTS_MOL,
    )

    logf_o2, y_o2 = module.compute_true_fo2_sensitivity(system_o2, solvent_species_name)
    logf_h2s, y_h2s = module.compute_true_fh2s_sensitivity(
        system_h2s, solvent_species_name
    )

    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 4.8))

    valid_o2 = np.isfinite(y_o2) & (y_o2 > 0.0)
    valid_h2s = np.isfinite(y_h2s) & (y_h2s > 0.0)

    axes[0].plot(logf_o2[valid_o2], y_o2[valid_o2], color="tab:red", linewidth=2.0)
    axes[0].set_yscale("log")
    axes[0].set_xlabel("log10(fO2 / bar)")
    axes[0].set_ylabel("Total dissolved Zn molality (mol/kg-H2O)")
    axes[0].set_title("Redox sweep")
    axes[0].grid(True, which="both", alpha=0.25, linestyle="--")

    axes[1].plot(
        logf_h2s[valid_h2s], y_h2s[valid_h2s], color="tab:green", linewidth=2.0
    )
    axes[1].set_yscale("log")
    axes[1].set_xlabel("log10(fH2S / bar)")
    axes[1].set_ylabel("Total dissolved Zn molality (mol/kg-H2O)")
    axes[1].set_title("Sulfur fugacity sweep")
    axes[1].grid(True, which="both", alpha=0.25, linestyle="--")

    fig.suptitle(
        "Zn Solubility vs True Fugacity Constraints\n"
        f"T={module.SENS_TEMPERATURE_C:.0f} C, P={module.SENS_PRESSURE_KBAR:.1f} kbar",
        fontsize=11,
    )
    fig.tight_layout()
    fig.savefig(OUTPUT_REDOX_SULFUR_FIGURE, dpi=250)
    plt.close(fig)

    return int(valid_o2.sum()), int(valid_h2s.sum())


def main():
    module = load_tutorial_module(TUTORIAL_PATH)

    try:
        module.Warnings.disable(906)
    except Exception:
        pass

    module.validate_user_inputs()

    # Ensure multi-mineral competition is active for phase classification.
    module.USE_COMPETING_ZN_MINERALS = True

    system = module.build_tutorial_system()

    specs = module.EquilibriumSpecs(system)
    specs.temperature()
    specs.pressure()
    specs.lgActivity(ACTIVITY_SPECIES_X)
    specs.lgActivity(ACTIVITY_SPECIES_Y)

    solver = module.make_equilibrium_solver(system, specs)
    conditions = module.EquilibriumConditions(specs)

    pressure_bar = PRESSURE_KBAR * 1000.0
    conditions.temperature(float(TEMPERATURE_C), "celsius")
    conditions.pressure(float(pressure_bar), "bar")

    categories = list(module.selected_mineral_names()) + [
        "AqueousOnly",
        "NoConvergence",
    ]
    category_to_index = {name: idx for idx, name in enumerate(categories)}

    evaluate, cache, eval_stats = make_label_evaluator(
        module,
        system,
        solver,
        conditions,
        category_to_index,
    )

    x_edges = np.linspace(LOG_A_X_MIN, LOG_A_X_MAX, COARSE_NX + 1)
    y_edges = np.linspace(LOG_A_Y_MIN, LOG_A_Y_MAX, COARSE_NY + 1)

    leaf_cells = []
    for j in range(COARSE_NY):
        y0 = y_edges[j]
        y1 = y_edges[j + 1]
        for i in range(COARSE_NX):
            x0 = x_edges[i]
            x1 = x_edges[i + 1]
            refine_cell(
                x0,
                x1,
                y0,
                y1,
                depth=0,
                max_depth=MAX_REFINEMENT_LEVEL,
                evaluate=evaluate,
                leaf_cells=leaf_cells,
            )

    cmap = ListedColormap(plt.cm.tab20(np.linspace(0.0, 1.0, len(categories))))

    fig, ax = plt.subplots(figsize=(10, 7.5))

    _plot_phase_rectangles(ax, leaf_cells, cmap)

    ax.set_xlabel(f"log10(a({ACTIVITY_SPECIES_X}))")
    ax.set_ylabel(f"log10(a({ACTIVITY_SPECIES_Y}))")
    ax.set_title(
        "Zn Mineral Phase Fields in Activity-Activity Space\n"
        f"T={TEMPERATURE_C:.0f} C, P={PRESSURE_KBAR:.1f} kbar"
    )
    ax.set_xlim(LOG_A_X_MIN, LOG_A_X_MAX)
    ax.set_ylim(LOG_A_Y_MIN, LOG_A_Y_MAX)
    ax.grid(True, alpha=0.25, linestyle="--")

    present_indices = sorted({int(cell[4]) for cell in leaf_cells})
    present_labels = [categories[idx] for idx in present_indices]
    leaf_label_counts = {name: 0 for name in categories}
    for _, _, _, _, label_idx, _ in leaf_cells:
        leaf_label_counts[categories[int(label_idx)]] += 1

    sample_label_counts = {name: 0 for name in categories}
    for idx in cache.values():
        sample_label_counts[categories[int(idx)]] += 1
    handles = [
        plt.Line2D(
            [0],
            [0],
            marker="s",
            linestyle="",
            markersize=9,
            markerfacecolor=cmap(idx),
            markeredgecolor="black",
            label=label,
        )
        for idx, label in zip(present_indices, present_labels)
    ]

    ax.legend(
        handles=handles, loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False
    )

    fig.tight_layout()
    fig.savefig(OUTPUT_FIGURE, dpi=250)
    plt.close(fig)

    print("Generated activity-activity phase diagram:")
    print(OUTPUT_FIGURE)
    print(
        "Adaptive mesh stats: "
        f"coarse={COARSE_NX}x{COARSE_NY}, max_level={MAX_REFINEMENT_LEVEL}, "
        f"leaf_cells={len(leaf_cells)}, equilibrium_solves={len(cache)}"
    )
    print("Leaf-cell label counts:")
    for name, count in leaf_label_counts.items():
        if count > 0:
            print(f"  {name}: {count}")
    print("Sample-point label counts:")
    for name, count in sample_label_counts.items():
        if count > 0:
            print(f"  {name}: {count}")
    if USE_STRICT_DOMAIN_CLIP:
        print(
            "Domain-clipped sample points (non-phase diagnostic): "
            f"{eval_stats['domain_clipped_points']}"
        )
        print(
            "Domain-clipped points with inferred mineral label: "
            f"{eval_stats['domain_clipped_inferred_points']}"
        )
        print(
            "Domain-clipped unresolved points: "
            f"{eval_stats['domain_clipped_unresolved_points']}"
        )
    print(
        "Inferred-from-nearby sample points (non-converged exact point): "
        f"{eval_stats['inferred_from_shift_points']}"
    )
    print(
        "True NoConvergence sample points after inference: "
        f"{eval_stats['true_noconvergence_points']}"
    )

    print("Computing regular-grid stability+solubility data...")
    x_vals, y_vals, label_grid, zn_grid, grid_stats = _solve_regular_grid(
        module,
        categories,
        category_to_index,
    )

    print("Saving contour/overlay/section figures...")
    contour_ok = _save_solubility_contour_plot(x_vals, y_vals, zn_grid)
    overlay_ok = _save_overlay_plot(
        x_vals, y_vals, label_grid, zn_grid, categories, cmap
    )
    sections_ok = _save_solubility_sections_plot(x_vals, y_vals, zn_grid)
    print("Computing redox/sulfur sweep figure...")
    valid_o2, valid_h2s = _save_redox_sulfur_sweeps(module)

    print(
        "Regular-grid stats: "
        f"points={grid_stats['grid_points']}, "
        f"converged={grid_stats['grid_converged']}, "
        f"inferred={grid_stats['grid_inferred']}, "
        f"true_noconvergence={grid_stats['grid_true_noconvergence']}, "
        f"domain_clipped={grid_stats['grid_domain_clipped']}"
    )
    if contour_ok:
        print("Generated solubility contour diagram:")
        print(OUTPUT_SOLUBILITY_CONTOUR_FIGURE)
    if overlay_ok:
        print("Generated phase+solubility overlay diagram:")
        print(OUTPUT_OVERLAY_FIGURE)
    if sections_ok:
        print("Generated solubility section diagram:")
        print(OUTPUT_CROSS_SECTION_FIGURE)
    print("Generated redox/sulfur sweep diagram:")
    print(OUTPUT_REDOX_SULFUR_FIGURE)
    print(
        "Fugacity sweep convergence: "
        f"fO2={valid_o2}/{len(module.LOG_FO2_RANGE)}, "
        f"fH2S={valid_h2s}/{len(module.LOG_FH2S_RANGE)}"
    )


if __name__ == "__main__":
    main()
