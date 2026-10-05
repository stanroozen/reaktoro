"""Matplotlib visualizations for MAP/Laplace/MCMC/grid posterior diagnostics.

This module is optional: importing it requires matplotlib, but no other
thermoinvert module depends on it. Every function takes an already-computed
result object (no recomputation) and returns the created ``Figure``.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Ellipse

from .bias import GroupedObservationBiasSummary
from .posterior import GridPosteriorResult, LaplaceApproximation, MCMCResult
from .sensitivity import PosteriorResolutionResult
from .validation import (
    BracketPosteriorPredictiveResult,
    LOOCVResult,
    PosteriorPredictiveResult,
)


def _covariance_ellipse(
    center: tuple[float, float],
    covariance: np.ndarray,
    sigma_level: float,
) -> Ellipse:
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    order = np.argsort(eigenvalues)[::-1]
    eigenvalues, eigenvectors = eigenvalues[order], eigenvectors[:, order]
    angle = float(np.degrees(np.arctan2(eigenvectors[1, 0], eigenvectors[0, 0])))
    width, height = 2.0 * sigma_level * np.sqrt(np.clip(eigenvalues, 0.0, None))
    return Ellipse(
        xy=center,
        width=float(width),
        height=float(height),
        angle=angle,
        edgecolor="tab:red",
        facecolor="none",
        linewidth=1.5,
        label=f"{sigma_level:g}\u03c3 Laplace covariance",
    )


def plot_mcmc_trace(mcmc_result: MCMCResult) -> plt.Figure:
    """Per-parameter trace plot with a running mean, for basic mixing checks."""
    n = len(mcmc_result.parameter_names)
    fig, axes = plt.subplots(n, 1, squeeze=False, figsize=(8.0, 2.5 * n))
    steps = np.arange(mcmc_result.samples.shape[0])
    for index, name in enumerate(mcmc_result.parameter_names):
        ax = axes[index, 0]
        values = mcmc_result.samples[:, index]
        ax.plot(steps, values, color="tab:blue", linewidth=0.5, alpha=0.7)
        running_mean = np.cumsum(values) / np.arange(1, len(values) + 1)
        ax.plot(
            steps, running_mean, color="tab:orange", linewidth=1.5, label="running mean"
        )
        ax.set_ylabel(name)
        ax.legend(loc="upper right", fontsize=8)
    axes[-1, 0].set_xlabel("post-burn-in sample")
    fig.tight_layout()
    return fig


def plot_posterior_scatter(
    mcmc_result: MCMCResult,
    laplace: LaplaceApproximation | None = None,
    sigma_level: float = 1.0,
) -> plt.Figure:
    """Scatter (2 parameters) or histogram (1 parameter) of posterior samples,
    optionally overlaid with a Laplace covariance ellipse."""
    n = len(mcmc_result.parameter_names)
    if n == 1:
        fig, ax = plt.subplots()
        ax.hist(mcmc_result.samples[:, 0], bins=40, color="tab:blue", alpha=0.7)
        ax.set_xlabel(mcmc_result.parameter_names[0])
        ax.set_ylabel("count")
        return fig
    if n != 2:
        raise ValueError("plot_posterior_scatter supports one or two parameters")

    fig, ax = plt.subplots()
    ax.scatter(
        mcmc_result.samples[:, 0],
        mcmc_result.samples[:, 1],
        s=4,
        alpha=0.25,
        color="tab:blue",
        label="posterior samples",
    )
    ax.set_xlabel(mcmc_result.parameter_names[0])
    ax.set_ylabel(mcmc_result.parameter_names[1])
    if laplace is not None and laplace.is_positive_definite:
        ellipse = _covariance_ellipse(
            (laplace.map_point[0], laplace.map_point[1]),
            laplace.covariance,
            sigma_level,
        )
        ax.add_patch(ellipse)
        ax.scatter(
            *laplace.map_point, color="tab:red", marker="x", s=60, label="MAP point"
        )
    ax.legend(loc="best", fontsize=8)
    return fig


def plot_grid_posterior(grid_result: GridPosteriorResult) -> plt.Figure:
    """1D density-with-credible-band, or 2D filled-contour, posterior plot."""
    n = len(grid_result.parameter_names)
    if n == 1:
        name = grid_result.parameter_names[0]
        axis = grid_result.axes[0]
        mean = grid_result.marginal_mean[name]
        std = grid_result.marginal_std[name]
        fig, ax = plt.subplots()
        ax.plot(axis, grid_result.posterior_density, color="tab:blue")
        ax.axvspan(
            mean - 1.96 * std,
            mean + 1.96 * std,
            color="tab:blue",
            alpha=0.15,
            label="~95% band",
        )
        ax.axvline(mean, color="tab:orange", linestyle="--", label="mean")
        ax.set_xlabel(name)
        ax.set_ylabel("posterior density")
        ax.legend(loc="best", fontsize=8)
        return fig
    if n == 2:
        x_axis, y_axis = grid_result.axes
        fig, ax = plt.subplots()
        contour = ax.contourf(
            x_axis, y_axis, grid_result.posterior_density.T, levels=20, cmap="viridis"
        )
        fig.colorbar(contour, ax=ax, label="posterior density")
        ax.set_xlabel(grid_result.parameter_names[0])
        ax.set_ylabel(grid_result.parameter_names[1])
        return fig
    raise ValueError("plot_grid_posterior supports one or two parameters")


def plot_posterior_predictive_check(
    ppc_result: PosteriorPredictiveResult,
) -> plt.Figure:
    """Observed value vs. predictive mean/95% interval, per observable."""
    names = ppc_result.observable_names
    y = np.arange(len(names))
    fig, ax = plt.subplots(figsize=(6.0, 0.4 * len(names) + 1.0))
    for index, name in enumerate(names):
        low, high = ppc_result.credible_interval_95[name]
        ax.plot([low, high], [index, index], color="tab:blue", linewidth=4, alpha=0.4)
        ax.scatter(
            [ppc_result.predictive_mean[index]],
            [index],
            color="tab:blue",
            marker="o",
            label="predictive mean" if index == 0 else None,
        )
        color = "tab:green" if ppc_result.within_95[name] else "tab:red"
        ax.scatter(
            [ppc_result.observed_values[index]],
            [index],
            color=color,
            marker="x",
            s=60,
            label="observed" if index == 0 else None,
        )
    ax.set_yticks(y)
    ax.set_yticklabels(names)
    ax.set_xlabel("value")
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    return fig


def plot_loocv_misfits(loocv_result: LOOCVResult) -> plt.Figure:
    """Per-experiment held-out misfit bar chart with the mean overlaid."""
    ids = [fold.held_out_experiment_id for fold in loocv_result.folds]
    misfits = [fold.misfit for fold in loocv_result.folds]
    spread = [
        float(np.std(fold.held_out_misfits, ddof=1))
        if len(fold.held_out_misfits) > 1
        else 0.0
        for fold in loocv_result.folds
    ]
    x = np.arange(len(ids))
    fig, ax = plt.subplots(figsize=(max(4.0, 0.6 * len(ids)), 4.0))
    ax.bar(
        x,
        misfits,
        yerr=spread if any(spread) else None,
        capsize=4.0,
        color="tab:blue",
    )
    ax.axhline(
        loocv_result.mean_misfit,
        color="tab:orange",
        linestyle="--",
        label="mean misfit",
    )
    ax.set_xticks(x)
    ax.set_xticklabels(ids, rotation=45, ha="right")
    ax.set_ylabel("held-out misfit")
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    return fig


def plot_grouped_observation_bias_summary(
    summary: GroupedObservationBiasSummary,
) -> plt.Figure:
    """Plot group-by-observation-key bias estimates as a labeled heatmap."""
    groups = tuple(summary.group_biases)
    keys = summary.observation_keys
    if not groups or not keys:
        raise ValueError("bias summary must contain groups and observation keys")
    values = np.asarray(
        [
            [summary.group_biases[group].get(key, np.nan) for key in keys]
            for group in groups
        ],
        dtype=float,
    )
    fig, ax = plt.subplots(
        figsize=(max(5.0, 1.2 * len(keys)), max(3.0, 0.6 * len(groups)))
    )
    image = ax.imshow(values, aspect="auto", cmap="coolwarm")
    fig.colorbar(image, ax=ax, label="observation bias")
    ax.set_xticks(np.arange(len(keys)))
    ax.set_xticklabels(keys, rotation=45, ha="right")
    ax.set_yticks(np.arange(len(groups)))
    ax.set_yticklabels(groups)
    ax.set_xlabel("observation key")
    ax.set_ylabel(summary.metadata_key)
    ax.set_title("Grouped observation bias")
    fig.tight_layout()
    return fig


def plot_posterior_resolution_summary(
    result: PosteriorResolutionResult,
) -> plt.Figure:
    """Plot mean local resolution and its elementwise posterior interval width."""
    names = result.parameter_names
    expected_shape = (len(names), len(names))
    mean = np.asarray(result.mean_matrix, dtype=float)
    lower = np.asarray(result.lower_matrix, dtype=float)
    upper = np.asarray(result.upper_matrix, dtype=float)
    if not names or mean.shape != expected_shape:
        raise ValueError("resolution matrices must match non-empty parameter names")
    if lower.shape != expected_shape or upper.shape != expected_shape:
        raise ValueError("resolution interval matrices must match parameter names")
    if not all(np.all(np.isfinite(matrix)) for matrix in (mean, lower, upper)):
        raise ValueError("resolution matrices must contain finite values")
    if np.any(lower > upper):
        raise ValueError("lower resolution bounds cannot exceed upper bounds")

    interval_width = upper - lower
    color_limit = float(np.max(np.abs(mean)))
    if color_limit == 0.0:
        color_limit = 1.0
    fig, axes = plt.subplots(1, 2, figsize=(max(8.0, 1.2 * len(names)), 4.5))
    mean_image = axes[0].imshow(
        mean,
        aspect="auto",
        cmap="coolwarm",
        norm=TwoSlopeNorm(vmin=-color_limit, vcenter=0.0, vmax=color_limit),
    )
    width_image = axes[1].imshow(interval_width, aspect="auto", cmap="viridis")
    axes[0].set_title("Mean local resolution")
    axes[1].set_title("Elementwise interval width")
    for ax in axes:
        ax.set_xticks(np.arange(len(names)))
        ax.set_xticklabels(names, rotation=45, ha="right")
        ax.set_yticks(np.arange(len(names)))
        ax.set_yticklabels(names)
        ax.set_xlabel("parameter")
    axes[0].set_ylabel("parameter")
    fig.colorbar(mean_image, ax=axes[0], label="mean resolution")
    fig.colorbar(width_image, ax=axes[1], label="upper - lower")
    fig.suptitle(f"Posterior-wide resolution ({result.sample_count} samples)")
    fig.tight_layout()
    return fig


def plot_bracket_posterior_predictive(
    result: BracketPosteriorPredictiveResult,
) -> plt.Figure:
    """Plot posterior bracket satisfaction and misfit summary."""
    fig, axes = plt.subplots(1, 2, figsize=(8.0, 3.5))
    axes[0].bar(
        ["satisfied", "violated"],
        [result.satisfied_fraction, 1.0 - result.satisfied_fraction],
        color=["tab:green", "tab:red"],
    )
    axes[0].set_ylim(0.0, 1.0)
    axes[0].set_ylabel("posterior fraction")
    axes[0].set_title("Bracket coverage")

    axes[1].bar(
        ["mean", "2.5%", "97.5%"],
        [
            result.mean_misfit,
            result.credible_interval_95[0],
            result.credible_interval_95[1],
        ],
        color="tab:blue",
    )
    axes[1].set_ylabel("bracket misfit")
    axes[1].set_title("Posterior misfit")
    fig.suptitle(
        f"{len(result.bracket_ids)} reaction brackets; n={result.samples_evaluated}"
    )
    fig.tight_layout()
    return fig


def plot_bracket_misfit_distribution(
    result: BracketPosteriorPredictiveResult,
) -> plt.Figure:
    """Plot the posterior distribution of aggregate bracket misfit."""
    if not result.misfits:
        raise ValueError("Bracket predictive result has no individual misfits")
    fig, ax = plt.subplots(figsize=(6.0, 4.0))
    ax.hist(result.misfits, bins="auto", color="tab:blue", alpha=0.75)
    ax.axvline(0.0, color="tab:green", linestyle="--", label="satisfied")
    ax.axvline(
        result.mean_misfit,
        color="tab:orange",
        linestyle="--",
        label="mean misfit",
    )
    ax.set_xlabel("aggregate bracket misfit")
    ax.set_ylabel("posterior sample count")
    ax.set_title(f"Bracket PPC; satisfied fraction={result.satisfied_fraction:.3f}")
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    return fig


def save_figure(fig: plt.Figure, path: str | Path) -> Path:
    """Save a figure as PDF (vector) or PNG (300 DPI raster) by extension."""
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    dpi = 300 if output.suffix.lower() == ".png" else None
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    return output
