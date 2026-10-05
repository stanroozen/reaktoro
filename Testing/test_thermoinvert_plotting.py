"""Tests for the optional matplotlib visualization module."""

import sys
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")

sys.path.insert(0, str(Path(__file__).parents[1]))

from thermoinvert.plotting import (
    plot_bracket_posterior_predictive,
    plot_bracket_misfit_distribution,
    plot_grid_posterior,
    plot_grouped_observation_bias_summary,
    plot_loocv_misfits,
    plot_mcmc_trace,
    plot_posterior_predictive_check,
    plot_posterior_resolution_summary,
    plot_posterior_scatter,
    save_figure,
)
from thermoinvert.bias import GroupedObservationBiasSummary
from thermoinvert.posterior import GridPosteriorResult, LaplaceApproximation, MCMCResult
from thermoinvert.sensitivity import PosteriorResolutionResult
from thermoinvert.validation import (
    BracketPosteriorPredictiveResult,
    LOOCVFold,
    LOOCVResult,
    PosteriorPredictiveResult,
)


def _fake_mcmc_result(n_parameters: int) -> MCMCResult:
    rng = np.random.default_rng(0)
    samples = rng.normal(size=(500, n_parameters))
    names = tuple(f"p{i}" for i in range(n_parameters))
    return MCMCResult(
        parameter_names=names,
        samples=samples,
        log_posterior=-np.sum(samples**2, axis=1),
        acceptance_rate=0.3,
        posterior_mean={
            name: float(np.mean(samples[:, i])) for i, name in enumerate(names)
        },
        posterior_std={
            name: float(np.std(samples[:, i])) for i, name in enumerate(names)
        },
        credible_interval_95={name: (-2.0, 2.0) for name in names},
    )


def test_plot_mcmc_trace_returns_figure_per_parameter():
    fig = plot_mcmc_trace(_fake_mcmc_result(2))
    assert len(fig.axes) == 2
    matplotlib.pyplot.close(fig)


def test_plot_posterior_scatter_handles_one_and_two_parameters():
    fig1 = plot_posterior_scatter(_fake_mcmc_result(1))
    assert len(fig1.axes) == 1
    matplotlib.pyplot.close(fig1)

    laplace = LaplaceApproximation(
        parameter_names=("p0", "p1"),
        map_point=np.array([0.0, 0.0]),
        hessian=np.eye(2),
        covariance=np.eye(2),
        correlation=np.eye(2),
        is_positive_definite=True,
    )
    fig2 = plot_posterior_scatter(_fake_mcmc_result(2), laplace=laplace)
    assert len(fig2.axes) == 1
    matplotlib.pyplot.close(fig2)


def test_plot_grid_posterior_handles_one_and_two_dimensions():
    axis = np.linspace(-1.0, 1.0, 21)
    density_1d = np.exp(-(axis**2))
    density_1d /= np.sum(density_1d) * (axis[1] - axis[0])
    grid_1d = GridPosteriorResult(
        parameter_names=("p0",),
        axes=(axis,),
        log_posterior=-(axis**2),
        posterior_density=density_1d,
        marginal_mean={"p0": 0.0},
        marginal_std={"p0": 0.5},
    )
    fig1 = plot_grid_posterior(grid_1d)
    matplotlib.pyplot.close(fig1)

    density_2d = np.outer(density_1d, density_1d)
    grid_2d = GridPosteriorResult(
        parameter_names=("p0", "p1"),
        axes=(axis, axis),
        log_posterior=np.zeros((21, 21)),
        posterior_density=density_2d,
        marginal_mean={"p0": 0.0, "p1": 0.0},
        marginal_std={"p0": 0.5, "p1": 0.5},
    )
    fig2 = plot_grid_posterior(grid_2d)
    matplotlib.pyplot.close(fig2)


def test_plot_posterior_predictive_check_marks_coverage():
    ppc = PosteriorPredictiveResult(
        observable_names=("phase_present:A", "phase_fraction:A"),
        observed_values=np.array([1.0, 1.0]),
        predictive_mean=np.array([1.0, 0.95]),
        predictive_std=np.array([0.0, 0.05]),
        credible_interval_95={
            "phase_present:A": (1.0, 1.0),
            "phase_fraction:A": (0.85, 1.05),
        },
        within_95={"phase_present:A": True, "phase_fraction:A": True},
    )
    fig = plot_posterior_predictive_check(ppc)
    assert len(fig.axes) == 1
    matplotlib.pyplot.close(fig)


def test_plot_loocv_misfits_bars_each_fold():
    result = LOOCVResult(
        folds=(
            LOOCVFold("exp-1", {"A": 100.0}, 0.5),
            LOOCVFold("exp-2", {"A": 110.0}, 0.7),
        ),
        mean_misfit=0.6,
    )
    fig = plot_loocv_misfits(result)
    assert len(fig.axes) == 1
    matplotlib.pyplot.close(fig)


def test_plot_loocv_misfits_supports_grouped_within_group_spread():
    result = LOOCVResult(
        folds=(
            LOOCVFold("study=one", {"A": 1.0}, 2.0, ("one-1", "one-2"), (1.0, 3.0)),
            LOOCVFold("study=two", {"A": 1.0}, 1.0, ("two-1", "two-2"), (1.0, 1.0)),
        ),
        mean_misfit=1.5,
        group_metadata_key="study",
    )
    fig = plot_loocv_misfits(result)
    assert len(fig.axes) == 1
    matplotlib.pyplot.close(fig)


def test_plot_grouped_observation_bias_summary_returns_figure():
    summary = GroupedObservationBiasSummary(
        "laboratory",
        ("composition:A:x", "composition:A:y"),
        {"lab-a": {"composition:A:x": 0.1, "composition:A:y": -0.1}},
        0.0,
        0.1,
        (-0.5, 0.5),
        1.0,
        0.2,
    )
    fig = plot_grouped_observation_bias_summary(summary)
    assert len(fig.axes) == 2
    matplotlib.pyplot.close(fig)


def test_plot_posterior_resolution_summary_shows_mean_and_interval_width():
    result = PosteriorResolutionResult(
        ("p0", "p1"),
        100,
        np.asarray([[0.7, -0.1], [0.2, 0.6]]),
        np.asarray([[0.4, -0.3], [0.0, 0.3]]),
        np.asarray([[0.9, 0.1], [0.4, 0.8]]),
        np.asarray([0.7, 0.6]),
        np.asarray([0.4, 0.3]),
        np.asarray([0.9, 0.8]),
    )
    fig = plot_posterior_resolution_summary(result)

    assert len(fig.axes) == 4
    assert [ax.get_title() for ax in fig.axes[:2]] == [
        "Mean local resolution",
        "Elementwise interval width",
    ]
    matplotlib.pyplot.close(fig)


def test_plot_posterior_resolution_summary_rejects_invalid_bounds():
    result = PosteriorResolutionResult(
        ("p0",),
        10,
        np.asarray([[0.5]]),
        np.asarray([[0.8]]),
        np.asarray([[0.2]]),
        np.asarray([0.5]),
        np.asarray([0.8]),
        np.asarray([0.2]),
    )
    try:
        plot_posterior_resolution_summary(result)
    except ValueError as error:
        assert "lower" in str(error)
    else:
        raise AssertionError("Expected ValueError for reversed resolution bounds")


def test_plot_bracket_posterior_predictive_returns_summary_figure():
    result = BracketPosteriorPredictiveResult(
        bracket_ids=("and-ky", "pressure-boundary"),
        samples_evaluated=100,
        satisfied_fraction=0.82,
        mean_misfit=18.0,
        credible_interval_95=(0.0, 100.0),
    )
    fig = plot_bracket_posterior_predictive(result)
    assert len(fig.axes) == 2
    matplotlib.pyplot.close(fig)


def test_plot_bracket_misfit_distribution_returns_histogram():
    result = BracketPosteriorPredictiveResult(
        bracket_ids=("and-ky",),
        samples_evaluated=5,
        satisfied_fraction=0.6,
        mean_misfit=20.0,
        credible_interval_95=(0.0, 100.0),
        misfits=(0.0, 0.0, 0.0, 100.0, 0.0),
    )
    fig = plot_bracket_misfit_distribution(result)
    assert len(fig.axes) == 1
    matplotlib.pyplot.close(fig)


def test_save_figure_writes_png_and_pdf(tmp_path):
    fig = plot_mcmc_trace(_fake_mcmc_result(1))
    png_path = save_figure(fig, tmp_path / "trace.png")
    pdf_path = save_figure(fig, tmp_path / "trace.pdf")
    assert png_path.is_file()
    assert pdf_path.is_file()
    matplotlib.pyplot.close(fig)
