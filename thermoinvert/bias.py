"""Profiled experiment and laboratory observation-bias utilities."""

from dataclasses import dataclass, replace
from typing import Sequence

import numpy as np

from .inversion import InversionProblem, MAPResult


@dataclass(frozen=True)
class ObservationBiasProfileResult:
    """MAP profile of one additive bias per metadata group."""

    biases: dict[str, float]
    objective: float
    corrections_j_mol: dict[str, float]
    bias_prior_penalty: float
    metadata_key: str
    observation_key: str
    records: tuple[dict[str, object], ...]
    bias_bounds: tuple[float, float] = (-1.0, 1.0)
    bias_prior_mean: float = 0.0
    bias_prior_sigma: float = 0.1


@dataclass(frozen=True)
class ObservationBiasSummary:
    """Compact summary of one profiled observation bias across groups."""

    metadata_key: str
    observation_key: str
    group_biases: dict[str, float]
    bias_prior_mean: float
    bias_prior_sigma: float
    bias_bounds: tuple[float, float]
    objective: float
    bias_prior_penalty: float


@dataclass(frozen=True)
class GroupedObservationBiasReport:
    """Compact summary across multiple observation keys and metadata groups."""

    metadata_key: str
    observation_keys: tuple[str, ...]
    group_biases: dict[str, dict[str, float]]
    bias_prior_mean: float
    bias_prior_sigma: float
    bias_bounds: tuple[float, float]
    objective: float
    bias_prior_penalty: float


@dataclass(frozen=True)
class GroupedObservationBiasSummary:
    """Table-ready summary for a grouped observation-bias report."""

    metadata_key: str
    observation_keys: tuple[str, ...]
    group_biases: dict[str, dict[str, float]]
    bias_prior_mean: float
    bias_prior_sigma: float
    bias_bounds: tuple[float, float]
    objective: float
    bias_prior_penalty: float


def profile_observation_bias_by_group(
    problem: InversionProblem,
    metadata_key: str,
    observation_key: str,
    bias_bounds: tuple[float, float] = (-1.0, 1.0),
    bias_prior_mean: float = 0.0,
    bias_prior_sigma: float = 0.1,
    bias_iterations: int = 8,
    optimize_starts: int = 4,
    optimize_iterations: int = 60,
) -> ObservationBiasProfileResult:
    """Estimate one additive bias per metadata group by nested MAP refits.

    A positive bias means a reported continuous observation is high by that
    amount. The profiler uses a Gaussian prior for each bias and preserves all
    parent inversion covariance and model-error settings.
    """
    lower, upper = bias_bounds
    if not metadata_key or not observation_key:
        raise ValueError("metadata_key and observation_key must be non-empty")
    if not np.isfinite(lower) or not np.isfinite(upper) or lower >= upper:
        raise ValueError("bias_bounds must be finite and increasing")
    if (
        not np.isfinite(bias_prior_mean)
        or not np.isfinite(bias_prior_sigma)
        or bias_prior_sigma <= 0.0
    ):
        raise ValueError("bias prior must have finite mean and positive sigma")
    if bias_iterations < 1:
        raise ValueError("bias_iterations must be positive")
    groups = {
        experiment.metadata.get(metadata_key) for experiment in problem.experiments
    }
    if not groups or any(not isinstance(group, str) or not group for group in groups):
        raise ValueError("all experiments need non-empty string bias groups")
    for experiment in problem.experiments:
        observed_keys = {
            *(f"fraction:{phase}" for phase in (experiment.phase_fractions or {})),
            *(
                f"composition:{phase}:{component}"
                for phase, values in (experiment.phase_compositions or {}).items()
                for component in values
            ),
        }
        if observation_key not in observed_keys:
            raise ValueError(
                f"Experiment {experiment.id!r} lacks observation {observation_key!r}"
            )
    biases = {group: float(np.clip(bias_prior_mean, lower, upper)) for group in groups}
    golden = (1.0 + np.sqrt(5.0)) / 2.0
    records: list[dict[str, object]] = []
    best_fit: MAPResult | None = None
    best_objective = float("inf")

    def fit_at(candidate: dict[str, float]) -> tuple[MAPResult, float]:
        experiments = []
        for experiment in problem.experiments:
            observation_bias = dict(experiment.observation_bias or {})
            observation_bias[observation_key] = candidate[
                experiment.metadata[metadata_key]
            ]
            experiments.append(replace(experiment, observation_bias=observation_bias))
        profiled = InversionProblem(
            problem.forward_model,
            experiments,
            problem.parameters,
            brackets=problem.brackets,
            missing_phase_penalty=problem.missing_phase_penalty,
            extra_phase_penalty=problem.extra_phase_penalty,
            marginalize_condition_uncertainty=problem.marginalize_condition_uncertainty,
            condition_quadrature_order=problem.condition_quadrature_order,
            objective_criterion=problem.objective_criterion,
            misfit_form=problem.misfit_form,
            prior_covariance=problem.prior_covariance,
            theoretical_covariance_scale=problem.theoretical_covariance_scale,
            theoretical_covariance_scale_metadata_key=problem.theoretical_covariance_scale_metadata_key,
            theoretical_covariance_scales=problem.theoretical_covariance_scales,
        )
        fit = (
            profiled.optimize_one_parameter(
                starts=optimize_starts, iterations=optimize_iterations
            )
            if len(problem.parameters) == 1
            else profiled.optimize_multi_parameter(
                starts=optimize_starts, iterations=optimize_iterations
            )
        )
        penalty = sum(
            0.5 * ((value - bias_prior_mean) / bias_prior_sigma) ** 2
            for value in candidate.values()
        )
        return fit, float(fit.objective + penalty)

    for group in sorted(biases):
        left, right = lower, upper
        while right - left > (upper - lower) / (golden**bias_iterations):
            first = right - (right - left) / golden
            second = left + (right - left) / golden
            first_candidate, second_candidate = dict(biases), dict(biases)
            first_candidate[group], second_candidate[group] = first, second
            first_fit, first_objective = fit_at(first_candidate)
            second_fit, second_objective = fit_at(second_candidate)
            records.extend(
                [
                    {"group": group, "bias": first, "objective": first_objective},
                    {"group": group, "bias": second, "objective": second_objective},
                ]
            )
            if first_objective < second_objective:
                right, biases = second, first_candidate
                chosen_fit, chosen_objective = first_fit, first_objective
            else:
                left, biases = first, second_candidate
                chosen_fit, chosen_objective = second_fit, second_objective
            if chosen_objective < best_objective:
                best_fit, best_objective = chosen_fit, chosen_objective
    best_fit, best_objective = fit_at(biases)
    penalty = sum(
        0.5 * ((value - bias_prior_mean) / bias_prior_sigma) ** 2
        for value in biases.values()
    )
    return ObservationBiasProfileResult(
        biases,
        best_objective,
        best_fit.corrections_j_mol,
        float(penalty),
        metadata_key,
        observation_key,
        tuple(records),
        (lower, upper),
        bias_prior_mean,
        bias_prior_sigma,
    )


def summarize_group_observation_bias(
    result: ObservationBiasProfileResult,
) -> ObservationBiasSummary:
    """Return a compact table-ready summary of profiled group biases."""
    return ObservationBiasSummary(
        result.metadata_key,
        result.observation_key,
        dict(result.biases),
        result.bias_prior_mean,
        result.bias_prior_sigma,
        result.bias_bounds,
        result.objective,
        result.bias_prior_penalty,
    )


def summarize_grouped_observation_biases(
    profiles: dict[str, ObservationBiasProfileResult],
) -> GroupedObservationBiasReport:
    """Summarize multiple profiled observation keys in one grouped report."""
    if not profiles:
        raise ValueError(
            "profiles must contain at least one ObservationBiasProfileResult"
        )
    metadata_key = next(iter(profiles.values())).metadata_key
    observation_keys = tuple(profile.observation_key for profile in profiles.values())
    group_biases = {
        group: {
            key: profile.biases[group]
            for key, profile in profiles.items()
            if group in profile.biases
        }
        for group in sorted(
            {group for profile in profiles.values() for group in profile.biases}
        )
    }
    first = next(iter(profiles.values()))
    objective = sum(profile.objective for profile in profiles.values())
    prior_penalty = sum(profile.bias_prior_penalty for profile in profiles.values())
    return GroupedObservationBiasReport(
        metadata_key,
        observation_keys,
        group_biases,
        first.bias_prior_mean,
        first.bias_prior_sigma,
        first.bias_bounds,
        objective,
        prior_penalty,
    )


def summarize_grouped_bias_report(
    report: GroupedObservationBiasReport,
) -> GroupedObservationBiasSummary:
    """Return a table-ready grouped-bias summary object for display."""
    return GroupedObservationBiasSummary(
        report.metadata_key,
        report.observation_keys,
        dict(report.group_biases),
        report.bias_prior_mean,
        report.bias_prior_sigma,
        report.bias_bounds,
        report.objective,
        report.bias_prior_penalty,
    )


def profile_grouped_bias_summary(
    problem: InversionProblem,
    metadata_key: str,
    observation_keys: str | Sequence[str],
    bias_bounds: tuple[float, float] = (-1.0, 1.0),
    bias_prior_mean: float = 0.0,
    bias_prior_sigma: float = 0.1,
    bias_iterations: int = 8,
    optimize_starts: int = 4,
    optimize_iterations: int = 60,
) -> GroupedObservationBiasSummary:
    """Profile grouped observation biases and return a compact summary object."""
    profiles = {
        key: profile_observation_bias_by_group(
            problem,
            metadata_key,
            key,
            bias_bounds=bias_bounds,
            bias_prior_mean=bias_prior_mean,
            bias_prior_sigma=bias_prior_sigma,
            bias_iterations=bias_iterations,
            optimize_starts=optimize_starts,
            optimize_iterations=optimize_iterations,
        )
        for key in (
            (observation_keys,)
            if isinstance(observation_keys, str)
            else tuple(observation_keys)
        )
    }
    report = summarize_grouped_observation_biases(profiles)
    return summarize_grouped_bias_report(report)


def profile_observation_biases_by_group(
    problem: InversionProblem,
    metadata_key: str,
    observation_keys: str | Sequence[str],
    bias_bounds: tuple[float, float] = (-1.0, 1.0),
    bias_prior_mean: float = 0.0,
    bias_prior_sigma: float = 0.1,
    bias_iterations: int = 8,
    optimize_starts: int = 4,
    optimize_iterations: int = 60,
) -> dict[str, dict[str, float]]:
    """Estimate one additive bias per metadata group for each observation key.

    This is a small but useful extension of the single-key profile: each key is
    profiled independently with the same shared bounds and prior, then the group
    results are assembled into a dictionary keyed by metadata group and
    observation key.
    """
    if isinstance(observation_keys, str):
        keys = (observation_keys,)
    else:
        keys = tuple(observation_keys)
    if not keys or any(not isinstance(key, str) or not key for key in keys):
        raise ValueError("observation_keys must contain non-empty string keys")
    per_key_profiles = {
        key: profile_observation_bias_by_group(
            problem,
            metadata_key,
            key,
            bias_bounds=bias_bounds,
            bias_prior_mean=bias_prior_mean,
            bias_prior_sigma=bias_prior_sigma,
            bias_iterations=bias_iterations,
            optimize_starts=optimize_starts,
            optimize_iterations=optimize_iterations,
        )
        for key in keys
    }
    groups = sorted(
        {group for profile in per_key_profiles.values() for group in profile.biases}
    )
    return {
        group: {key: per_key_profiles[key].biases[group] for key in keys}
        for group in groups
    }


def apply_group_observation_biases(
    problem: InversionProblem,
    metadata_key: str,
    group_biases: dict[str, dict[str, float]],
) -> InversionProblem:
    """Apply a precomputed group-bias map to every experiment in a problem.

    The map is keyed by metadata group and observation key, e.g.
    {"lab-a": {"composition:A:x": 0.05}}. Any existing observation bias on the
    same key is overwritten for the matched group.
    """
    if not isinstance(metadata_key, str) or not metadata_key:
        raise ValueError("metadata_key must be a non-empty string")
    if not isinstance(group_biases, dict) or not group_biases:
        raise ValueError("group_biases must be a non-empty dictionary")
    experiments = []
    for experiment in problem.experiments:
        group = experiment.metadata.get(metadata_key)
        if group not in group_biases:
            raise ValueError(f"No group-bias map for {metadata_key}={group!r}")
        bias_map = dict(experiment.observation_bias or {})
        for key, value in group_biases[group].items():
            if not isinstance(key, str) or not key:
                raise ValueError("group bias keys must be non-empty strings")
            if not isinstance(value, (int, float)) or not np.isfinite(value):
                raise ValueError(f"Bias for {key!r} must be finite")
            bias_map[key] = float(value)
        experiments.append(replace(experiment, observation_bias=bias_map))
    return InversionProblem(
        problem.forward_model,
        experiments,
        problem.parameters,
        brackets=problem.brackets,
        missing_phase_penalty=problem.missing_phase_penalty,
        extra_phase_penalty=problem.extra_phase_penalty,
        marginalize_condition_uncertainty=problem.marginalize_condition_uncertainty,
        condition_quadrature_order=problem.condition_quadrature_order,
        objective_criterion=problem.objective_criterion,
        misfit_form=problem.misfit_form,
        prior_covariance=problem.prior_covariance,
        theoretical_covariance_scale=problem.theoretical_covariance_scale,
        theoretical_covariance_scale_metadata_key=problem.theoretical_covariance_scale_metadata_key,
        theoretical_covariance_scales=problem.theoretical_covariance_scales,
    )
