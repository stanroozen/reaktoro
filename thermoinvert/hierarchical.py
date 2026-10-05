"""Hierarchical nuisance-uncertainty utilities."""

from dataclasses import dataclass
from dataclasses import replace
from typing import Mapping, Sequence

import numpy as np

from .inversion import InversionProblem, MAPResult


@dataclass(frozen=True)
class HierarchicalMCMCResult:
    """Joint MCMC samples for inversion corrections and a global model-error scale."""

    parameter_names: tuple[str, ...]
    parameter_units: tuple[str, ...]
    samples: np.ndarray
    log_posterior: np.ndarray
    acceptance_rate: float
    posterior_mean: dict[str, float]
    posterior_std: dict[str, float]
    credible_interval_95: dict[str, tuple[float, float]]
    scale_name: str
    scale_prior_mean: float
    scale_prior_sigma: float
    scale_bounds: tuple[float, float]
    seed: int
    burn_in: int


@dataclass(frozen=True)
class GroupedHierarchicalMCMCResult:
    """Joint MCMC samples for corrections and metadata-group model scales."""

    parameter_names: tuple[str, ...]
    parameter_units: tuple[str, ...]
    samples: np.ndarray
    log_posterior: np.ndarray
    acceptance_rate: float
    posterior_mean: dict[str, float]
    posterior_std: dict[str, float]
    credible_interval_95: dict[str, tuple[float, float]]
    metadata_key: str
    scale_names: tuple[str, ...]
    scale_prior_means: dict[str, float]
    scale_prior_sigmas: dict[str, float]
    scale_bounds: tuple[float, float]
    seed: int
    burn_in: int


@dataclass(frozen=True)
class GroupedObservationBiasMCMCResult:
    """Joint MCMC samples for corrections and one bias per metadata group."""

    parameter_names: tuple[str, ...]
    parameter_units: tuple[str, ...]
    samples: np.ndarray
    log_posterior: np.ndarray
    acceptance_rate: float
    posterior_mean: dict[str, float]
    posterior_std: dict[str, float]
    credible_interval_95: dict[str, tuple[float, float]]
    metadata_key: str
    observation_key: str
    bias_names: tuple[str, ...]
    bias_prior_mean: float
    bias_prior_sigma: float
    bias_bounds: tuple[float, float]
    seed: int
    burn_in: int
    observation_keys: tuple[str, ...] = ()


@dataclass(frozen=True)
class CombinedGroupedUncertaintyMCMCResult:
    """Joint MCMC result for corrections, group model scales, and biases."""

    parameter_names: tuple[str, ...]
    parameter_units: tuple[str, ...]
    samples: np.ndarray
    log_posterior: np.ndarray
    acceptance_rate: float
    posterior_mean: dict[str, float]
    posterior_std: dict[str, float]
    credible_interval_95: dict[str, tuple[float, float]]
    metadata_key: str
    observation_key: str
    groups: tuple[str, ...]
    scale_bounds: tuple[float, float]
    scale_prior_mean: float
    scale_prior_sigma: float
    bias_bounds: tuple[float, float]
    bias_prior_mean: float
    bias_prior_sigma: float
    seed: int
    burn_in: int
    observation_keys: tuple[str, ...] = ()
    bias_labels: tuple[tuple[str, str], ...] = ()


def mcmc_sample_combined_grouped_uncertainty(
    problem: InversionProblem,
    metadata_key: str,
    observation_key: str | Sequence[str],
    initial_corrections_j_mol: Sequence[float],
    initial_scales: Mapping[str, float],
    initial_biases: Mapping[str, float] | Mapping[tuple[str, str], float],
    n_samples: int = 5000,
    burn_in: int = 500,
    proposal_sigma_j_mol: Sequence[float] | None = None,
    proposal_sigma_scale: float = 0.1,
    proposal_sigma_bias: float = 0.01,
    scale_bounds: tuple[float, float] = (0.0, 5.0),
    scale_prior_mean: float = 1.0,
    scale_prior_sigma: float = 1.0,
    bias_bounds: tuple[float, float] = (-1.0, 1.0),
    bias_prior_mean: float = 0.0,
    bias_prior_sigma: float = 0.1,
    seed: int = 0,
) -> CombinedGroupedUncertaintyMCMCResult:
    """Sample model-error scales and observation biases jointly by MH."""
    if problem.theoretical_covariance_scale_metadata_key:
        raise ValueError("combined grouped sampling expects an unprofiled problem")
    observation_keys = (
        (observation_key,)
        if isinstance(observation_key, str)
        else tuple(observation_key)
    )
    if (
        not observation_keys
        or len(set(observation_keys)) != len(observation_keys)
        or any(not isinstance(key, str) or not key for key in observation_keys)
    ):
        raise ValueError("observation_key(s) must be unique non-empty strings")
    groups = sorted(
        {experiment.metadata.get(metadata_key) for experiment in problem.experiments}
    )
    if (
        not metadata_key
        or not observation_key
        or not groups
        or any(not isinstance(group, str) or not group for group in groups)
    ):
        raise ValueError(
            "metadata_key, observation_key, and group labels must be non-empty"
        )
    if set(initial_scales) != set(groups):
        raise ValueError("initial scales must match metadata groups")
    if isinstance(observation_key, str):
        if set(initial_biases) != set(groups):
            raise ValueError("initial biases must match metadata groups")
        bias_labels = tuple(groups)
        bias_values = [initial_biases[group] for group in groups]
    else:
        expected_labels = tuple(
            (key, group) for key in observation_keys for group in groups
        )
        if set(initial_biases) != set(expected_labels):
            raise ValueError(
                "initial biases must match observation-key and group pairs"
            )
        bias_labels = expected_labels
        bias_values = [initial_biases[label] for label in expected_labels]
    scale_lower, scale_upper = scale_bounds
    bias_lower, bias_upper = bias_bounds
    if (
        not np.isfinite(scale_lower)
        or not np.isfinite(scale_upper)
        or scale_lower < 0.0
        or scale_lower >= scale_upper
    ):
        raise ValueError("scale_bounds must be finite, non-negative, and increasing")
    if (
        not np.isfinite(bias_lower)
        or not np.isfinite(bias_upper)
        or bias_lower >= bias_upper
    ):
        raise ValueError("bias_bounds must be finite and increasing")
    if (
        not np.isfinite(scale_prior_mean)
        or scale_prior_mean < 0.0
        or not np.isfinite(scale_prior_sigma)
        or scale_prior_sigma <= 0.0
    ):
        raise ValueError(
            "scale prior must be finite with non-negative mean and positive sigma"
        )
    if (
        not np.isfinite(bias_prior_mean)
        or not np.isfinite(bias_prior_sigma)
        or bias_prior_sigma <= 0.0
    ):
        raise ValueError("bias prior must be finite with positive sigma")
    if (
        not np.isfinite(proposal_sigma_scale)
        or proposal_sigma_scale <= 0.0
        or not np.isfinite(proposal_sigma_bias)
        or proposal_sigma_bias <= 0.0
    ):
        raise ValueError("nuisance proposal sigmas must be finite and positive")
    if len(initial_corrections_j_mol) != len(problem.parameters):
        raise ValueError("initial corrections do not match parameter count")
    for experiment in problem.experiments:
        available = {
            *(f"fraction:{phase}" for phase in (experiment.phase_fractions or {})),
            *(
                f"composition:{phase}:{component}"
                for phase, values in (experiment.phase_compositions or {}).items()
                for component in values
            ),
        }
        if any(key not in available for key in observation_keys):
            raise ValueError(
                f"Experiment {experiment.id!r} lacks a requested observation"
            )
    if proposal_sigma_j_mol is None:
        correction_sigmas = np.asarray(
            [
                max(abs(parameter.prior_sigma) * 0.1, 1.0e-6)
                for parameter in problem.parameters
            ]
        )
    else:
        if len(proposal_sigma_j_mol) != len(problem.parameters):
            raise ValueError("proposal_sigma_j_mol does not match parameter count")
        correction_sigmas = np.asarray(proposal_sigma_j_mol, dtype=float)
    if not np.all(np.isfinite(correction_sigmas)) or np.any(correction_sigmas <= 0.0):
        raise ValueError("proposal_sigma_j_mol must be finite and positive")
    initial = np.asarray(
        [
            *initial_corrections_j_mol,
            *(initial_scales[group] for group in groups),
            *bias_values,
        ],
        dtype=float,
    )
    bounds = np.asarray(
        [parameter.bounds() for parameter in problem.parameters]
        + [scale_bounds for _ in groups]
        + [bias_bounds for _ in bias_labels],
        dtype=float,
    )
    if (
        not np.all(np.isfinite(initial))
        or np.any(initial < bounds[:, 0])
        or np.any(initial > bounds[:, 1])
    ):
        raise ValueError(
            "initial combined nuisance point must be finite and within bounds"
        )
    proposal_sigmas = np.asarray(
        [
            *correction_sigmas,
            *([proposal_sigma_scale] * len(groups)),
            *([proposal_sigma_bias] * len(bias_labels)),
        ]
    )

    def objective(values: np.ndarray) -> float:
        n_corrections = len(problem.parameters)
        scales = {
            group: float(values[n_corrections + index])
            for index, group in enumerate(groups)
        }
        biases = {
            label: float(values[n_corrections + len(groups) + index])
            for index, label in enumerate(bias_labels)
        }
        experiments = []
        for experiment in problem.experiments:
            observation_bias = dict(experiment.observation_bias or {})
            group = experiment.metadata[metadata_key]
            for key in observation_keys:
                label = group if isinstance(observation_key, str) else (key, group)
                observation_bias[key] = biases[label]
            experiments.append(replace(experiment, observation_bias=observation_bias))
        nuisance_problem = InversionProblem(
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
            theoretical_covariance_scale_metadata_key=metadata_key,
            theoretical_covariance_scales=scales,
        )
        penalty = sum(
            0.5 * ((value - scale_prior_mean) / scale_prior_sigma) ** 2
            for value in scales.values()
        )
        penalty += sum(
            0.5 * ((value - bias_prior_mean) / bias_prior_sigma) ** 2
            for value in biases.values()
        )
        return float(nuisance_problem.objective(values[:n_corrections]) + penalty)

    if n_samples <= 0 or burn_in < 0 or burn_in >= n_samples:
        raise ValueError("n_samples must be positive and exceed burn_in")
    rng = np.random.default_rng(seed)
    current = initial.copy()
    current_objective = objective(current)
    samples = np.empty((n_samples, len(initial)))
    log_posterior = np.empty(n_samples)
    accepted = 0
    for step in range(n_samples):
        proposal = current + proposal_sigmas * rng.normal(size=len(initial))
        proposal_objective = float("inf")
        if np.all(proposal >= bounds[:, 0]) and np.all(proposal <= bounds[:, 1]):
            proposal_objective = objective(proposal)
        if np.log(rng.uniform()) < current_objective - proposal_objective:
            current, current_objective = proposal, proposal_objective
            accepted += 1
        samples[step] = current
        log_posterior[step] = -current_objective
    kept = samples[burn_in:]
    names = (
        tuple(parameter.name for parameter in problem.parameters)
        + tuple(f"scale_C_T[{group}]" for group in groups)
        + tuple(
            f"bias[{label}]"
            if not isinstance(observation_key, str)
            else f"bias[{label}]"
            for label in bias_labels
        )
    )
    units = (
        tuple(getattr(parameter, "unit", "J/mol") for parameter in problem.parameters)
        + tuple("dimensionless" for _ in groups)
        + tuple("native observation unit" for _ in groups)
    )
    return CombinedGroupedUncertaintyMCMCResult(
        names,
        units,
        kept,
        log_posterior[burn_in:],
        accepted / n_samples,
        {name: float(np.mean(kept[:, index])) for index, name in enumerate(names)},
        {name: float(np.std(kept[:, index])) for index, name in enumerate(names)},
        {
            name: (
                float(np.quantile(kept[:, index], 0.025)),
                float(np.quantile(kept[:, index], 0.975)),
            )
            for index, name in enumerate(names)
        },
        metadata_key,
        observation_keys[0],
        tuple(groups),
        (float(scale_lower), float(scale_upper)),
        float(scale_prior_mean),
        float(scale_prior_sigma),
        (float(bias_lower), float(bias_upper)),
        float(bias_prior_mean),
        float(bias_prior_sigma),
        seed,
        burn_in,
        observation_keys,
        tuple((key, group) for key in observation_keys for group in groups),
    )


def mcmc_sample_group_observation_biases(
    problem: InversionProblem,
    metadata_key: str,
    observation_key: str,
    initial_corrections_j_mol: Sequence[float],
    initial_biases: Mapping[str, float],
    n_samples: int = 5000,
    burn_in: int = 500,
    proposal_sigma_j_mol: Sequence[float] | None = None,
    proposal_sigma_bias: float = 0.01,
    bias_bounds: tuple[float, float] = (-1.0, 1.0),
    bias_prior_mean: float = 0.0,
    bias_prior_sigma: float = 0.1,
    seed: int = 0,
) -> GroupedObservationBiasMCMCResult:
    """Sample corrections and one additive observation bias per metadata group."""
    if not metadata_key or not observation_key:
        raise ValueError("metadata_key and observation_key must be non-empty")
    groups = sorted(
        {experiment.metadata.get(metadata_key) for experiment in problem.experiments}
    )
    return mcmc_sample_group_observation_biases_multi_key(
        problem,
        metadata_key,
        (observation_key,),
        initial_corrections_j_mol,
        {(observation_key, group): value for group, value in initial_biases.items()},
        n_samples=n_samples,
        burn_in=burn_in,
        proposal_sigma_j_mol=proposal_sigma_j_mol,
        proposal_sigma_bias=proposal_sigma_bias,
        bias_bounds=bias_bounds,
        bias_prior_mean=bias_prior_mean,
        bias_prior_sigma=bias_prior_sigma,
        seed=seed,
    )


def mcmc_sample_group_observation_biases_multi_key(
    problem: InversionProblem,
    metadata_key: str,
    observation_keys: Sequence[str],
    initial_corrections_j_mol: Sequence[float],
    initial_biases: Mapping[tuple[str, str], float],
    n_samples: int = 5000,
    burn_in: int = 500,
    proposal_sigma_j_mol: Sequence[float] | None = None,
    proposal_sigma_bias: float = 0.01,
    bias_bounds: tuple[float, float] = (-1.0, 1.0),
    bias_prior_mean: float = 0.0,
    bias_prior_sigma: float = 0.1,
    seed: int = 0,
) -> GroupedObservationBiasMCMCResult:
    """Sample corrections and biases for several observation keys jointly."""
    keys = tuple(observation_keys)
    if not keys or len(set(keys)) != len(keys) or any(not key for key in keys):
        raise ValueError("observation_keys must be unique non-empty strings")
    groups = sorted(
        {experiment.metadata.get(metadata_key) for experiment in problem.experiments}
    )
    if not groups or any(not isinstance(group, str) or not group for group in groups):
        raise ValueError("all experiments need non-empty string bias groups")
    labels = tuple((key, group) for key in keys for group in groups)
    if set(initial_biases) != set(labels):
        raise ValueError(
            "initial_biases must match observation-key and metadata groups"
        )
    for experiment in problem.experiments:
        available = {
            *(f"fraction:{phase}" for phase in (experiment.phase_fractions or {})),
            *(
                f"composition:{phase}:{component}"
                for phase, values in (experiment.phase_compositions or {}).items()
                for component in values
            ),
        }
        if any(key not in available for key in keys):
            raise ValueError(
                f"Experiment {experiment.id!r} lacks a requested observation key"
            )
    if len(initial_corrections_j_mol) != len(problem.parameters):
        raise ValueError("initial corrections do not match parameter count")
    lower, upper = bias_bounds
    if not np.isfinite(lower) or not np.isfinite(upper) or lower >= upper:
        raise ValueError("bias_bounds must be finite and increasing")
    if (
        not np.isfinite(bias_prior_mean)
        or not np.isfinite(bias_prior_sigma)
        or bias_prior_sigma <= 0.0
    ):
        raise ValueError("bias prior must have finite mean and positive sigma")
    if n_samples <= 0 or burn_in < 0 or burn_in >= n_samples:
        raise ValueError("n_samples must be positive and exceed burn_in")
    if not np.isfinite(proposal_sigma_bias) or proposal_sigma_bias <= 0.0:
        raise ValueError("proposal_sigma_bias must be finite and positive")
    if proposal_sigma_j_mol is None:
        correction_sigmas = np.asarray(
            [
                max(abs(parameter.prior_sigma) * 0.1, 1.0e-6)
                for parameter in problem.parameters
            ]
        )
    else:
        if len(proposal_sigma_j_mol) != len(problem.parameters):
            raise ValueError("proposal_sigma_j_mol does not match parameter count")
        correction_sigmas = np.asarray(proposal_sigma_j_mol, dtype=float)
    if not np.all(np.isfinite(correction_sigmas)) or np.any(correction_sigmas <= 0.0):
        raise ValueError("proposal_sigma_j_mol must be finite and positive")
    initial = np.asarray(
        [*initial_corrections_j_mol, *(initial_biases[label] for label in labels)],
        dtype=float,
    )
    bounds = np.asarray(
        [parameter.bounds() for parameter in problem.parameters]
        + [bias_bounds for _ in labels],
        dtype=float,
    )
    if (
        not np.all(np.isfinite(initial))
        or np.any(initial < bounds[:, 0])
        or np.any(initial > bounds[:, 1])
    ):
        raise ValueError(
            "initial multi-key bias point must be finite and within bounds"
        )
    proposal_sigmas = np.asarray(
        [*correction_sigmas, *([proposal_sigma_bias] * len(labels))]
    )

    def objective(values: np.ndarray) -> float:
        biases = {
            label: float(values[len(problem.parameters) + index])
            for index, label in enumerate(labels)
        }
        experiments = []
        for experiment in problem.experiments:
            observation_bias = dict(experiment.observation_bias or {})
            group = experiment.metadata[metadata_key]
            for key in keys:
                observation_bias[key] = biases[(key, group)]
            experiments.append(replace(experiment, observation_bias=observation_bias))
        biased_problem = InversionProblem(
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
        penalty = sum(
            0.5 * ((value - bias_prior_mean) / bias_prior_sigma) ** 2
            for value in biases.values()
        )
        return float(
            biased_problem.objective(values[: len(problem.parameters)]) + penalty
        )

    rng = np.random.default_rng(seed)
    current = initial.copy()
    current_objective = objective(current)
    samples = np.empty((n_samples, len(initial)))
    log_posterior = np.empty(n_samples)
    accepted = 0
    for step in range(n_samples):
        proposal = current + proposal_sigmas * rng.normal(size=len(initial))
        proposal_objective = float("inf")
        if np.all(proposal >= bounds[:, 0]) and np.all(proposal <= bounds[:, 1]):
            proposal_objective = objective(proposal)
        if np.log(rng.uniform()) < current_objective - proposal_objective:
            current, current_objective = proposal, proposal_objective
            accepted += 1
        samples[step] = current
        log_posterior[step] = -current_objective
    kept = samples[burn_in:]
    names = tuple(parameter.name for parameter in problem.parameters) + tuple(
        f"bias[{key}|{group}]" for key, group in labels
    )
    units = tuple(
        getattr(parameter, "unit", "J/mol") for parameter in problem.parameters
    ) + tuple("native observation unit" for _ in labels)
    result = GroupedObservationBiasMCMCResult(
        names,
        units,
        kept,
        log_posterior[burn_in:],
        accepted / n_samples,
        {name: float(np.mean(kept[:, index])) for index, name in enumerate(names)},
        {name: float(np.std(kept[:, index])) for index, name in enumerate(names)},
        {
            name: (
                float(np.quantile(kept[:, index], 0.025)),
                float(np.quantile(kept[:, index], 0.975)),
            )
            for index, name in enumerate(names)
        },
        metadata_key,
        keys[0],
        tuple(f"{key}|{group}" for key, group in labels),
        float(bias_prior_mean),
        float(bias_prior_sigma),
        (float(lower), float(upper)),
        seed,
        burn_in,
        keys,
    )
    return result
    if not groups or any(not isinstance(group, str) or not group for group in groups):
        raise ValueError("all experiments need non-empty string bias groups")
    if set(initial_biases) != set(groups):
        raise ValueError("initial_biases must match metadata groups")
    lower, upper = bias_bounds
    if not np.isfinite(lower) or not np.isfinite(upper) or lower >= upper:
        raise ValueError("bias_bounds must be finite and increasing")
    if (
        not np.isfinite(bias_prior_mean)
        or not np.isfinite(bias_prior_sigma)
        or bias_prior_sigma <= 0.0
    ):
        raise ValueError("bias prior must have finite mean and positive sigma")
    if not np.isfinite(proposal_sigma_bias) or proposal_sigma_bias <= 0.0:
        raise ValueError("proposal_sigma_bias must be finite and positive")
    if len(initial_corrections_j_mol) != len(problem.parameters):
        raise ValueError("initial corrections do not match parameter count")
    for experiment in problem.experiments:
        keys = {
            *(f"fraction:{phase}" for phase in (experiment.phase_fractions or {})),
            *(
                f"composition:{phase}:{component}"
                for phase, values in (experiment.phase_compositions or {}).items()
                for component in values
            ),
        }
        if observation_key not in keys:
            raise ValueError(
                f"Experiment {experiment.id!r} lacks observation {observation_key!r}"
            )
    initial = np.asarray(
        [*initial_corrections_j_mol, *(initial_biases[group] for group in groups)],
        dtype=float,
    )
    bounds = np.asarray(
        [parameter.bounds() for parameter in problem.parameters]
        + [bias_bounds for _ in groups],
        dtype=float,
    )
    if (
        not np.all(np.isfinite(initial))
        or np.any(initial < bounds[:, 0])
        or np.any(initial > bounds[:, 1])
    ):
        raise ValueError("initial bias point must be finite and within bounds")
    if proposal_sigma_j_mol is None:
        correction_sigmas = np.asarray(
            [
                max(abs(parameter.prior_sigma) * 0.1, 1.0e-6)
                for parameter in problem.parameters
            ]
        )
    else:
        if len(proposal_sigma_j_mol) != len(problem.parameters):
            raise ValueError("proposal_sigma_j_mol does not match parameter count")
        correction_sigmas = np.asarray(proposal_sigma_j_mol, dtype=float)
    if not np.all(np.isfinite(correction_sigmas)) or np.any(correction_sigmas <= 0.0):
        raise ValueError("proposal_sigma_j_mol must be finite and positive")
    proposal_sigmas = np.asarray(
        [*correction_sigmas, *([proposal_sigma_bias] * len(groups))]
    )

    def objective(values: np.ndarray) -> float:
        biases = {
            group: float(values[len(problem.parameters) + index])
            for index, group in enumerate(groups)
        }
        experiments = []
        for experiment in problem.experiments:
            observation_bias = dict(experiment.observation_bias or {})
            observation_bias[observation_key] = biases[
                experiment.metadata[metadata_key]
            ]
            experiments.append(replace(experiment, observation_bias=observation_bias))
        biased_problem = InversionProblem(
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
        penalty = sum(
            0.5 * ((value - bias_prior_mean) / bias_prior_sigma) ** 2
            for value in biases.values()
        )
        return float(
            biased_problem.objective(values[: len(problem.parameters)]) + penalty
        )

    if n_samples <= 0 or burn_in < 0 or burn_in >= n_samples:
        raise ValueError("n_samples must be positive and exceed burn_in")
    rng = np.random.default_rng(seed)
    current = initial.copy()
    current_objective = objective(current)
    samples = np.empty((n_samples, len(initial)))
    log_posterior = np.empty(n_samples)
    accepted = 0
    for step in range(n_samples):
        proposal = current + proposal_sigmas * rng.normal(size=len(initial))
        proposal_objective = float("inf")
        if np.all(proposal >= bounds[:, 0]) and np.all(proposal <= bounds[:, 1]):
            proposal_objective = objective(proposal)
        if np.log(rng.uniform()) < current_objective - proposal_objective:
            current, current_objective = proposal, proposal_objective
            accepted += 1
        samples[step] = current
        log_posterior[step] = -current_objective
    kept = samples[burn_in:]
    names = tuple(parameter.name for parameter in problem.parameters) + tuple(
        f"bias[{group}]" for group in groups
    )
    units = tuple(
        getattr(parameter, "unit", "J/mol") for parameter in problem.parameters
    ) + tuple("native observation unit" for _ in groups)
    return GroupedObservationBiasMCMCResult(
        names,
        units,
        kept,
        log_posterior[burn_in:],
        accepted / n_samples,
        {name: float(np.mean(kept[:, index])) for index, name in enumerate(names)},
        {name: float(np.std(kept[:, index])) for index, name in enumerate(names)},
        {
            name: (
                float(np.quantile(kept[:, index], 0.025)),
                float(np.quantile(kept[:, index], 0.975)),
            )
            for index, name in enumerate(names)
        },
        metadata_key,
        observation_key,
        tuple(groups),
        float(bias_prior_mean),
        float(bias_prior_sigma),
        (float(lower), float(upper)),
        seed,
        burn_in,
        observation_keys,
    )


def mcmc_sample_theoretical_covariance_group_scales(
    problem: InversionProblem,
    metadata_key: str,
    initial_corrections_j_mol: Sequence[float],
    initial_scales: Mapping[str, float],
    n_samples: int = 5000,
    burn_in: int = 500,
    proposal_sigma_j_mol: Sequence[float] | None = None,
    proposal_sigma_scale: float = 0.1,
    scale_bounds: tuple[float, float] = (0.0, 5.0),
    scale_prior_mean: float = 1.0,
    scale_prior_sigma: float = 1.0,
    seed: int = 0,
) -> GroupedHierarchicalMCMCResult:
    """Sample corrections and one model-error scale per metadata group jointly."""
    if not metadata_key:
        raise ValueError("metadata_key must be non-empty")
    if problem.theoretical_covariance_scale_metadata_key:
        raise ValueError("group scale sampling expects an unprofiled problem")
    groups = sorted(
        {
            experiment.metadata.get(metadata_key)
            for experiment in problem.experiments
            if experiment.theoretical_covariance is not None
        }
    )
    if not groups or any(not isinstance(group, str) or not group for group in groups):
        raise ValueError(
            "theoretical-covariance experiments need non-empty string groups"
        )
    if set(initial_scales) != set(groups):
        raise ValueError("initial_scales must match theoretical-covariance groups")
    if n_samples <= 0 or burn_in < 0 or burn_in >= n_samples:
        raise ValueError("n_samples must be positive and exceed burn_in")
    lower, upper = scale_bounds
    if (
        not np.isfinite(lower)
        or not np.isfinite(upper)
        or lower < 0.0
        or lower >= upper
    ):
        raise ValueError("scale_bounds must be finite, non-negative, and increasing")
    if not np.isfinite(scale_prior_mean) or scale_prior_mean < 0.0:
        raise ValueError("scale_prior_mean must be finite and non-negative")
    if not np.isfinite(scale_prior_sigma) or scale_prior_sigma <= 0.0:
        raise ValueError("scale_prior_sigma must be finite and positive")
    if not np.isfinite(proposal_sigma_scale) or proposal_sigma_scale <= 0.0:
        raise ValueError("proposal_sigma_scale must be finite and positive")
    if len(initial_corrections_j_mol) != len(problem.parameters):
        raise ValueError("initial corrections do not match parameter count")
    initial = np.asarray(
        [*initial_corrections_j_mol, *(initial_scales[group] for group in groups)],
        dtype=float,
    )
    bounds = np.asarray(
        [parameter.bounds() for parameter in problem.parameters]
        + [scale_bounds for _ in groups],
        dtype=float,
    )
    if (
        not np.all(np.isfinite(initial))
        or np.any(initial < bounds[:, 0])
        or np.any(initial > bounds[:, 1])
    ):
        raise ValueError(
            "initial grouped hierarchical point must be finite and within bounds"
        )
    if proposal_sigma_j_mol is None:
        correction_sigmas = np.asarray(
            [
                max(abs(parameter.prior_sigma) * 0.1, 1.0e-6)
                for parameter in problem.parameters
            ]
        )
    else:
        if len(proposal_sigma_j_mol) != len(problem.parameters):
            raise ValueError("proposal_sigma_j_mol does not match parameter count")
        correction_sigmas = np.asarray(proposal_sigma_j_mol, dtype=float)
    if not np.all(np.isfinite(correction_sigmas)) or np.any(correction_sigmas <= 0.0):
        raise ValueError("proposal_sigma_j_mol must be finite and positive")
    proposal_sigmas = np.asarray(
        [*correction_sigmas, *([proposal_sigma_scale] * len(groups))]
    )

    def objective(values: np.ndarray) -> float:
        scales = {
            group: float(values[len(problem.parameters) + index])
            for index, group in enumerate(groups)
        }
        scaled_problem = InversionProblem(
            problem.forward_model,
            problem.experiments,
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
            theoretical_covariance_scale_metadata_key=metadata_key,
            theoretical_covariance_scales=scales,
        )
        penalty = sum(
            0.5 * ((value - scale_prior_mean) / scale_prior_sigma) ** 2
            for value in scales.values()
        )
        return float(
            scaled_problem.objective(values[: len(problem.parameters)]) + penalty
        )

    rng = np.random.default_rng(seed)
    current = initial.copy()
    current_objective = objective(current)
    samples = np.empty((n_samples, len(initial)))
    log_posterior = np.empty(n_samples)
    accepted = 0
    for step in range(n_samples):
        proposal = current + proposal_sigmas * rng.normal(size=len(initial))
        proposal_objective = float("inf")
        if np.all(proposal >= bounds[:, 0]) and np.all(proposal <= bounds[:, 1]):
            proposal_objective = objective(proposal)
        if np.log(rng.uniform()) < current_objective - proposal_objective:
            current, current_objective = proposal, proposal_objective
            accepted += 1
        samples[step] = current
        log_posterior[step] = -current_objective
    kept = samples[burn_in:]
    names = tuple(parameter.name for parameter in problem.parameters) + tuple(
        f"scale_C_T[{group}]" for group in groups
    )
    units = tuple(
        getattr(parameter, "unit", "J/mol") for parameter in problem.parameters
    ) + tuple("dimensionless" for _ in groups)
    return GroupedHierarchicalMCMCResult(
        names,
        units,
        kept,
        log_posterior[burn_in:],
        accepted / n_samples,
        {name: float(np.mean(kept[:, index])) for index, name in enumerate(names)},
        {name: float(np.std(kept[:, index])) for index, name in enumerate(names)},
        {
            name: (
                float(np.quantile(kept[:, index], 0.025)),
                float(np.quantile(kept[:, index], 0.975)),
            )
            for index, name in enumerate(names)
        },
        metadata_key,
        tuple(groups),
        {group: float(scale_prior_mean) for group in groups},
        {group: float(scale_prior_sigma) for group in groups},
        (float(lower), float(upper)),
        seed,
        burn_in,
    )


@dataclass(frozen=True)
class TheoreticalScaleProfileResult:
    """Profile-MAP estimate of a global theoretical-error scale."""

    scale: float
    objective: float
    corrections_j_mol: dict[str, float]
    scale_prior_penalty: float
    records: tuple[dict[str, object], ...]


@dataclass(frozen=True)
class TheoreticalGroupScaleProfileResult:
    """Profile-MAP estimates of metadata-group theoretical-error scales."""

    scales: dict[str, float]
    objective: float
    corrections_j_mol: dict[str, float]
    scale_prior_penalty: float
    metadata_key: str
    records: tuple[dict[str, object], ...]


def mcmc_sample_theoretical_covariance_scale(
    problem: InversionProblem,
    initial_corrections_j_mol: Sequence[float],
    initial_scale: float = 1.0,
    n_samples: int = 5000,
    burn_in: int = 500,
    proposal_sigma_j_mol: Sequence[float] | None = None,
    proposal_sigma_scale: float = 0.1,
    scale_bounds: tuple[float, float] = (0.0, 5.0),
    scale_prior_mean: float = 1.0,
    scale_prior_sigma: float = 1.0,
    seed: int = 0,
) -> HierarchicalMCMCResult:
    """Sample corrections and one global ``scale_C_T`` jointly by MH."""
    if len(initial_corrections_j_mol) != len(problem.parameters):
        raise ValueError("initial corrections do not match parameter count")
    if problem.theoretical_covariance_scale_metadata_key:
        raise ValueError(
            "global scale sampling cannot be combined with metadata-specific scales"
        )
    if n_samples <= 0 or burn_in < 0 or burn_in >= n_samples:
        raise ValueError("n_samples must be positive and exceed burn_in")
    lower_scale, upper_scale = scale_bounds
    if (
        not np.isfinite(lower_scale)
        or not np.isfinite(upper_scale)
        or lower_scale < 0.0
        or lower_scale >= upper_scale
    ):
        raise ValueError("scale_bounds must be finite, non-negative, and increasing")
    if not np.isfinite(scale_prior_mean) or scale_prior_mean < 0.0:
        raise ValueError("scale_prior_mean must be finite and non-negative")
    if not np.isfinite(scale_prior_sigma) or scale_prior_sigma <= 0.0:
        raise ValueError("scale_prior_sigma must be finite and positive")
    if not np.isfinite(proposal_sigma_scale) or proposal_sigma_scale <= 0.0:
        raise ValueError("proposal_sigma_scale must be finite and positive")
    initial = np.asarray([*initial_corrections_j_mol, initial_scale], dtype=float)
    bounds = np.asarray(
        [*[(parameter.bounds()) for parameter in problem.parameters], scale_bounds],
        dtype=float,
    )
    if (
        not np.all(np.isfinite(initial))
        or np.any(initial < bounds[:, 0])
        or np.any(initial > bounds[:, 1])
    ):
        raise ValueError("initial hierarchical point must be finite and within bounds")
    if proposal_sigma_j_mol is None:
        correction_sigmas = np.asarray(
            [
                max(abs(parameter.prior_sigma) * 0.1, 1.0e-6)
                for parameter in problem.parameters
            ]
        )
    else:
        if len(proposal_sigma_j_mol) != len(problem.parameters):
            raise ValueError("proposal_sigma_j_mol does not match parameter count")
        correction_sigmas = np.asarray(proposal_sigma_j_mol, dtype=float)
        if not np.all(np.isfinite(correction_sigmas)) or np.any(
            correction_sigmas <= 0.0
        ):
            raise ValueError("proposal_sigma_j_mol must be finite and positive")
    proposal_sigmas = np.asarray([*correction_sigmas, proposal_sigma_scale])

    def objective(values: np.ndarray) -> float:
        scale = float(values[-1])
        scaled_problem = InversionProblem(
            problem.forward_model,
            problem.experiments,
            problem.parameters,
            brackets=problem.brackets,
            missing_phase_penalty=problem.missing_phase_penalty,
            extra_phase_penalty=problem.extra_phase_penalty,
            marginalize_condition_uncertainty=problem.marginalize_condition_uncertainty,
            condition_quadrature_order=problem.condition_quadrature_order,
            objective_criterion=problem.objective_criterion,
            misfit_form=problem.misfit_form,
            prior_covariance=problem.prior_covariance,
            theoretical_covariance_scale=scale,
            theoretical_covariance_scale_metadata_key=problem.theoretical_covariance_scale_metadata_key,
            theoretical_covariance_scales=problem.theoretical_covariance_scales,
        )
        correction_objective = scaled_problem.objective(values[:-1])
        scale_penalty = 0.5 * ((scale - scale_prior_mean) / scale_prior_sigma) ** 2
        return float(correction_objective + scale_penalty)

    rng = np.random.default_rng(seed)
    current = initial.copy()
    current_objective = objective(current)
    samples = np.empty((n_samples, len(initial)))
    log_posterior = np.empty(n_samples)
    accepted = 0
    for step in range(n_samples):
        proposal = current + proposal_sigmas * rng.normal(size=len(initial))
        proposal_objective = float("inf")
        if np.all(proposal >= bounds[:, 0]) and np.all(proposal <= bounds[:, 1]):
            proposal_objective = objective(proposal)
        if np.log(rng.uniform()) < current_objective - proposal_objective:
            current, current_objective = proposal, proposal_objective
            accepted += 1
        samples[step] = current
        log_posterior[step] = -current_objective

    kept = samples[burn_in:]
    kept_log_posterior = log_posterior[burn_in:]
    names = tuple(parameter.name for parameter in problem.parameters) + ("scale_C_T",)
    units = tuple(
        getattr(parameter, "unit", "J/mol") for parameter in problem.parameters
    ) + ("dimensionless",)
    return HierarchicalMCMCResult(
        names,
        units,
        kept,
        kept_log_posterior,
        accepted / n_samples,
        {name: float(np.mean(kept[:, index])) for index, name in enumerate(names)},
        {name: float(np.std(kept[:, index])) for index, name in enumerate(names)},
        {
            name: (
                float(np.quantile(kept[:, index], 0.025)),
                float(np.quantile(kept[:, index], 0.975)),
            )
            for index, name in enumerate(names)
        },
        "scale_C_T",
        float(scale_prior_mean),
        float(scale_prior_sigma),
        (float(lower_scale), float(upper_scale)),
        seed,
        burn_in,
    )


def profile_theoretical_covariance_scale(
    problem: InversionProblem,
    scale_bounds: tuple[float, float] = (0.0, 5.0),
    scale_prior_mean: float = 1.0,
    scale_prior_sigma: float = 1.0,
    scale_iterations: int = 25,
    optimize_starts: int = 4,
    optimize_iterations: int = 60,
) -> TheoreticalScaleProfileResult:
    """Estimate a global model-error scale by nested bounded MAP refits.

    The profiled objective is the base inversion objective with
    ``C_D + scale**2 * C_T`` plus a Gaussian prior on the non-negative scale.
    This is a MAP nuisance estimate, not a posterior sample of the scale.
    """
    lower, upper = scale_bounds
    if (
        not np.isfinite(lower)
        or not np.isfinite(upper)
        or lower < 0.0
        or lower >= upper
    ):
        raise ValueError("scale_bounds must be finite, non-negative, and increasing")
    if not np.isfinite(scale_prior_mean) or scale_prior_mean < 0.0:
        raise ValueError("scale_prior_mean must be finite and non-negative")
    if not np.isfinite(scale_prior_sigma) or scale_prior_sigma <= 0.0:
        raise ValueError("scale_prior_sigma must be finite and positive")
    if scale_iterations < 1:
        raise ValueError("scale_iterations must be positive")
    if not any(
        experiment.theoretical_covariance is not None
        for experiment in problem.experiments
    ):
        raise ValueError("theoretical covariance is required for scale profiling")

    def fit_at(scale: float) -> tuple[MAPResult, float]:
        profiled = InversionProblem(
            problem.forward_model,
            problem.experiments,
            problem.parameters,
            brackets=problem.brackets,
            missing_phase_penalty=problem.missing_phase_penalty,
            extra_phase_penalty=problem.extra_phase_penalty,
            marginalize_condition_uncertainty=problem.marginalize_condition_uncertainty,
            condition_quadrature_order=problem.condition_quadrature_order,
            objective_criterion=problem.objective_criterion,
            misfit_form=problem.misfit_form,
            prior_covariance=problem.prior_covariance,
            theoretical_covariance_scale=scale,
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
        penalty = 0.5 * ((scale - scale_prior_mean) / scale_prior_sigma) ** 2
        return fit, float(fit.objective + penalty)

    golden = (1.0 + np.sqrt(5.0)) / 2.0
    left, right = lower, upper
    records: list[dict[str, object]] = []
    best_fit: MAPResult | None = None
    best_scale = left
    best_objective = float("inf")
    for _ in range(scale_iterations):
        first = right - (right - left) / golden
        second = left + (right - left) / golden
        first_fit, first_objective = fit_at(first)
        second_fit, second_objective = fit_at(second)
        records.extend(
            [
                {"scale": first, "objective": first_objective},
                {"scale": second, "objective": second_objective},
            ]
        )
        for scale, fit, objective in (
            (first, first_fit, first_objective),
            (second, second_fit, second_objective),
        ):
            if objective < best_objective:
                best_scale, best_fit, best_objective = scale, fit, objective
        if first_objective < second_objective:
            right = second
        else:
            left = first

    if best_fit is None:  # pragma: no cover - bounded loop always evaluates twice
        raise RuntimeError("theoretical scale profile produced no fit")
    penalty = 0.5 * ((best_scale - scale_prior_mean) / scale_prior_sigma) ** 2
    return TheoreticalScaleProfileResult(
        best_scale,
        best_objective,
        best_fit.corrections_j_mol,
        float(penalty),
        tuple(records),
    )


def profile_theoretical_covariance_group_scales(
    problem: InversionProblem,
    metadata_key: str,
    scale_bounds: tuple[float, float] = (0.0, 5.0),
    scale_prior_mean: float = 1.0,
    scale_prior_sigma: float = 1.0,
    scale_iterations: int = 8,
    optimize_starts: int = 4,
    optimize_iterations: int = 60,
) -> TheoreticalGroupScaleProfileResult:
    """Estimate one theoretical-error scale per metadata group by MAP profiling."""
    lower, upper = scale_bounds
    if (
        not metadata_key
        or lower < 0.0
        or not np.isfinite(lower)
        or not np.isfinite(upper)
        or lower >= upper
    ):
        raise ValueError("metadata key and scale bounds are invalid")
    if scale_prior_mean < 0.0 or not np.isfinite(scale_prior_mean):
        raise ValueError("scale_prior_mean must be finite and non-negative")
    if scale_prior_sigma <= 0.0 or not np.isfinite(scale_prior_sigma):
        raise ValueError("scale_prior_sigma must be finite and positive")
    if scale_iterations < 1:
        raise ValueError("scale_iterations must be positive")
    groups = {
        experiment.metadata.get(metadata_key)
        for experiment in problem.experiments
        if experiment.theoretical_covariance is not None
    }
    if not groups or any(not isinstance(group, str) or not group for group in groups):
        raise ValueError(
            "theoretical-covariance experiments need non-empty string groups"
        )
    scales = {group: float(np.clip(scale_prior_mean, lower, upper)) for group in groups}
    golden = (1.0 + np.sqrt(5.0)) / 2.0
    records: list[dict[str, object]] = []
    best_fit: MAPResult | None = None
    best_objective = float("inf")

    def fit_at(candidate_scales: dict[str, float]) -> tuple[MAPResult, float]:
        profiled = InversionProblem(
            problem.forward_model,
            problem.experiments,
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
            theoretical_covariance_scale_metadata_key=metadata_key,
            theoretical_covariance_scales=candidate_scales,
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
            0.5 * ((value - scale_prior_mean) / scale_prior_sigma) ** 2
            for value in candidate_scales.values()
        )
        return fit, float(fit.objective + penalty)

    for _ in range(scale_iterations):
        for group in sorted(scales):
            left, right = lower, upper
            while right - left > (upper - lower) / (golden**scale_iterations):
                first = right - (right - left) / golden
                second = left + (right - left) / golden
                first_scales, second_scales = dict(scales), dict(scales)
                first_scales[group], second_scales[group] = first, second
                first_fit, first_objective = fit_at(first_scales)
                second_fit, second_objective = fit_at(second_scales)
                records.extend(
                    [
                        {"group": group, "scale": first, "objective": first_objective},
                        {
                            "group": group,
                            "scale": second,
                            "objective": second_objective,
                        },
                    ]
                )
                candidate, fit, objective = (
                    (first_scales, first_fit, first_objective)
                    if first_objective < second_objective
                    else (second_scales, second_fit, second_objective)
                )
                if objective < best_objective:
                    scales, best_fit, best_objective = candidate, fit, objective
                if first_objective < second_objective:
                    right = second
                else:
                    left = first
    if best_fit is None:
        raise RuntimeError("theoretical group-scale profile produced no fit")
    penalty = sum(
        0.5 * ((value - scale_prior_mean) / scale_prior_sigma) ** 2
        for value in scales.values()
    )
    return TheoreticalGroupScaleProfileResult(
        scales,
        best_objective,
        best_fit.corrections_j_mol,
        float(penalty),
        metadata_key,
        tuple(records),
    )
