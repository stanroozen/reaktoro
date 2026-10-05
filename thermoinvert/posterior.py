"""Posterior diagnostics: Laplace approximation, grid evaluation, and MCMC sampling."""

from dataclasses import dataclass
from typing import Sequence

import numpy as np

from .inversion import PosteriorProblem

UNCERTAINTY_SOURCES = ("all", "data", "prior")


def _component_objective(problem: PosteriorProblem, uncertainty_source: str):
    """Return an objective callable restricted to one component of
    `InversionProblem.objective_components` (data misfit, prior penalty, or
    both). This isolates each source's contribution to the inferred
    uncertainty, in the spirit of mc_fit's analytical-vs-thermodynamic
    decomposition -- though here the split is by objective term, not by
    perturbing raw inputs, so results are not directly comparable to mc_fit's.
    """
    if uncertainty_source not in UNCERTAINTY_SOURCES:
        raise ValueError(f"uncertainty_source must be one of {UNCERTAINTY_SOURCES}")

    def objective(values: Sequence[float]) -> float:
        if uncertainty_source == "data":
            data, _ = problem.objective_components(values)
            return data
        if uncertainty_source == "prior":
            _, prior = problem.objective_components(values)
            return prior
        return problem.objective(values)  # includes any configured Bayes penalty

    return objective


@dataclass(frozen=True)
class LaplaceApproximation:
    """Local Gaussian approximation of the posterior around a MAP point.

    This is a numerical-Hessian diagnostic, not a posterior sample. Sparse
    phase-equilibrium posteriors may be asymmetric, ridge-like, bounded, or
    multimodal; treat this result as a first-order summary only, per the
    guardrails in the implementation design document.
    """

    parameter_names: tuple[str, ...]
    map_point: np.ndarray
    hessian: np.ndarray
    covariance: np.ndarray | None
    correlation: np.ndarray | None
    is_positive_definite: bool
    uncertainty_source: str = "all"
    boundary_limited: bool = False
    boundary_parameters: tuple[str, ...] = ()
    parameter_units: tuple[str, ...] = ()


def laplace_covariance(
    problem: PosteriorProblem,
    map_point: Sequence[float],
    step_fraction: float = 1.0e-2,
    uncertainty_source: str = "all",
) -> LaplaceApproximation:
    """Estimate posterior covariance from the numerical Hessian of the objective."""
    if len(map_point) != len(problem.parameters):
        raise ValueError("map_point does not match parameter count")
    if step_fraction <= 0.0:
        raise ValueError("step_fraction must be positive")
    objective_fn = _component_objective(problem, uncertainty_source)

    x0 = np.asarray(map_point, dtype=float)
    n = len(x0)
    bounds = np.asarray([parameter.bounds() for parameter in problem.parameters])
    if np.any(x0 < bounds[:, 0]) or np.any(x0 > bounds[:, 1]):
        raise ValueError("map_point must lie within parameter bounds")
    requested_steps = np.asarray(
        [
            max(abs(parameter.prior_sigma) * step_fraction, 1.0e-6)
            for parameter in problem.parameters
        ]
    )
    clearance = np.minimum(x0 - bounds[:, 0], bounds[:, 1] - x0)
    boundary_mask = clearance <= 0.0
    boundary_parameters = tuple(
        parameter.name
        for parameter, limited in zip(problem.parameters, boundary_mask)
        if limited
    )
    if boundary_parameters:
        return LaplaceApproximation(
            tuple(parameter.name for parameter in problem.parameters),
            x0,
            np.full((n, n), np.nan),
            None,
            None,
            False,
            uncertainty_source,
            True,
            boundary_parameters,
            tuple(
                getattr(parameter, "unit", "J/mol") for parameter in problem.parameters
            ),
        )
    steps = np.minimum(requested_steps, clearance * 0.5)

    def evaluate(offsets: dict[int, float]) -> float:
        point = x0.copy()
        for index, offset in offsets.items():
            point[index] += offset
        return objective_fn(point)

    phi0 = evaluate({})
    hessian = np.zeros((n, n))
    for i in range(n):
        phi_plus = evaluate({i: steps[i]})
        phi_minus = evaluate({i: -steps[i]})
        hessian[i, i] = (phi_plus - 2.0 * phi0 + phi_minus) / steps[i] ** 2
    for i in range(n):
        for j in range(i + 1, n):
            pp = evaluate({i: steps[i], j: steps[j]})
            pm = evaluate({i: steps[i], j: -steps[j]})
            mp = evaluate({i: -steps[i], j: steps[j]})
            mm = evaluate({i: -steps[i], j: -steps[j]})
            value = (pp - pm - mp + mm) / (4.0 * steps[i] * steps[j])
            hessian[i, j] = hessian[j, i] = value

    eigenvalues = np.linalg.eigvalsh(hessian)
    is_positive_definite = bool(np.all(eigenvalues > 0.0))
    covariance = None
    correlation = None
    if is_positive_definite:
        covariance = np.linalg.inv(hessian)
        sigma = np.sqrt(np.diag(covariance))
        correlation = covariance / np.outer(sigma, sigma)

    return LaplaceApproximation(
        tuple(parameter.name for parameter in problem.parameters),
        x0,
        hessian,
        covariance,
        correlation,
        is_positive_definite,
        uncertainty_source,
        False,
        (),
        tuple(getattr(parameter, "unit", "J/mol") for parameter in problem.parameters),
    )


@dataclass(frozen=True)
class GridPosteriorResult:
    """Brute-force posterior density evaluated on a regular parameter grid.

    Intended only as a verification tool for one or two parameters, per the
    design document's guidance that grid evaluation is impractical beyond a
    handful of dimensions but useful for sanity-checking MAP/Laplace results.
    """

    parameter_names: tuple[str, ...]
    axes: tuple[np.ndarray, ...]
    log_posterior: np.ndarray
    posterior_density: np.ndarray
    marginal_mean: dict[str, float]
    marginal_std: dict[str, float]


def grid_posterior(
    problem: PosteriorProblem,
    bounds: Sequence[tuple[float, float]],
    resolution: int = 25,
) -> GridPosteriorResult:
    """Evaluate a normalized posterior density on a regular grid of 1-2 parameters."""
    if len(bounds) != len(problem.parameters):
        raise ValueError("bounds must match parameter count")
    if len(bounds) not in (1, 2):
        raise ValueError("grid_posterior only supports one or two parameters")
    if resolution < 2:
        raise ValueError("resolution must be at least 2")
    declared_bounds = [parameter.bounds() for parameter in problem.parameters]
    for index, ((lower, upper), (declared_lower, declared_upper)) in enumerate(
        zip(bounds, declared_bounds)
    ):
        if lower >= upper:
            raise ValueError(f"grid bounds must increase for parameter {index}")
        if lower < declared_lower or upper > declared_upper:
            raise ValueError(
                f"grid bounds for parameter {index} must lie within declared parameter bounds"
            )

    axes = tuple(np.linspace(low, high, resolution) for low, high in bounds)
    if len(axes) == 1:
        phi = np.asarray([problem.objective([value]) for value in axes[0]])
    else:
        phi = np.asarray(
            [[problem.objective([a, b]) for b in axes[1]] for a in axes[0]]
        )

    if not np.all(np.isfinite(phi)):
        raise ValueError("grid contains no finite posterior objective everywhere")
    log_posterior = -(phi - np.min(phi))
    cell_volume = np.prod([axis[1] - axis[0] for axis in axes])
    density = np.exp(log_posterior)
    normalization = np.sum(density) * cell_volume
    posterior_density = density / normalization

    marginal_mean: dict[str, float] = {}
    marginal_std: dict[str, float] = {}
    for index, parameter in enumerate(problem.parameters):
        if len(axes) == 1:
            marginal = posterior_density
        else:
            marginal = np.sum(posterior_density, axis=1 - index) * (
                axes[1 - index][1] - axes[1 - index][0]
            )
        axis = axes[index]
        weights = marginal / (np.sum(marginal) if np.sum(marginal) > 0.0 else 1.0)
        mean = float(np.sum(axis * weights))
        variance = float(np.sum(weights * (axis - mean) ** 2))
        marginal_mean[parameter.name] = mean
        marginal_std[parameter.name] = float(np.sqrt(max(variance, 0.0)))

    return GridPosteriorResult(
        tuple(parameter.name for parameter in problem.parameters),
        axes,
        log_posterior,
        posterior_density,
        marginal_mean,
        marginal_std,
    )


@dataclass(frozen=True)
class MCMCResult:
    """Metropolis-Hastings posterior samples for a MAP inversion problem.

    This is a plain random-walk Metropolis sampler with a fixed Gaussian
    proposal (diagonal by default; optionally full covariance; no NUTS/HMC,
    scipy/emcee dependency). It is a
    practical baseline for small parameter counts; convergence must still be
    checked (e.g. acceptance rate, trace behaviour) before trusting the
    reported summaries, per the guardrails in the design document.
    """

    parameter_names: tuple[str, ...]
    samples: np.ndarray
    log_posterior: np.ndarray
    acceptance_rate: float
    posterior_mean: dict[str, float]
    posterior_std: dict[str, float]
    credible_interval_95: dict[str, tuple[float, float]]
    uncertainty_source: str = "all"
    proposal_covariance_j_mol2: np.ndarray | None = None
    seed: int | None = None
    burn_in: int = 0
    parameter_units: tuple[str, ...] = ()

    @property
    def proposal_covariance(self) -> np.ndarray | None:
        """Unit-neutral proposal covariance; legacy field name is retained."""
        return self.proposal_covariance_j_mol2


@dataclass(frozen=True)
class ParallelTemperingResult:
    """Cold-chain posterior and mode diagnostics from parallel tempering."""

    parameter_names: tuple[str, ...]
    parameter_units: tuple[str, ...]
    temperatures: tuple[float, ...]
    samples: np.ndarray
    log_posterior: np.ndarray
    mode_labels: np.ndarray
    mode_occupancy: dict[int, float]
    mode_means: dict[int, dict[str, float]]
    mode_best_objective: dict[int, float]
    acceptance_rates: tuple[float, ...]
    swap_acceptance_rates: tuple[float, ...]
    crossed_modes: bool
    seed: int
    burn_in: int
    proposal_sigma: tuple[float, ...]
    temperature_mode_occupancy: dict[float, dict[int, float]] | None = None
    mode_transition_count: int = 0
    mode_transition_rate: float = 0.0
    all_temperature_samples: np.ndarray | None = None
    temperature_mode_labels: np.ndarray | None = None
    all_temperature_objectives: np.ndarray | None = None
    uniform_reference_objectives: np.ndarray | None = None


@dataclass(frozen=True)
class ThermodynamicIntegrationResult:
    """Thermodynamic-integration estimate of a bounded target normalizer."""

    inverse_temperatures: tuple[float, ...]
    mean_objective: tuple[float, ...]
    log_normalizer_relative_to_uniform_bounds: float
    uniform_reference_sample_count: int
    monte_carlo_standard_error: float
    monte_carlo_interval_95: tuple[float, float]
    bootstrap_replicates: int
    block_length: int
    seed: int
    quadratic_ladder_log_normalizer: float | None
    quadrature_discrepancy: float | None
    maximum_inverse_temperature_step: float


@dataclass(frozen=True)
class GlobalModeAnalysis:
    """Shared-basin diagnostics across all retained tempering chains."""

    mode_labels: np.ndarray
    mode_probability: dict[int, float]
    temperature_mode_occupancy: dict[float, dict[int, float]]
    mode_means: dict[int, dict[str, float]]
    mode_best_objective: dict[int, float]
    temperature_mode_labels: np.ndarray
    mode_probability_by_separation: dict[float, dict[int, float]]
    clustering_mode_counts: dict[float, int]
    mode_probability_sensitivity: dict[int, tuple[float, float]]


def _cluster_mode_labels(points: np.ndarray, threshold: float) -> np.ndarray:
    labels = np.full(len(points), -1, dtype=int)
    centers: list[np.ndarray] = []
    for index, point in enumerate(points):
        distances = [float(np.linalg.norm(point - center)) for center in centers]
        if not distances or min(distances) > threshold:
            labels[index] = len(centers)
            centers.append(point.copy())
        else:
            labels[index] = int(np.argmin(distances))
            centers[labels[index]] = np.mean(
                points[: index + 1][labels[: index + 1] == labels[index]], axis=0
            )
    return labels


def analyze_parallel_tempering_modes(
    result: ParallelTemperingResult,
    mode_separation_fraction: float = 0.05,
    mode_separation_factors: Sequence[float] = (0.75, 1.0, 1.25),
) -> GlobalModeAnalysis:
    """Assign stable basin identities using every retained temperature chain.

    Mode identities are clustered jointly across all temperatures. Posterior
    mode probabilities are estimated from cold-chain occupancy, not from
    Bayesian evidence; hotter chains are used only to identify shared basins.
    """
    if mode_separation_fraction <= 0.0 or not np.isfinite(mode_separation_fraction):
        raise ValueError("mode_separation_fraction must be finite and positive")
    if not mode_separation_factors or any(
        not np.isfinite(factor) or factor <= 0.0 for factor in mode_separation_factors
    ):
        raise ValueError("mode_separation_factors must be finite and positive")
    if not any(np.isclose(factor, 1.0) for factor in mode_separation_factors):
        raise ValueError("mode_separation_factors must include 1.0")
    points_by_temperature = result.all_temperature_samples
    if points_by_temperature is None:
        points_by_temperature = result.samples[:, np.newaxis, :]
    if points_by_temperature.ndim != 3:
        raise ValueError("all_temperature_samples must be a 3D array")
    if points_by_temperature.shape[1] != len(result.temperatures):
        raise ValueError("temperature sample count must match temperatures")
    flat = points_by_temperature.reshape(-1, points_by_temperature.shape[-1])
    scale = np.ptp(flat, axis=0)
    scale = np.where(scale > 0.0, scale, 1.0)
    labels = _cluster_mode_labels(
        (flat / scale),
        mode_separation_fraction * np.sqrt(flat.shape[-1]),
    ).reshape(points_by_temperature.shape[:2])
    modes = sorted(set(int(value) for value in labels.ravel()))
    cold_labels = labels[:, 0]
    probability = {mode: float(np.mean(cold_labels == mode)) for mode in modes}
    occupancy = {
        temperature: {
            mode: float(np.mean(labels[:, index] == mode))
            for mode in modes
            if np.any(labels[:, index] == mode)
        }
        for index, temperature in enumerate(result.temperatures)
    }
    cold_samples = points_by_temperature[:, 0, :]
    mode_means = {
        mode: {
            name: float(np.mean(cold_samples[cold_labels == mode, index]))
            for index, name in enumerate(result.parameter_names)
        }
        for mode in modes
        if np.any(cold_labels == mode)
    }
    mode_best_objective = {
        mode: float(np.min(-result.log_posterior[cold_labels == mode]))
        for mode in mode_means
    }
    mode_probability_by_separation: dict[float, dict[int, float]] = {}
    clustering_mode_counts: dict[float, int] = {}
    for factor in mode_separation_factors:
        key = float(mode_separation_fraction * factor)
        alternative_labels = _cluster_mode_labels(
            flat / scale,
            mode_separation_fraction * factor * np.sqrt(flat.shape[-1]),
        ).reshape(points_by_temperature.shape[:2])
        alternative_cold = alternative_labels[:, 0]
        alternative_modes = sorted(
            set(int(value) for value in alternative_labels.ravel())
        )
        mapped_counts = {mode: 0 for mode in modes}
        for alternative_mode in alternative_modes:
            selected = alternative_cold == alternative_mode
            if not np.any(selected):
                continue
            overlap = np.bincount(
                cold_labels[selected], minlength=max(modes, default=-1) + 1
            )
            mapped_mode = int(np.argmax(overlap))
            mapped_counts[mapped_mode] += int(np.count_nonzero(selected))
        total_cold = max(len(alternative_cold), 1)
        mode_probability_by_separation[key] = {
            mode: float(count / total_cold) for mode, count in mapped_counts.items()
        }
        clustering_mode_counts[key] = len(alternative_modes)
    mode_probability_sensitivity = {
        mode: (
            float(
                min(values[mode] for values in mode_probability_by_separation.values())
            ),
            float(
                max(values[mode] for values in mode_probability_by_separation.values())
            ),
        )
        for mode in modes
    }
    return GlobalModeAnalysis(
        cold_labels,
        probability,
        occupancy,
        mode_means,
        mode_best_objective,
        labels,
        mode_probability_by_separation,
        clustering_mode_counts,
        mode_probability_sensitivity,
    )


def _quadratic_ladder_integral(x: np.ndarray, y: np.ndarray) -> float:
    """Integrate a sampled curve with quadratic panels over irregular nodes."""
    total = 0.0
    index = 0
    while index + 2 < len(x):
        left_width = x[index + 1] - x[index]
        right_width = x[index + 2] - x[index + 1]
        width = left_width + right_width
        total += (
            width
            / 6.0
            * (
                (2.0 - right_width / left_width) * y[index]
                + width**2 / (left_width * right_width) * y[index + 1]
                + (2.0 - left_width / right_width) * y[index + 2]
            )
        )
        index += 2
    if index + 1 < len(x):
        total += (x[index + 1] - x[index]) * (y[index] + y[index + 1]) / 2.0
    return float(total)


def thermodynamic_integration(
    result: ParallelTemperingResult,
    bootstrap_replicates: int = 500,
    block_length: int | None = None,
    seed: int = 0,
) -> ThermodynamicIntegrationResult:
    """Estimate log E_uniform[exp(-objective)] using the sampled temperature ladder.

    This is a thermodynamic-integration estimate for the objective-defined
    target relative to uniform measure on the declared bounded parameter box.
    It is not a Bayes factor: likelihood and prior normalization constants are
    not added by the inversion objective.
    """
    objectives = result.all_temperature_objectives
    reference = result.uniform_reference_objectives
    if objectives is None or reference is None or len(reference) < 2:
        raise ValueError(
            "thermodynamic integration requires retained temperature objectives "
            "and at least two uniform-reference objectives"
        )
    if (
        objectives.ndim != 2
        or objectives.shape[0] < 2
        or objectives.shape[1] != len(result.temperatures)
    ):
        raise ValueError(
            "temperature objective matrix must have at least two rows and match "
            "the temperature ladder"
        )
    if not np.all(np.isfinite(objectives)) or not np.all(np.isfinite(reference)):
        raise ValueError("thermodynamic integration objectives must be finite")
    if bootstrap_replicates < 2:
        raise ValueError("bootstrap_replicates must be at least two")
    if block_length is None:
        block_length = max(1, int(np.ceil(np.sqrt(len(objectives)))))
    if block_length < 1 or block_length > len(objectives):
        raise ValueError("block_length must be between one and retained sample count")
    inverse_temperatures = np.asarray(
        [0.0, *(1.0 / value for value in result.temperatures[::-1])], dtype=float
    )
    mean_objective = np.asarray(
        [float(np.mean(reference)), *np.mean(objectives, axis=0)[::-1]], dtype=float
    )
    log_normalizer = -float(np.trapezoid(mean_objective, inverse_temperatures))
    quadratic_log_normalizer = None
    quadrature_discrepancy = None
    if len(inverse_temperatures) >= 3:
        quadratic_log_normalizer = -_quadratic_ladder_integral(
            inverse_temperatures, mean_objective
        )
        quadrature_discrepancy = abs(quadratic_log_normalizer - log_normalizer)
    rng = np.random.default_rng(seed)
    sample_count = len(objectives)
    n_blocks = int(np.ceil(sample_count / block_length))
    bootstrap_estimates = np.empty(bootstrap_replicates, dtype=float)
    for replicate in range(bootstrap_replicates):
        starts = rng.integers(0, sample_count, size=n_blocks)
        indices = np.concatenate(
            [(start + np.arange(block_length)) % sample_count for start in starts]
        )[:sample_count]
        chain_means = np.mean(objectives[indices], axis=0)[::-1]
        reference_mean = float(
            np.mean(rng.choice(reference, size=len(reference), replace=True))
        )
        replicate_curve = np.concatenate(([reference_mean], chain_means))
        bootstrap_estimates[replicate] = -float(
            np.trapezoid(replicate_curve, inverse_temperatures)
        )
    standard_error = float(np.std(bootstrap_estimates, ddof=1))
    interval = tuple(
        float(value) for value in np.quantile(bootstrap_estimates, [0.025, 0.975])
    )
    return ThermodynamicIntegrationResult(
        tuple(float(value) for value in inverse_temperatures),
        tuple(float(value) for value in mean_objective),
        log_normalizer,
        int(len(reference)),
        standard_error,
        interval,
        bootstrap_replicates,
        block_length,
        seed,
        quadratic_log_normalizer,
        quadrature_discrepancy,
        float(np.max(np.diff(inverse_temperatures))),
    )


def parallel_tempering_sample(
    problem: PosteriorProblem,
    initial_points: Sequence[Sequence[float]],
    temperatures: Sequence[float] | None = None,
    n_samples: int = 5000,
    burn_in: int = 500,
    proposal_sigma: Sequence[float] | None = None,
    swap_interval: int = 1,
    mode_separation_fraction: float = 0.05,
    seed: int = 0,
    uncertainty_source: str = "all",
    uniform_reference_samples: int = 128,
) -> ParallelTemperingResult:
    """Explore bounded, potentially multimodal posteriors with temperature swaps."""
    if len(initial_points) < 2:
        raise ValueError("parallel tempering requires at least two chains")
    if n_samples <= 0 or burn_in < 0 or burn_in >= n_samples:
        raise ValueError("n_samples must be positive and exceed burn_in")
    if swap_interval <= 0:
        raise ValueError("swap_interval must be positive")
    if mode_separation_fraction <= 0.0 or not np.isfinite(mode_separation_fraction):
        raise ValueError("mode_separation_fraction must be finite and positive")
    if uniform_reference_samples < 2:
        raise ValueError("uniform_reference_samples must be at least 2")
    objective_fn = _component_objective(problem, uncertainty_source)
    n_parameters = len(problem.parameters)
    bounds = np.asarray(
        [parameter.bounds() for parameter in problem.parameters], dtype=float
    )
    if temperatures is None:
        temperatures = tuple(float(2.0**index) for index in range(len(initial_points)))
    temperatures = tuple(float(value) for value in temperatures)
    if len(temperatures) != len(initial_points) or temperatures[0] != 1.0:
        raise ValueError("temperatures must match chains and start with 1.0")
    if any(not np.isfinite(value) or value < 1.0 for value in temperatures):
        raise ValueError("temperatures must be finite and at least 1")
    if any(left >= right for left, right in zip(temperatures, temperatures[1:])):
        raise ValueError("temperatures must be strictly increasing")
    states = np.asarray(initial_points, dtype=float)
    if states.shape != (len(initial_points), n_parameters):
        raise ValueError("initial_points must match chain and parameter counts")
    if (
        not np.all(np.isfinite(states))
        or np.any(states < bounds[:, 0])
        or np.any(states > bounds[:, 1])
    ):
        raise ValueError("initial_points must be finite and within parameter bounds")
    if proposal_sigma is None:
        proposal = np.asarray(
            [
                max(abs(parameter.prior_sigma) * 0.1, 1.0e-6)
                for parameter in problem.parameters
            ]
        )
    else:
        if len(proposal_sigma) != n_parameters:
            raise ValueError("proposal_sigma must match parameter count")
        proposal = np.asarray(proposal_sigma, dtype=float)
    if not np.all(np.isfinite(proposal)) or np.any(proposal <= 0.0):
        raise ValueError("proposal_sigma must be finite and positive")
    objectives = np.asarray([objective_fn(state) for state in states], dtype=float)
    if not np.all(np.isfinite(objectives)):
        raise ValueError("initial_points must have finite posterior objectives")

    rng = np.random.default_rng(seed)
    chain_accepts = np.zeros(len(states), dtype=int)
    swap_attempts = np.zeros(len(states) - 1, dtype=int)
    swap_accepts = np.zeros(len(states) - 1, dtype=int)
    all_chain_samples = np.empty((n_samples, len(states), n_parameters))
    all_chain_objectives = np.empty((n_samples, len(states)))
    cold_samples = np.empty((n_samples, n_parameters))
    cold_log_posterior = np.empty(n_samples)
    for step in range(n_samples):
        for chain_index, temperature in enumerate(temperatures):
            candidate = states[chain_index] + proposal * rng.normal(size=n_parameters)
            if np.any(candidate < bounds[:, 0]) or np.any(candidate > bounds[:, 1]):
                continue
            candidate_objective = objective_fn(candidate)
            if (
                np.log(rng.uniform())
                < (objectives[chain_index] - candidate_objective) / temperature
            ):
                states[chain_index] = candidate
                objectives[chain_index] = candidate_objective
                chain_accepts[chain_index] += 1
        if (step + 1) % swap_interval == 0:
            for index in range(len(states) - 1):
                swap_attempts[index] += 1
                log_ratio = (
                    1.0 / temperatures[index] - 1.0 / temperatures[index + 1]
                ) * (objectives[index] - objectives[index + 1])
                if np.log(rng.uniform()) < log_ratio:
                    states[[index, index + 1]] = states[[index + 1, index]]
                    objectives[index], objectives[index + 1] = (
                        objectives[index + 1],
                        objectives[index],
                    )
                    swap_accepts[index] += 1
        cold_samples[step] = states[0]
        cold_log_posterior[step] = -objectives[0]
        all_chain_samples[step] = states
        all_chain_objectives[step] = objectives

    samples = cold_samples[burn_in:]
    log_posterior = cold_log_posterior[burn_in:]
    scaled = (samples - bounds[:, 0]) / (bounds[:, 1] - bounds[:, 0])
    threshold = mode_separation_fraction * np.sqrt(n_parameters)
    labels = _cluster_mode_labels(scaled, threshold)
    mode_transition_count = int(np.count_nonzero(labels[1:] != labels[:-1]))
    mode_transition_rate = float(mode_transition_count / max(len(labels) - 1, 1))
    mode_occupancy = {
        mode: float(np.mean(labels == mode)) for mode in sorted(set(labels))
    }
    names = tuple(parameter.name for parameter in problem.parameters)
    mode_means = {
        mode: {
            name: float(np.mean(samples[labels == mode, index]))
            for index, name in enumerate(names)
        }
        for mode in mode_occupancy
    }
    mode_best_objective = {
        mode: float(np.min(-log_posterior[labels == mode])) for mode in mode_occupancy
    }
    post_burn_all = all_chain_samples[burn_in:]
    post_burn_objectives = all_chain_objectives[burn_in:]
    temperature_mode_occupancy = {
        temperature: {
            mode: float(np.mean(chain_labels == mode))
            for mode in sorted(set(chain_labels))
        }
        for temperature, chain_labels in zip(
            temperatures,
            (
                _cluster_mode_labels(
                    (post_burn_all[:, index, :] - bounds[:, 0])
                    / (bounds[:, 1] - bounds[:, 0]),
                    threshold,
                )
                for index in range(len(temperatures))
            ),
        )
    }
    global_mode_analysis = analyze_parallel_tempering_modes(
        ParallelTemperingResult(
            names,
            tuple(
                getattr(parameter, "unit", "J/mol") for parameter in problem.parameters
            ),
            temperatures,
            samples,
            log_posterior,
            labels,
            mode_occupancy,
            mode_means,
            mode_best_objective,
            tuple(chain_accepts / n_samples),
            tuple(
                np.divide(
                    swap_accepts,
                    swap_attempts,
                    out=np.zeros_like(swap_accepts, dtype=float),
                    where=swap_attempts > 0,
                )
            ),
            len(set(labels)) > 1,
            seed,
            burn_in,
            tuple(float(value) for value in proposal),
            temperature_mode_occupancy,
            mode_transition_count,
            mode_transition_rate,
            post_burn_all,
            None,
            post_burn_objectives,
        ),
        mode_separation_fraction,
    )
    return ParallelTemperingResult(
        names,
        tuple(getattr(parameter, "unit", "J/mol") for parameter in problem.parameters),
        temperatures,
        samples,
        log_posterior,
        labels,
        mode_occupancy,
        mode_means,
        mode_best_objective,
        tuple(chain_accepts / n_samples),
        tuple(
            np.divide(
                swap_accepts,
                swap_attempts,
                out=np.zeros_like(swap_accepts, dtype=float),
                where=swap_attempts > 0,
            )
        ),
        len(set(labels)) > 1,
        seed,
        burn_in,
        tuple(float(value) for value in proposal),
        temperature_mode_occupancy,
        mode_transition_count,
        mode_transition_rate,
        post_burn_all,
        global_mode_analysis.temperature_mode_labels,
        post_burn_objectives,
        np.asarray(
            [
                objective_fn(
                    bounds[:, 0]
                    + rng.random(n_parameters) * (bounds[:, 1] - bounds[:, 0])
                )
                for _ in range(uniform_reference_samples)
            ],
            dtype=float,
        ),
    )


def mcmc_sample(
    problem: PosteriorProblem,
    initial_point: Sequence[float],
    n_samples: int = 5000,
    burn_in: int = 500,
    proposal_sigma_j_mol: Sequence[float] | None = None,
    proposal_covariance_j_mol2: np.ndarray | None = None,
    seed: int = 0,
    uncertainty_source: str = "all",
    proposal_sigma: Sequence[float] | None = None,
    proposal_covariance: np.ndarray | None = None,
) -> MCMCResult:
    """Draw posterior samples with a random-walk Metropolis-Hastings sampler."""
    if len(initial_point) != len(problem.parameters):
        raise ValueError("initial_point does not match parameter count")
    if n_samples <= 0:
        raise ValueError("n_samples must be positive")
    if burn_in < 0 or burn_in >= n_samples:
        raise ValueError("burn_in must be non-negative and smaller than n_samples")
    bounds = np.asarray([parameter.bounds() for parameter in problem.parameters])
    initial = np.asarray(initial_point, dtype=float)
    if not np.all(np.isfinite(initial)):
        raise ValueError("initial_point must contain finite values")
    if np.any(initial < bounds[:, 0]) or np.any(initial > bounds[:, 1]):
        raise ValueError("initial_point must lie within parameter bounds")
    objective_fn = _component_objective(problem, uncertainty_source)

    if proposal_sigma is not None:
        if proposal_sigma_j_mol is not None:
            raise ValueError("Specify proposal_sigma or proposal_sigma_j_mol, not both")
        proposal_sigma_j_mol = proposal_sigma
    if proposal_covariance is not None:
        if proposal_covariance_j_mol2 is not None:
            raise ValueError(
                "Specify proposal_covariance or proposal_covariance_j_mol2, not both"
            )
        proposal_covariance_j_mol2 = proposal_covariance
    if proposal_sigma_j_mol is not None and proposal_covariance_j_mol2 is not None:
        raise ValueError("Specify proposal covariance or proposal sigmas, not both")
    if proposal_covariance_j_mol2 is not None:
        proposal_covariance = np.asarray(proposal_covariance_j_mol2, dtype=float)
        expected_shape = (len(problem.parameters), len(problem.parameters))
        if proposal_covariance.shape != expected_shape:
            raise ValueError("proposal_covariance_j_mol2 must match parameter count")
        if not np.all(np.isfinite(proposal_covariance)):
            raise ValueError("proposal_covariance_j_mol2 must contain finite values")
        proposal_covariance = 0.5 * (proposal_covariance + proposal_covariance.T)
        if np.any(np.linalg.eigvalsh(proposal_covariance) <= 0.0):
            raise ValueError("proposal_covariance_j_mol2 must be positive definite")
        proposal_factor = np.linalg.cholesky(proposal_covariance)
    else:
        if proposal_sigma_j_mol is None:
            proposal_sigma = np.asarray(
                [
                    max(abs(parameter.prior_sigma) * 0.1, 1.0e-6)
                    for parameter in problem.parameters
                ]
            )
        else:
            if len(proposal_sigma_j_mol) != len(problem.parameters):
                raise ValueError("proposal_sigma_j_mol does not match parameter count")
            proposal_sigma = np.asarray(proposal_sigma_j_mol, dtype=float)
            if not np.all(np.isfinite(proposal_sigma)):
                raise ValueError("proposal_sigma_j_mol entries must be finite")
            if np.any(proposal_sigma <= 0.0):
                raise ValueError("proposal_sigma_j_mol entries must be positive")
        proposal_covariance = np.diag(proposal_sigma**2)
        proposal_factor = np.diag(proposal_sigma)

    rng = np.random.default_rng(seed)
    n = len(initial_point)
    current = initial
    current_phi = objective_fn(current)
    if not np.isfinite(current_phi):
        raise ValueError("initial_point has a non-finite posterior objective")

    samples = np.empty((n_samples, n))
    log_posterior = np.empty(n_samples)
    accepted = 0
    for step in range(n_samples):
        proposal = current + proposal_factor @ rng.normal(size=n)
        proposal_phi = objective_fn(proposal)
        if np.log(rng.uniform()) < current_phi - proposal_phi:
            current, current_phi = proposal, proposal_phi
            accepted += 1
        samples[step] = current
        log_posterior[step] = -current_phi

    kept = samples[burn_in:]
    kept_log_posterior = log_posterior[burn_in:]
    names = tuple(parameter.name for parameter in problem.parameters)
    posterior_mean, posterior_std, credible_interval_95 = _summarize_samples(
        names, kept
    )

    return MCMCResult(
        names,
        kept,
        kept_log_posterior,
        accepted / n_samples,
        posterior_mean,
        posterior_std,
        credible_interval_95,
        uncertainty_source,
        proposal_covariance,
        seed,
        burn_in,
        tuple(getattr(parameter, "unit", "J/mol") for parameter in problem.parameters),
    )


def _summarize_samples(
    parameter_names: Sequence[str], samples: np.ndarray
) -> tuple[dict[str, float], dict[str, float], dict[str, tuple[float, float]]]:
    posterior_mean: dict[str, float] = {}
    posterior_std: dict[str, float] = {}
    credible_interval_95: dict[str, tuple[float, float]] = {}
    for index, name in enumerate(parameter_names):
        column = samples[:, index]
        posterior_mean[name] = float(np.mean(column))
        posterior_std[name] = float(np.std(column, ddof=1)) if len(column) > 1 else 0.0
        low, high = np.percentile(column, [2.5, 97.5])
        credible_interval_95[name] = (float(low), float(high))
    return posterior_mean, posterior_std, credible_interval_95


def filter_outliers_mcmc(
    mcmc_result: MCMCResult,
    z_threshold: float = 3.0,
    max_iterations: int = 5,
) -> MCMCResult:
    """Iteratively drop samples whose objective is a one-sided outlier.

    Mirrors mc_fit's robust one-sided z-score filter on the objective
    function: samples with objective (``-log_posterior``) more than
    ``z_threshold`` standard deviations above the current mean are dropped,
    then the mean/std are recomputed on the remaining samples and the process
    repeats until stable or ``max_iterations`` is reached. Only the
    high-objective tail is filtered (a low objective is never "too good").
    """
    if z_threshold <= 0.0:
        raise ValueError("z_threshold must be positive")
    if max_iterations < 1:
        raise ValueError("max_iterations must be at least one")

    objective = -mcmc_result.log_posterior
    keep = np.ones(len(objective), dtype=bool)
    for _ in range(max_iterations):
        kept_objective = objective[keep]
        if len(kept_objective) < 2:
            break
        mean, std = np.mean(kept_objective), np.std(kept_objective)
        if std <= 0.0:
            break
        new_keep = keep & (objective <= mean + z_threshold * std)
        if np.array_equal(new_keep, keep):
            break
        keep = new_keep

    filtered_samples = mcmc_result.samples[keep]
    filtered_log_posterior = mcmc_result.log_posterior[keep]
    posterior_mean, posterior_std, credible_interval_95 = _summarize_samples(
        mcmc_result.parameter_names, filtered_samples
    )
    return MCMCResult(
        mcmc_result.parameter_names,
        filtered_samples,
        filtered_log_posterior,
        mcmc_result.acceptance_rate,
        posterior_mean,
        posterior_std,
        credible_interval_95,
        mcmc_result.uncertainty_source,
        mcmc_result.proposal_covariance_j_mol2,
        mcmc_result.seed,
        mcmc_result.burn_in,
        mcmc_result.parameter_units,
    )


def filter_by_acceptance_threshold(
    mcmc_result: MCMCResult,
    max_objective: float,
) -> MCMCResult:
    """Keep only samples with objective (``-log_posterior``) at or below a
    threshold, mirroring mc_fit's ``oktol`` model-acceptance cutoff."""
    objective = -mcmc_result.log_posterior
    keep = objective <= max_objective
    if not np.any(keep):
        raise ValueError("No samples satisfy the acceptance threshold")

    filtered_samples = mcmc_result.samples[keep]
    filtered_log_posterior = mcmc_result.log_posterior[keep]
    posterior_mean, posterior_std, credible_interval_95 = _summarize_samples(
        mcmc_result.parameter_names, filtered_samples
    )
    return MCMCResult(
        mcmc_result.parameter_names,
        filtered_samples,
        filtered_log_posterior,
        mcmc_result.acceptance_rate,
        posterior_mean,
        posterior_std,
        credible_interval_95,
        mcmc_result.uncertainty_source,
        mcmc_result.proposal_covariance_j_mol2,
        mcmc_result.seed,
        mcmc_result.burn_in,
        mcmc_result.parameter_units,
    )


def effective_sample_size(samples_1d: np.ndarray) -> float:
    """Autocorrelation-based effective sample size for one parameter's chain.

    Uses Geyer's initial positive sequence estimator: sums consecutive pairs
    of the autocorrelation function until a pair sum turns negative.
    """
    n = len(samples_1d)
    if n < 2:
        return float(n)
    centered = np.asarray(samples_1d, dtype=float) - np.mean(samples_1d)
    variance = np.var(centered)
    if variance <= 0.0:
        return float(n)
    autocovariance = np.correlate(centered, centered, mode="full")[n - 1 :] / n
    acf = autocovariance / variance

    rho_sum = 0.0
    lag = 1
    while lag + 1 < len(acf):
        pair_sum = acf[lag] + acf[lag + 1]
        if pair_sum < 0.0:
            break
        rho_sum += pair_sum
        lag += 2
    ess = n / (1.0 + 2.0 * rho_sum)
    return float(np.clip(ess, 1.0, n))


def potential_scale_reduction(chains: Sequence[np.ndarray]) -> float:
    """Gelman-Rubin R-hat for several equal-length chains of one parameter.

    Values near 1.0 indicate the chains agree; values above roughly 1.1
    indicate the chains have not converged to the same distribution.
    """
    stacked = np.asarray(chains, dtype=float)
    if stacked.ndim != 2:
        raise ValueError("chains must be a 2D array of shape (n_chains, n_samples)")
    m, n = stacked.shape
    if m < 2:
        raise ValueError("At least two chains are required for R-hat")
    if n < 2:
        raise ValueError("Each chain needs at least two samples")

    chain_means = np.mean(stacked, axis=1)
    chain_variances = np.var(stacked, axis=1, ddof=1)
    within_chain_variance = np.mean(chain_variances)
    between_chain_variance = n * np.var(chain_means, ddof=1)
    pooled_variance = ((n - 1) / n) * within_chain_variance + between_chain_variance / n
    if within_chain_variance <= 0.0:
        return float("inf") if pooled_variance > 0.0 else 1.0
    return float(np.sqrt(pooled_variance / within_chain_variance))


@dataclass(frozen=True)
class MultiChainMCMCResult:
    """Several independent MCMC chains plus their convergence diagnostics."""

    parameter_names: tuple[str, ...]
    chains: tuple[MCMCResult, ...]
    r_hat: dict[str, float]
    effective_sample_size: dict[str, float]
    initial_points: tuple[tuple[float, ...], ...] = ()
    seeds: tuple[int, ...] = ()
    parameter_units: tuple[str, ...] = ()


def mcmc_sample_multi_chain(
    problem: PosteriorProblem,
    initial_points: Sequence[Sequence[float]],
    n_samples: int = 5000,
    burn_in: int = 500,
    proposal_sigma_j_mol: Sequence[float] | None = None,
    proposal_covariance_j_mol2: np.ndarray | None = None,
    seeds: Sequence[int] | None = None,
    uncertainty_source: str = "all",
    proposal_sigma: Sequence[float] | None = None,
    proposal_covariance: np.ndarray | None = None,
) -> MultiChainMCMCResult:
    """Run several independent chains from different starting points and
    report Gelman-Rubin R-hat and total effective sample size per parameter.

    At least two chains are required; convergence cannot be assessed from a
    single chain. Each starting point should ideally be over-dispersed
    relative to the expected posterior for R-hat to be a meaningful check.
    """
    if len(initial_points) < 2:
        raise ValueError("At least two chains are required for convergence diagnostics")
    if seeds is None:
        seeds = list(range(len(initial_points)))
    if len(seeds) != len(initial_points):
        raise ValueError("seeds must match initial_points count")

    chains = tuple(
        mcmc_sample(
            problem,
            point,
            n_samples=n_samples,
            burn_in=burn_in,
            proposal_sigma_j_mol=proposal_sigma_j_mol,
            proposal_covariance_j_mol2=proposal_covariance_j_mol2,
            seed=seed,
            uncertainty_source=uncertainty_source,
            proposal_sigma=proposal_sigma,
            proposal_covariance=proposal_covariance,
        )
        for point, seed in zip(initial_points, seeds)
    )
    names = chains[0].parameter_names
    parameter_units = chains[0].parameter_units
    if any(chain.parameter_units != parameter_units for chain in chains[1:]):
        raise ValueError("all MCMC chains must use the same parameter units")
    r_hat: dict[str, float] = {}
    ess: dict[str, float] = {}
    for index, name in enumerate(names):
        stacked = np.stack([chain.samples[:, index] for chain in chains])
        r_hat[name] = potential_scale_reduction(stacked)
        ess[name] = float(
            sum(effective_sample_size(chain.samples[:, index]) for chain in chains)
        )
    return MultiChainMCMCResult(
        names,
        chains,
        r_hat,
        ess,
        tuple(tuple(float(value) for value in point) for point in initial_points),
        tuple(int(seed) for seed in seeds),
        parameter_units,
    )
