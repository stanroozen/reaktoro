"""Combined identifiability diagnostic: prior-vs-posterior shrinkage plus
structural (SVD) weak-direction cross-referencing.

Uncertainty width alone (a credible interval or posterior std) does not tell
you whether a parameter is data-constrained: a narrow posterior can simply
reflect a narrow prior. This module compares a "data+prior" result against a
"prior-only" result (both produced by `posterior.laplace_covariance` or
`posterior.mcmc_sample` with `uncertainty_source`) to compute how much the
data actually shrank each parameter's uncertainty, and cross-references
`sensitivity.identify_sensitivity`'s structurally weak directions so a
parameter that only looks constrained because of prior curvature is flagged
either way.
"""

from dataclasses import dataclass
from typing import Sequence

import numpy as np

from .inversion import InversionProblem, MAPResult
from .posterior import GridPosteriorResult, LaplaceApproximation, MCMCResult
from .sensitivity import IdentifiabilityResult

STATUSES = ("data_constrained", "ambiguous", "prior_dominated")


@dataclass(frozen=True)
class ParameterIdentifiability:
    """Shrinkage-based and structural identifiability summary for one parameter."""

    name: str
    prior_variance: float
    posterior_variance: float
    shrinkage: float
    status: str
    structurally_weak: bool


@dataclass(frozen=True)
class IdentifiabilitySummary:
    """Per-parameter identifiability summary for an `InversionProblem`."""

    parameters: tuple[ParameterIdentifiability, ...]

    def __getitem__(self, name: str) -> ParameterIdentifiability:
        for parameter in self.parameters:
            if parameter.name == name:
                return parameter
        raise KeyError(name)


@dataclass(frozen=True)
class MAPBasin:
    """One cluster of cold-start MAP endpoints in normalized bound space."""

    representative_j_mol: dict[str, float]
    objective: float
    start_count: int
    near_optimal: bool


@dataclass(frozen=True)
class GlobalityDiagnostic:
    """Evidence from multi-start MAP searches, not a global-resolution proof."""

    basins: tuple[MAPBasin, ...]
    near_optimal_basin_count: int
    consistent_cold_starts: bool


@dataclass(frozen=True)
class ParameterPosteriorShape:
    """Sample-based marginal shape measures for one posterior parameter."""

    name: str
    skewness: float
    excess_kurtosis: float
    standardized_mean_median_difference: float
    non_gaussian_indicated: bool


@dataclass(frozen=True)
class PosteriorShapeDiagnostic:
    """Heuristic evidence for non-Gaussian marginal posterior structure."""

    parameters: tuple[ParameterPosteriorShape, ...]

    @property
    def non_gaussian_indicated(self) -> bool:
        return any(parameter.non_gaussian_indicated for parameter in self.parameters)


@dataclass(frozen=True)
class GridPosteriorMode:
    """A local maximum of a one- or two-parameter grid posterior."""

    coordinates: dict[str, float]
    relative_density: float


@dataclass(frozen=True)
class GridModeDiagnostic:
    """Discrete mode summary from an explicitly evaluated grid posterior."""

    modes: tuple[GridPosteriorMode, ...]

    @property
    def multimodal(self) -> bool:
        return len(self.modes) > 1


@dataclass(frozen=True)
class GridInformationGain:
    """Numerical KL information gain over an evaluated bounded posterior grid."""

    nats: float
    grid_covers_parameter_bounds: bool


def diagnose_grid_modes(
    result: GridPosteriorResult, minimum_relative_density: float = 0.05
) -> GridModeDiagnostic:
    """Detect separated local maxima in a one- or two-dimensional grid posterior.

    A reported mode must be strictly higher than at least one neighboring grid
    cell and no lower than any neighbor, avoiding repeated reports across a
    flat plateau. Only modes at least ``minimum_relative_density`` of the
    global mode are retained. The diagnostic is limited to the provided grid
    bounds and resolution; it cannot establish modes outside that domain.
    """
    if not 0.0 < minimum_relative_density <= 1.0:
        raise ValueError("minimum_relative_density must be in (0, 1]")
    density = np.asarray(result.posterior_density, dtype=float)
    if density.ndim not in (1, 2) or density.shape != tuple(
        len(axis) for axis in result.axes
    ):
        raise ValueError("grid posterior density must match one or two axes")
    global_density = float(np.max(density))
    if not np.isfinite(global_density) or global_density <= 0.0:
        raise ValueError("grid posterior density must have a positive maximum")

    modes: list[GridPosteriorMode] = []
    for index in np.ndindex(density.shape):
        value = float(density[index])
        if value < minimum_relative_density * global_density:
            continue
        neighbor_values: list[float] = []
        for offset in np.ndindex(*(3,) * density.ndim):
            shift = tuple(component - 1 for component in offset)
            if not any(shift):
                continue
            neighbor = tuple(position + delta for position, delta in zip(index, shift))
            if all(
                0 <= position < size for position, size in zip(neighbor, density.shape)
            ):
                neighbor_values.append(float(density[neighbor]))
        if (
            neighbor_values
            and value >= max(neighbor_values)
            and value > min(neighbor_values)
        ):
            modes.append(
                GridPosteriorMode(
                    {
                        name: float(axis[position])
                        for name, axis, position in zip(
                            result.parameter_names, result.axes, index
                        )
                    },
                    value / global_density,
                )
            )
    return GridModeDiagnostic(tuple(modes))


def grid_information_gain(
    problem: InversionProblem, result: GridPosteriorResult
) -> GridInformationGain:
    """Estimate KL(posterior || prior) directly from a one- or two-D grid.

    Both distributions are normalized over the supplied grid domain. The KL is
    therefore a full truncated-prior result only when every grid axis spans the
    corresponding declared parameter bounds; otherwise it is conditional on the
    selected grid window.
    """
    parameter_names = tuple(parameter.name for parameter in problem.parameters)
    if result.parameter_names != parameter_names:
        raise ValueError("grid parameter names must match the problem order")
    if len(result.axes) not in (1, 2):
        raise ValueError("grid information gain supports one or two parameters")
    density = np.asarray(result.posterior_density, dtype=float)
    if density.shape != tuple(len(axis) for axis in result.axes):
        raise ValueError("grid posterior density must match grid axes")
    if np.any(density < 0.0) or not np.all(np.isfinite(density)):
        raise ValueError("grid posterior density must be finite and non-negative")
    spacings = []
    for axis in result.axes:
        values = np.asarray(axis, dtype=float)
        if values.ndim != 1 or len(values) < 2:
            raise ValueError("each grid axis needs at least two values")
        differences = np.diff(values)
        if np.any(differences <= 0.0) or not np.allclose(differences, differences[0]):
            raise ValueError("grid axes must be uniformly increasing")
        spacings.append(float(differences[0]))

    mesh = np.meshgrid(*result.axes, indexing="ij")
    points = np.column_stack([coordinate.ravel() for coordinate in mesh])
    mean = np.asarray(
        [parameter.prior_mean for parameter in problem.parameters], dtype=float
    )
    deviations = points - mean
    precision = np.linalg.inv(problem.prior_covariance)
    log_prior = -0.5 * np.einsum("ij,jk,ik->i", deviations, precision, deviations)
    prior_mass = np.exp(log_prior - np.max(log_prior))
    prior_mass /= np.sum(prior_mass)
    posterior_mass = density.ravel() * np.prod(spacings)
    posterior_mass /= np.sum(posterior_mass)
    positive = posterior_mass > 0.0
    nats = float(
        np.sum(
            posterior_mass[positive]
            * np.log(posterior_mass[positive] / prior_mass[positive])
        )
    )
    covers_bounds = all(
        np.isclose(axis[0], parameter.bounds()[0])
        and np.isclose(axis[-1], parameter.bounds()[1])
        for axis, parameter in zip(result.axes, problem.parameters)
    )
    return GridInformationGain(nats, covers_bounds)


def diagnose_posterior_shape(
    result: MCMCResult,
    skewness_threshold: float = 0.5,
    excess_kurtosis_threshold: float = 1.0,
    mean_median_threshold: float = 0.25,
) -> PosteriorShapeDiagnostic:
    """Flag asymmetric or heavy/light-tailed MCMC marginal distributions.

    The thresholds are screening criteria, not hypothesis tests. A positive
    indication means that Laplace/Gaussian summaries are questionable and
    samples, quantiles, grids, and separated-chain checks should be examined.
    It cannot detect every multimodal or jointly non-Gaussian posterior.
    """
    if (
        skewness_threshold < 0.0
        or excess_kurtosis_threshold < 0.0
        or mean_median_threshold < 0.0
    ):
        raise ValueError("posterior-shape thresholds must be non-negative")
    samples = np.asarray(result.samples, dtype=float)
    if samples.ndim != 2 or samples.shape[1] != len(result.parameter_names):
        raise ValueError("MCMC samples must match parameter names")
    if samples.shape[0] < 4:
        raise ValueError("At least four MCMC samples are required")

    parameters: list[ParameterPosteriorShape] = []
    for index, name in enumerate(result.parameter_names):
        values = samples[:, index]
        mean = float(np.mean(values))
        standard_deviation = float(np.std(values, ddof=1))
        if standard_deviation <= 0.0:
            skewness = excess_kurtosis = mean_median_difference = 0.0
        else:
            standardized = (values - mean) / standard_deviation
            skewness = float(np.mean(standardized**3))
            excess_kurtosis = float(np.mean(standardized**4) - 3.0)
            mean_median_difference = (
                abs(mean - float(np.median(values))) / standard_deviation
            )
        indicated = (
            abs(skewness) >= skewness_threshold
            or abs(excess_kurtosis) >= excess_kurtosis_threshold
            or mean_median_difference >= mean_median_threshold
        )
        parameters.append(
            ParameterPosteriorShape(
                name,
                skewness,
                excess_kurtosis,
                mean_median_difference,
                indicated,
            )
        )
    return PosteriorShapeDiagnostic(tuple(parameters))


def diagnose_map_globality(
    problem: InversionProblem,
    map_result: MAPResult,
    separation_fraction: float = 0.05,
    objective_relative_tolerance: float = 1.0e-3,
    objective_absolute_tolerance: float = 1.0e-8,
) -> GlobalityDiagnostic:
    """Cluster cold-start MAP endpoints and flag competing near-optimal basins.

    Endpoint separation is measured after dividing each parameter difference by
    its bounded search width. A single near-optimal basin is evidence that the
    cold starts agree, but it does not establish a global resolution operator.
    """
    if separation_fraction <= 0.0:
        raise ValueError("separation_fraction must be positive")
    if objective_relative_tolerance < 0.0 or objective_absolute_tolerance < 0.0:
        raise ValueError("objective tolerances must be non-negative")
    if not map_result.starts:
        raise ValueError("map_result must contain at least one start record")

    names = tuple(parameter.key for parameter in problem.parameters)
    widths = np.asarray(
        [
            upper - lower
            for lower, upper in (parameter.bounds() for parameter in problem.parameters)
        ],
        dtype=float,
    )
    records = sorted(map_result.starts, key=lambda record: float(record["objective"]))
    representatives: list[np.ndarray] = []
    objectives: list[float] = []
    counts: list[int] = []
    for record in records:
        endpoint = np.asarray(record["map_j_mol"], dtype=float)
        if endpoint.shape != widths.shape:
            raise ValueError("MAP endpoint does not match problem parameter count")
        objective = float(record["objective"])
        distances = [
            np.linalg.norm((endpoint - point) / widths) for point in representatives
        ]
        matching = next(
            (
                index
                for index, distance in enumerate(distances)
                if distance <= separation_fraction
            ),
            None,
        )
        if matching is None:
            representatives.append(endpoint)
            objectives.append(objective)
            counts.append(1)
        else:
            counts[matching] += 1

    best_objective = min(objectives)
    objective_tolerance = max(
        objective_absolute_tolerance,
        objective_relative_tolerance * max(abs(best_objective), 1.0),
    )
    basins = tuple(
        MAPBasin(
            dict(zip(names, point.tolist())),
            objective,
            count,
            objective <= best_objective + objective_tolerance,
        )
        for point, objective, count in zip(representatives, objectives, counts)
    )
    near_optimal_basin_count = sum(basin.near_optimal for basin in basins)
    return GlobalityDiagnostic(
        basins,
        near_optimal_basin_count,
        near_optimal_basin_count == 1,
    )


def _variance_by_name(
    result: LaplaceApproximation | MCMCResult,
) -> dict[str, float | None]:
    if isinstance(result, LaplaceApproximation):
        if result.covariance is None:
            return {name: None for name in result.parameter_names}
        diagonal = np.diag(result.covariance)
        return dict(zip(result.parameter_names, diagonal.tolist()))
    if isinstance(result, MCMCResult):
        return {name: std**2 for name, std in result.posterior_std.items()}
    raise TypeError("result must be a LaplaceApproximation or MCMCResult")


def identifiability_summary(
    problem: InversionProblem,
    all_result: LaplaceApproximation | MCMCResult,
    prior_result: LaplaceApproximation | MCMCResult,
    sensitivity: IdentifiabilityResult | None = None,
    data_constrained_threshold: float = 0.5,
    prior_dominated_threshold: float = 0.1,
    weak_component_threshold: float = 0.5,
) -> IdentifiabilitySummary:
    """Combine prior-vs-posterior shrinkage with structural weak-direction flags.

    ``all_result`` and ``prior_result`` should come from the same diagnostic
    (both `laplace_covariance` or both `mcmc_sample`) run with
    ``uncertainty_source="all"`` and ``uncertainty_source="prior"``
    respectively, at a comparable point. ``shrinkage = 1 - posterior_var /
    prior_var``: near 0 means the data barely narrowed the prior belief
    (prior-dominated / unconstrained by this data); near 1 means the data
    dominates (well constrained). If ``sensitivity`` (from
    `sensitivity.identify_sensitivity`) is given, a parameter with a large
    component in any weak singular direction is flagged
    ``structurally_weak=True`` regardless of its shrinkage, since that
    reflects a forward-model limitation no amount of data can fix.
    """
    if not 0.0 <= prior_dominated_threshold < data_constrained_threshold <= 1.0:
        raise ValueError(
            "Require 0 <= prior_dominated_threshold < data_constrained_threshold <= 1"
        )

    all_variance = _variance_by_name(all_result)
    prior_variance = _variance_by_name(prior_result)

    weak_names: set[str] = set()
    if sensitivity is not None:
        for direction in sensitivity.weak_directions:
            for name, component in direction.items():
                if abs(component) >= weak_component_threshold:
                    weak_names.add(name)

    parameters = []
    for parameter in problem.parameters:
        name = parameter.name
        prior_var = prior_variance.get(name)
        posterior_var = all_variance.get(name)
        if prior_var is None or posterior_var is None or prior_var <= 0.0:
            shrinkage = float("nan")
            status = "ambiguous"
        else:
            shrinkage = 1.0 - posterior_var / prior_var
            if shrinkage >= data_constrained_threshold:
                status = "data_constrained"
            elif shrinkage <= prior_dominated_threshold:
                status = "prior_dominated"
            else:
                status = "ambiguous"
        parameters.append(
            ParameterIdentifiability(
                name,
                float(prior_var) if prior_var is not None else float("nan"),
                float(posterior_var) if posterior_var is not None else float("nan"),
                float(shrinkage),
                status,
                name in weak_names,
            )
        )
    return IdentifiabilitySummary(tuple(parameters))


def _mean_by_name(result: LaplaceApproximation | MCMCResult) -> dict[str, float]:
    if isinstance(result, LaplaceApproximation):
        return dict(zip(result.parameter_names, np.asarray(result.map_point).tolist()))
    if isinstance(result, MCMCResult):
        return dict(result.posterior_mean)
    raise TypeError("result must be a LaplaceApproximation or MCMCResult")


def gaussian_kl_divergence(
    posterior_mean: float,
    posterior_variance: float,
    prior_mean: float,
    prior_variance: float,
) -> float:
    """KL(posterior || prior) in nats for two univariate Gaussians.

    This is the standard "information gain" measure: how many nats of
    information the data added relative to the prior belief. Zero means the
    posterior is identical to the prior (no information gained).
    """
    if prior_variance <= 0.0 or posterior_variance <= 0.0:
        raise ValueError("Variances must be positive")
    return 0.5 * (
        np.log(prior_variance / posterior_variance)
        + (posterior_variance + (posterior_mean - prior_mean) ** 2) / prior_variance
        - 1.0
    )


def gaussian_multivariate_kl_divergence(
    posterior_mean: Sequence[float],
    posterior_covariance: np.ndarray,
    prior_mean: Sequence[float],
    prior_covariance: np.ndarray,
) -> float:
    """Return KL(posterior || prior) in nats for multivariate Gaussians.

    Unlike the per-parameter KL summary, this includes posterior correlations
    and covariance-volume changes, so it measures information in the joint
    parameter space rather than summing independent marginal approximations.
    """
    posterior_mean = np.asarray(posterior_mean, dtype=float)
    prior_mean = np.asarray(prior_mean, dtype=float)
    posterior_covariance = np.asarray(posterior_covariance, dtype=float)
    prior_covariance = np.asarray(prior_covariance, dtype=float)
    if posterior_mean.ndim != 1 or prior_mean.ndim != 1:
        raise ValueError("Gaussian means must be one-dimensional")
    if posterior_mean.shape != prior_mean.shape:
        raise ValueError("Gaussian means must have matching shapes")
    if not np.all(np.isfinite(posterior_mean)) or not np.all(np.isfinite(prior_mean)):
        raise ValueError("Gaussian means must contain finite values")
    dimension = len(posterior_mean)
    if posterior_covariance.shape != (dimension, dimension):
        raise ValueError("posterior covariance shape does not match mean")
    if prior_covariance.shape != (dimension, dimension):
        raise ValueError("prior covariance shape does not match mean")
    if not np.all(np.isfinite(posterior_covariance)) or not np.all(
        np.isfinite(prior_covariance)
    ):
        raise ValueError("Gaussian covariance matrices must contain finite values")

    posterior_covariance = 0.5 * (posterior_covariance + posterior_covariance.T)
    prior_covariance = 0.5 * (prior_covariance + prior_covariance.T)
    prior_sign, prior_logdet = np.linalg.slogdet(prior_covariance)
    posterior_sign, posterior_logdet = np.linalg.slogdet(posterior_covariance)
    if prior_sign <= 0.0 or posterior_sign <= 0.0:
        raise ValueError("Gaussian covariance matrices must be positive definite")
    delta = posterior_mean - prior_mean
    precision_prior = np.linalg.inv(prior_covariance)
    return float(
        0.5
        * (
            prior_logdet
            - posterior_logdet
            - dimension
            + np.trace(precision_prior @ posterior_covariance)
            + delta @ precision_prior @ delta
        )
    )


def _result_mean_covariance(
    result: LaplaceApproximation | MCMCResult,
) -> tuple[np.ndarray, np.ndarray]:
    if isinstance(result, LaplaceApproximation):
        if result.covariance is None:
            raise ValueError("Laplace result has no positive-definite covariance")
        return np.asarray(result.map_point, dtype=float), np.asarray(result.covariance)
    if isinstance(result, MCMCResult):
        if len(result.samples) < 2:
            raise ValueError("MCMC result needs at least two samples")
        return np.mean(result.samples, axis=0), np.cov(result.samples, rowvar=False)
    raise TypeError("result must be a LaplaceApproximation or MCMCResult")


def multivariate_information_gain(
    problem: InversionProblem,
    posterior_result: LaplaceApproximation | MCMCResult,
) -> float:
    """Return joint KL(posterior || the problem's Gaussian prior) in nats."""
    posterior_mean, posterior_covariance = _result_mean_covariance(posterior_result)
    prior_mean = np.asarray(
        [parameter.prior_mean for parameter in problem.parameters], dtype=float
    )
    return gaussian_multivariate_kl_divergence(
        posterior_mean,
        posterior_covariance,
        prior_mean,
        problem.prior_covariance,
    )


def posterior_correlation_matrix(
    result: LaplaceApproximation | MCMCResult,
) -> np.ndarray:
    """Return the full posterior correlation matrix."""
    _, covariance = _result_mean_covariance(result)
    standard_deviation = np.sqrt(np.diag(covariance))
    if np.any(standard_deviation <= 0.0):
        raise ValueError("Posterior covariance has a non-positive variance")
    return covariance / np.outer(standard_deviation, standard_deviation)


def information_gain_summary(
    problem: InversionProblem,
    all_result: LaplaceApproximation | MCMCResult,
) -> dict[str, float]:
    """Per-parameter KL(posterior || prior) in nats, using each parameter's
    exact Gaussian prior (`prior_mean_j_mol`/`prior_sigma_j_mol`) directly --
    no separate "prior-only" diagnostic run is needed, since the prior term
    in `InversionProblem.objective` is already exactly Gaussian by
    construction. NaN is reported where the posterior variance is
    unavailable (e.g. a non-positive-definite Laplace Hessian).
    """
    means = _mean_by_name(all_result)
    variances = _variance_by_name(all_result)
    gains: dict[str, float] = {}
    for parameter in problem.parameters:
        name = parameter.name
        posterior_mean = means.get(name)
        posterior_variance = variances.get(name)
        if (
            posterior_mean is None
            or posterior_variance is None
            or posterior_variance <= 0.0
        ):
            gains[name] = float("nan")
            continue
        gains[name] = float(
            gaussian_kl_divergence(
                posterior_mean,
                posterior_variance,
                parameter.prior_mean,
                parameter.prior_sigma**2,
            )
        )
    return gains
