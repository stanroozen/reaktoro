"""Numerical sensitivity and SVD-based identifiability diagnostics."""

import csv
import json
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np

from .data import DirichletConcentrationParameter, Experiment, ThermodynamicParameter
from .forward_model import EquilibriumResult, ReaktoroForwardModel
from .inversion import InversionProblem
from .posterior import MCMCResult, MultiChainMCMCResult
from .pt_inversion import PTInversionProblem


@dataclass(frozen=True)
class IdentifiabilityResult:
    """SVD summary for a numerical forward-model sensitivity matrix."""

    singular_values: np.ndarray
    rank: int
    condition_number: float
    parameter_names: tuple[str, ...]
    weak_directions: tuple[dict[str, float], ...]


@dataclass(frozen=True)
class SensitivityResult:
    """Finite-difference sensitivity matrix and its identifiability summary."""

    observable_names: tuple[str, ...]
    parameter_names: tuple[str, ...]
    values: np.ndarray
    jacobian: np.ndarray
    identifiability: IdentifiabilityResult


@dataclass(frozen=True)
class LocalResolutionResult:
    """Local linear-Gaussian parameter resolution around one sensitivity point."""

    parameter_names: tuple[str, ...]
    matrix: np.ndarray
    diagonal: np.ndarray
    posterior_covariance: np.ndarray


@dataclass(frozen=True)
class PosteriorResolutionResult:
    """Elementwise posterior summaries of local resolution matrices."""

    parameter_names: tuple[str, ...]
    sample_count: int
    mean_matrix: np.ndarray
    lower_matrix: np.ndarray
    upper_matrix: np.ndarray
    mean_diagonal: np.ndarray
    lower_diagonal: np.ndarray
    upper_diagonal: np.ndarray
    selected_chain_indices: tuple[int, ...] = ()
    selected_draw_indices: tuple[int, ...] = ()
    selection_strategy: str = "unspecified"
    observable_names: tuple[str, ...] = ()
    credible_interval: tuple[float, float] = (0.05, 0.95)
    step_fraction: float | None = None
    observation_covariance_source: str = "unspecified"
    theoretical_covariance_source: str = "unspecified"
    observation_covariance: np.ndarray | None = None
    theoretical_covariance: np.ndarray | None = None
    prior_covariance: np.ndarray | None = None


def local_resolution_matrix(
    jacobian: np.ndarray,
    parameter_names: Sequence[str],
    prior_covariance: np.ndarray,
    observation_covariance: np.ndarray | None = None,
    theoretical_covariance: np.ndarray | None = None,
) -> LocalResolutionResult:
    """Return the Tarantola local resolution matrix for a linearized problem.

    For a Jacobian $G$, analytical/model observation covariance $C_D$, and
    parameter-prior covariance $C_M$, the resolution is
    ``(G.T @ C_D^-1 @ G + C_M^-1)^-1 @ G.T @ C_D^-1 @ G``. It describes only
    the local differentiable Gaussian approximation and is not valid as a
    global diagnostic across phase-boundary discontinuities.
    """
    matrix = np.asarray(jacobian, dtype=float)
    if matrix.ndim != 2:
        raise ValueError("jacobian must be two-dimensional")
    n_observations, n_parameters = matrix.shape
    if n_parameters != len(parameter_names):
        raise ValueError("parameter_names must match jacobian columns")
    prior = np.asarray(prior_covariance, dtype=float)
    if prior.shape != (n_parameters, n_parameters):
        raise ValueError("prior_covariance must match jacobian columns")
    if not np.all(np.isfinite(prior)):
        raise ValueError("prior_covariance must contain finite values")
    prior = 0.5 * (prior + prior.T)
    if np.any(np.linalg.eigvalsh(prior) <= 0.0):
        raise ValueError("prior_covariance must be positive definite")
    observations = np.eye(n_observations)
    if observation_covariance is not None:
        observations = np.asarray(observation_covariance, dtype=float)
        if observations.shape != (n_observations, n_observations):
            raise ValueError("observation_covariance must match jacobian rows")
        if not np.all(np.isfinite(observations)):
            raise ValueError("observation_covariance must contain finite values")
    if theoretical_covariance is not None:
        theoretical = np.asarray(theoretical_covariance, dtype=float)
        if theoretical.shape != (n_observations, n_observations):
            raise ValueError("theoretical_covariance must match jacobian rows")
        if not np.all(np.isfinite(theoretical)):
            raise ValueError("theoretical_covariance must contain finite values")
        observations = observations + theoretical
    observations = 0.5 * (observations + observations.T)
    if np.any(np.linalg.eigvalsh(observations) <= 0.0):
        raise ValueError("combined observation covariance must be positive definite")

    data_information = matrix.T @ np.linalg.solve(observations, matrix)
    posterior_covariance = np.linalg.inv(data_information + np.linalg.inv(prior))
    resolution = posterior_covariance @ data_information
    return LocalResolutionResult(
        tuple(parameter_names),
        resolution,
        np.diag(resolution).copy(),
        posterior_covariance,
    )


def posterior_averaged_resolution(
    jacobians: np.ndarray,
    parameter_names: Sequence[str],
    prior_covariance: np.ndarray,
    observation_covariance: np.ndarray | None = None,
    theoretical_covariance: np.ndarray | None = None,
    credible_interval: tuple[float, float] = (0.05, 0.95),
    observable_names: Sequence[str] = (),
) -> PosteriorResolutionResult:
    """Summarize local resolution over Jacobians evaluated across a posterior.

    ``jacobians`` has shape ``(samples, observations, parameters)``. Bounds
    are elementwise quantiles and are descriptive, not simultaneous intervals.
    """
    samples = np.asarray(jacobians, dtype=float)
    if samples.ndim != 3 or samples.shape[0] == 0:
        raise ValueError(
            "jacobians must have shape (samples, observations, parameters)"
        )
    if not np.all(np.isfinite(samples)):
        raise ValueError("jacobians must contain finite values")
    if samples.shape[2] != len(parameter_names):
        raise ValueError("parameter_names must match jacobian columns")
    if observable_names and (
        len(observable_names) != samples.shape[1]
        or any(not isinstance(name, str) or not name for name in observable_names)
        or len(set(observable_names)) != len(observable_names)
    ):
        raise ValueError("observable_names must uniquely match jacobian rows")
    lower, upper = credible_interval
    if not (0.0 <= lower < upper <= 1.0):
        raise ValueError("credible_interval must satisfy 0 <= lower < upper <= 1")

    matrices = np.stack(
        [
            local_resolution_matrix(
                jacobian,
                parameter_names,
                prior_covariance,
                observation_covariance,
                theoretical_covariance,
            ).matrix
            for jacobian in samples
        ]
    )
    lower_matrix, upper_matrix = np.quantile(matrices, [lower, upper], axis=0)
    diagonals = np.diagonal(matrices, axis1=1, axis2=2)
    lower_diagonal, upper_diagonal = np.quantile(diagonals, [lower, upper], axis=0)
    return PosteriorResolutionResult(
        tuple(parameter_names),
        samples.shape[0],
        np.mean(matrices, axis=0),
        lower_matrix,
        upper_matrix,
        np.mean(diagonals, axis=0),
        lower_diagonal,
        upper_diagonal,
        observable_names=tuple(observable_names),
        credible_interval=(float(lower), float(upper)),
        observation_covariance=(
            np.asarray(observation_covariance, dtype=float).copy()
            if observation_covariance is not None
            else None
        ),
        theoretical_covariance=(
            np.asarray(theoretical_covariance, dtype=float).copy()
            if theoretical_covariance is not None
            else None
        ),
        prior_covariance=np.asarray(prior_covariance, dtype=float).copy(),
    )


def _resolution_observation_model(
    problem: InversionProblem,
    explicit_observation_covariance: np.ndarray | None,
    explicit_theoretical_covariance: np.ndarray | None,
) -> tuple[list[tuple[Experiment, str]], np.ndarray, np.ndarray]:
    specs: list[tuple[Experiment, str]] = []
    all_keys_by_id: dict[str, list[str]] = {}
    for experiment in problem.experiments:
        keys = [
            *(
                f"fraction:{phase}"
                for phase in sorted(experiment.phase_fractions or {})
            ),
            *(
                f"composition:{phase}:{component}"
                for phase in sorted(experiment.phase_compositions or {})
                for component in sorted((experiment.phase_compositions or {})[phase])
            ),
        ]
        all_keys_by_id[experiment.id] = keys
        if not keys:
            continue
        include_all = any(
            value is not None
            for value in (
                explicit_observation_covariance,
                explicit_theoretical_covariance,
                experiment.observation_covariance,
                experiment.theoretical_covariance,
            )
        )
        active = keys if include_all else []
        if not include_all:
            for key in keys:
                if key.startswith("fraction:"):
                    sigma = (experiment.phase_fraction_sigma or {}).get(key[9:])
                else:
                    _, phase, component = key.split(":", 2)
                    sigma = (
                        (experiment.phase_composition_sigma or {})
                        .get(phase, {})
                        .get(component)
                    )
                if sigma is not None and np.isfinite(sigma) and sigma > 0.0:
                    active.append(key)
        specs.extend((experiment, key) for key in active)

    size = len(specs)
    if size == 0:
        raise ValueError(
            "posterior resolution requires continuous observations with uncertainty"
        )
    analytical = np.zeros((size, size))
    theoretical = np.zeros((size, size))
    offset = 0
    for experiment in problem.experiments:
        keys = all_keys_by_id[experiment.id]
        active = [key for owner, key in specs if owner.id == experiment.id]
        count = len(active)
        if not count:
            continue
        if explicit_observation_covariance is not None:
            pass
        elif experiment.observation_covariance is not None:
            covariance_keys = experiment.observation_covariance_keys
            if set(covariance_keys) != set(keys) or len(covariance_keys) != len(keys):
                raise ValueError(
                    "observation_covariance_keys must identify every continuous observation"
                )
            covariance = np.asarray(experiment.observation_covariance, dtype=float)
            if covariance.shape != (len(keys), len(keys)) or not np.all(
                np.isfinite(covariance)
            ):
                raise ValueError(
                    "observation_covariance must be finite and match its keys"
                )
            covariance = 0.5 * (covariance + covariance.T)
            if np.any(np.linalg.eigvalsh(covariance) <= 0.0):
                raise ValueError("observation_covariance must be positive definite")
            order = [covariance_keys.index(key) for key in active]
            analytical[offset : offset + count, offset : offset + count] = covariance[
                np.ix_(order, order)
            ]
        else:
            sigmas: list[float] = []
            for key in active:
                if key.startswith("fraction:"):
                    sigma = (experiment.phase_fraction_sigma or {}).get(key[9:])
                else:
                    _, phase, component = key.split(":", 2)
                    sigma = (
                        (experiment.phase_composition_sigma or {})
                        .get(phase, {})
                        .get(component)
                    )
                if sigma is None or not np.isfinite(sigma) or sigma <= 0.0:
                    raise ValueError(
                        "continuous observations require positive analytical sigmas"
                    )
                sigmas.append(float(sigma))
            analytical[offset : offset + count, offset : offset + count] = np.diag(
                np.square(sigmas)
            )

        if (
            explicit_theoretical_covariance is None
            and experiment.theoretical_covariance is not None
        ):
            covariance_keys = experiment.theoretical_covariance_keys
            if set(covariance_keys) != set(keys) or len(covariance_keys) != len(keys):
                raise ValueError(
                    "theoretical_covariance_keys must identify every continuous observation"
                )
            covariance = np.asarray(experiment.theoretical_covariance, dtype=float)
            if covariance.shape != (len(keys), len(keys)) or not np.all(
                np.isfinite(covariance)
            ):
                raise ValueError(
                    "theoretical_covariance must be finite and match its keys"
                )
            covariance = 0.5 * (covariance + covariance.T)
            if np.any(np.linalg.eigvalsh(covariance) <= 0.0):
                raise ValueError("theoretical_covariance must be positive definite")
            order = [covariance_keys.index(key) for key in active]
            if problem.theoretical_covariance_scale_metadata_key:
                group = experiment.metadata.get(
                    problem.theoretical_covariance_scale_metadata_key
                )
                scale = problem.theoretical_covariance_scales[group]
            else:
                scale = problem.theoretical_covariance_scale
            theoretical[offset : offset + count, offset : offset + count] = (
                scale**2 * covariance[np.ix_(order, order)]
            )
        offset += count

    observation = (
        np.asarray(explicit_observation_covariance, dtype=float)
        if explicit_observation_covariance is not None
        else analytical
    )
    model_error = (
        np.asarray(explicit_theoretical_covariance, dtype=float)
        if explicit_theoretical_covariance is not None
        else theoretical
    )
    return specs, observation, model_error


def posterior_resolution_from_mcmc(
    problem: InversionProblem,
    posterior: MCMCResult | MultiChainMCMCResult,
    observation_covariance: np.ndarray | None = None,
    theoretical_covariance: np.ndarray | None = None,
    step_fraction: float = 1.0e-3,
    max_samples: int = 100,
    credible_interval: tuple[float, float] = (0.05, 0.95),
) -> PosteriorResolutionResult:
    """Compute posterior-wide local resolution from sampled parameter vectors.

    The function evaluates experiment observables only; bracket inequalities
    are discontinuous and are therefore excluded. Draws are selected at
    evenly spaced indices when ``max_samples`` is smaller than the posterior.
    For multiple chains, the selection is balanced across chains and requires
    ``max_samples`` to be at least the chain count.
    The observable layout must remain constant under finite differences and
    across selected draws, otherwise local Jacobians are not comparable.
    When covariance matrices are omitted, continuous-observation covariance
    is assembled from experiment sigmas, keyed covariance, and model scales.
    """
    if problem.brackets:
        raise ValueError("posterior resolution requires experiments, not brackets")
    if not problem.experiments:
        raise ValueError("posterior resolution requires at least one experiment")
    if not np.isfinite(step_fraction) or step_fraction <= 0.0:
        raise ValueError("step_fraction must be finite and positive")
    if (
        isinstance(max_samples, bool)
        or not isinstance(max_samples, int)
        or max_samples <= 0
    ):
        raise ValueError("max_samples must be a positive integer")
    parameter_names = tuple(parameter.name for parameter in problem.parameters)
    parameter_units = tuple(
        getattr(parameter, "unit", "J/mol") for parameter in problem.parameters
    )
    bounds = np.asarray([parameter.bounds() for parameter in problem.parameters])
    if posterior.parameter_names != parameter_names:
        raise ValueError("posterior parameter names must match problem parameter order")

    def validate_units(units: tuple[str, ...], source: str) -> None:
        if units and (len(units) != len(parameter_units) or units != parameter_units):
            raise ValueError(
                f"{source} parameter units must match problem parameter order"
            )

    validate_units(posterior.parameter_units, "posterior")
    if isinstance(posterior, MultiChainMCMCResult):
        if len(posterior.chains) < 2:
            raise ValueError("multi-chain posterior must contain at least two chains")
        if max_samples < len(posterior.chains):
            raise ValueError("max_samples must be at least the number of chains")
        chain_samples = [
            np.asarray(chain.samples, dtype=float) for chain in posterior.chains
        ]
        for chain, samples in zip(posterior.chains, chain_samples):
            if chain.parameter_names != parameter_names:
                raise ValueError(
                    "all chain parameter names must match problem parameter order"
                )
            validate_units(chain.parameter_units, "chain")
            if (
                samples.ndim != 2
                or samples.shape[0] == 0
                or samples.shape[1] != len(parameter_names)
            ):
                raise ValueError("each chain must contain a non-empty parameter matrix")
            if not np.all(np.isfinite(samples)):
                raise ValueError("posterior samples must contain finite values")
            if np.any(samples < bounds[:, 0]) or np.any(samples > bounds[:, 1]):
                raise ValueError("posterior samples must lie within parameter bounds")

        quotas = np.ones(len(chain_samples), dtype=int)
        remaining = max_samples - len(chain_samples)
        while remaining:
            advanced = False
            for index, samples in enumerate(chain_samples):
                if quotas[index] < len(samples):
                    quotas[index] += 1
                    remaining -= 1
                    advanced = True
                    if remaining == 0:
                        break
            if not advanced:
                break
        chain_indices = tuple(
            chain_index
            for chain_index, quota in enumerate(quotas)
            for _ in range(quota)
        )
        draw_indices = tuple(
            int(draw_index)
            for samples, quota in zip(chain_samples, quotas)
            for draw_index in np.linspace(0, len(samples) - 1, quota, dtype=int)
        )
        selected = np.stack(
            [
                chain_samples[chain_index][draw_index]
                for chain_index, draw_index in zip(chain_indices, draw_indices)
            ]
        )
        selection_strategy = "evenly_spaced_balanced_chains"
    else:
        samples = np.asarray(posterior.samples, dtype=float)
        if (
            samples.ndim != 2
            or samples.shape[0] == 0
            or samples.shape[1] != len(parameter_names)
        ):
            raise ValueError("posterior samples must be a non-empty parameter matrix")
        if not np.all(np.isfinite(samples)):
            raise ValueError("posterior samples must contain finite values")
        if np.any(samples < bounds[:, 0]) or np.any(samples > bounds[:, 1]):
            raise ValueError("posterior samples must lie within parameter bounds")
        draw_indices = tuple(
            int(index)
            for index in np.linspace(
                0, len(samples) - 1, min(len(samples), max_samples), dtype=int
            )
        )
        chain_indices = tuple(0 for _ in draw_indices)
        selected = samples[np.asarray(draw_indices, dtype=int)]
        selection_strategy = "evenly_spaced_single_chain"
    observation_covariance_was_explicit = observation_covariance is not None
    theoretical_covariance_was_explicit = theoretical_covariance is not None
    observable_specs, observation_covariance, theoretical_covariance = (
        _resolution_observation_model(
            problem, observation_covariance, theoretical_covariance
        )
    )

    def evaluate(
        values: np.ndarray,
        expected_names: tuple[str, ...] | None = None,
        expected_assemblages: tuple[tuple[str, ...], ...] | None = None,
    ):
        corrections = problem.corrections_for_values(values)
        predictions: dict[str, EquilibriumResult] = {}
        assemblages: list[tuple[str, ...]] = []
        for experiment in problem.experiments:
            prediction = problem.forward_model.predict(experiment, corrections)
            if not prediction.valid:
                raise ValueError(
                    f"forward model failed for experiment {experiment.id!r} "
                    "during posterior resolution"
                )
            predictions[experiment.id] = prediction
            assemblages.append(tuple(sorted(prediction.stable_phases)))
        names: list[str] = []
        observables: list[float] = []
        for experiment, key in observable_specs:
            prediction = predictions[experiment.id]
            if key.startswith("fraction:"):
                value = prediction.phase_fractions.get(key[9:], 0.0)
            else:
                _, phase, component = key.split(":", 2)
                value = prediction.phase_compositions.get(phase, {}).get(component, 0.0)
            names.append(f"{experiment.id}:{key}")
            observables.append(float(value))
        result_names = tuple(names)
        if expected_names is not None and result_names != expected_names:
            raise ValueError(
                "observable layout changed across posterior draws or finite differences"
            )
        result_assemblages = tuple(assemblages)
        if (
            expected_assemblages is not None
            and result_assemblages != expected_assemblages
        ):
            raise ValueError(
                "stable-phase assemblage changed across posterior draws or "
                "finite differences; local resolution is not valid across phase boundaries"
            )
        return result_names, np.asarray(observables, dtype=float), result_assemblages

    jacobians: list[np.ndarray] = []
    observable_names: tuple[str, ...] | None = None
    reference_assemblages: tuple[tuple[str, ...], ...] | None = None
    for sample in selected:
        observable_names, baseline, reference_assemblages = evaluate(
            sample, observable_names, reference_assemblages
        )
        columns: list[np.ndarray] = []
        for index, parameter in enumerate(problem.parameters):
            step = max(abs(parameter.prior_sigma) * step_fraction, 1.0e-9)
            lower_room = sample[index] - bounds[index, 0]
            upper_room = bounds[index, 1] - sample[index]
            if lower_room >= step and upper_room >= step:
                plus = sample.copy()
                minus = sample.copy()
                plus[index] += step
                minus[index] -= step
                _, plus_values, _ = evaluate(
                    plus, observable_names, reference_assemblages
                )
                _, minus_values, _ = evaluate(
                    minus, observable_names, reference_assemblages
                )
                column = (plus_values - minus_values) / (2.0 * step)
            elif upper_room > 0.0:
                step = min(step, upper_room)
                plus = sample.copy()
                plus[index] += step
                _, plus_values, _ = evaluate(
                    plus, observable_names, reference_assemblages
                )
                column = (plus_values - baseline) / step
            elif lower_room > 0.0:
                step = min(step, lower_room)
                minus = sample.copy()
                minus[index] -= step
                _, minus_values, _ = evaluate(
                    minus, observable_names, reference_assemblages
                )
                column = (baseline - minus_values) / step
            else:
                raise ValueError(
                    f"cannot perturb parameter {parameter.name!r} within bounds"
                )
            columns.append(column)
        jacobians.append(np.column_stack(columns))

    result = posterior_averaged_resolution(
        np.stack(jacobians),
        parameter_names,
        problem.prior_covariance,
        observation_covariance,
        theoretical_covariance,
        credible_interval,
        observable_names or (),
    )
    return replace(
        result,
        selected_chain_indices=chain_indices,
        selected_draw_indices=draw_indices,
        selection_strategy=selection_strategy,
        step_fraction=float(step_fraction),
        observation_covariance_source=(
            "explicit" if observation_covariance_was_explicit else "problem_assembled"
        ),
        theoretical_covariance_source=(
            "explicit" if theoretical_covariance_was_explicit else "problem_assembled"
        ),
        prior_covariance=np.asarray(problem.prior_covariance, dtype=float).copy(),
    )


def prior_scaled_information(
    jacobian: np.ndarray,
    parameters: Sequence[ThermodynamicParameter],
    prior_covariance: np.ndarray | None = None,
) -> np.ndarray:
    """Return Fisher-like information in prior-whitened coordinates.

    With no covariance, columns are scaled by independent prior sigmas. With a
    correlated covariance $C_M=L L^T$, the Jacobian is transformed as $J L$,
    preserving prior correlations instead of replacing them with marginal
    variances.
    """
    if np.asarray(jacobian).shape[1] != len(parameters):
        raise ValueError("parameters must match jacobian columns")
    if prior_covariance is None:
        scales = np.asarray(
            [parameter.prior_sigma for parameter in parameters], dtype=float
        )
        if np.any(scales <= 0.0):
            raise ValueError("Prior sigmas must be positive")
        transform = np.diag(scales)
    else:
        covariance = np.asarray(prior_covariance, dtype=float)
        expected = (len(parameters), len(parameters))
        if covariance.shape != expected:
            raise ValueError("prior_covariance must match parameter count")
        if not np.all(np.isfinite(covariance)):
            raise ValueError("prior_covariance must contain finite values")
        covariance = 0.5 * (covariance + covariance.T)
        try:
            transform = np.linalg.cholesky(covariance)
        except np.linalg.LinAlgError as error:
            raise ValueError("prior_covariance must be positive definite") from error
    scaled = np.asarray(jacobian, dtype=float) @ transform
    return scaled.T @ scaled


def write_sensitivity_report(
    result: SensitivityResult,
    parameters: Sequence[ThermodynamicParameter],
    output_dir: str | Path,
    prior_covariance: np.ndarray | None = None,
) -> None:
    """Write Jacobian, singular values, weak directions, and information matrices."""
    if tuple(parameter.name for parameter in parameters) != result.parameter_names:
        raise ValueError("parameters must match the sensitivity result")
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    with (output / "sensitivity_matrix.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.writer(handle)
        writer.writerow(["observable"] + list(result.parameter_names))
        writer.writerows(
            [name, *row.tolist()]
            for name, row in zip(result.observable_names, result.jacobian)
        )
    np.savetxt(
        output / "singular_values.csv",
        result.identifiability.singular_values,
        delimiter=",",
    )
    np.savetxt(
        output / "prior_scaled_information.csv",
        prior_scaled_information(result.jacobian, parameters, prior_covariance),
        delimiter=",",
    )
    report = {
        "parameter_names": list(result.parameter_names),
        "observable_names": list(result.observable_names),
        "prior_covariance_j_mol2": (
            np.asarray(prior_covariance, dtype=float).tolist()
            if prior_covariance is not None
            else None
        ),
        "rank": result.identifiability.rank,
        "condition_number": result.identifiability.condition_number,
        "singular_values": result.identifiability.singular_values.tolist(),
        "weak_directions": list(result.identifiability.weak_directions),
    }
    with (output / "identifiability.json").open("w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2)
    np.savetxt(output / "values.csv", result.values, delimiter=",")


def read_sensitivity_report(output_dir: str | Path) -> SensitivityResult:
    """Reload a sensitivity report written by `write_sensitivity_report`."""
    output = Path(output_dir)
    with (output / "identifiability.json").open(encoding="utf-8") as handle:
        report = json.load(handle)
    observable_names = tuple(report["observable_names"])
    parameter_names = tuple(report["parameter_names"])
    rows: list[list[float]] = []
    with (output / "sensitivity_matrix.csv").open(
        newline="", encoding="utf-8"
    ) as handle:
        reader = csv.reader(handle)
        header = next(reader)
        if tuple(header[1:]) != parameter_names:
            raise ValueError("sensitivity report parameter names do not match header")
        for row in reader:
            if len(row) != len(parameter_names) + 1:
                raise ValueError("invalid sensitivity matrix row")
            if row[0] != observable_names[len(rows)]:
                raise ValueError("sensitivity report observable order does not match")
            rows.append([float(value) for value in row[1:]])
    if len(rows) != len(observable_names):
        raise ValueError("sensitivity report observable count does not match")
    values = np.atleast_1d(np.loadtxt(output / "values.csv", delimiter=","))
    jacobian = np.asarray(rows, dtype=float)
    singular_values = np.atleast_1d(
        np.loadtxt(output / "singular_values.csv", delimiter=",")
    )
    identifiability = IdentifiabilityResult(
        singular_values,
        int(report["rank"]),
        float(report["condition_number"]),
        parameter_names,
        tuple(
            {name: float(value) for name, value in direction.items()}
            for direction in report["weak_directions"]
        ),
    )
    if values.shape != (len(observable_names),):
        raise ValueError("sensitivity report values do not match observables")
    return SensitivityResult(
        observable_names, parameter_names, values, jacobian, identifiability
    )


def result_observables(
    experiment: Experiment, prediction: EquilibriumResult
) -> tuple[tuple[str, ...], np.ndarray]:
    """Flatten supported prediction values into a stable, named observable vector."""
    names: list[str] = []
    values: list[float] = []
    for phase in sorted(
        set(experiment.observed_phases) | set(prediction.stable_phases)
    ):
        names.append(f"phase_present:{phase}")
        values.append(float(phase in prediction.stable_phases))
    if experiment.phase_fractions:
        for phase in sorted(experiment.phase_fractions):
            names.append(f"phase_fraction:{phase}")
            values.append(float(prediction.phase_fractions.get(phase, 0.0)))
    if experiment.phase_compositions:
        for phase in sorted(experiment.phase_compositions):
            calculated = prediction.phase_compositions.get(phase, {})
            for component in sorted(experiment.phase_compositions[phase]):
                names.append(f"composition:{phase}:{component}")
                values.append(float(calculated.get(component, 0.0)))
    if not values:
        raise ValueError(f"Experiment {experiment.id!r} has no supported observables")
    return tuple(names), np.asarray(values, dtype=float)


def numerical_sensitivity(
    forward_model: ReaktoroForwardModel,
    experiments: Sequence[Experiment],
    parameters: Sequence[ThermodynamicParameter],
    step_fraction: float = 1.0e-3,
) -> SensitivityResult:
    """Calculate a central finite-difference Jacobian around zero corrections."""
    if not experiments:
        raise ValueError("At least one experiment is required")
    if not parameters:
        raise ValueError("At least one parameter is required")
    if step_fraction <= 0.0:
        raise ValueError("step_fraction must be positive")

    baseline_names: list[str] = []
    baseline_values: list[float] = []
    plus_values: list[np.ndarray] = []
    minus_values: list[np.ndarray] = []
    for experiment in experiments:
        prediction = forward_model.predict(experiment, {})
        names, values = result_observables(experiment, prediction)
        baseline_names.extend(f"{experiment.id}:{name}" for name in names)
        baseline_values.extend(values)

    for parameter in parameters:
        step = max(abs(parameter.prior_sigma) * step_fraction, 1.0e-9)
        if isinstance(parameter, DirichletConcentrationParameter):
            baseline = np.asarray(baseline_values, dtype=float)
            plus_values.append(baseline)
            minus_values.append(baseline)
            continue
        plus = {parameter.key: parameter.prior_mean + step}
        minus = {parameter.key: parameter.prior_mean - step}
        plus_values.append(
            np.concatenate(
                [
                    result_observables(
                        experiment, forward_model.predict(experiment, plus)
                    )[1]
                    for experiment in experiments
                ]
            )
        )
        minus_values.append(
            np.concatenate(
                [
                    result_observables(
                        experiment, forward_model.predict(experiment, minus)
                    )[1]
                    for experiment in experiments
                ]
            )
        )

    baseline = np.asarray(baseline_values, dtype=float)
    jacobian = np.column_stack(
        [
            (plus - minus)
            / (2.0 * max(abs(parameter.prior_sigma) * step_fraction, 1.0e-9))
            for parameter, plus, minus in zip(parameters, plus_values, minus_values)
        ]
    )
    identifiability = identify_sensitivity(
        jacobian, [parameter.name for parameter in parameters]
    )
    return SensitivityResult(
        tuple(baseline_names),
        tuple(parameter.name for parameter in parameters),
        baseline,
        jacobian,
        identifiability,
    )


def numerical_pt_sensitivity(
    problem: PTInversionProblem,
    point: Sequence[float],
    step_fraction: float = 1.0e-3,
) -> SensitivityResult:
    """Central-difference P-T sensitivity of mineral-composition observables.

    This is the P-T counterpart to `numerical_sensitivity`: it perturbs the
    candidate pressure and temperature rather than database corrections.
    It applies only to composition/fraction observations. Pure-phase bracket
    inequalities are discontinuous by construction and should be diagnosed
    with the posterior/Laplace geometry instead of a finite-difference
    observable Jacobian.
    """
    if problem.brackets:
        raise ValueError(
            "numerical_pt_sensitivity requires composition/fraction observations, not brackets"
        )
    if len(point) != 2:
        raise ValueError("point must contain [pressure_bar, temperature_c]")
    if step_fraction <= 0.0:
        raise ValueError("step_fraction must be positive")

    center = np.asarray(point, dtype=float)
    baseline_experiment = problem._candidate_experiment(center)
    baseline_prediction = problem.forward_model.predict(
        baseline_experiment, problem.corrections
    )
    names, baseline = result_observables(baseline_experiment, baseline_prediction)
    steps = np.asarray(
        [
            max(abs(parameter.prior_sigma) * step_fraction, 1.0e-9)
            for parameter in problem.parameters
        ]
    )
    columns: list[np.ndarray] = []
    for index, step in enumerate(steps):
        plus_point = center.copy()
        minus_point = center.copy()
        plus_point[index] += step
        minus_point[index] -= step
        plus_experiment = problem._candidate_experiment(plus_point)
        minus_experiment = problem._candidate_experiment(minus_point)
        plus = result_observables(
            plus_experiment,
            problem.forward_model.predict(plus_experiment, problem.corrections),
        )[1]
        minus = result_observables(
            minus_experiment,
            problem.forward_model.predict(minus_experiment, problem.corrections),
        )[1]
        columns.append((plus - minus) / (2.0 * step))

    jacobian = np.column_stack(columns)
    identifiability = identify_sensitivity(
        jacobian, [parameter.name for parameter in problem.parameters]
    )
    return SensitivityResult(
        tuple(f"{problem.experiment.id}:{name}" for name in names),
        tuple(parameter.name for parameter in problem.parameters),
        baseline,
        jacobian,
        identifiability,
    )


def identify_sensitivity(
    jacobian: np.ndarray,
    parameter_names: Sequence[str],
    relative_tolerance: float = 1.0e-8,
) -> IdentifiabilityResult:
    """Compute rank and weak parameter directions from a sensitivity matrix."""
    matrix = np.asarray(jacobian, dtype=float)
    if matrix.ndim != 2:
        raise ValueError("jacobian must be two-dimensional")
    if matrix.shape[1] != len(parameter_names):
        raise ValueError("parameter_names must match jacobian columns")
    # full_matrices=True is required so V includes the null-space rows for
    # underdetermined systems (fewer observables than parameters) -- exactly
    # the sparse-data regime this diagnostic exists to catch. The economy SVD
    # (full_matrices=False) silently drops those rows, hiding genuinely
    # unconstrained parameter combinations.
    _, singular_values, vt = np.linalg.svd(matrix, full_matrices=True)
    scale = singular_values[0] if singular_values.size else 0.0
    threshold = relative_tolerance * scale
    rank = int(np.count_nonzero(singular_values > threshold)) if scale else 0
    condition_number = (
        float(singular_values[0] / singular_values[rank - 1])
        if rank and singular_values[rank - 1] > 0.0
        else float("inf")
    )
    weak_directions = tuple(
        {name: float(coefficient) for name, coefficient in zip(parameter_names, vector)}
        for vector in vt[rank:]
    )
    return IdentifiabilityResult(
        singular_values,
        rank,
        condition_number,
        tuple(parameter_names),
        weak_directions,
    )
