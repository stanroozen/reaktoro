"""Hierarchical group-level pressure and temperature calibration models."""

from collections.abc import Mapping, Sequence
from dataclasses import replace
from math import erf, log, sqrt

import numpy as np

from .data import (
    DirichletConcentrationParameter,
    Experiment,
    Parameter,
    ThermodynamicParameter,
)
from .forward_model import ReaktoroForwardModel
from .inversion import MAPResult
from .likelihood import MISFIT_FORMS, experiment_misfit
from .posterior import MultiChainMCMCResult, mcmc_sample_multi_chain


class GroupedPTCalibrationProblem:
    """Infer group P-T offsets and shared half-normal scale hyperparameters.

    The parameter order is all pressure offsets, all temperature offsets,
    pressure offset scale, and temperature offset scale. Group offsets have
    zero-centered Normal distributions conditional on their sampled scale;
    each scale has a fixed half-Normal hyperprior.
    """

    def __init__(
        self,
        forward_model: ReaktoroForwardModel,
        experiments: Sequence[Experiment],
        metadata_key: str,
        pressure_scale_hyperprior: float,
        temperature_scale_hyperprior: float,
        corrections: Mapping[str, float] | None = None,
        missing_phase_penalty: float = 100.0,
        extra_phase_penalty: float = 10.0,
        misfit_form: str = "weighted_lsq",
        theoretical_covariance_scale: float = 1.0,
        thermodynamic_parameters: Sequence[ThermodynamicParameter] = (),
        thermodynamic_prior_covariance: np.ndarray | None = None,
        dirichlet_concentration_parameters: Sequence[
            DirichletConcentrationParameter
        ] = (),
        observation_bias_keys: Sequence[str] = (),
        observation_bias_prior_sigma: float = 0.1,
        observation_bias_bounds: tuple[float, float] = (-1.0, 1.0),
        model_error_scale_prior_sigma: float | None = None,
        model_error_scale_bounds: tuple[float, float] = (1.0e-3, 5.0),
    ) -> None:
        if not metadata_key:
            raise ValueError("metadata_key must be non-empty")
        if not experiments:
            raise ValueError("at least one experiment is required")
        if any(
            not np.isfinite(value) or value <= 0.0
            for value in (pressure_scale_hyperprior, temperature_scale_hyperprior)
        ):
            raise ValueError(
                "calibration scale hyperpriors must be finite and positive"
            )
        if misfit_form not in MISFIT_FORMS:
            raise ValueError(f"misfit_form must be one of {MISFIT_FORMS}")
        if (
            not np.isfinite(theoretical_covariance_scale)
            or theoretical_covariance_scale < 0.0
        ):
            raise ValueError(
                "theoretical_covariance_scale must be finite and non-negative"
            )
        labels: list[str] = []
        for experiment in experiments:
            label = experiment.metadata.get(metadata_key)
            if not isinstance(label, str) or not label:
                raise ValueError(
                    f"Experiment {experiment.id!r} needs a non-empty {metadata_key!r} group"
                )
            labels.append(label)
        self.forward_model = forward_model
        self.experiments = tuple(experiments)
        self.metadata_key = metadata_key
        self.groups = tuple(sorted(set(labels)))
        self.pressure_scale_hyperprior = float(pressure_scale_hyperprior)
        self.temperature_scale_hyperprior = float(temperature_scale_hyperprior)
        self.corrections = dict(corrections or {})
        self.missing_phase_penalty = float(missing_phase_penalty)
        self.extra_phase_penalty = float(extra_phase_penalty)
        self.misfit_form = misfit_form
        self.theoretical_covariance_scale = float(theoretical_covariance_scale)
        bias_keys = tuple(observation_bias_keys)
        if len(set(bias_keys)) != len(bias_keys) or any(
            not isinstance(key, str) or not key for key in bias_keys
        ):
            raise ValueError("observation_bias_keys must be unique non-empty strings")
        if (
            not np.isfinite(observation_bias_prior_sigma)
            or observation_bias_prior_sigma <= 0.0
        ):
            raise ValueError("observation_bias_prior_sigma must be finite and positive")
        bias_lower, bias_upper = observation_bias_bounds
        if (
            not np.isfinite(bias_lower)
            or not np.isfinite(bias_upper)
            or bias_lower >= bias_upper
        ):
            raise ValueError("observation_bias_bounds must be finite and increasing")
        for experiment in experiments:
            available_keys = {
                *(f"fraction:{phase}" for phase in (experiment.phase_fractions or {})),
                *(
                    f"composition:{phase}:{component}"
                    for phase, composition in (
                        experiment.phase_compositions or {}
                    ).items()
                    for component in composition
                ),
            }
            if any(key not in available_keys for key in bias_keys):
                raise ValueError(
                    f"Experiment {experiment.id!r} lacks a requested observation-bias key"
                )
        if model_error_scale_prior_sigma is not None:
            scale_lower, scale_upper = model_error_scale_bounds
            if (
                not np.isfinite(model_error_scale_prior_sigma)
                or model_error_scale_prior_sigma <= 0.0
                or not np.isfinite(scale_lower)
                or not np.isfinite(scale_upper)
                or scale_lower <= 0.0
                or scale_lower >= scale_upper
            ):
                raise ValueError("model-error scale prior and bounds are invalid")
            if not any(
                experiment.theoretical_covariance is not None
                for experiment in experiments
            ):
                raise ValueError(
                    "model-error scale inference requires theoretical covariance data"
                )
        self.observation_bias_keys = bias_keys
        self.observation_bias_prior_sigma = float(observation_bias_prior_sigma)
        self.observation_bias_bounds = (float(bias_lower), float(bias_upper))
        self.model_error_scale_prior_sigma = model_error_scale_prior_sigma
        self.model_error_scale_bounds = tuple(
            float(value) for value in model_error_scale_bounds
        )
        self.model_error_scale_groups = (
            tuple(
                group
                for group in self.groups
                if any(
                    experiment.metadata[metadata_key] == group
                    and experiment.theoretical_covariance is not None
                    for experiment in experiments
                )
            )
            if model_error_scale_prior_sigma is not None
            else ()
        )
        self.observation_bias_labels = tuple(
            (key, group) for key in self.observation_bias_keys for group in self.groups
        )
        self.thermodynamic_parameters = tuple(thermodynamic_parameters)
        self.dirichlet_concentration_parameters = tuple(
            dirichlet_concentration_parameters
        )
        if self.dirichlet_concentration_parameters and misfit_form != "dirichlet_nll":
            raise ValueError(
                "Dirichlet concentration parameters require misfit_form='dirichlet_nll'"
            )
        concentration_phases = [
            parameter.phase for parameter in self.dirichlet_concentration_parameters
        ]
        if len(set(concentration_phases)) != len(concentration_phases):
            raise ValueError("Dirichlet concentration phases must be unique")
        observed_composition_phases = {
            phase
            for experiment in experiments
            for phase in (experiment.phase_compositions or {})
        }
        if any(
            phase not in observed_composition_phases for phase in concentration_phases
        ):
            raise ValueError(
                "each Dirichlet concentration parameter must target an observed phase composition"
            )
        thermo_keys = [parameter.key for parameter in self.thermodynamic_parameters]
        if len(set(thermo_keys)) != len(thermo_keys):
            raise ValueError("thermodynamic correction keys must be unique")
        if set(thermo_keys).intersection(self.corrections):
            raise ValueError(
                "inferred correction keys cannot also be fixed corrections"
            )
        if thermodynamic_prior_covariance is None:
            self.thermodynamic_prior_covariance = np.diag(
                [
                    parameter.prior_sigma**2
                    for parameter in self.thermodynamic_parameters
                ]
            )
        else:
            covariance = np.asarray(thermodynamic_prior_covariance, dtype=float)
            count = len(self.thermodynamic_parameters)
            if covariance.shape != (count, count) or not np.all(
                np.isfinite(covariance)
            ):
                raise ValueError(
                    "thermodynamic_prior_covariance must be finite and match parameter count"
                )
            covariance = 0.5 * (covariance + covariance.T)
            if count and np.any(np.linalg.eigvalsh(covariance) <= 0.0):
                raise ValueError(
                    "thermodynamic_prior_covariance must be positive definite"
                )
            if not count and covariance.size:
                raise ValueError("thermodynamic_prior_covariance requires parameters")
            self.thermodynamic_prior_covariance = covariance
        self._thermodynamic_parameter_count = len(self.thermodynamic_parameters)
        self._dirichlet_parameter_count = len(self.dirichlet_concentration_parameters)
        self._pressure_bias_start = (
            self._thermodynamic_parameter_count + self._dirichlet_parameter_count
        )
        self._observation_bias_start = self._pressure_bias_start + 2 * len(self.groups)
        self._model_error_scale_start = self._observation_bias_start + len(
            self.observation_bias_labels
        )
        self._pressure_scale_index = self._model_error_scale_start + len(
            self.model_error_scale_groups
        )
        self.parameters = [
            *self.thermodynamic_parameters,
            *self.dirichlet_concentration_parameters,
            *(
                Parameter(
                    f"pressure_bias[{group}]",
                    0.0,
                    pressure_scale_hyperprior,
                    -5.0 * pressure_scale_hyperprior,
                    5.0 * pressure_scale_hyperprior,
                    "bar",
                )
                for group in self.groups
            ),
            *(
                Parameter(
                    f"temperature_bias[{group}]",
                    0.0,
                    temperature_scale_hyperprior,
                    -5.0 * temperature_scale_hyperprior,
                    5.0 * temperature_scale_hyperprior,
                    "degC",
                )
                for group in self.groups
            ),
            *(
                Parameter(
                    f"bias[{key}|{group}]",
                    0.0,
                    observation_bias_prior_sigma,
                    bias_lower,
                    bias_upper,
                    "observation unit",
                )
                for key, group in self.observation_bias_labels
            ),
            *(
                Parameter(
                    f"scale_C_T[{group}]",
                    min(
                        max(
                            float(model_error_scale_prior_sigma),
                            self.model_error_scale_bounds[0],
                        ),
                        self.model_error_scale_bounds[1],
                    ),
                    float(model_error_scale_prior_sigma),
                    self.model_error_scale_bounds[0],
                    self.model_error_scale_bounds[1],
                    "1",
                )
                for group in self.model_error_scale_groups
            ),
            Parameter(
                "pressure_bias_scale",
                pressure_scale_hyperprior,
                pressure_scale_hyperprior,
                max(pressure_scale_hyperprior * 1.0e-3, 1.0e-8),
                5.0 * pressure_scale_hyperprior,
                "bar",
            ),
            Parameter(
                "temperature_bias_scale",
                temperature_scale_hyperprior,
                temperature_scale_hyperprior,
                max(temperature_scale_hyperprior * 1.0e-3, 1.0e-8),
                5.0 * temperature_scale_hyperprior,
                "degC",
            ),
        ]

    @property
    def initial_point(self) -> list[float]:
        return (
            [parameter.prior_mean for parameter in self.thermodynamic_parameters]
            + [
                parameter.prior_mean
                for parameter in self.dirichlet_concentration_parameters
            ]
            + [0.0] * (2 * len(self.groups))
            + [0.0] * len(self.observation_bias_labels)
            + [
                parameter.prior_mean
                for parameter in self.parameters[self._model_error_scale_start : -2]
            ]
            + [
                self.pressure_scale_hyperprior,
                self.temperature_scale_hyperprior,
            ]
        )

    def sample_unseen_group_offsets(
        self,
        pressure_scale: float,
        temperature_scale: float,
        rng: np.random.Generator,
        groups: Sequence[str] | None = None,
    ) -> dict[str, tuple[float, float]]:
        """Draw bounded predictive offsets for groups absent from the fit."""
        if (
            not np.isfinite(pressure_scale)
            or not np.isfinite(temperature_scale)
            or pressure_scale <= 0.0
            or temperature_scale <= 0.0
        ):
            raise ValueError(
                "predictive calibration scales must be finite and positive"
            )
        requested_groups = tuple(groups or ())
        if any(not isinstance(group, str) or not group for group in requested_groups):
            raise ValueError("predictive group labels must be non-empty strings")
        pressure = self._sample_truncated_normal(
            pressure_scale,
            5.0 * self.pressure_scale_hyperprior,
            len(requested_groups),
            rng,
        )
        temperature = self._sample_truncated_normal(
            temperature_scale,
            5.0 * self.temperature_scale_hyperprior,
            len(requested_groups),
            rng,
        )
        return {
            group: (float(pressure[index]), float(temperature[index]))
            for index, group in enumerate(requested_groups)
        }

    def corrections_for_values(self, values: Sequence[float]) -> dict[str, float]:
        """Combine fixed corrections with inferred thermodynamic parameters."""
        count = len(self.thermodynamic_parameters)
        if len(values) != len(self.parameters):
            raise ValueError("parameter values do not match grouped calibration model")
        return {
            **self.corrections,
            **{
                parameter.key: float(value)
                for parameter, value in zip(
                    self.thermodynamic_parameters, values[:count]
                )
            },
        }

    def dirichlet_concentrations_for_values(
        self, values: Sequence[float]
    ) -> dict[str, float]:
        if len(values) != len(self.parameters):
            raise ValueError("parameter values do not match grouped calibration model")
        start = self._thermodynamic_parameter_count
        return {
            parameter.phase: float(value)
            for parameter, value in zip(
                self.dirichlet_concentration_parameters,
                values[start : start + self._dirichlet_parameter_count],
            )
        }

    def observation_biases_for_values(
        self, values: Sequence[float]
    ) -> dict[tuple[str, str], float]:
        if len(values) != len(self.parameters):
            raise ValueError("parameter values do not match grouped calibration model")
        return {
            label: float(values[self._observation_bias_start + index])
            for index, label in enumerate(self.observation_bias_labels)
        }

    def theoretical_scales_for_values(
        self, values: Sequence[float]
    ) -> dict[str, float]:
        if len(values) != len(self.parameters):
            raise ValueError("parameter values do not match grouped calibration model")
        return {
            group: float(values[self._model_error_scale_start + index])
            for index, group in enumerate(self.model_error_scale_groups)
        }

    def sample_unseen_group_nuisances(
        self, groups: Sequence[str], rng: np.random.Generator
    ) -> tuple[dict[tuple[str, str], float], dict[str, float]]:
        """Draw bias and model-error nuisances from their priors for new groups."""
        requested = tuple(groups)
        if any(group not in self.groups for group in requested):
            raise ValueError("unseen groups must be known calibration groups")
        biases: dict[tuple[str, str], float] = {}
        if self.observation_bias_keys:
            lower, upper = self.observation_bias_bounds
            for key in self.observation_bias_keys:
                for group in requested:
                    biases[(key, group)] = self._sample_truncated_normal_interval(
                        self.observation_bias_prior_sigma,
                        lower,
                        upper,
                        1,
                        rng,
                    )[0]
        scales: dict[str, float] = {}
        if self.model_error_scale_prior_sigma is not None:
            lower, upper = self.model_error_scale_bounds
            eligible_groups = set(self.model_error_scale_groups)
            for group in requested:
                if group not in eligible_groups:
                    continue
                for _attempt in range(1000):
                    scale = abs(
                        float(rng.normal(0.0, self.model_error_scale_prior_sigma))
                    )
                    if lower <= scale <= upper:
                        scales[group] = scale
                        break
                else:
                    raise ValueError("could not draw a bounded model-error scale")
        return biases, scales

    @staticmethod
    def _sample_truncated_normal(
        scale: float,
        bound: float,
        count: int,
        rng: np.random.Generator,
    ) -> np.ndarray:
        draws = np.empty(count, dtype=float)
        pending = np.arange(count)
        while len(pending):
            proposed = rng.normal(0.0, scale, size=len(pending))
            accepted = np.abs(proposed) <= bound
            draws[pending[accepted]] = proposed[accepted]
            pending = pending[~accepted]
        return draws

    @staticmethod
    def _sample_truncated_normal_interval(
        sigma: float,
        lower: float,
        upper: float,
        count: int,
        rng: np.random.Generator,
    ) -> np.ndarray:
        draws = np.empty(count, dtype=float)
        pending = np.arange(count)
        while len(pending):
            proposed = rng.normal(0.0, sigma, size=len(pending))
            accepted = (proposed >= lower) & (proposed <= upper)
            draws[pending[accepted]] = proposed[accepted]
            pending = pending[~accepted]
        return draws

    def objective_components(self, values: Sequence[float]) -> tuple[float, float]:
        thermo_count = len(self.thermodynamic_parameters)
        expected = len(self.parameters)
        if len(values) != expected:
            raise ValueError(f"expected {expected} grouped calibration parameters")
        numeric = np.asarray(values, dtype=float)
        bounds = np.asarray([parameter.bounds() for parameter in self.parameters])
        if (
            not np.all(np.isfinite(numeric))
            or np.any(numeric < bounds[:, 0])
            or np.any(numeric > bounds[:, 1])
        ):
            return float("inf"), float("inf")

        group_index = {group: index for index, group in enumerate(self.groups)}
        corrections = self.corrections_for_values(numeric)
        dirichlet_concentrations = self.dirichlet_concentrations_for_values(numeric)
        pressure_biases = numeric[
            self._pressure_bias_start : self._pressure_bias_start + len(self.groups)
        ]
        temperature_start = self._pressure_bias_start + len(self.groups)
        temperature_biases = numeric[
            temperature_start : temperature_start + len(self.groups)
        ]
        observation_biases = self.observation_biases_for_values(numeric)
        theoretical_scales = self.theoretical_scales_for_values(numeric)
        pressure_scale, temperature_scale = numeric[-2:]
        prior = 0.0
        if thermo_count:
            thermo_deviation = numeric[:thermo_count] - np.asarray(
                [parameter.prior_mean for parameter in self.thermodynamic_parameters]
            )
            prior += 0.5 * float(
                thermo_deviation
                @ np.linalg.solve(self.thermodynamic_prior_covariance, thermo_deviation)
            )
        for parameter, value in zip(
            self.dirichlet_concentration_parameters,
            numeric[thermo_count : thermo_count + self._dirichlet_parameter_count],
        ):
            prior += 0.5 * ((value - parameter.prior_mean) / parameter.prior_sigma) ** 2
        data = 0.0
        for experiment in self.experiments:
            group = experiment.metadata[self.metadata_key]
            index = group_index[group]
            applied_biases = dict(experiment.observation_bias or {})
            for key in self.observation_bias_keys:
                label = (key, group)
                applied_biases[key] = (
                    applied_biases.get(key, 0.0) + observation_biases[label]
                )
            candidate = replace(
                experiment,
                pressure_bar=experiment.pressure_bar + pressure_biases[index],
                temperature_c=experiment.temperature_c + temperature_biases[index],
                observation_bias=applied_biases or None,
            )
            covariance_scale = theoretical_scales.get(
                group, self.theoretical_covariance_scale
            )
            data += experiment_misfit(
                candidate,
                self.forward_model.predict(candidate, corrections),
                self.missing_phase_penalty,
                self.extra_phase_penalty,
                self.misfit_form,
                covariance_scale,
                dirichlet_concentrations,
            )

        bias_values = numeric[
            self._observation_bias_start : self._model_error_scale_start
        ]
        if len(bias_values):
            prior += 0.5 * float(
                np.sum((bias_values / self.observation_bias_prior_sigma) ** 2)
            )
        if self.model_error_scale_groups:
            scale_values = numeric[
                self._model_error_scale_start : self._pressure_scale_index
            ]
            prior += 0.5 * float(
                np.sum((scale_values / float(self.model_error_scale_prior_sigma)) ** 2)
            )
        pressure_standardized = pressure_biases / pressure_scale
        temperature_standardized = temperature_biases / temperature_scale
        prior += 0.5 * float(pressure_standardized @ pressure_standardized)
        prior += len(self.groups) * np.log(
            pressure_scale / self.pressure_scale_hyperprior
        )
        pressure_mass = erf(
            5.0 * self.pressure_scale_hyperprior / (sqrt(2.0) * pressure_scale)
        )
        prior += len(self.groups) * log(pressure_mass)
        prior += 0.5 * float(temperature_standardized @ temperature_standardized)
        prior += len(self.groups) * np.log(
            temperature_scale / self.temperature_scale_hyperprior
        )
        temperature_mass = erf(
            5.0 * self.temperature_scale_hyperprior / (sqrt(2.0) * temperature_scale)
        )
        prior += len(self.groups) * log(temperature_mass)
        prior += 0.5 * (pressure_scale / self.pressure_scale_hyperprior) ** 2
        prior += 0.5 * (temperature_scale / self.temperature_scale_hyperprior) ** 2
        return float(data), float(prior)

    def objective(self, values: Sequence[float]) -> float:
        data, prior = self.objective_components(values)
        return data + prior

    def optimize(
        self,
        starts: int = 8,
        iterations: int = 120,
        seed: int = 0,
        central_point: Sequence[float] | None = None,
        hot_start_spread_fraction: float = 0.1,
    ) -> MAPResult:
        """Find a bounded MAP estimate using derivative-free pattern search."""
        if starts < 1 or iterations < 1:
            raise ValueError("starts and iterations must be positive")
        if not 0.0 < hot_start_spread_fraction <= 1.0:
            raise ValueError("hot_start_spread_fraction must be in (0, 1]")
        bounds = np.asarray([parameter.bounds() for parameter in self.parameters])
        prior_center = np.asarray(
            [parameter.prior_mean for parameter in self.parameters], dtype=float
        )
        if central_point is not None:
            if len(central_point) != len(self.parameters):
                raise ValueError("central_point must match parameter count")
            center = np.asarray(central_point, dtype=float)
            if (
                not np.all(np.isfinite(center))
                or np.any(center < bounds[:, 0])
                or np.any(center > bounds[:, 1])
            ):
                raise ValueError(
                    "central_point must be finite and within parameter bounds"
                )
        else:
            center = np.clip(prior_center, bounds[:, 0], bounds[:, 1])
        rng = np.random.default_rng(seed)
        span = bounds[:, 1] - bounds[:, 0]
        start_values = [center]
        start_values.extend(
            np.clip(
                center
                + (2.0 * rng.random(len(center)) - 1.0)
                * hot_start_spread_fraction
                * span,
                bounds[:, 0],
                bounds[:, 1],
            )
            for _ in range(starts - 1)
        )

        records: list[dict[str, object]] = []
        best_values = center.copy()
        best_objective = self.objective(best_values.tolist())
        for start in start_values:
            values = np.asarray(start, dtype=float)
            objective = self.objective(values.tolist())
            steps = span / 4.0
            for _ in range(iterations):
                improved = False
                for index in range(len(values)):
                    for direction in (-1.0, 1.0):
                        candidate = values.copy()
                        candidate[index] = np.clip(
                            candidate[index] + direction * steps[index],
                            bounds[index, 0],
                            bounds[index, 1],
                        )
                        candidate_objective = self.objective(candidate.tolist())
                        if candidate_objective < objective:
                            values, objective = candidate, candidate_objective
                            improved = True
                if not improved:
                    steps *= 0.5
                    if np.max(steps) <= 1.0e-8:
                        break
            records.append(
                {
                    "start": start.tolist(),
                    "map": values.tolist(),
                    "objective": float(objective),
                }
            )
            if objective < best_objective:
                best_values, best_objective = values.copy(), objective
        return MAPResult(
            {
                parameter.name: float(value)
                for parameter, value in zip(self.parameters, best_values)
            },
            float(best_objective),
            records,
        )


def mcmc_sample_grouped_pt_calibration(
    problem: GroupedPTCalibrationProblem,
    n_chains: int = 4,
    n_samples: int = 5000,
    burn_in: int = 500,
    initial_points: Sequence[Sequence[float]] | None = None,
    seeds: Sequence[int] | None = None,
    seed: int = 0,
    proposal_sigma: Sequence[float] | None = None,
) -> MultiChainMCMCResult:
    """Sample grouped calibration parameters and report multi-chain diagnostics.

    If starting points are omitted, generate reproducible in-bounds starts by
    perturbing the prior-center vector. Explicit starts are passed through to
    the common sampler, which validates their dimensions and bounds.
    """
    if initial_points is None:
        if n_chains < 2:
            raise ValueError("n_chains must be at least two")
        rng = np.random.default_rng(seed)
        center = np.asarray(problem.initial_point, dtype=float)
        bounds = np.asarray([parameter.bounds() for parameter in problem.parameters])
        scale = np.asarray([parameter.prior_sigma for parameter in problem.parameters])
        generated: list[list[float]] = []
        for _ in range(n_chains):
            point = center.copy()
            thermo_count = len(problem.thermodynamic_parameters)
            if thermo_count:
                thermo_mean = center[:thermo_count]
                thermo_factor = np.linalg.cholesky(
                    problem.thermodynamic_prior_covariance
                )
                for _attempt in range(1000):
                    candidate = thermo_mean + thermo_factor @ rng.normal(
                        size=thermo_count
                    )
                    if np.all(candidate >= bounds[:thermo_count, 0]) and np.all(
                        candidate <= bounds[:thermo_count, 1]
                    ):
                        point[:thermo_count] = candidate
                        break
                else:
                    raise ValueError(
                        "could not draw an in-bounds correlated thermodynamic start"
                    )
            for index in range(thermo_count, len(point)):
                for _attempt in range(100):
                    candidate = center[index] + rng.normal(0.0, 0.5 * scale[index])
                    if bounds[index, 0] <= candidate <= bounds[index, 1]:
                        point[index] = candidate
                        break
            if not np.isfinite(problem.objective(point.tolist())):
                raise ValueError("could not generate finite grouped calibration start")
            generated.append(point.tolist())
        initial_points = generated
    if len(initial_points) < 2:
        raise ValueError("at least two initial_points are required")
    if seeds is None:
        seeds = [seed + index for index in range(len(initial_points))]
    if len(seeds) != len(initial_points):
        raise ValueError("seeds must match initial_points count")
    return mcmc_sample_multi_chain(
        problem,
        initial_points,
        n_samples=n_samples,
        burn_in=burn_in,
        seeds=seeds,
        proposal_sigma=proposal_sigma,
    )
