"""Deterministic MAP inversion for selected thermodynamic corrections."""

from dataclasses import dataclass, replace
from math import inf
from typing import Any, Mapping, Protocol, Sequence

import numpy as np

from .data import (
    DirichletConcentrationParameter,
    EquilibriumBracket,
    Experiment,
    Parameter,
    ThermodynamicParameter,
)
from .forward_model import ReaktoroForwardModel
from .likelihood import MISFIT_FORMS, bracket_misfit, experiment_misfit


@dataclass
class MAPResult:
    corrections_j_mol: dict[str, float]
    objective: float
    starts: list[dict[str, object]]
    parameter_values: dict[str, float] | None = None


class PosteriorProblem(Protocol):
    """Structural contract shared by database, P-T, and joint inversions."""

    parameters: Any

    def objective_components(
        self, *args: Any, **kwargs: Any
    ) -> tuple[float, float]: ...

    def objective(self, *args: Any, **kwargs: Any) -> float: ...


class InversionProblem:
    """Small-dimensional, prior-regularized inversion around an existing database."""

    def __init__(
        self,
        forward_model: ReaktoroForwardModel,
        experiments: Sequence[Experiment],
        parameters: Sequence[
            ThermodynamicParameter | Parameter | DirichletConcentrationParameter
        ],
        brackets: Sequence[EquilibriumBracket] = (),
        missing_phase_penalty: float = 100.0,
        extra_phase_penalty: float = 10.0,
        marginalize_condition_uncertainty: bool = False,
        condition_quadrature_order: int = 3,
        objective_criterion: str = "map",
        misfit_form: str = "weighted_lsq",
        prior_covariance_j_mol2: np.ndarray | None = None,
        prior_covariance: np.ndarray | None = None,
        theoretical_covariance_scale: float = 1.0,
        theoretical_covariance_scale_metadata_key: str | None = None,
        theoretical_covariance_scales: Mapping[str, float] | None = None,
    ):
        if not parameters:
            raise ValueError("At least one inversion parameter is required")
        self.forward_model = forward_model
        self.experiments = list(experiments)
        self.parameters = list(parameters)
        self.brackets = list(brackets)
        keys = [parameter.key for parameter in self.parameters]
        if len(set(keys)) != len(keys):
            raise ValueError(
                "Each inversion parameter must target a unique (phase, field) pair"
            )
        self.missing_phase_penalty = missing_phase_penalty
        self.extra_phase_penalty = extra_phase_penalty
        if condition_quadrature_order < 1:
            raise ValueError("condition_quadrature_order must be positive")
        self.marginalize_condition_uncertainty = marginalize_condition_uncertainty
        self.condition_quadrature_order = condition_quadrature_order
        if objective_criterion not in ("map", "bayes"):
            raise ValueError("objective_criterion must be 'map' or 'bayes'")
        self.objective_criterion = objective_criterion
        if misfit_form not in MISFIT_FORMS:
            raise ValueError(f"misfit_form must be one of {MISFIT_FORMS}")
        self.misfit_form = misfit_form
        concentration_parameters = [
            parameter
            for parameter in self.parameters
            if isinstance(parameter, DirichletConcentrationParameter)
        ]
        if concentration_parameters and misfit_form != "dirichlet_nll":
            raise ValueError(
                "DirichletConcentrationParameter requires misfit_form='dirichlet_nll'"
            )
        observed_composition_phases = {
            phase
            for experiment in self.experiments
            for phase in (experiment.phase_compositions or {})
        }
        if any(
            parameter.phase not in observed_composition_phases
            for parameter in concentration_parameters
        ):
            raise ValueError(
                "each Dirichlet concentration parameter must target an observed phase composition"
            )
        if (
            not np.isfinite(theoretical_covariance_scale)
            or theoretical_covariance_scale < 0.0
        ):
            raise ValueError(
                "theoretical_covariance_scale must be finite and non-negative"
            )
        self.theoretical_covariance_scale = float(theoretical_covariance_scale)
        if (
            theoretical_covariance_scale_metadata_key
            and not theoretical_covariance_scales
        ):
            raise ValueError(
                "theoretical_covariance_scales are required with a metadata key"
            )
        if theoretical_covariance_scales is not None:
            if any(
                not isinstance(key, str)
                or not key
                or not np.isfinite(value)
                or value < 0.0
                for key, value in theoretical_covariance_scales.items()
            ):
                raise ValueError(
                    "theoretical covariance group scales must be finite and non-negative"
                )
        self.theoretical_covariance_scale_metadata_key = (
            theoretical_covariance_scale_metadata_key
        )
        self.theoretical_covariance_scales = dict(theoretical_covariance_scales or {})
        if theoretical_covariance_scale_metadata_key:
            missing_groups = {
                experiment.metadata.get(theoretical_covariance_scale_metadata_key)
                for experiment in self.experiments
                if experiment.theoretical_covariance is not None
            } - self.theoretical_covariance_scales.keys()
            if missing_groups:
                raise ValueError(
                    f"No theoretical covariance scale for groups {missing_groups}"
                )
        if prior_covariance is not None and prior_covariance_j_mol2 is not None:
            raise ValueError(
                "Specify prior_covariance or prior_covariance_j_mol2, not both"
            )
        covariance = (
            prior_covariance
            if prior_covariance is not None
            else prior_covariance_j_mol2
        )
        self.prior_covariance = self._prepare_prior_covariance(covariance)
        self.prior_covariance_j_mol2 = self.prior_covariance

    def _prepare_prior_covariance(self, covariance: np.ndarray | None) -> np.ndarray:
        if covariance is None:
            return np.diag([parameter.prior_sigma**2 for parameter in self.parameters])
        matrix = np.asarray(covariance, dtype=float)
        expected = len(self.parameters)
        if matrix.shape != (expected, expected):
            raise ValueError("prior_covariance must match parameter count")
        if not np.all(np.isfinite(matrix)):
            raise ValueError("prior_covariance must contain finite values")
        matrix = 0.5 * (matrix + matrix.T)
        eigenvalues = np.linalg.eigvalsh(matrix)
        if np.any(eigenvalues <= 0.0):
            raise ValueError("prior_covariance must be positive definite")
        return matrix

    def _condition_nodes(
        self, experiment: Experiment
    ) -> list[tuple[Experiment, float]]:
        """Return deterministic P-T nuisance nodes and normalized weights."""
        if not self.marginalize_condition_uncertainty:
            return [(experiment, 1.0)]
        if not experiment.pressure_sigma_bar and not experiment.temperature_sigma_c:
            return [(experiment, 1.0)]
        nodes, weights = np.polynomial.hermite.hermgauss(
            self.condition_quadrature_order
        )
        nodes = nodes / np.sqrt(2.0)
        weights = weights / np.sqrt(np.pi)
        pressure_sigma = experiment.pressure_sigma_bar or 0.0
        temperature_sigma = experiment.temperature_sigma_c or 0.0
        output: list[tuple[Experiment, float]] = []
        for pressure_node, pressure_weight in zip(nodes, weights):
            for temperature_node, temperature_weight in zip(nodes, weights):
                output.append(
                    (
                        replace(
                            experiment,
                            pressure_bar=experiment.pressure_bar
                            + pressure_sigma * pressure_node,
                            temperature_c=experiment.temperature_c
                            + temperature_sigma * temperature_node,
                        ),
                        float(pressure_weight * temperature_weight),
                    )
                )
        return output

    def _experiment_misfit_with_uncertainty(
        self,
        experiment: Experiment,
        corrections: dict[str, float],
        dirichlet_concentrations: dict[str, float],
    ) -> float:
        scale = self.theoretical_covariance_scale
        if self.theoretical_covariance_scale_metadata_key:
            group = experiment.metadata.get(
                self.theoretical_covariance_scale_metadata_key
            )
            if group not in self.theoretical_covariance_scales:
                raise ValueError(f"No theoretical covariance scale for group {group!r}")
            scale = self.theoretical_covariance_scales[group]
        return sum(
            weight
            * experiment_misfit(
                node,
                self.forward_model.predict(node, corrections),
                self.missing_phase_penalty,
                self.extra_phase_penalty,
                self.misfit_form,
                scale,
                dirichlet_concentrations,
            )
            for node, weight in self._condition_nodes(experiment)
        )

    def corrections_for_values(self, values: Sequence[float]) -> dict[str, float]:
        """Return only database-backed corrections from the full parameter vector."""
        if len(values) != len(self.parameters):
            raise ValueError("parameter values do not match inversion parameter count")
        return {
            parameter.key: float(value)
            for parameter, value in zip(self.parameters, values)
            if not isinstance(parameter, DirichletConcentrationParameter)
        }

    def dirichlet_concentrations_for_values(
        self, values: Sequence[float]
    ) -> dict[str, float]:
        """Return fitted per-phase concentrations from a parameter vector."""
        if len(values) != len(self.parameters):
            raise ValueError("parameter values do not match inversion parameter count")
        return {
            parameter.phase: float(value)
            for parameter, value in zip(self.parameters, values)
            if isinstance(parameter, DirichletConcentrationParameter)
        }

    def objective_components(
        self, values_j_mol: Sequence[float]
    ) -> tuple[float, float]:
        """Return (data_misfit, prior_penalty) separately, so callers can
        evaluate an uncertainty-source-restricted objective (data only,
        prior only, or both -- see `posterior.uncertainty_source`)."""
        if len(values_j_mol) != len(self.parameters):
            raise ValueError("Correction vector does not match parameter count")
        values = np.asarray(values_j_mol, dtype=float)
        bounds = np.asarray([parameter.bounds() for parameter in self.parameters])
        if np.any(values < bounds[:, 0]) or np.any(values > bounds[:, 1]):
            return inf, inf
        corrections = self.corrections_for_values(values)
        dirichlet_concentrations = self.dirichlet_concentrations_for_values(values)
        prior = 0.0
        for parameter in self.parameters:
            if parameter.prior_sigma <= 0.0:
                raise ValueError(f"Prior sigma must be positive: {parameter.name}")
        deviations = np.asarray(
            [
                value - parameter.prior_mean
                for parameter, value in zip(self.parameters, values)
            ],
            dtype=float,
        )
        prior = 0.5 * deviations @ np.linalg.solve(self.prior_covariance, deviations)
        data = 0.0
        for experiment in self.experiments:
            data += self._experiment_misfit_with_uncertainty(
                experiment, corrections, dirichlet_concentrations
            )
        for bracket in self.brackets:
            lower_nodes = self._condition_nodes(bracket.lower)
            upper_nodes = self._condition_nodes(bracket.upper)
            data += sum(
                lower_weight
                * upper_weight
                * bracket_misfit(
                    bracket,
                    self.forward_model.predict(lower_node, corrections),
                    self.forward_model.predict(upper_node, corrections),
                )
                for lower_node, lower_weight in lower_nodes
                for upper_node, upper_weight in upper_nodes
            )
        return data, prior

    def bayes_penalty(self) -> float:
        """BIC-style model-complexity penalty: 0.5 * n_parameters * ln(n_observations).

        ``n_observations`` is approximated as the number of experiments plus
        brackets (each treated as one data point), not the finer-grained count
        of individual residual terms. This mirrors mc_fit's ln(Bayes)
        criterion in spirit (penalize extra free parameters relative to the
        amount of independent data) but is a coarser approximation.
        """
        n_observations = max(len(self.experiments) + len(self.brackets), 2)
        return 0.5 * len(self.parameters) * float(np.log(n_observations))

    def objective(self, values_j_mol: Sequence[float]) -> float:
        data, prior = self.objective_components(values_j_mol)
        total = data + prior
        if self.objective_criterion == "bayes":
            total += self.bayes_penalty()
        return total

    def optimize_multi_parameter(
        self,
        starts: int = 8,
        iterations: int = 100,
        seed: int = 0,
        central_point: Sequence[float] | None = None,
        hot_start_spread_fraction: float = 0.1,
    ) -> MAPResult:
        """Run deterministic bounded coordinate-pattern searches from multiple starts.

        This derivative-free method is intended for the small, nonsmooth parameter
        spaces produced by phase-equilibrium calculations. It is not a substitute
        for posterior sampling and records all starts for later diagnostics.

        By default (``central_point=None``) this is a "cold start": one search
        begins at the prior mean and the rest are drawn uniformly across each
        parameter's full bounds. Passing ``central_point`` (e.g. a previous
        MAP result) switches to a "hot start": every search instead begins
        near that point, perturbed by up to ``hot_start_spread_fraction`` of
        each parameter's bound width, mirroring mc_fit's hot/cold start modes.
        """
        if starts < 1 or iterations < 1:
            raise ValueError("starts and iterations must be positive")
        if not 0.0 < hot_start_spread_fraction <= 1.0:
            raise ValueError("hot_start_spread_fraction must be in (0, 1]")
        bounds = np.asarray([parameter.bounds() for parameter in self.parameters])
        centers = np.asarray(
            [parameter.prior_mean for parameter in self.parameters], dtype=float
        )
        rng = np.random.default_rng(seed)
        if central_point is not None:
            if len(central_point) != len(self.parameters):
                raise ValueError("central_point does not match parameter count")
            center = np.clip(
                np.asarray(central_point, dtype=float), bounds[:, 0], bounds[:, 1]
            )
            spread = hot_start_spread_fraction * (bounds[:, 1] - bounds[:, 0])
            start_values = [center]
            start_values.extend(
                np.clip(
                    center + (rng.random(len(self.parameters)) * 2.0 - 1.0) * spread,
                    bounds[:, 0],
                    bounds[:, 1],
                )
                for _ in range(starts - 1)
            )
        else:
            start_values = [np.clip(centers, bounds[:, 0], bounds[:, 1])]
            start_values.extend(
                bounds[:, 0]
                + rng.random(len(self.parameters)) * (bounds[:, 1] - bounds[:, 0])
                for _ in range(starts - 1)
            )

        records: list[dict[str, object]] = []
        best_values = centers.copy()
        best_objective = self.objective(best_values)
        for start in start_values:
            values = np.asarray(start, dtype=float)
            objective = self.objective(values)
            steps = (bounds[:, 1] - bounds[:, 0]) / 4.0
            evaluations = 1
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
                        candidate_objective = self.objective(candidate)
                        evaluations += 1
                        if candidate_objective < objective:
                            values, objective = candidate, candidate_objective
                            improved = True
                if not improved:
                    steps *= 0.5
                    if np.max(steps) <= 1.0e-6:
                        break
            record = {
                "start_j_mol": start.tolist(),
                "map_j_mol": values.tolist(),
                "objective": float(objective),
                "evaluations": evaluations,
            }
            records.append(record)
            if objective < best_objective:
                best_values, best_objective = values.copy(), objective
        return MAPResult(
            self.corrections_for_values(best_values),
            float(best_objective),
            records,
            {
                parameter.name: float(value)
                for parameter, value in zip(self.parameters, best_values)
            },
        )

    def optimize_one_parameter(
        self,
        parameter_index: int = 0,
        starts: int = 5,
        iterations: int = 60,
    ) -> MAPResult:
        """Bounded golden-section searches from prior-spread starting intervals."""
        if len(self.parameters) != 1:
            raise ValueError("optimize_one_parameter requires exactly one parameter")
        parameter = self.parameters[parameter_index]
        lower, upper = parameter.bounds()
        start_values = np.linspace(lower, upper, max(starts, 1))
        records: list[dict[str, object]] = []
        best_value, best_objective = parameter.prior_mean, inf
        phi = (1.0 + np.sqrt(5.0)) / 2.0
        for start in start_values:
            left, right = lower, upper
            for _ in range(iterations):
                first = right - (right - left) / phi
                second = left + (right - left) / phi
                if self.objective([first]) < self.objective([second]):
                    right = second
                else:
                    left = first
            value = (left + right) / 2.0
            objective = self.objective([value])
            records.append(
                {
                    "start_j_mol": float(start),
                    "map_j_mol": value,
                    "objective": objective,
                }
            )
            if objective < best_objective:
                best_value, best_objective = value, objective
        return MAPResult(
            self.corrections_for_values([best_value]),
            best_objective,
            records,
            {parameter.name: best_value},
        )

    def sample_theoretical_covariance_scale(self, initial_corrections, **kwargs):
        """Sample corrections and a global model-error scale for this problem."""
        from .hierarchical import mcmc_sample_theoretical_covariance_scale

        return mcmc_sample_theoretical_covariance_scale(
            self, initial_corrections, **kwargs
        )

    def sample_theoretical_covariance_group_scales(
        self, metadata_key: str, initial_corrections, initial_scales, **kwargs
    ):
        """Sample corrections and one model-error scale per metadata group."""
        from .hierarchical import mcmc_sample_theoretical_covariance_group_scales

        return mcmc_sample_theoretical_covariance_group_scales(
            self, metadata_key, initial_corrections, initial_scales, **kwargs
        )

    def sample_combined_grouped_uncertainty(
        self,
        metadata_key: str,
        observation_key,
        initial_corrections,
        initial_scales,
        initial_biases,
        **kwargs,
    ):
        """Sample corrections, grouped model scales, and grouped biases jointly."""
        from .hierarchical import mcmc_sample_combined_grouped_uncertainty

        return mcmc_sample_combined_grouped_uncertainty(
            self,
            metadata_key,
            observation_key,
            initial_corrections,
            initial_scales,
            initial_biases,
            **kwargs,
        )

    def sample_group_observation_biases(
        self,
        metadata_key: str,
        observation_key: str,
        initial_corrections,
        initial_biases,
        **kwargs,
    ):
        """Sample corrections and one continuous-observation bias per group."""
        from .hierarchical import mcmc_sample_group_observation_biases

        return mcmc_sample_group_observation_biases(
            self,
            metadata_key,
            observation_key,
            initial_corrections,
            initial_biases,
            **kwargs,
        )

    def sample_group_observation_biases_multi_key(
        self,
        metadata_key: str,
        observation_keys,
        initial_corrections,
        initial_biases,
        **kwargs,
    ):
        """Sample corrections and grouped biases for several observations."""
        from .hierarchical import mcmc_sample_group_observation_biases_multi_key

        return mcmc_sample_group_observation_biases_multi_key(
            self,
            metadata_key,
            observation_keys,
            initial_corrections,
            initial_biases,
            **kwargs,
        )
