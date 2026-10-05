"""Invert pressure and temperature from a sample's observed mineral phases
and compositions (classic geothermobarometry), reusing the same likelihood,
prior, and posterior-diagnostic machinery as `InversionProblem`.

`InversionProblem` perturbs *database parameters* at fixed P-T (e.g. an
enthalpy correction). `PTInversionProblem` is the complementary case: it
perturbs the *experiment's own P-T* at a fixed database, scoring candidates
against the experiment's observed phases/compositions via the same
`likelihood.experiment_misfit`. Because it exposes the same duck-typed
surface (`.parameters`, `.objective(values)`, `.objective_components(values)`),
`posterior.laplace_covariance`/`mcmc_sample`/`grid_posterior` and
`identifiability.identifiability_summary`/`information_gain_summary` all work
on it unchanged. Two utilities do NOT: `sensitivity.numerical_sensitivity`
and `validation.posterior_predictive_check` both build corrections dicts
keyed by `parameter.key` (a database species/field), which has no meaning
here -- structural identifiability should instead be read from the Laplace
covariance's correlation matrix (see `laplace_covariance`).
"""

from collections.abc import Mapping, Sequence
from dataclasses import replace as dataclass_replace

import numpy as np

from .data import (
    EquilibriumBracket,
    Experiment,
    Parameter,
    ThermodynamicParameter,
)
from .forward_model import ReaktoroForwardModel
from .inversion import MAPResult
from .likelihood import MISFIT_FORMS, bracket_misfit, experiment_misfit


class PTInversionProblem:
    """MAP/posterior inversion of one experiment's pressure and temperature."""

    def __init__(
        self,
        forward_model: ReaktoroForwardModel,
        experiment: Experiment,
        pressure_prior: ThermodynamicParameter | Parameter,
        temperature_prior: ThermodynamicParameter | Parameter,
        bracket: EquilibriumBracket | None = None,
        brackets: Sequence[EquilibriumBracket] = (),
        corrections_j_mol: Mapping[str, float] | None = None,
        missing_phase_penalty: float = 100.0,
        extra_phase_penalty: float = 10.0,
        misfit_form: str = "weighted_lsq",
        prior_covariance_j_mol2: np.ndarray | None = None,
        prior_covariance: np.ndarray | None = None,
        theoretical_covariance_scale: float = 1.0,
    ):
        if misfit_form not in MISFIT_FORMS:
            raise ValueError(f"misfit_form must be one of {MISFIT_FORMS}")
        self.forward_model = forward_model
        self.experiment = experiment
        if bracket is not None and brackets:
            raise ValueError("Use either bracket or brackets, not both")
        self.brackets = list(brackets) or ([bracket] if bracket is not None else [])
        self.bracket = (
            bracket
            if bracket is not None
            else (self.brackets[0] if len(self.brackets) == 1 else None)
        )
        # Order is fixed: [pressure, temperature]; posterior.py only relies on
        # positional correspondence between `self.parameters` and value vectors.
        self.parameters = [pressure_prior, temperature_prior]
        if prior_covariance is not None and prior_covariance_j_mol2 is not None:
            raise ValueError(
                "Specify prior_covariance or prior_covariance_j_mol2, not both"
            )
        covariance_input = (
            prior_covariance
            if prior_covariance is not None
            else prior_covariance_j_mol2
        )
        if covariance_input is None:
            self.prior_covariance_j_mol2 = np.diag(
                [parameter.prior_sigma**2 for parameter in self.parameters]
            )
        else:
            covariance = np.asarray(covariance_input, dtype=float)
            if covariance.shape != (2, 2) or not np.all(np.isfinite(covariance)):
                raise ValueError("prior_covariance_j_mol2 must be a finite 2x2 matrix")
            covariance = 0.5 * (covariance + covariance.T)
            if np.any(np.linalg.eigvalsh(covariance) <= 0.0):
                raise ValueError("prior_covariance_j_mol2 must be positive definite")
            self.prior_covariance_j_mol2 = covariance
        self.prior_covariance = self.prior_covariance_j_mol2
        self.corrections = dict(corrections_j_mol or {})
        self.missing_phase_penalty = missing_phase_penalty
        self.extra_phase_penalty = extra_phase_penalty
        self.misfit_form = misfit_form
        if (
            not np.isfinite(theoretical_covariance_scale)
            or theoretical_covariance_scale < 0.0
        ):
            raise ValueError(
                "theoretical_covariance_scale must be finite and non-negative"
            )
        self.theoretical_covariance_scale = float(theoretical_covariance_scale)

    def _candidate_experiment(self, values: Sequence[float]) -> Experiment:
        if len(values) != 2:
            raise ValueError("Expected exactly [pressure_bar, temperature_c]")
        pressure_bar, temperature_c = values
        return dataclass_replace(
            self.experiment,
            pressure_bar=float(pressure_bar),
            temperature_c=float(temperature_c),
        )

    def _candidate_bracket(self, values: Sequence[float]) -> EquilibriumBracket:
        """Translate a bracket's endpoint geometry around a candidate P-T center."""
        if self.bracket is None:
            raise ValueError("No equilibrium bracket was supplied")
        if len(values) != 2:
            raise ValueError("Expected exactly [pressure_bar, temperature_c]")
        center_pressure, center_temperature = values
        reference_pressure = (
            self.bracket.lower.pressure_bar + self.bracket.upper.pressure_bar
        ) / 2.0
        reference_temperature = (
            self.bracket.lower.temperature_c + self.bracket.upper.temperature_c
        ) / 2.0

        def translate(endpoint: Experiment) -> Experiment:
            return dataclass_replace(
                endpoint,
                pressure_bar=float(
                    center_pressure + endpoint.pressure_bar - reference_pressure
                ),
                temperature_c=float(
                    center_temperature + endpoint.temperature_c - reference_temperature
                ),
            )

        return EquilibriumBracket(
            self.bracket.id,
            translate(self.bracket.lower),
            translate(self.bracket.upper),
            self.bracket.lower_required_phases,
            self.bracket.upper_required_phases,
            self.bracket.lower_forbidden_phases,
            self.bracket.upper_forbidden_phases,
        )

    def _translate_bracket(
        self, bracket: EquilibriumBracket, values: Sequence[float]
    ) -> EquilibriumBracket:
        """Translate one bracket's endpoint geometry around a candidate P-T center."""
        center_pressure, center_temperature = values
        reference_pressure = (
            bracket.lower.pressure_bar + bracket.upper.pressure_bar
        ) / 2.0
        reference_temperature = (
            bracket.lower.temperature_c + bracket.upper.temperature_c
        ) / 2.0

        def translate(endpoint: Experiment) -> Experiment:
            return dataclass_replace(
                endpoint,
                pressure_bar=float(
                    center_pressure + endpoint.pressure_bar - reference_pressure
                ),
                temperature_c=float(
                    center_temperature + endpoint.temperature_c - reference_temperature
                ),
            )

        return EquilibriumBracket(
            bracket.id,
            translate(bracket.lower),
            translate(bracket.upper),
            bracket.lower_required_phases,
            bracket.upper_required_phases,
            bracket.lower_forbidden_phases,
            bracket.upper_forbidden_phases,
        )

    def objective_components(self, values: Sequence[float]) -> tuple[float, float]:
        if len(values) != 2:
            raise ValueError("Expected exactly [pressure_bar, temperature_c]")
        numeric_values = np.asarray(values, dtype=float)
        if not np.all(np.isfinite(numeric_values)):
            return float("inf"), float("inf")
        bounds = np.asarray([parameter.bounds() for parameter in self.parameters])
        if np.any(numeric_values < bounds[:, 0]) or np.any(
            numeric_values > bounds[:, 1]
        ):
            return float("inf"), float("inf")
        values_list: Sequence[float] = numeric_values.tolist()
        prior = 0.0
        for parameter, value in zip(self.parameters, values_list):
            if parameter.prior_sigma <= 0.0:
                raise ValueError(f"Prior sigma must be positive: {parameter.name}")
        deviations = numeric_values - np.asarray(
            [parameter.prior_mean for parameter in self.parameters]
        )
        prior = 0.5 * deviations @ np.linalg.solve(self.prior_covariance, deviations)
        if self.brackets:
            data = 0.0
            for bracket in self.brackets:
                candidate_bracket = self._translate_bracket(bracket, values)
                data += bracket_misfit(
                    candidate_bracket,
                    self.forward_model.predict(
                        candidate_bracket.lower, self.corrections
                    ),
                    self.forward_model.predict(
                        candidate_bracket.upper, self.corrections
                    ),
                    self.missing_phase_penalty,
                )
        else:
            candidate = self._candidate_experiment(values)
            prediction = self.forward_model.predict(candidate, self.corrections)
            data = experiment_misfit(
                candidate,
                prediction,
                self.missing_phase_penalty,
                self.extra_phase_penalty,
                self.misfit_form,
                self.theoretical_covariance_scale,
            )
        return float(data), float(prior)

    def objective(self, values: Sequence[float]) -> float:
        data, prior = self.objective_components(values)
        return data + prior

    def optimize(
        self,
        starts: int = 8,
        iterations: int = 100,
        seed: int = 0,
        central_point: Sequence[float] | None = None,
        hot_start_spread_fraction: float = 0.1,
    ) -> MAPResult:
        """Bounded coordinate-pattern search over [pressure_bar, temperature_c].

        Same derivative-free algorithm as `InversionProblem.optimize_multi_parameter`,
        kept as a separate small implementation here rather than shared, since the
        result is keyed by ``"pressure_bar"``/``"temperature_c"`` instead of a
        database species key.
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
        if central_point is None:
            start_values = [np.clip(centers, bounds[:, 0], bounds[:, 1])]
            start_values.extend(
                bounds[:, 0] + rng.random(2) * (bounds[:, 1] - bounds[:, 0])
                for _ in range(starts - 1)
            )
        else:
            if len(central_point) != 2:
                raise ValueError(
                    "central_point must contain [pressure_bar, temperature_c]"
                )
            center = np.clip(
                np.asarray(central_point, dtype=float), bounds[:, 0], bounds[:, 1]
            )
            spread = hot_start_spread_fraction * (bounds[:, 1] - bounds[:, 0])
            start_values = [center]
            start_values.extend(
                np.clip(
                    center + (rng.random(2) * 2.0 - 1.0) * spread,
                    bounds[:, 0],
                    bounds[:, 1],
                )
                for _ in range(starts - 1)
            )

        records: list[dict[str, object]] = []
        best_values = centers.copy()
        best_objective = self.objective(best_values.tolist())
        for start in start_values:
            values = np.asarray(start, dtype=float)
            objective = self.objective(values.tolist())
            steps = (bounds[:, 1] - bounds[:, 0]) / 4.0
            for _ in range(iterations):
                improved = False
                for index in range(2):
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
                    if np.max(steps) <= 1.0e-6:
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
                "pressure_bar": float(best_values[0]),
                "temperature_c": float(best_values[1]),
            },
            float(best_objective),
            records,
        )


class PTCalibrationBiasProblem:
    """P-T inversion with inferred pressure and temperature calibration offsets."""

    def __init__(
        self,
        base_problem: PTInversionProblem,
        pressure_bias_prior: Parameter,
        temperature_bias_prior: Parameter,
        prior_covariance: np.ndarray | None = None,
    ) -> None:
        if pressure_bias_prior.unit != "bar":
            raise ValueError("pressure bias prior must use unit 'bar'")
        if temperature_bias_prior.unit not in {"degC", "C", "celsius"}:
            raise ValueError("temperature bias prior must use a Celsius unit")
        self.base_problem = base_problem
        self.forward_model = base_problem.forward_model
        self.experiment = base_problem.experiment
        self.brackets = base_problem.brackets
        self.parameters = [
            *base_problem.parameters,
            pressure_bias_prior,
            temperature_bias_prior,
        ]
        if len({parameter.name for parameter in self.parameters}) != 4:
            raise ValueError("P-T calibration parameter names must be unique")
        if prior_covariance is None:
            covariance = np.diag(
                [parameter.prior_sigma**2 for parameter in self.parameters]
            )
            covariance[:2, :2] = base_problem.prior_covariance
            self.prior_covariance = covariance
        else:
            covariance = np.asarray(prior_covariance, dtype=float)
            if covariance.shape != (4, 4) or not np.all(np.isfinite(covariance)):
                raise ValueError("prior_covariance must be a finite 4x4 matrix")
            covariance = 0.5 * (covariance + covariance.T)
            if np.any(np.linalg.eigvalsh(covariance) <= 0.0):
                raise ValueError("prior_covariance must be positive definite")
            self.prior_covariance = covariance

    def objective_components(self, values: Sequence[float]) -> tuple[float, float]:
        if len(values) != 4:
            raise ValueError(
                "Expected [pressure, temperature, pressure_bias, temperature_bias]"
            )
        numeric = np.asarray(values, dtype=float)
        bounds = np.asarray([parameter.bounds() for parameter in self.parameters])
        if (
            not np.all(np.isfinite(numeric))
            or np.any(numeric < bounds[:, 0])
            or np.any(numeric > bounds[:, 1])
        ):
            return float("inf"), float("inf")
        corrected = numeric[:2] + numeric[2:]
        if self.brackets:
            data = 0.0
            for bracket in self.brackets:
                candidate_bracket = self.base_problem._translate_bracket(
                    bracket, corrected
                )
                data += bracket_misfit(
                    candidate_bracket,
                    self.forward_model.predict(
                        candidate_bracket.lower, self.base_problem.corrections
                    ),
                    self.forward_model.predict(
                        candidate_bracket.upper, self.base_problem.corrections
                    ),
                    self.base_problem.missing_phase_penalty,
                )
        else:
            candidate = self.base_problem._candidate_experiment(corrected)
            data = experiment_misfit(
                candidate,
                self.forward_model.predict(candidate, self.base_problem.corrections),
                self.base_problem.missing_phase_penalty,
                self.base_problem.extra_phase_penalty,
                self.base_problem.misfit_form,
                self.base_problem.theoretical_covariance_scale,
            )
        means = np.asarray([parameter.prior_mean for parameter in self.parameters])
        deviation = numeric - means
        prior = 0.5 * deviation @ np.linalg.solve(self.prior_covariance, deviation)
        return float(data), float(prior)

    def objective(self, values: Sequence[float]) -> float:
        data, prior = self.objective_components(values)
        return data + prior

    def optimize(
        self, starts: int = 8, iterations: int = 120, seed: int = 0
    ) -> MAPResult:
        """Run bounded coordinate-pattern MAP search over P-T and both biases."""
        if starts < 1 or iterations < 1:
            raise ValueError("starts and iterations must be positive")
        bounds = np.asarray([parameter.bounds() for parameter in self.parameters])
        centers = np.asarray([parameter.prior_mean for parameter in self.parameters])
        rng = np.random.default_rng(seed)
        start_values = [np.clip(centers, bounds[:, 0], bounds[:, 1])]
        start_values.extend(
            bounds[:, 0] + rng.random(4) * (bounds[:, 1] - bounds[:, 0])
            for _ in range(starts - 1)
        )
        best_values = start_values[0].copy()
        best_objective = self.objective(best_values.tolist())
        records: list[dict[str, object]] = []
        for start in start_values:
            values = start.copy()
            objective = self.objective(values.tolist())
            steps = (bounds[:, 1] - bounds[:, 0]) / 4.0
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
                    if np.max(steps) <= 1.0e-6:
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
