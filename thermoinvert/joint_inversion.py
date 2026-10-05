"""Joint thermodynamic-correction and pressure-temperature inversion."""

from collections.abc import Sequence
from dataclasses import replace
from math import inf

import numpy as np

from .data import Experiment, Parameter, ThermodynamicParameter
from .forward_model import ReaktoroForwardModel
from .inversion import MAPResult
from .likelihood import MISFIT_FORMS, experiment_misfit


class JointInversionProblem:
    """Infer database corrections and one experiment's pressure/temperature jointly."""

    def __init__(
        self,
        forward_model: ReaktoroForwardModel,
        experiment: Experiment | Sequence[Experiment],
        thermodynamic_parameters: Sequence[ThermodynamicParameter],
        pressure_parameter: Parameter | Sequence[Parameter],
        temperature_parameter: Parameter | Sequence[Parameter],
        missing_phase_penalty: float = 100.0,
        extra_phase_penalty: float = 10.0,
        misfit_form: str = "weighted_lsq",
        prior_covariance: np.ndarray | None = None,
        theoretical_covariance_scale: float = 1.0,
    ):
        if not thermodynamic_parameters:
            raise ValueError("At least one thermodynamic parameter is required")
        if misfit_form not in MISFIT_FORMS:
            raise ValueError(f"misfit_form must be one of {MISFIT_FORMS}")
        experiments = (
            (experiment,) if isinstance(experiment, Experiment) else tuple(experiment)
        )
        if not experiments:
            raise ValueError("At least one experiment is required")
        if len({item.id for item in experiments}) != len(experiments):
            raise ValueError("joint experiment IDs must be unique")
        pressure_parameters = (
            (pressure_parameter,)
            if isinstance(pressure_parameter, Parameter)
            else tuple(pressure_parameter)
        )
        temperature_parameters = (
            (temperature_parameter,)
            if isinstance(temperature_parameter, Parameter)
            else tuple(temperature_parameter)
        )
        if len(pressure_parameters) != len(experiments) or len(
            temperature_parameters
        ) != len(experiments):
            raise ValueError(
                "provide one pressure and temperature parameter per experiment"
            )
        if any(parameter.unit != "bar" for parameter in pressure_parameters):
            raise ValueError("pressure parameters must use unit 'bar'")
        if any(
            parameter.unit not in {"degC", "C", "celsius"}
            for parameter in temperature_parameters
        ):
            raise ValueError("temperature parameters must use a Celsius unit")
        parameters: list[Parameter | ThermodynamicParameter] = [
            *thermodynamic_parameters
        ]
        for pressure, temperature in zip(pressure_parameters, temperature_parameters):
            parameters.extend((pressure, temperature))
        names = [parameter.name for parameter in parameters]
        if len(set(names)) != len(names):
            raise ValueError("joint parameter names must be unique")
        keys = [parameter.key for parameter in thermodynamic_parameters]
        if len(set(keys)) != len(keys):
            raise ValueError("thermodynamic correction keys must be unique")
        if (
            not np.isfinite(theoretical_covariance_scale)
            or theoretical_covariance_scale < 0.0
        ):
            raise ValueError(
                "theoretical_covariance_scale must be finite and non-negative"
            )
        self.forward_model = forward_model
        self.experiments = experiments
        self.experiment = experiments[0] if len(experiments) == 1 else experiments
        self.thermodynamic_parameters = list(thermodynamic_parameters)
        self.pressure_parameters = pressure_parameters
        self.temperature_parameters = temperature_parameters
        self.pressure_parameter = (
            pressure_parameters[0] if len(experiments) == 1 else pressure_parameters
        )
        self.temperature_parameter = (
            temperature_parameters[0]
            if len(experiments) == 1
            else temperature_parameters
        )
        self.parameters = parameters
        self.missing_phase_penalty = float(missing_phase_penalty)
        self.extra_phase_penalty = float(extra_phase_penalty)
        self.misfit_form = misfit_form
        self.theoretical_covariance_scale = float(theoretical_covariance_scale)
        self.prior_covariance = self._prepare_prior_covariance(prior_covariance)
        self.prior_covariance_j_mol2 = self.prior_covariance

    def _prepare_prior_covariance(self, covariance: np.ndarray | None) -> np.ndarray:
        if covariance is None:
            return np.diag([parameter.prior_sigma**2 for parameter in self.parameters])
        matrix = np.asarray(covariance, dtype=float)
        expected = len(self.parameters)
        if matrix.shape != (expected, expected) or not np.all(np.isfinite(matrix)):
            raise ValueError("prior_covariance must be a finite square matrix")
        matrix = 0.5 * (matrix + matrix.T)
        if np.any(np.linalg.eigvalsh(matrix) <= 0.0):
            raise ValueError("prior_covariance must be positive definite")
        return matrix

    def _candidate_experiment(
        self, experiment: Experiment, pressure: float, temperature: float
    ) -> Experiment:
        return replace(
            experiment,
            pressure_bar=float(pressure),
            temperature_c=float(temperature),
        )

    def objective_components(self, values: Sequence[float]) -> tuple[float, float]:
        if len(values) != len(self.parameters):
            raise ValueError("joint parameter vector does not match parameter count")
        numeric = np.asarray(values, dtype=float)
        bounds = np.asarray([parameter.bounds() for parameter in self.parameters])
        if (
            not np.all(np.isfinite(numeric))
            or np.any(numeric < bounds[:, 0])
            or np.any(numeric > bounds[:, 1])
        ):
            return inf, inf
        corrections = {
            parameter.key: float(value)
            for parameter, value in zip(
                self.thermodynamic_parameters,
                numeric[: len(self.thermodynamic_parameters)],
            )
        }
        data = 0.0
        offset = len(self.thermodynamic_parameters)
        for index, experiment in enumerate(self.experiments):
            pressure = numeric[offset + 2 * index]
            temperature = numeric[offset + 2 * index + 1]
            candidate = self._candidate_experiment(experiment, pressure, temperature)
            prediction = self.forward_model.predict(candidate, corrections)
            data += experiment_misfit(
                candidate,
                prediction,
                self.missing_phase_penalty,
                self.extra_phase_penalty,
                self.misfit_form,
                self.theoretical_covariance_scale,
            )
        means = np.asarray([parameter.prior_mean for parameter in self.parameters])
        deviation = numeric - means
        prior = 0.5 * deviation @ np.linalg.solve(self.prior_covariance, deviation)
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
    ) -> MAPResult:
        """Run bounded coordinate-pattern searches over all joint parameters."""
        if starts < 1 or iterations < 1:
            raise ValueError("starts and iterations must be positive")
        bounds = np.asarray([parameter.bounds() for parameter in self.parameters])
        centers = np.asarray(
            [parameter.prior_mean for parameter in self.parameters], dtype=float
        )
        rng = np.random.default_rng(seed)
        if central_point is not None:
            if len(central_point) != len(self.parameters):
                raise ValueError("central_point does not match parameter count")
            centers = np.asarray(central_point, dtype=float)
            if np.any(centers < bounds[:, 0]) or np.any(centers > bounds[:, 1]):
                raise ValueError("central_point must lie within parameter bounds")
        start_values = [centers]
        start_values.extend(
            bounds[:, 0]
            + rng.random(len(self.parameters)) * (bounds[:, 1] - bounds[:, 0])
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
                steps *= 0.5
                if not improved and np.max(steps) < 1.0e-8:
                    break
            records.append(
                {
                    "start": start.tolist(),
                    "map": values.tolist(),
                    "objective": objective,
                }
            )
            if objective < best_objective:
                best_values, best_objective = values, objective
        correction_names = [
            parameter.key for parameter in self.thermodynamic_parameters
        ]
        corrections = {
            key: float(value)
            for key, value in zip(
                correction_names, best_values[: len(correction_names)]
            )
        }
        if len(self.experiments) == 1:
            corrections["pressure_bar"] = float(best_values[-2])
            corrections["temperature_c"] = float(best_values[-1])
        else:
            offset = len(self.thermodynamic_parameters)
            for index, experiment in enumerate(self.experiments):
                corrections[f"pressure_bar:{experiment.id}"] = float(
                    best_values[offset + 2 * index]
                )
                corrections[f"temperature_c:{experiment.id}"] = float(
                    best_values[offset + 2 * index + 1]
                )
        return MAPResult(corrections, float(best_objective), records)
