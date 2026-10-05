"""Input-space uncertainty propagation for MAP thermodynamic inversion.

This module complements `posterior.uncertainty_source`: that feature splits
curvature by objective term (data versus prior), whereas this module perturbs
the actual reported observations and/or thermodynamic nuisance corrections,
then reruns MAP. It is the closer analogue of mc_fit's `invunc` perturbation
modes, while remaining explicit about the distinction.
"""

from dataclasses import dataclass, replace
import json
from pathlib import Path
from typing import Mapping

import numpy as np

from .data import EquilibriumBracket, Experiment
from .inversion import InversionProblem, MAPResult

PERTURBATION_SOURCES = ("analytical", "thermodynamic", "all")


@dataclass(frozen=True)
class InputPerturbationConfig:
    """Configuration for repeated input perturbation and MAP refitting."""

    n_iterations: int = 100
    source: str = "analytical"
    seed: int = 0
    analytical_scale: float = 1.0
    thermodynamic_sigma_j_mol: Mapping[str, float] | None = None
    thermodynamic_covariance_j_mol2: np.ndarray | None = None
    thermodynamic_covariance_keys: tuple[str, ...] | None = None
    optimize_starts: int = 4
    optimize_iterations: int = 60


@dataclass(frozen=True)
class InputPerturbationResult:
    """MAP corrections obtained from repeated perturbed input datasets."""

    parameter_names: tuple[str, ...]
    source: str
    corrections_samples: np.ndarray
    objectives: np.ndarray
    correction_mean: dict[str, float]
    correction_std: dict[str, float]
    success_count: int
    thermodynamic_covariance_keys: tuple[str, ...] = ()
    thermodynamic_covariance_diagnostics: dict[str, float | bool | int] | None = None


def load_packed_covariance_json(
    path: str | Path, scale: float = 1.0
) -> tuple[tuple[str, ...], np.ndarray]:
    """Load ``Entities``/``PackedUpperTriangle`` covariance JSON.

    Set ``scale=1000`` for Holland-Powell covariance exports stored in
    (kJ/mol)^2 when thermoinvert parameters are expressed in J/mol.
    """
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    keys = tuple(payload["Entities"])
    packed = payload["PackedUpperTriangle"]
    expected = len(keys) * (len(keys) + 1) // 2
    if len(packed) != expected:
        raise ValueError(f"Invalid packed covariance length: expected {expected}")
    covariance = np.zeros((len(keys), len(keys)), dtype=float)
    position = 0
    for row in range(len(keys)):
        for column in range(row, len(keys)):
            value = float(packed[position]) * scale
            covariance[row, column] = covariance[column, row] = value
            position += 1
    return keys, covariance


def covariance_diagnostics(
    covariance: np.ndarray,
) -> dict[str, float | bool | int]:
    """Return PSD, eigenvalue, effective-rank, and condition diagnostics."""
    matrix = np.asarray(covariance, dtype=float)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("covariance must be a square matrix")
    if not np.all(np.isfinite(matrix)):
        raise ValueError("covariance must contain only finite values")
    symmetric = 0.5 * (matrix + matrix.T)
    eigenvalues = np.linalg.eigvalsh(symmetric)
    scale = max(float(np.max(eigenvalues)) if len(eigenvalues) else 0.0, 1.0)
    tolerance = max(1.0e-20, 1.0e-12 * scale)
    positive = eigenvalues[eigenvalues > tolerance]
    return {
        "is_psd": bool(np.min(eigenvalues) >= -tolerance) if len(eigenvalues) else True,
        "min_eigenvalue": float(np.min(eigenvalues)) if len(eigenvalues) else 0.0,
        "max_eigenvalue": float(np.max(eigenvalues)) if len(eigenvalues) else 0.0,
        "effective_rank": int(len(positive)),
        "condition_number": float(np.max(positive) / np.min(positive))
        if len(positive)
        else float("inf"),
    }


def _sample_covariance(covariance: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Draw one correlated normal vector from a positive-semidefinite matrix."""
    covariance = np.asarray(covariance, dtype=float)
    if covariance.ndim != 2 or covariance.shape[0] != covariance.shape[1]:
        raise ValueError("covariance must be a square matrix")
    if not np.all(np.isfinite(covariance)):
        raise ValueError("covariance must contain only finite values")
    covariance = 0.5 * (covariance + covariance.T)
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    tolerance = max(1.0e-15, 1.0e-12 * max(float(np.max(eigenvalues)), 1.0))
    if np.any(eigenvalues < -tolerance):
        raise ValueError("covariance must be positive semidefinite")
    factor = eigenvectors @ np.diag(np.sqrt(np.clip(eigenvalues, 0.0, None)))
    return factor @ rng.standard_normal(covariance.shape[0])


def _perturb_continuous_observations(
    experiment: Experiment, rng: np.random.Generator, scale: float
) -> tuple[dict[str, float] | None, dict[str, dict[str, float]] | None]:
    """Perturb measured fractions/compositions using analytical uncertainty only."""
    fractions = dict(experiment.phase_fractions or {})
    compositions = {
        phase: dict(values)
        for phase, values in (experiment.phase_compositions or {}).items()
    }
    if experiment.observation_covariance is not None:
        keys = experiment.observation_covariance_keys
        covariance = np.asarray(experiment.observation_covariance, dtype=float)
        if covariance.shape != (len(keys), len(keys)):
            raise ValueError("observation_covariance shape must match its keys")
        offsets = _sample_covariance(covariance * scale**2, rng)
        for key, offset in zip(keys, offsets):
            kind, phase, *component = key.split(":")
            if kind == "fraction" and not component and phase in fractions:
                fractions[phase] += float(offset)
            elif (
                kind == "composition"
                and len(component) == 1
                and phase in compositions
                and component[0] in compositions[phase]
            ):
                compositions[phase][component[0]] += float(offset)
            else:
                raise ValueError(
                    "observation_covariance_keys must identify measured fractions or compositions"
                )
    else:
        if experiment.phase_fraction_sigma:
            fractions = {
                phase: value
                + rng.normal(
                    0.0, scale * experiment.phase_fraction_sigma.get(phase, 0.0)
                )
                for phase, value in fractions.items()
            }
        if experiment.phase_composition_sigma:
            compositions = {
                phase: {
                    component: value
                    + rng.normal(
                        0.0,
                        scale
                        * experiment.phase_composition_sigma.get(phase, {}).get(
                            component, 0.0
                        ),
                    )
                    for component, value in composition.items()
                }
                for phase, composition in compositions.items()
            }
    return fractions or None, compositions or None


def _perturb_experiment(
    experiment: Experiment, rng: np.random.Generator, scale: float
) -> Experiment:
    pressure = experiment.pressure_bar
    temperature = experiment.temperature_c
    if experiment.pressure_sigma_bar:
        pressure += rng.normal(0.0, scale * experiment.pressure_sigma_bar)
    if experiment.temperature_sigma_c:
        temperature += rng.normal(0.0, scale * experiment.temperature_sigma_c)

    phase_fractions, phase_compositions = _perturb_continuous_observations(
        experiment, rng, scale
    )

    return replace(
        experiment,
        pressure_bar=float(pressure),
        temperature_c=float(temperature),
        phase_fractions=phase_fractions,
        phase_compositions=phase_compositions,
    )


def _perturb_bracket(
    bracket: EquilibriumBracket, rng: np.random.Generator, scale: float
) -> EquilibriumBracket:
    return replace(
        bracket,
        lower=_perturb_experiment(bracket.lower, rng, scale),
        upper=_perturb_experiment(bracket.upper, rng, scale),
    )


class _ThermodynamicNoiseForwardModel:
    def __init__(self, forward_model, offsets: Mapping[str, float]):
        self.forward_model = forward_model
        self.offsets = dict(offsets)

    def predict(self, experiment, corrections=None):
        merged = dict(self.offsets)
        for key, value in (corrections or {}).items():
            merged[key] = merged.get(key, 0.0) + float(value)
        return self.forward_model.predict(experiment, merged)


def _fit(problem: InversionProblem, config: InputPerturbationConfig) -> MAPResult:
    if len(problem.parameters) == 1:
        return problem.optimize_one_parameter(
            starts=config.optimize_starts, iterations=config.optimize_iterations
        )
    return problem.optimize_multi_parameter(
        starts=config.optimize_starts, iterations=config.optimize_iterations
    )


def run_input_perturbation(
    problem: InversionProblem,
    config: InputPerturbationConfig = InputPerturbationConfig(),
) -> InputPerturbationResult:
    """Perturb input sources repeatedly and refit the MAP model.

    Analytical perturbations use each experiment's declared P-T uncertainty and
    either its joint analytical observation covariance or independent phase-
    fraction/composition sigmas. Thermodynamic perturbations use the supplied
    `thermodynamic_sigma_j_mol` mapping of correction keys to nuisance standard
    deviations; those offsets are injected into every forward-model call for
    one replicate. The target inversion parameters remain the fitted values.
    """
    if config.n_iterations < 1:
        raise ValueError("n_iterations must be positive")
    if config.source not in PERTURBATION_SOURCES:
        raise ValueError(f"source must be one of {PERTURBATION_SOURCES}")
    if config.analytical_scale < 0.0:
        raise ValueError("analytical_scale must be non-negative")
    thermo_sigma = dict(config.thermodynamic_sigma_j_mol or {})
    covariance = config.thermodynamic_covariance_j_mol2
    covariance_keys = config.thermodynamic_covariance_keys
    covariance_info = None
    if covariance is not None:
        covariance = np.asarray(covariance, dtype=float)
        if covariance.ndim != 2 or covariance.shape[0] != covariance.shape[1]:
            raise ValueError("thermodynamic covariance must be square")
        if covariance_keys is None or len(covariance_keys) != covariance.shape[0]:
            raise ValueError("thermodynamic_covariance_keys must match covariance size")
        if thermo_sigma:
            raise ValueError("Specify thermodynamic covariance or sigmas, not both")
        covariance_info = covariance_diagnostics(covariance)
    if (
        config.source in ("thermodynamic", "all")
        and covariance is None
        and not thermo_sigma
    ):
        raise ValueError(
            "thermodynamic sigma or covariance is required for thermodynamic perturbations"
        )
    if any(value < 0.0 for value in thermo_sigma.values()):
        raise ValueError("thermodynamic sigmas must be non-negative")

    rng = np.random.default_rng(config.seed)
    names = tuple(parameter.name for parameter in problem.parameters)
    samples = np.full((config.n_iterations, len(names)), np.nan)
    objectives = np.full(config.n_iterations, np.inf)
    success_count = 0

    for iteration in range(config.n_iterations):
        source_analytical = config.source in ("analytical", "all")
        experiments = [
            _perturb_experiment(experiment, rng, config.analytical_scale)
            if source_analytical
            else experiment
            for experiment in problem.experiments
        ]
        brackets = [
            _perturb_bracket(bracket, rng, config.analytical_scale)
            if source_analytical
            else bracket
            for bracket in problem.brackets
        ]
        if config.source in ("thermodynamic", "all"):
            if covariance is not None:
                offsets = dict(
                    zip(covariance_keys, _sample_covariance(covariance, rng))
                )
            else:
                offsets = {
                    key: float(rng.normal(0.0, sigma))
                    for key, sigma in thermo_sigma.items()
                }
        else:
            offsets = {}
        forward_model = (
            _ThermodynamicNoiseForwardModel(problem.forward_model, offsets)
            if offsets
            else problem.forward_model
        )
        replicate = InversionProblem(
            forward_model,
            experiments,
            problem.parameters,
            brackets=brackets,
            missing_phase_penalty=problem.missing_phase_penalty,
            extra_phase_penalty=problem.extra_phase_penalty,
            marginalize_condition_uncertainty=problem.marginalize_condition_uncertainty,
            condition_quadrature_order=problem.condition_quadrature_order,
            prior_covariance=problem.prior_covariance,
            objective_criterion=problem.objective_criterion,
            misfit_form=problem.misfit_form,
            theoretical_covariance_scale=problem.theoretical_covariance_scale,
            theoretical_covariance_scale_metadata_key=problem.theoretical_covariance_scale_metadata_key,
            theoretical_covariance_scales=problem.theoretical_covariance_scales,
        )
        try:
            result = _fit(replicate, config)
        except (RuntimeError, ValueError, FloatingPointError):
            continue
        if result.parameter_values is None:
            raise RuntimeError("MAP result is missing fitted parameter values")
        samples[iteration] = [
            result.parameter_values[parameter.name] for parameter in problem.parameters
        ]
        objectives[iteration] = result.objective
        success_count += 1

    finite = np.isfinite(objectives)
    if not np.any(finite):
        raise RuntimeError("No input-perturbation MAP replicates converged")
    mean = np.nanmean(samples, axis=0)
    std = (
        np.nanstd(samples, axis=0, ddof=1)
        if success_count > 1
        else np.zeros(len(names))
    )
    return InputPerturbationResult(
        parameter_names=names,
        source=config.source,
        corrections_samples=samples,
        objectives=objectives,
        correction_mean={name: float(mean[index]) for index, name in enumerate(names)},
        correction_std={name: float(std[index]) for index, name in enumerate(names)},
        success_count=success_count,
        thermodynamic_covariance_keys=tuple(covariance_keys or ()),
        thermodynamic_covariance_diagnostics=covariance_info,
    )
