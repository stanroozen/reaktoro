"""Leave-one-out cross-validation and posterior-predictive checks."""

from dataclasses import dataclass, replace
from pathlib import Path
from typing import Sequence

import numpy as np

from .data import DirichletConcentrationParameter, EquilibriumBracket, Experiment
from .forward_model import database_file_provenance
from .inversion import InversionProblem, MAPResult
from .likelihood import bracket_misfit, experiment_log_likelihood, experiment_misfit
from .posterior import MCMCResult
from .pt_calibration import (
    GroupedPTCalibrationProblem,
    mcmc_sample_grouped_pt_calibration,
)
from .pt_inversion import PTInversionProblem


def _fitted_dirichlet_concentrations(
    problem: InversionProblem, fit: MAPResult
) -> dict[str, float]:
    parameters = [
        parameter
        for parameter in problem.parameters
        if isinstance(parameter, DirichletConcentrationParameter)
    ]
    if not parameters:
        return {}
    if fit.parameter_values is None:
        raise ValueError("MAP result is missing fitted likelihood parameter values")
    return {
        parameter.phase: fit.parameter_values[parameter.name]
        for parameter in parameters
    }


def _held_out_log_score(
    problem: InversionProblem,
    experiment: Experiment,
    prediction,
    fit: MAPResult,
    phase_error_probability: float,
) -> float:
    scale = problem.theoretical_covariance_scale
    if problem.theoretical_covariance_scale_metadata_key:
        group = experiment.metadata.get(
            problem.theoretical_covariance_scale_metadata_key
        )
        if group not in problem.theoretical_covariance_scales:
            raise ValueError(f"No theoretical covariance scale for group {group!r}")
        scale = problem.theoretical_covariance_scales[group]
    return experiment_log_likelihood(
        experiment,
        prediction,
        problem.misfit_form,
        scale,
        _fitted_dirichlet_concentrations(problem, fit),
        phase_error_probability,
    )


def _has_normalized_predictive_likelihood(problem: InversionProblem) -> bool:
    return problem.misfit_form in {
        "weighted_lsq",
        "logit_lsq",
        "logistic_normal_nll",
        "dirichlet_nll",
    }


@dataclass(frozen=True)
class LOOCVFold:
    """One leave-one-out fold: refit corrections and held-out misfit."""

    held_out_experiment_id: str
    corrections_j_mol: dict[str, float]
    misfit: float
    held_out_experiment_ids: tuple[str, ...] = ()
    held_out_misfits: tuple[float, ...] = ()
    held_out_log_predictive_densities: tuple[float, ...] = ()


@dataclass(frozen=True)
class LOOCVResult:
    """Leave-one-out cross-validation summary across all experiments."""

    folds: tuple[LOOCVFold, ...]
    mean_misfit: float
    group_metadata_key: str | None = None
    mean_log_predictive_density: float | None = None


def _log_mean_exp(log_values: Sequence[float]) -> float:
    values = np.asarray(log_values, dtype=float)
    if values.size == 0:
        raise ValueError("at least one log likelihood is required")
    maximum = float(np.max(values))
    if not np.isfinite(maximum):
        return maximum
    return maximum + float(np.log(np.mean(np.exp(values - maximum))))


def validate_database_provenance(
    experiments: Sequence[Experiment],
    path_metadata_key: str = "database_path",
    hash_metadata_key: str = "database_sha256",
) -> tuple[dict[str, object], ...]:
    """Validate and return database hashes recorded in experiment metadata."""
    if not path_metadata_key or not hash_metadata_key:
        raise ValueError("database provenance metadata keys must be non-empty")
    records: list[dict[str, object]] = []
    for experiment in experiments:
        path = experiment.metadata.get(path_metadata_key)
        recorded_hash = experiment.metadata.get(hash_metadata_key)
        if not isinstance(path, str) or not path:
            raise ValueError(
                f"Experiment {experiment.id!r} lacks database path provenance"
            )
        provenance = database_file_provenance(path)
        if recorded_hash is not None and recorded_hash != provenance["database_sha256"]:
            raise ValueError(
                f"Experiment {experiment.id!r} has a mismatched database SHA-256"
            )
        records.append(provenance)
    return tuple(records)


def _training_brackets_for_groups(
    problem: InversionProblem,
    metadata_key: str,
    held_out_groups: set[object],
) -> list[EquilibriumBracket]:
    experiment_groups = {
        experiment.id: experiment.metadata[metadata_key]
        for experiment in problem.experiments
    }
    training_ids = {
        experiment_id
        for experiment_id, group in experiment_groups.items()
        if group not in held_out_groups
    }
    training_brackets: list[EquilibriumBracket] = []
    for bracket in problem.brackets:
        endpoint_groups: list[object] = []
        for endpoint in (bracket.lower, bracket.upper):
            group = experiment_groups.get(endpoint.id)
            if group is None:
                group = endpoint.metadata.get(metadata_key)
            if group is None:
                raise ValueError(
                    f"Bracket {bracket.id!r} endpoint {endpoint.id!r} lacks "
                    f"metadata {metadata_key!r}; cannot prevent validation leakage"
                )
            endpoint_groups.append(group)
        if all(group not in held_out_groups for group in endpoint_groups):
            if all(
                endpoint.id in training_ids
                for endpoint in (bracket.lower, bracket.upper)
            ):
                training_brackets.append(bracket)
    return training_brackets


def leave_one_out_cross_validation(
    problem: InversionProblem,
    optimize_starts: int = 4,
    optimize_iterations: int = 60,
    phase_error_probability: float = 0.01,
) -> LOOCVResult:
    """Refit corrections excluding each experiment and score it out-of-sample.

    Each fold drops one experiment, re-runs MAP optimization on the rest
    (using the same parameters, brackets, and penalties as `problem`), then
    scores the held-out experiment with the refit corrections. This checks
    whether the model generalizes rather than memorizing individual points.
    """
    if len(problem.experiments) < 2:
        raise ValueError("Cross-validation requires at least two experiments")

    folds: list[LOOCVFold] = []
    for index, held_out in enumerate(problem.experiments):
        training_experiments = [
            experiment
            for position, experiment in enumerate(problem.experiments)
            if position != index
        ]
        fold_problem = InversionProblem(
            problem.forward_model,
            training_experiments,
            problem.parameters,
            brackets=problem.brackets,
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
        if len(problem.parameters) == 1:
            fit = fold_problem.optimize_one_parameter(
                starts=optimize_starts, iterations=optimize_iterations
            )
        else:
            fit = fold_problem.optimize_multi_parameter(
                starts=optimize_starts, iterations=optimize_iterations
            )
        prediction = problem.forward_model.predict(held_out, fit.corrections_j_mol)
        concentrations = _fitted_dirichlet_concentrations(problem, fit)
        misfit = experiment_misfit(
            held_out,
            prediction,
            problem.missing_phase_penalty,
            problem.extra_phase_penalty,
            problem.misfit_form,
            problem.theoretical_covariance_scale,
            concentrations,
        )
        log_scores = (
            (
                _held_out_log_score(
                    problem, held_out, prediction, fit, phase_error_probability
                ),
            )
            if _has_normalized_predictive_likelihood(problem)
            else ()
        )
        folds.append(
            LOOCVFold(
                held_out.id,
                fit.corrections_j_mol,
                misfit,
                (held_out.id,),
                (misfit,),
                log_scores,
            )
        )

    finite_misfits = [fold.misfit for fold in folds if np.isfinite(fold.misfit)]
    mean_misfit = float(np.mean(finite_misfits)) if finite_misfits else float("inf")
    log_scores = [
        value for fold in folds for value in fold.held_out_log_predictive_densities
    ]
    mean_log_score = float(np.mean(log_scores)) if log_scores else None
    return LOOCVResult(tuple(folds), mean_misfit, None, mean_log_score)


def leave_one_group_out_cross_validation(
    problem: InversionProblem,
    metadata_key: str,
    optimize_starts: int = 4,
    optimize_iterations: int = 60,
    phase_error_probability: float = 0.01,
) -> LOOCVResult:
    """Hold out all experiments sharing a metadata group and refit.

    Use this for study, laboratory, or phase-assemblage groups. Each returned
    fold contains the group's mean held-out misfit, so it can be plotted with
    the existing LOOCV tools. At least two groups are required.
    """
    if not metadata_key:
        raise ValueError("metadata_key must be non-empty")
    if len(problem.experiments) < 2:
        raise ValueError("Grouped cross-validation requires at least two experiments")
    groups: dict[object, list[Experiment]] = {}
    for experiment in problem.experiments:
        if metadata_key not in experiment.metadata:
            raise ValueError(
                f"Experiment {experiment.id!r} lacks metadata {metadata_key!r}"
            )
        value = experiment.metadata[metadata_key]
        try:
            groups.setdefault(value, []).append(experiment)
        except TypeError as error:
            raise ValueError(
                f"Metadata group {metadata_key!r} must be hashable"
            ) from error
    if len(groups) < 2:
        raise ValueError("Grouped cross-validation requires at least two groups")

    folds: list[LOOCVFold] = []
    for group, held_out_experiments in groups.items():
        held_out_ids = {experiment.id for experiment in held_out_experiments}
        training_experiments = [
            experiment
            for experiment in problem.experiments
            if experiment.id not in held_out_ids
        ]
        training_brackets = _training_brackets_for_groups(
            problem, metadata_key, {group}
        )
        fold_problem = InversionProblem(
            problem.forward_model,
            training_experiments,
            problem.parameters,
            brackets=training_brackets,
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
        if len(problem.parameters) == 1:
            fit = fold_problem.optimize_one_parameter(
                starts=optimize_starts, iterations=optimize_iterations
            )
        else:
            fit = fold_problem.optimize_multi_parameter(
                starts=optimize_starts, iterations=optimize_iterations
            )
        predictions = [
            problem.forward_model.predict(experiment, fit.corrections_j_mol)
            for experiment in held_out_experiments
        ]
        concentrations = _fitted_dirichlet_concentrations(problem, fit)
        misfits = [
            experiment_misfit(
                experiment,
                prediction,
                problem.missing_phase_penalty,
                problem.extra_phase_penalty,
                problem.misfit_form,
                problem.theoretical_covariance_scale,
                concentrations,
            )
            for experiment, prediction in zip(held_out_experiments, predictions)
        ]
        log_scores = (
            tuple(
                _held_out_log_score(
                    problem, experiment, prediction, fit, phase_error_probability
                )
                for experiment, prediction in zip(held_out_experiments, predictions)
            )
            if _has_normalized_predictive_likelihood(problem)
            else ()
        )
        finite_misfits = [misfit for misfit in misfits if np.isfinite(misfit)]
        group_misfit = (
            float(np.mean(finite_misfits)) if finite_misfits else float("inf")
        )
        folds.append(
            LOOCVFold(
                f"{metadata_key}={group}",
                fit.corrections_j_mol,
                group_misfit,
                tuple(experiment.id for experiment in held_out_experiments),
                tuple(float(misfit) for misfit in misfits),
                log_scores,
            )
        )

    finite_misfits = [fold.misfit for fold in folds if np.isfinite(fold.misfit)]
    log_scores = [
        value for fold in folds for value in fold.held_out_log_predictive_densities
    ]
    return LOOCVResult(
        tuple(folds),
        float(np.mean(finite_misfits)) if finite_misfits else float("inf"),
        metadata_key,
        float(np.mean(log_scores)) if log_scores else None,
    )


def repeated_group_holdout_cross_validation(
    problem: InversionProblem,
    metadata_key: str,
    repeats: int = 5,
    test_fraction: float = 0.25,
    seed: int = 0,
    optimize_starts: int = 4,
    optimize_iterations: int = 60,
    phase_error_probability: float = 0.01,
) -> LOOCVResult:
    """Evaluate repeated random held-out group partitions.

    Every split holds out complete metadata groups, ensuring that study or
    laboratory effects cannot leak from training into validation observations.
    """
    if not metadata_key:
        raise ValueError("metadata_key must be non-empty")
    if repeats < 1 or not 0.0 < test_fraction < 1.0:
        raise ValueError(
            "repeats must be positive and test_fraction must lie in (0, 1)"
        )
    groups: dict[object, list[Experiment]] = {}
    for experiment in problem.experiments:
        if metadata_key not in experiment.metadata:
            raise ValueError(
                f"Experiment {experiment.id!r} lacks metadata {metadata_key!r}"
            )
        try:
            groups.setdefault(experiment.metadata[metadata_key], []).append(experiment)
        except TypeError as error:
            raise ValueError(
                f"Metadata group {metadata_key!r} must be hashable"
            ) from error
    group_values = list(groups)
    if len(group_values) < 2:
        raise ValueError("Repeated grouped validation requires at least two groups")
    holdout_count = min(
        len(group_values) - 1,
        max(1, int(round(len(group_values) * test_fraction))),
    )
    rng = np.random.default_rng(seed)
    folds: list[LOOCVFold] = []
    for repeat in range(repeats):
        held_out_groups = set(
            rng.choice(group_values, size=holdout_count, replace=False)
        )
        held_out_experiments = [
            experiment for group in held_out_groups for experiment in groups[group]
        ]
        training_experiments = [
            experiment
            for group, members in groups.items()
            if group not in held_out_groups
            for experiment in members
        ]
        training_brackets = _training_brackets_for_groups(
            problem, metadata_key, held_out_groups
        )
        fold_problem = InversionProblem(
            problem.forward_model,
            training_experiments,
            problem.parameters,
            brackets=training_brackets,
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
        fit = (
            fold_problem.optimize_one_parameter(
                starts=optimize_starts, iterations=optimize_iterations
            )
            if len(problem.parameters) == 1
            else fold_problem.optimize_multi_parameter(
                starts=optimize_starts, iterations=optimize_iterations
            )
        )
        predictions = [
            problem.forward_model.predict(experiment, fit.corrections_j_mol)
            for experiment in held_out_experiments
        ]
        concentrations = _fitted_dirichlet_concentrations(problem, fit)
        misfits = [
            experiment_misfit(
                experiment,
                prediction,
                problem.missing_phase_penalty,
                problem.extra_phase_penalty,
                problem.misfit_form,
                problem.theoretical_covariance_scale,
                concentrations,
            )
            for experiment, prediction in zip(held_out_experiments, predictions)
        ]
        log_scores = (
            tuple(
                _held_out_log_score(
                    problem, experiment, prediction, fit, phase_error_probability
                )
                for experiment, prediction in zip(held_out_experiments, predictions)
            )
            if _has_normalized_predictive_likelihood(problem)
            else ()
        )
        finite = [misfit for misfit in misfits if np.isfinite(misfit)]
        folds.append(
            LOOCVFold(
                f"repeat={repeat}",
                fit.corrections_j_mol,
                float(np.mean(finite)) if finite else float("inf"),
                tuple(experiment.id for experiment in held_out_experiments),
                tuple(float(misfit) for misfit in misfits),
                log_scores,
            )
        )
    finite = [fold.misfit for fold in folds if np.isfinite(fold.misfit)]
    log_scores = [
        value for fold in folds for value in fold.held_out_log_predictive_densities
    ]
    return LOOCVResult(
        tuple(folds),
        float(np.mean(finite)) if finite else float("inf"),
        metadata_key,
        float(np.mean(log_scores)) if log_scores else None,
    )


def repeated_grouped_pt_calibration_validation(
    problem: GroupedPTCalibrationProblem,
    repeats: int = 5,
    test_fraction: float = 0.25,
    n_samples: int = 2000,
    burn_in: int = 500,
    predictive_draws: int = 100,
    seed: int = 0,
    n_chains: int = 4,
    phase_error_probability: float = 0.01,
) -> LOOCVResult:
    """Validate hierarchical P-T calibration on unseen metadata groups.

    Each repeat refits the hyperposterior using training groups, then draws
    shared pressure/temperature offsets for each held-out group from posterior
    scale draws. Returned held-out misfits are averages over those predictive
    draws; they are not plug-in scores at zero group bias.
    """
    if repeats < 1 or not 0.0 < test_fraction < 1.0:
        raise ValueError(
            "repeats must be positive and test_fraction must lie in (0, 1)"
        )
    if n_samples < 2 or not 0 <= burn_in < n_samples:
        raise ValueError("n_samples must exceed burn_in and be at least two")
    if predictive_draws < 1:
        raise ValueError("predictive_draws must be positive")
    if n_chains < 2:
        raise ValueError("n_chains must be at least two")
    if (
        not np.isfinite(phase_error_probability)
        or not 0.0 < phase_error_probability < 0.5
    ):
        raise ValueError("phase_error_probability must lie strictly between 0 and 0.5")
    grouped: dict[str, list[Experiment]] = {}
    for experiment in problem.experiments:
        group = experiment.metadata[problem.metadata_key]
        grouped.setdefault(group, []).append(experiment)
    group_values = list(problem.groups)
    if len(group_values) < 2:
        raise ValueError("grouped P-T validation requires at least two groups")
    holdout_count = min(
        len(group_values) - 1,
        max(1, int(round(len(group_values) * test_fraction))),
    )
    rng = np.random.default_rng(seed)
    folds: list[LOOCVFold] = []
    for repeat in range(repeats):
        held_out_groups = set(
            rng.choice(group_values, size=holdout_count, replace=False).tolist()
        )
        training = [
            experiment
            for group, experiments in grouped.items()
            if group not in held_out_groups
            for experiment in experiments
        ]
        held_out = [
            experiment for group in held_out_groups for experiment in grouped[group]
        ]
        training_problem = GroupedPTCalibrationProblem(
            problem.forward_model,
            training,
            problem.metadata_key,
            problem.pressure_scale_hyperprior,
            problem.temperature_scale_hyperprior,
            corrections=problem.corrections,
            missing_phase_penalty=problem.missing_phase_penalty,
            extra_phase_penalty=problem.extra_phase_penalty,
            misfit_form=problem.misfit_form,
            theoretical_covariance_scale=problem.theoretical_covariance_scale,
            thermodynamic_parameters=problem.thermodynamic_parameters,
            thermodynamic_prior_covariance=problem.thermodynamic_prior_covariance,
            dirichlet_concentration_parameters=problem.dirichlet_concentration_parameters,
            observation_bias_keys=problem.observation_bias_keys,
            observation_bias_prior_sigma=problem.observation_bias_prior_sigma,
            observation_bias_bounds=problem.observation_bias_bounds,
            model_error_scale_prior_sigma=problem.model_error_scale_prior_sigma,
            model_error_scale_bounds=problem.model_error_scale_bounds,
        )
        posterior = mcmc_sample_grouped_pt_calibration(
            training_problem,
            n_chains=n_chains,
            n_samples=n_samples,
            burn_in=burn_in,
            seed=seed + repeat,
        )
        posterior_samples = np.concatenate(
            [chain.samples for chain in posterior.chains], axis=0
        )
        draw_count = min(predictive_draws, len(posterior_samples))
        sample_indices = rng.choice(
            len(posterior_samples), size=draw_count, replace=False
        )
        predictive_misfits: dict[str, list[float]] = {item.id: [] for item in held_out}
        predictive_log_likelihoods: dict[str, list[float]] = {
            item.id: [] for item in held_out
        }
        for sample_index in sample_indices:
            posterior_values = posterior_samples[sample_index]
            scale_pressure, scale_temperature = posterior_values[-2:]
            correction_map = training_problem.corrections_for_values(posterior_values)
            biases = problem.sample_unseen_group_offsets(
                float(scale_pressure),
                float(scale_temperature),
                rng,
                tuple(sorted(held_out_groups)),
            )
            observation_biases, theoretical_scales = (
                problem.sample_unseen_group_nuisances(
                    tuple(sorted(held_out_groups)), rng
                )
            )
            for experiment in held_out:
                group = experiment.metadata[problem.metadata_key]
                pressure_bias, temperature_bias = biases[group]
                fixed_and_drawn_biases = dict(experiment.observation_bias or {})
                for key in problem.observation_bias_keys:
                    fixed_and_drawn_biases[key] = (
                        fixed_and_drawn_biases.get(key, 0.0)
                        + observation_biases[(key, group)]
                    )
                candidate = replace(
                    experiment,
                    pressure_bar=experiment.pressure_bar + pressure_bias,
                    temperature_c=experiment.temperature_c + temperature_bias,
                    observation_bias=fixed_and_drawn_biases or None,
                )
                covariance_scale = theoretical_scales.get(
                    group, problem.theoretical_covariance_scale
                )
                prediction = problem.forward_model.predict(candidate, correction_map)
                predictive_misfits[experiment.id].append(
                    experiment_misfit(
                        candidate,
                        prediction,
                        problem.missing_phase_penalty,
                        problem.extra_phase_penalty,
                        problem.misfit_form,
                        covariance_scale,
                        training_problem.dirichlet_concentrations_for_values(
                            posterior_values
                        ),
                    )
                )
                if _has_normalized_predictive_likelihood(problem):
                    predictive_log_likelihoods[experiment.id].append(
                        experiment_log_likelihood(
                            candidate,
                            prediction,
                            problem.misfit_form,
                            covariance_scale,
                            training_problem.dirichlet_concentrations_for_values(
                                posterior_values
                            ),
                            phase_error_probability,
                        )
                    )
        held_out_misfits = tuple(
            float(np.mean(predictive_misfits[item.id])) for item in held_out
        )
        finite = [value for value in held_out_misfits if np.isfinite(value)]
        held_out_log_scores = (
            tuple(
                _log_mean_exp(predictive_log_likelihoods[item.id]) for item in held_out
            )
            if _has_normalized_predictive_likelihood(problem)
            else ()
        )
        posterior_parameter_means = {
            name: float(np.mean(posterior_samples[:, index]))
            for index, name in enumerate(posterior.parameter_names)
        }
        folds.append(
            LOOCVFold(
                f"repeat={repeat}",
                posterior_parameter_means,
                float(np.mean(finite)) if finite else float("inf"),
                tuple(item.id for item in held_out),
                held_out_misfits,
                held_out_log_scores,
            )
        )
    finite_scores = [fold.misfit for fold in folds if np.isfinite(fold.misfit)]
    predictive_scores = [
        score
        for fold in folds
        for score in fold.held_out_log_predictive_densities
    ]
    return LOOCVResult(
        tuple(folds),
        float(np.mean(finite_scores)) if finite_scores else float("inf"),
        problem.metadata_key,
        float(np.mean(predictive_scores)) if predictive_scores else None,
    )


def leave_one_phase_assemblage_out_cross_validation(
    problem: InversionProblem,
    optimize_starts: int = 4,
    optimize_iterations: int = 60,
    phase_error_probability: float = 0.01,
) -> LOOCVResult:
    """Hold out each observed phase-assemblage group in turn."""
    if not problem.experiments:
        raise ValueError("Phase-assemblage validation requires experiments")
    metadata_key = "__thermoinvert_phase_assemblage__"
    grouped_experiments = [
        replace(
            experiment,
            metadata={
                **experiment.metadata,
                metadata_key: "|".join(
                    (
                        *sorted(experiment.observed_phases),
                        *(f"!{phase}" for phase in sorted(experiment.absent_phases)),
                    )
                ),
            },
        )
        for experiment in problem.experiments
    ]
    grouped_problem = InversionProblem(
        problem.forward_model,
        grouped_experiments,
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
    return leave_one_group_out_cross_validation(
        grouped_problem,
        metadata_key,
        optimize_starts=optimize_starts,
        optimize_iterations=optimize_iterations,
        phase_error_probability=phase_error_probability,
    )


def leave_one_database_family_out_cross_validation(
    problem: InversionProblem,
    metadata_key: str = "database_family",
    optimize_starts: int = 4,
    optimize_iterations: int = 60,
    phase_error_probability: float = 0.01,
) -> LOOCVResult:
    """Hold out each database-family group in turn using experiment metadata."""
    if not metadata_key:
        raise ValueError("metadata_key must be non-empty")
    grouped_experiments: list[Experiment] = []
    for experiment in problem.experiments:
        family = experiment.metadata.get(metadata_key)
        if family is None and "database_path" in experiment.metadata:
            provenance = database_file_provenance(experiment.metadata["database_path"])
            family = provenance["database_sha256"]
        if not isinstance(family, str) or not family:
            raise ValueError(
                f"Experiment {experiment.id!r} needs a non-empty database-family metadata value"
            )
        grouped_experiments.append(
            replace(
                experiment,
                metadata={**experiment.metadata, metadata_key: family},
            )
        )
    grouped_problem = InversionProblem(
        problem.forward_model,
        grouped_experiments,
        problem.parameters,
        brackets=problem.brackets,
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
    return leave_one_group_out_cross_validation(
        grouped_problem,
        metadata_key,
        optimize_starts=optimize_starts,
        optimize_iterations=optimize_iterations,
        phase_error_probability=phase_error_probability,
    )


def leave_one_out_of_domain_pt_cross_validation(
    problem: InversionProblem,
    pressure_bounds: tuple[float, float],
    temperature_bounds: tuple[float, float],
    optimize_starts: int = 4,
    optimize_iterations: int = 60,
    phase_error_probability: float = 0.01,
) -> LOOCVResult:
    """Train inside a P-T rectangle and score experiments outside it."""
    pressure_low, pressure_high = pressure_bounds
    temperature_low, temperature_high = temperature_bounds
    if not (
        np.isfinite(pressure_low)
        and np.isfinite(pressure_high)
        and pressure_low < pressure_high
        and np.isfinite(temperature_low)
        and np.isfinite(temperature_high)
        and temperature_low < temperature_high
    ):
        raise ValueError(
            "pressure and temperature bounds must be finite and increasing"
        )
    held_out = [
        experiment
        for experiment in problem.experiments
        if not (
            pressure_low <= experiment.pressure_bar <= pressure_high
            and temperature_low <= experiment.temperature_c <= temperature_high
        )
    ]
    training = [
        experiment for experiment in problem.experiments if experiment not in held_out
    ]
    if not held_out or not training:
        raise ValueError(
            "out-of-domain validation requires non-empty training and held-out sets"
        )
    fold_problem = InversionProblem(
        problem.forward_model,
        training,
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
        fold_problem.optimize_one_parameter(
            starts=optimize_starts, iterations=optimize_iterations
        )
        if len(problem.parameters) == 1
        else fold_problem.optimize_multi_parameter(
            starts=optimize_starts, iterations=optimize_iterations
        )
    )
    predictions = [
        problem.forward_model.predict(experiment, fit.corrections_j_mol)
        for experiment in held_out
    ]
    concentrations = _fitted_dirichlet_concentrations(problem, fit)
    misfits = [
        experiment_misfit(
            experiment,
            prediction,
            problem.missing_phase_penalty,
            problem.extra_phase_penalty,
            problem.misfit_form,
            problem.theoretical_covariance_scale,
            concentrations,
        )
        for experiment, prediction in zip(held_out, predictions)
    ]
    log_scores = (
        tuple(
            _held_out_log_score(
                problem, experiment, prediction, fit, phase_error_probability
            )
            for experiment, prediction in zip(held_out, predictions)
        )
        if _has_normalized_predictive_likelihood(problem)
        else ()
    )
    finite = [misfit for misfit in misfits if np.isfinite(misfit)]
    return LOOCVResult(
        (
            LOOCVFold(
                "out_of_domain_pt",
                fit.corrections_j_mol,
                float(np.mean(finite)) if finite else float("inf"),
                tuple(experiment.id for experiment in held_out),
                tuple(float(misfit) for misfit in misfits),
                log_scores,
            ),
        ),
        float(np.mean(finite)) if finite else float("inf"),
        "out_of_domain_pt",
        float(np.mean(log_scores)) if log_scores else None,
    )


def _observed_observables(experiment: Experiment) -> tuple[tuple[str, ...], np.ndarray]:
    names: list[str] = []
    values: list[float] = []
    for phase in sorted(experiment.observed_phases):
        names.append(f"phase_present:{phase}")
        values.append(1.0)
    if experiment.phase_fractions:
        for phase in sorted(experiment.phase_fractions):
            names.append(f"phase_fraction:{phase}")
            values.append(experiment.phase_fractions[phase])
    if experiment.phase_compositions:
        for phase in sorted(experiment.phase_compositions):
            for component in sorted(experiment.phase_compositions[phase]):
                names.append(f"composition:{phase}:{component}")
                values.append(experiment.phase_compositions[phase][component])
    if not values:
        raise ValueError(f"Experiment {experiment.id!r} has no supported observables")
    return tuple(names), np.asarray(values, dtype=float)


def _predicted_observables(names: Sequence[str], prediction) -> np.ndarray:
    values: list[float] = []
    for name in names:
        kind, _, rest = name.partition(":")
        if kind == "phase_present":
            values.append(1.0 if rest in prediction.stable_phases else 0.0)
        elif kind == "phase_fraction":
            values.append(prediction.phase_fractions.get(rest, 0.0))
        else:
            phase, _, component = rest.partition(":")
            values.append(
                prediction.phase_compositions.get(phase, {}).get(component, 0.0)
            )
    return np.asarray(values, dtype=float)


def _observation_noise(
    experiment: Experiment,
    names: Sequence[str],
    rng: np.random.Generator,
) -> np.ndarray:
    """Draw one continuous-observation noise vector in observable-name order."""
    continuous_indices = [
        index
        for index, name in enumerate(names)
        if not name.startswith("phase_present:")
    ]
    if not continuous_indices:
        return np.zeros(len(names))
    keys = []
    sigmas = []
    for index in continuous_indices:
        kind, _, rest = names[index].partition(":")
        key = f"fraction:{rest}" if kind == "phase_fraction" else f"composition:{rest}"
        keys.append(key)
        if kind == "phase_fraction":
            sigma = (
                experiment.phase_fraction_sigma.get(rest)
                if experiment.phase_fraction_sigma
                else None
            )
        else:
            phase, _, component = rest.partition(":")
            sigma = (
                experiment.phase_composition_sigma.get(phase, {}).get(component)
                if experiment.phase_composition_sigma
                else None
            )
        sigmas.append(sigma)

    covariance = np.diag(
        [
            float(sigma) ** 2 if sigma is not None and sigma > 0.0 else 0.0
            for sigma in sigmas
        ]
    )
    if experiment.observation_covariance is not None:
        declared_keys = experiment.observation_covariance_keys
        declared = np.asarray(experiment.observation_covariance, dtype=float)
        if set(declared_keys) != set(keys) or len(declared_keys) != len(keys):
            raise ValueError(
                "observation_covariance_keys must identify every continuous observable"
            )
        if declared.shape != (len(declared_keys), len(declared_keys)):
            raise ValueError("observation_covariance shape must match its keys")
        indices = [declared_keys.index(key) for key in keys]
        covariance = declared[np.ix_(indices, indices)]
    if experiment.theoretical_covariance is not None:
        declared_keys = experiment.theoretical_covariance_keys
        theoretical = np.asarray(experiment.theoretical_covariance, dtype=float)
        if set(declared_keys) != set(keys) or len(declared_keys) != len(keys):
            raise ValueError(
                "theoretical_covariance_keys must identify every continuous observable"
            )
        if theoretical.shape != (len(declared_keys), len(declared_keys)):
            raise ValueError("theoretical_covariance shape must match its keys")
        indices = [declared_keys.index(key) for key in keys]
        covariance = covariance + theoretical[np.ix_(indices, indices)]
    covariance = 0.5 * (covariance + covariance.T)
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    tolerance = max(1.0e-15, 1.0e-12 * max(float(np.max(eigenvalues)), 1.0))
    if np.any(eigenvalues < -tolerance):
        raise ValueError(
            "predictive observation covariance must be positive semidefinite"
        )
    noise = np.zeros(len(names))
    factor = eigenvectors @ np.diag(np.sqrt(np.clip(eigenvalues, 0.0, None)))
    noise[continuous_indices] = factor @ rng.standard_normal(len(keys))
    return noise


def _posterior_indices(
    total_samples: int, n_predictive_samples: int, selection: str, seed: int
) -> np.ndarray:
    if selection not in {"even", "random"}:
        raise ValueError("sample_selection must be 'even' or 'random'")
    n = min(n_predictive_samples, total_samples)
    if n < 1:
        raise ValueError("mcmc_result contains no samples")
    if selection == "even":
        return np.linspace(0, total_samples - 1, n).astype(int)
    return np.random.default_rng(seed).choice(total_samples, n, replace=False)


@dataclass(frozen=True)
class PosteriorPredictiveResult:
    """Posterior predictive distribution for one experiment's observables."""

    observable_names: tuple[str, ...]
    observed_values: np.ndarray
    predictive_mean: np.ndarray
    predictive_std: np.ndarray
    credible_interval_95: dict[str, tuple[float, float]]
    within_95: dict[str, bool]
    include_observation_uncertainty: bool = False
    seed: int | None = None
    sample_selection: str = "even"


@dataclass(frozen=True)
class BracketPosteriorPredictiveResult:
    """Posterior predictive satisfaction of pure-phase reaction brackets."""

    bracket_ids: tuple[str, ...]
    samples_evaluated: int
    satisfied_fraction: float
    mean_misfit: float
    credible_interval_95: tuple[float, float]
    misfits: tuple[float, ...] = ()
    seed: int | None = None
    sample_selection: str = "even"


def posterior_predictive_check_pt_brackets(
    problem: PTInversionProblem,
    mcmc_result: MCMCResult,
    n_predictive_samples: int = 200,
    sample_selection: str = "even",
    seed: int = 0,
) -> BracketPosteriorPredictiveResult:
    """Check whether posterior P-T samples reproduce bracket phase inequalities.

    Unlike composition PPC, this reports a satisfaction probability and
    bracket-misfit distribution because phase-presence inequalities are binary
    and discontinuous. A sample is satisfied only when every bracket has zero
    `bracket_misfit` at its translated endpoints.
    """
    if not problem.brackets:
        raise ValueError("posterior_predictive_check_pt_brackets requires brackets")
    if mcmc_result.samples.shape[1] != 2:
        raise ValueError("mcmc_result must contain pressure and temperature columns")
    total_samples = len(mcmc_result.samples)
    indices = _posterior_indices(
        total_samples, n_predictive_samples, sample_selection, seed
    )
    n = len(indices)

    misfits: list[float] = []
    satisfied = 0
    for sample_index in indices:
        values = mcmc_result.samples[sample_index]
        sample_misfit = 0.0
        for bracket in problem.brackets:
            candidate = problem._translate_bracket(bracket, values)
            sample_misfit += bracket_misfit(
                candidate,
                problem.forward_model.predict(candidate.lower, problem.corrections),
                problem.forward_model.predict(candidate.upper, problem.corrections),
                problem.missing_phase_penalty,
            )
        misfits.append(float(sample_misfit))
        if sample_misfit == 0.0:
            satisfied += 1

    return BracketPosteriorPredictiveResult(
        tuple(bracket.id for bracket in problem.brackets),
        n,
        satisfied / n,
        float(np.mean(misfits)),
        (float(np.percentile(misfits, 2.5)), float(np.percentile(misfits, 97.5))),
        tuple(misfits),
        seed,
        sample_selection,
    )


def posterior_predictive_check(
    problem: InversionProblem,
    experiment: Experiment,
    mcmc_result: MCMCResult,
    n_predictive_samples: int = 200,
    include_observation_uncertainty: bool = False,
    seed: int = 0,
    sample_selection: str = "even",
) -> PosteriorPredictiveResult:
    """Compare an experiment's observations to predictions over posterior samples."""
    names, observed = _observed_observables(experiment)
    total_samples = len(mcmc_result.samples)
    indices = _posterior_indices(
        total_samples, n_predictive_samples, sample_selection, seed
    )
    n = len(indices)

    predicted = np.empty((n, len(names)))
    rng = np.random.default_rng(seed)
    for row, sample_index in enumerate(indices):
        corrections = problem.corrections_for_values(mcmc_result.samples[sample_index])
        prediction = problem.forward_model.predict(experiment, corrections)
        predicted[row] = _predicted_observables(names, prediction)
        if include_observation_uncertainty:
            predicted[row] += _observation_noise(experiment, names, rng)

    predictive_mean = predicted.mean(axis=0)
    predictive_std = predicted.std(axis=0, ddof=1) if n > 1 else np.zeros(len(names))
    low = np.percentile(predicted, 2.5, axis=0)
    high = np.percentile(predicted, 97.5, axis=0)

    credible_interval_95 = {
        name: (float(low[i]), float(high[i])) for i, name in enumerate(names)
    }
    within_95 = {
        name: bool(low[i] <= observed[i] <= high[i]) for i, name in enumerate(names)
    }
    return PosteriorPredictiveResult(
        names,
        observed,
        predictive_mean,
        predictive_std,
        credible_interval_95,
        within_95,
        include_observation_uncertainty,
        seed,
        sample_selection,
    )


def posterior_predictive_check_pt(
    problem: PTInversionProblem,
    mcmc_result: MCMCResult,
    n_predictive_samples: int = 200,
    include_observation_uncertainty: bool = False,
    seed: int = 0,
    sample_selection: str = "even",
) -> PosteriorPredictiveResult:
    """Posterior predictive check for composition/fraction P-T inversion.

    Replays MCMC samples as candidate ``[pressure_bar, temperature_c]``
    conditions. Bracket-only P-T problems are excluded because their binary
    phase inequalities require a reaction-boundary-specific predictive check,
    not a comparison to one experiment's continuous observables.
    """
    if problem.brackets:
        raise ValueError(
            "posterior_predictive_check_pt requires composition/fraction observations, not brackets"
        )
    names, observed = _observed_observables(problem.experiment)
    total_samples = len(mcmc_result.samples)
    indices = _posterior_indices(
        total_samples, n_predictive_samples, sample_selection, seed
    )
    n = len(indices)
    if mcmc_result.samples.shape[1] != 2:
        raise ValueError("mcmc_result must contain pressure and temperature columns")
    predicted = np.empty((n, len(names)))
    rng = np.random.default_rng(seed)
    for row, sample_index in enumerate(indices):
        candidate = problem._candidate_experiment(mcmc_result.samples[sample_index])
        prediction = problem.forward_model.predict(candidate, problem.corrections)
        predicted[row] = _predicted_observables(names, prediction)
        if include_observation_uncertainty:
            predicted[row] += _observation_noise(problem.experiment, names, rng)

    predictive_mean = predicted.mean(axis=0)
    predictive_std = predicted.std(axis=0, ddof=1) if n > 1 else np.zeros(len(names))
    low = np.percentile(predicted, 2.5, axis=0)
    high = np.percentile(predicted, 97.5, axis=0)
    credible_interval_95 = {
        name: (float(low[i]), float(high[i])) for i, name in enumerate(names)
    }
    within_95 = {
        name: bool(low[i] <= observed[i] <= high[i]) for i, name in enumerate(names)
    }
    return PosteriorPredictiveResult(
        names,
        observed,
        predictive_mean,
        predictive_std,
        credible_interval_95,
        within_95,
        include_observation_uncertainty,
        seed,
        sample_selection,
    )
