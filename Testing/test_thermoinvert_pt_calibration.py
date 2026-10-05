"""Tests for hierarchical study/lab P-T calibration offsets."""

import sys
from dataclasses import replace
from math import erf, sqrt
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parents[1]))

from thermoinvert.data import (
    DirichletConcentrationParameter,
    Experiment,
    ThermodynamicParameter,
)
from thermoinvert.forward_model import EquilibriumResult
from thermoinvert.posterior import mcmc_sample
from thermoinvert.pt_calibration import (
    GroupedPTCalibrationProblem,
    mcmc_sample_grouped_pt_calibration,
)
import thermoinvert.validation as validation_module
from thermoinvert.archive import (
    read_grouped_pt_calibration_configuration,
    write_grouped_pt_calibration_configuration,
)
from thermoinvert.validation import repeated_grouped_pt_calibration_validation


class ConstantCalibrationForwardModel:
    def predict(self, experiment, corrections):
        return EquilibriumResult(
            True,
            ("phase",),
            {"phase": 1.0},
            {},
            {},
        )


def _experiment(identifier: str, group: str) -> Experiment:
    return Experiment(
        identifier,
        2500.0,
        600.0,
        {"A": 1.0},
        ("phase",),
        phase_fractions={"phase": 1.0},
        phase_fraction_sigma={"phase": 0.1},
        metadata={"study": group},
    )


def test_grouped_pt_calibration_has_shared_scale_hyperpriors_and_samples():
    problem = GroupedPTCalibrationProblem(
        ConstantCalibrationForwardModel(),
        [_experiment("one", "lab-a"), _experiment("two", "lab-b")],
        "study",
        pressure_scale_hyperprior=100.0,
        temperature_scale_hyperprior=10.0,
    )
    data, prior = problem.objective_components(problem.initial_point)

    assert len(problem.parameters) == 6
    assert problem.parameters[-2].name == "pressure_bias_scale"
    assert problem.parameters[-1].name == "temperature_bias_scale"
    assert np.isclose(data, 0.0)
    expected_prior = 1.0 + 4.0 * np.log(erf(5.0 / sqrt(2.0)))
    assert np.isclose(prior, expected_prior)

    posterior = mcmc_sample(
        problem,
        problem.initial_point,
        n_samples=30,
        burn_in=5,
        seed=14,
    )

    assert posterior.samples.shape == (25, 6)
    assert posterior.parameter_names[-2:] == (
        "pressure_bias_scale",
        "temperature_bias_scale",
    )
    assert posterior.parameter_units[-2:] == ("bar", "degC")

    offsets = problem.sample_unseen_group_offsets(
        500.0, 50.0, np.random.default_rng(9), ("new-lab", "other-lab")
    )
    assert set(offsets) == {"new-lab", "other-lab"}
    assert all(abs(value[0]) <= 500.0 for value in offsets.values())
    assert all(abs(value[1]) <= 50.0 for value in offsets.values())


def test_grouped_pt_calibration_configuration_round_trips(tmp_path):
    problem = GroupedPTCalibrationProblem(
        ConstantCalibrationForwardModel(),
        [_experiment("one", "lab-a"), _experiment("two", "lab-b")],
        "study",
        pressure_scale_hyperprior=100.0,
        temperature_scale_hyperprior=10.0,
    )

    write_grouped_pt_calibration_configuration(problem, tmp_path)
    restored = read_grouped_pt_calibration_configuration(tmp_path)

    assert restored["configuration"]["groups"] == ["lab-a", "lab-b"]
    assert restored["configuration"]["parameter_order"] == [
        "pressure_bias[lab-a]",
        "pressure_bias[lab-b]",
        "temperature_bias[lab-a]",
        "temperature_bias[lab-b]",
        "pressure_bias_scale",
        "temperature_bias_scale",
    ]
    assert len(restored["experiments"]) == 2

    inferred_problem = GroupedPTCalibrationProblem(
        ConstantCalibrationForwardModel(),
        [_experiment("one", "lab-a"), _experiment("two", "lab-b")],
        "study",
        pressure_scale_hyperprior=100.0,
        temperature_scale_hyperprior=10.0,
        thermodynamic_parameters=[
            ThermodynamicParameter("dH_phase", "phase", prior_sigma_j_mol=20.0)
        ],
        thermodynamic_prior_covariance=np.asarray([[400.0]]),
    )
    inferred_path = tmp_path / "with_thermo"
    write_grouped_pt_calibration_configuration(inferred_problem, inferred_path)
    inferred_archive = read_grouped_pt_calibration_configuration(inferred_path)
    assert inferred_archive["configuration"]["parameter_order"][0] == "dH_phase"
    assert inferred_archive["configuration"]["thermodynamic_prior_covariance"] == [
        [400.0]
    ]


def test_grouped_pt_calibration_jointly_infers_thermodynamic_corrections():
    problem = GroupedPTCalibrationProblem(
        ConstantCalibrationForwardModel(),
        [_experiment("one", "lab-a"), _experiment("two", "lab-b")],
        "study",
        pressure_scale_hyperprior=100.0,
        temperature_scale_hyperprior=10.0,
        thermodynamic_parameters=[
            ThermodynamicParameter("dH_phase", "phase", prior_sigma_j_mol=20.0)
        ],
        thermodynamic_prior_covariance=np.asarray([[400.0]]),
    )

    data, prior = problem.objective_components(problem.initial_point)
    posterior = mcmc_sample(
        problem,
        problem.initial_point,
        n_samples=18,
        burn_in=3,
        seed=17,
    )

    assert len(problem.parameters) == 7
    assert np.isclose(data, 0.0)
    assert np.isfinite(prior)
    assert posterior.parameter_names[0] == "dH_phase"
    assert posterior.samples.shape == (15, 7)
    assert problem.corrections_for_values(problem.initial_point)["phase"] == 0.0


def test_grouped_pt_calibration_map_is_bounded_and_improves_center():
    problem = GroupedPTCalibrationProblem(
        ConstantCalibrationForwardModel(),
        [_experiment("one", "lab-a"), _experiment("two", "lab-b")],
        "study",
        pressure_scale_hyperprior=100.0,
        temperature_scale_hyperprior=10.0,
    )
    center_objective = problem.objective(problem.initial_point)

    result = problem.optimize(starts=3, iterations=12, seed=8)

    values = np.asarray(
        [result.corrections_j_mol[parameter.name] for parameter in problem.parameters]
    )
    bounds = np.asarray([parameter.bounds() for parameter in problem.parameters])
    assert result.objective <= center_objective
    assert np.all(values >= bounds[:, 0])
    assert np.all(values <= bounds[:, 1])
    assert len(result.starts) == 3


def test_grouped_pt_calibration_multichain_reports_diagnostics():
    problem = GroupedPTCalibrationProblem(
        ConstantCalibrationForwardModel(),
        [_experiment("one", "lab-a"), _experiment("two", "lab-b")],
        "study",
        pressure_scale_hyperprior=100.0,
        temperature_scale_hyperprior=10.0,
    )

    result = mcmc_sample_grouped_pt_calibration(
        problem,
        n_chains=3,
        n_samples=24,
        burn_in=4,
        seed=31,
    )

    assert len(result.chains) == 3
    assert result.parameter_names == tuple(
        parameter.name for parameter in problem.parameters
    )
    assert set(result.r_hat) == set(result.parameter_names)
    assert all(value >= 1.0 for value in result.r_hat.values())
    assert all(value > 0.0 for value in result.effective_sample_size.values())
    assert len(set(result.seeds)) == 3


def test_grouped_calibration_jointly_samples_observation_bias_and_model_error_scale(
    tmp_path,
):
    experiments = [
        replace(
            _experiment("one", "lab-a"),
            observation_covariance=[[0.01]],
            observation_covariance_keys=("fraction:phase",),
            theoretical_covariance=[[0.02]],
            theoretical_covariance_keys=("fraction:phase",),
        ),
        replace(
            _experiment("two", "lab-b"),
            observation_covariance=[[0.01]],
            observation_covariance_keys=("fraction:phase",),
            theoretical_covariance=[[0.02]],
            theoretical_covariance_keys=("fraction:phase",),
        ),
    ]
    problem = GroupedPTCalibrationProblem(
        ConstantCalibrationForwardModel(),
        experiments,
        "study",
        pressure_scale_hyperprior=100.0,
        temperature_scale_hyperprior=10.0,
        observation_bias_keys=("fraction:phase",),
        observation_bias_prior_sigma=0.1,
        observation_bias_bounds=(-0.5, 0.5),
        model_error_scale_prior_sigma=0.5,
        model_error_scale_bounds=(0.01, 3.0),
    )

    data, prior = problem.objective_components(problem.initial_point)
    posterior = mcmc_sample(
        problem,
        problem.initial_point,
        n_samples=20,
        burn_in=3,
        seed=18,
    )

    assert len(problem.parameters) == 10
    assert problem.parameters[4].name == "bias[fraction:phase|lab-a]"
    assert problem.parameters[5].name == "bias[fraction:phase|lab-b]"
    assert problem.parameters[6].name == "scale_C_T[lab-a]"
    assert problem.parameters[7].name == "scale_C_T[lab-b]"
    assert np.isfinite(data) and np.isfinite(prior)
    assert posterior.samples.shape == (17, 10)
    write_grouped_pt_calibration_configuration(problem, tmp_path)
    archived = read_grouped_pt_calibration_configuration(tmp_path)["configuration"]
    assert archived["observation_bias_keys"] == ["fraction:phase"]
    assert archived["model_error_scale_groups"] == ["lab-a", "lab-b"]
    assert archived["parameter_order"] == [
        parameter.name for parameter in problem.parameters
    ]


def test_grouped_calibration_can_jointly_fit_dirichlet_concentration_and_bias(tmp_path):
    class CompositionCalibrationForwardModel:
        def predict(self, experiment, corrections):
            return EquilibriumResult(
                True,
                ("phase",),
                {"phase": 1.0},
                {"phase": {"Zn": 0.7, "Mg": 0.3}},
                {},
            )

    experiments = [
        Experiment(
            identifier,
            2400.0,
            580.0,
            {"Zn": 1.0},
            ("phase",),
            phase_compositions={"phase": {"Zn": 0.8, "Mg": 0.2}},
            metadata={"study": group},
        )
        for identifier, group in (("one", "lab-a"), ("two", "lab-b"))
    ]
    problem = GroupedPTCalibrationProblem(
        CompositionCalibrationForwardModel(),
        experiments,
        "study",
        pressure_scale_hyperprior=100.0,
        temperature_scale_hyperprior=10.0,
        misfit_form="dirichlet_nll",
        dirichlet_concentration_parameters=[
            DirichletConcentrationParameter(
                "phase",
                prior_mean=25.0,
                prior_sigma=8.0,
                lower_bound=0.5,
                upper_bound=120.0,
            )
        ],
        observation_bias_keys=("composition:phase:Zn",),
        observation_bias_prior_sigma=0.05,
        observation_bias_bounds=(-0.2, 0.2),
    )

    data, prior = problem.objective_components(problem.initial_point)
    posterior = mcmc_sample(
        problem,
        problem.initial_point,
        n_samples=16,
        burn_in=3,
        seed=26,
    )

    assert np.isfinite(data) and np.isfinite(prior)
    assert problem.parameters[0].name == "dirichlet_kappa[phase]"
    assert "bias[composition:phase:Zn|lab-a]" in posterior.parameter_names
    assert posterior.samples.shape[1] == len(problem.parameters)

    write_grouped_pt_calibration_configuration(problem, tmp_path)
    archived = read_grouped_pt_calibration_configuration(tmp_path)["configuration"]
    assert archived["dirichlet_concentration_parameter_count"] == 1
    assert archived["parameter_order"] == [
        parameter.name for parameter in problem.parameters
    ]


def test_grouped_multichain_starts_respect_correlated_thermodynamic_prior():
    problem = GroupedPTCalibrationProblem(
        ConstantCalibrationForwardModel(),
        [_experiment("one", "lab-a"), _experiment("two", "lab-b")],
        "study",
        pressure_scale_hyperprior=100.0,
        temperature_scale_hyperprior=10.0,
        thermodynamic_parameters=[
            ThermodynamicParameter("dH_one", "phase-one", prior_sigma_j_mol=10.0),
            ThermodynamicParameter("dH_two", "phase-two", prior_sigma_j_mol=20.0),
        ],
        thermodynamic_prior_covariance=np.asarray([[100.0, 150.0], [150.0, 400.0]]),
    )

    first = mcmc_sample_grouped_pt_calibration(
        problem, n_chains=4, n_samples=4, burn_in=0, seed=73
    )
    repeated = mcmc_sample_grouped_pt_calibration(
        problem, n_chains=4, n_samples=4, burn_in=0, seed=73
    )
    correlation_check = mcmc_sample_grouped_pt_calibration(
        problem, n_chains=32, n_samples=2, burn_in=0, seed=73
    )
    starts = np.asarray(first.initial_points)
    correlated_starts = np.asarray(correlation_check.initial_points)

    assert np.array_equal(starts, np.asarray(repeated.initial_points))
    assert np.corrcoef(correlated_starts[:, 0], correlated_starts[:, 1])[0, 1] > 0.4
    bounds = np.asarray([parameter.bounds() for parameter in problem.parameters])
    assert np.all(starts >= bounds[:, 0])
    assert np.all(starts <= bounds[:, 1])
    assert np.all(np.isfinite(starts))


def test_repeated_grouped_pt_validation_predicts_unseen_group_offsets():
    experiments = [
        replace(
            _experiment("one", "lab-a"),
            observation_covariance=[[0.01]],
            observation_covariance_keys=("fraction:phase",),
            theoretical_covariance=[[0.02]],
            theoretical_covariance_keys=("fraction:phase",),
        ),
        replace(
            _experiment("two", "lab-b"),
            observation_covariance=[[0.01]],
            observation_covariance_keys=("fraction:phase",),
            theoretical_covariance=[[0.02]],
            theoretical_covariance_keys=("fraction:phase",),
        ),
    ]
    problem = GroupedPTCalibrationProblem(
        ConstantCalibrationForwardModel(),
        experiments,
        "study",
        pressure_scale_hyperprior=100.0,
        temperature_scale_hyperprior=10.0,
        thermodynamic_parameters=[
            ThermodynamicParameter("dH_phase", "phase", prior_sigma_j_mol=20.0)
        ],
        observation_bias_keys=("fraction:phase",),
        observation_bias_prior_sigma=0.1,
        model_error_scale_prior_sigma=0.5,
        model_error_scale_bounds=(0.01, 3.0),
    )

    result = repeated_grouped_pt_calibration_validation(
        problem,
        repeats=2,
        test_fraction=0.5,
        n_samples=16,
        burn_in=4,
        predictive_draws=3,
        seed=22,
    )

    assert len(result.folds) == 2
    assert all(len(fold.held_out_experiment_ids) == 1 for fold in result.folds)
    assert all(np.isfinite(fold.misfit) for fold in result.folds)
    assert all(
        len(fold.held_out_log_predictive_densities) == 1 for fold in result.folds
    )
    assert all(
        np.isfinite(score)
        for fold in result.folds
        for score in fold.held_out_log_predictive_densities
    )
    assert np.isfinite(result.mean_log_predictive_density)
    assert all("pressure_bias_scale" in fold.corrections_j_mol for fold in result.folds)
    assert all("dH_phase" in fold.corrections_j_mol for fold in result.folds)
    assert all(
        "bias[fraction:phase|lab-a]" in fold.corrections_j_mol
        or "bias[fraction:phase|lab-b]" in fold.corrections_j_mol
        for fold in result.folds
    )


def test_repeated_grouped_pt_validation_preserves_impossible_log_scores(monkeypatch):
    problem = GroupedPTCalibrationProblem(
        ConstantCalibrationForwardModel(),
        [_experiment("one", "lab-a"), _experiment("two", "lab-b")],
        "study",
        pressure_scale_hyperprior=100.0,
        temperature_scale_hyperprior=10.0,
    )
    monkeypatch.setattr(
        validation_module, "experiment_log_likelihood", lambda *args, **kwargs: -np.inf
    )

    result = repeated_grouped_pt_calibration_validation(
        problem,
        repeats=1,
        test_fraction=0.5,
        n_samples=8,
        burn_in=2,
        predictive_draws=2,
        seed=22,
        n_chains=2,
    )

    assert result.mean_log_predictive_density == -np.inf
