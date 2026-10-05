"""Tests for PTInversionProblem: geothermobarometry (P-T from mineral
compositions) reusing the same likelihood/prior/posterior machinery as
InversionProblem.
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parents[1]))

from thermoinvert.data import (
    EquilibriumBracket,
    Experiment,
    Parameter,
    ScalarParameter,
    ThermodynamicParameter,
)
from thermoinvert.forward_model import EquilibriumResult
from thermoinvert.identifiability import (
    identifiability_summary,
    information_gain_summary,
)
from thermoinvert.likelihood import bracket_misfit
from thermoinvert.posterior import laplace_covariance, mcmc_sample
from thermoinvert.pt_inversion import PTCalibrationBiasProblem, PTInversionProblem
from thermoinvert.archive import (
    read_pt_inversion_configuration,
    write_pt_inversion_configuration,
)
import json
from thermoinvert.sensitivity import numerical_pt_sensitivity
from thermoinvert.validation import (
    posterior_predictive_check_pt,
    posterior_predictive_check_pt_brackets,
)

TRUE_PRESSURE_BAR = 3200.0
TRUE_TEMPERATURE_C = 650.0


class LinearPTForwardModel:
    """Two independent mineral-composition observables, one mostly pressure-
    sensitive and one mostly temperature-sensitive -- a stand-in for a real
    paired barometer + thermometer reaction. A single composition constraint
    generally cannot separate P from T (only their combination is
    identifiable), so this deliberately mirrors standard practice."""

    def predict(self, experiment, corrections):
        dp = experiment.pressure_bar - TRUE_PRESSURE_BAR
        dt = experiment.temperature_c - TRUE_TEMPERATURE_C
        return EquilibriumResult(
            True,
            ("Grt", "Bt"),
            {"Grt": 1.0, "Bt": 1.0},
            {
                "Grt": {"Mg": 0.5 + dp / 50000.0},
                "Bt": {"Fe": 0.4 + dt / 5000.0},
            },
            {},
        )


def _observed_experiment() -> Experiment:
    return Experiment(
        id="natural-sample",
        pressure_bar=2800.0,
        temperature_c=600.0,
        bulk_composition={},
        observed_phases=("Grt", "Bt"),
        phase_compositions={"Grt": {"Mg": 0.5}, "Bt": {"Fe": 0.4}},
        phase_composition_sigma={"Grt": {"Mg": 0.002}, "Bt": {"Fe": 0.002}},
    )


def _pt_problem(**kwargs) -> PTInversionProblem:
    pressure_prior = ThermodynamicParameter(
        "pressure_bar",
        "pressure",
        prior_mean_j_mol=3000.0,
        prior_sigma_j_mol=1000.0,
        lower_bound_j_mol=500.0,
        upper_bound_j_mol=8000.0,
    )
    temperature_prior = ThermodynamicParameter(
        "temperature_c",
        "temperature",
        prior_mean_j_mol=600.0,
        prior_sigma_j_mol=150.0,
        lower_bound_j_mol=300.0,
        upper_bound_j_mol=1000.0,
    )
    return PTInversionProblem(
        LinearPTForwardModel(),
        _observed_experiment(),
        pressure_prior,
        temperature_prior,
        **kwargs,
    )


def test_optimize_recovers_true_pressure_and_temperature():
    problem = _pt_problem()
    result = problem.optimize(starts=6, iterations=100)
    assert abs(result.corrections_j_mol["pressure_bar"] - TRUE_PRESSURE_BAR) < 50.0
    assert abs(result.corrections_j_mol["temperature_c"] - TRUE_TEMPERATURE_C) < 10.0


def test_pt_calibration_bias_problem_corrects_reported_conditions():
    problem = PTCalibrationBiasProblem(
        _pt_problem(),
        Parameter("pressure_bias", 0.0, 500.0, -1000.0, 1000.0, "bar"),
        Parameter("temperature_bias", 0.0, 50.0, -100.0, 100.0, "degC"),
    )

    data, prior = problem.objective_components([2800.0, 600.0, 400.0, 50.0])

    assert data < 1.0e-10
    assert prior > 0.0

    fit = problem.optimize(starts=2, iterations=80)
    corrected_pressure = (
        fit.corrections_j_mol["pressure_bar"] + fit.corrections_j_mol["pressure_bias"]
    )
    corrected_temperature = (
        fit.corrections_j_mol["temperature_c"]
        + fit.corrections_j_mol["temperature_bias"]
    )
    assert abs(corrected_pressure - TRUE_PRESSURE_BAR) < 50.0
    assert abs(corrected_temperature - TRUE_TEMPERATURE_C) < 10.0


def test_pt_calibration_bias_configuration_round_trips(tmp_path):
    problem = PTCalibrationBiasProblem(
        _pt_problem(prior_covariance=np.asarray([[40000.0, 140.0], [140.0, 25.0]])),
        Parameter("pressure_bias", 0.0, 50.0, -200.0, 200.0, "bar"),
        Parameter("temperature_bias", 0.0, 5.0, -20.0, 20.0, "degC"),
    )

    write_pt_inversion_configuration(problem, tmp_path)
    archived = read_pt_inversion_configuration(tmp_path)
    configuration = archived["configuration"]

    assert configuration["calibration_bias"] is True
    assert configuration["parameter_order"] == [
        "pressure_bar",
        "temperature_c",
        "pressure_bias",
        "temperature_bias",
    ]
    assert np.allclose(configuration["prior_covariance"], problem.prior_covariance)
    assert configuration["parameters"][2]["unit"] == "bar"
    assert configuration["parameters"][3]["unit"] == "degC"


def test_pt_optimize_supports_hot_start_and_validates_central_point():
    problem = _pt_problem()
    result = problem.optimize(
        starts=4,
        iterations=80,
        central_point=[TRUE_PRESSURE_BAR - 40.0, TRUE_TEMPERATURE_C - 5.0],
        hot_start_spread_fraction=0.05,
    )
    assert abs(result.corrections_j_mol["pressure_bar"] - TRUE_PRESSURE_BAR) < 100.0
    assert abs(result.corrections_j_mol["temperature_c"] - TRUE_TEMPERATURE_C) < 20.0

    try:
        problem.optimize(central_point=[TRUE_PRESSURE_BAR])
    except ValueError:
        pass
    else:
        raise AssertionError("Expected ValueError for a mismatched PT central_point")

    try:
        problem.optimize(hot_start_spread_fraction=0.0)
    except ValueError:
        pass
    else:
        raise AssertionError(
            "Expected ValueError for a non-positive PT hot-start spread"
        )


def test_laplace_covariance_runs_directly_on_pt_problem():
    problem = _pt_problem()
    map_point = [TRUE_PRESSURE_BAR, TRUE_TEMPERATURE_C]
    approximation = laplace_covariance(problem, map_point)
    assert approximation.parameter_names == ("pressure_bar", "temperature_c")
    assert approximation.hessian.shape == (2, 2)
    assert approximation.is_positive_definite


def test_mcmc_sample_runs_directly_on_pt_problem():
    problem = _pt_problem()
    result = mcmc_sample(
        problem,
        [TRUE_PRESSURE_BAR, TRUE_TEMPERATURE_C],
        n_samples=6000,
        burn_in=1000,
        proposal_sigma_j_mol=[70.0, 7.0],
        seed=0,
    )
    assert result.samples.shape[1] == 2
    assert result.samples.shape[0] == 5000
    assert result.parameter_units == ("J/mol", "J/mol")
    assert abs(result.posterior_mean["pressure_bar"] - TRUE_PRESSURE_BAR) < 100.0
    assert abs(result.posterior_mean["temperature_c"] - TRUE_TEMPERATURE_C) < 20.0


def test_posterior_predictive_check_runs_directly_on_pt_problem():
    problem = _pt_problem()
    samples = mcmc_sample(
        problem,
        [TRUE_PRESSURE_BAR, TRUE_TEMPERATURE_C],
        n_samples=4000,
        burn_in=500,
        proposal_sigma_j_mol=[70.0, 7.0],
        seed=2,
    )
    result = posterior_predictive_check_pt(problem, samples, n_predictive_samples=100)
    grt_index = result.observable_names.index("composition:Grt:Mg")
    bt_index = result.observable_names.index("composition:Bt:Fe")
    assert abs(result.predictive_mean[grt_index] - 0.5) < 0.01
    assert abs(result.predictive_mean[bt_index] - 0.4) < 0.01
    assert result.within_95["composition:Grt:Mg"]
    assert result.within_95["composition:Bt:Fe"]


def test_pt_objective_rejects_nonfinite_and_out_of_bounds_values():
    problem = _pt_problem()
    assert np.isinf(problem.objective([np.nan, TRUE_TEMPERATURE_C]))
    assert np.isinf(problem.objective([TRUE_PRESSURE_BAR + 1.0e6, TRUE_TEMPERATURE_C]))


def test_pt_inversion_supports_correlated_pressure_temperature_prior():
    problem = _pt_problem(
        prior_covariance_j_mol2=np.asarray([[40000.0, 140.0], [140.0, 25.0]])
    )
    values = np.asarray([TRUE_PRESSURE_BAR + 10.0, TRUE_TEMPERATURE_C + 2.0])
    means = np.asarray([parameter.prior_mean_j_mol for parameter in problem.parameters])
    deviations = values - means
    expected_prior = (
        0.5 * deviations @ np.linalg.solve(problem.prior_covariance_j_mol2, deviations)
    )

    assert np.isclose(problem.objective_components(values)[1], expected_prior)


def test_pt_inversion_accepts_native_unit_scalar_parameters():
    problem = PTInversionProblem(
        LinearPTForwardModel(),
        _observed_experiment(),
        ScalarParameter(
            "pressure_bar",
            prior_mean=3000.0,
            prior_sigma=1000.0,
            lower_bound=500.0,
            upper_bound=8000.0,
            unit="bar",
        ),
        ScalarParameter(
            "temperature_c",
            prior_mean=600.0,
            prior_sigma=150.0,
            lower_bound=300.0,
            upper_bound=1000.0,
            unit="degC",
        ),
    )

    result = problem.optimize(starts=3, iterations=40)
    posterior = mcmc_sample(
        problem,
        [3000.0, 600.0],
        n_samples=20,
        burn_in=2,
        proposal_sigma_j_mol=[10.0, 2.0],
    )

    assert result.corrections_j_mol["pressure_bar"] > 0.0
    assert result.corrections_j_mol["temperature_c"] > 0.0
    assert posterior.parameter_units == ("bar", "degC")


def test_pt_inversion_accepts_generic_parameter_base_and_archives_it(tmp_path):
    problem = PTInversionProblem(
        LinearPTForwardModel(),
        _observed_experiment(),
        Parameter("pressure_bar", 3000.0, 1000.0, 500.0, 8000.0, "bar"),
        Parameter("temperature_c", 600.0, 150.0, 300.0, 1000.0, "degC"),
    )
    result = problem.optimize(starts=2, iterations=10)

    write_pt_inversion_configuration(problem, tmp_path)
    archived = read_pt_inversion_configuration(tmp_path)

    assert result.corrections_j_mol["pressure_bar"] > 0.0
    assert archived["configuration"]["parameters"][0]["parameter_kind"] == "scalar"
    assert archived["configuration"]["parameters"][0]["unit"] == "bar"


def test_scalar_parameter_rejects_invalid_native_unit_prior():
    for kwargs in (
        {"prior_sigma": 0.0},
        {"prior_mean": np.nan},
        {"unit": ""},
        {"lower_bound": 2.0, "upper_bound": 1.0},
    ):
        try:
            ScalarParameter("pressure", **kwargs)
        except ValueError:
            pass
        else:
            raise AssertionError("Expected ValueError for invalid ScalarParameter")


def test_pt_inversion_configuration_archives_correlated_prior(tmp_path):
    problem = _pt_problem(
        prior_covariance_j_mol2=np.asarray([[40000.0, 140.0], [140.0, 25.0]])
    )

    write_pt_inversion_configuration(problem, tmp_path)

    configuration = json.loads(
        (tmp_path / "pt_inversion_configuration.json").read_text(encoding="utf-8")
    )
    assert configuration["parameter_order"] == ["pressure_bar", "temperature_c"]
    assert configuration["prior_covariance_j_mol2"] == [[40000.0, 140.0], [140.0, 25.0]]
    assert configuration["problem_mode"] == "composition"
    assert (tmp_path / "experiment.json").is_file()
    restored = read_pt_inversion_configuration(tmp_path)
    assert restored["configuration"]["prior_covariance_j_mol2"] == [
        [40000.0, 140.0],
        [140.0, 25.0],
    ]


def test_pt_scalar_configuration_archive_preserves_native_units(tmp_path):
    problem = PTInversionProblem(
        LinearPTForwardModel(),
        _observed_experiment(),
        ScalarParameter("pressure_bar", 3000.0, 1000.0, 500.0, 8000.0, "bar"),
        ScalarParameter("temperature_c", 600.0, 150.0, 300.0, 1000.0, "degC"),
    )

    write_pt_inversion_configuration(problem, tmp_path)
    restored = read_pt_inversion_configuration(tmp_path)

    assert restored["configuration"]["parameters"][0]["parameter_kind"] == "scalar"
    assert restored["configuration"]["parameters"][0]["unit"] == "bar"


def test_identifiability_and_information_gain_run_directly_on_pt_problem():
    problem = _pt_problem()
    map_point = [TRUE_PRESSURE_BAR, TRUE_TEMPERATURE_C]
    all_result = laplace_covariance(problem, map_point, uncertainty_source="all")
    prior_result = laplace_covariance(problem, map_point, uncertainty_source="prior")

    summary = identifiability_summary(problem, all_result, prior_result)
    assert summary["pressure_bar"].status in (
        "data_constrained",
        "ambiguous",
        "prior_dominated",
    )

    gains = information_gain_summary(problem, all_result)
    assert np.isfinite(gains["pressure_bar"])
    assert np.isfinite(gains["temperature_c"])


def test_bracket_mode_recovers_p_t_reaction_boundary():
    transition_pressure = 3200.0
    transition_temperature = 650.0

    class PurePhaseForwardModel:
        def predict(self, experiment, corrections):
            phase = (
                "and"
                if (
                    experiment.pressure_bar < transition_pressure
                    and experiment.temperature_c < transition_temperature
                )
                else "ky"
            )
            return EquilibriumResult(True, (phase,), {phase: 1.0}, {}, {})

    bracket = EquilibriumBracket(
        "and-ky-boundary",
        Experiment("lower", 3150.0, 640.0, {}, ("and",)),
        Experiment("upper", 3250.0, 660.0, {}, ("ky",)),
        lower_required_phases=("and",),
        upper_required_phases=("ky",),
        lower_forbidden_phases=("ky",),
        upper_forbidden_phases=("and",),
    )
    problem = PTInversionProblem(
        PurePhaseForwardModel(),
        bracket.lower,
        ThermodynamicParameter(
            "pressure_bar",
            "pressure",
            prior_mean_j_mol=3000.0,
            prior_sigma_j_mol=1000.0,
            lower_bound_j_mol=2500.0,
            upper_bound_j_mol=4000.0,
        ),
        ThermodynamicParameter(
            "temperature_c",
            "temperature",
            prior_mean_j_mol=600.0,
            prior_sigma_j_mol=100.0,
            lower_bound_j_mol=500.0,
            upper_bound_j_mol=800.0,
        ),
        bracket=bracket,
    )

    result = problem.optimize(starts=4, iterations=80)
    candidate = problem._candidate_bracket(
        [
            result.corrections_j_mol["pressure_bar"],
            result.corrections_j_mol["temperature_c"],
        ]
    )
    lower_prediction = problem.forward_model.predict(candidate.lower, {})
    upper_prediction = problem.forward_model.predict(candidate.upper, {})
    assert result.objective < 1.0
    assert bracket_misfit(bracket, lower_prediction, upper_prediction) == 0.0
    assert abs(result.corrections_j_mol["pressure_bar"] - transition_pressure) < 100.0


def test_multiple_brackets_recover_unique_pressure_and_temperature():
    target_pressure = 3200.0
    target_temperature = 650.0

    class TwoBoundaryForwardModel:
        def predict(self, experiment, corrections):
            reaction = experiment.metadata["reaction"]
            if reaction == "pressure":
                phase = "lowP" if experiment.pressure_bar < target_pressure else "highP"
            else:
                phase = (
                    "lowT" if experiment.temperature_c < target_temperature else "highT"
                )
            return EquilibriumResult(True, (phase,), {phase: 1.0}, {}, {})

    pressure_bracket = EquilibriumBracket(
        "pressure-boundary",
        Experiment(
            "p-low", 3150.0, 650.0, {}, ("lowP",), metadata={"reaction": "pressure"}
        ),
        Experiment(
            "p-high", 3250.0, 650.0, {}, ("highP",), metadata={"reaction": "pressure"}
        ),
        lower_required_phases=("lowP",),
        upper_required_phases=("highP",),
        lower_forbidden_phases=("highP",),
        upper_forbidden_phases=("lowP",),
    )
    temperature_bracket = EquilibriumBracket(
        "temperature-boundary",
        Experiment(
            "t-low", 3200.0, 640.0, {}, ("lowT",), metadata={"reaction": "temperature"}
        ),
        Experiment(
            "t-high",
            3200.0,
            660.0,
            {},
            ("highT",),
            metadata={"reaction": "temperature"},
        ),
        lower_required_phases=("lowT",),
        upper_required_phases=("highT",),
        lower_forbidden_phases=("highT",),
        upper_forbidden_phases=("lowT",),
    )
    problem = PTInversionProblem(
        TwoBoundaryForwardModel(),
        pressure_bracket.lower,
        ThermodynamicParameter(
            "pressure_bar",
            "pressure",
            prior_mean_j_mol=3200.0,
            prior_sigma_j_mol=1000.0,
            lower_bound_j_mol=2500.0,
            upper_bound_j_mol=4000.0,
        ),
        ThermodynamicParameter(
            "temperature_c",
            "temperature",
            prior_mean_j_mol=650.0,
            prior_sigma_j_mol=100.0,
            lower_bound_j_mol=500.0,
            upper_bound_j_mol=800.0,
        ),
        brackets=[pressure_bracket, temperature_bracket],
    )

    target_data, _ = problem.objective_components([target_pressure, target_temperature])
    off_target_data, _ = problem.objective_components([2800.0, 600.0])
    assert target_data == 0.0
    assert off_target_data > target_data

    result = problem.optimize(starts=4, iterations=80)
    assert abs(result.corrections_j_mol["pressure_bar"] - target_pressure) < 100.0
    assert abs(result.corrections_j_mol["temperature_c"] - target_temperature) < 20.0


def test_bracket_posterior_predictive_check_reports_satisfaction_probability():
    target_pressure = 3200.0
    target_temperature = 650.0

    class BoundaryForwardModel:
        def predict(self, experiment, corrections):
            reaction = experiment.metadata["reaction"]
            if reaction == "pressure":
                phase = "lowP" if experiment.pressure_bar < target_pressure else "highP"
            else:
                phase = (
                    "lowT" if experiment.temperature_c < target_temperature else "highT"
                )
            return EquilibriumResult(True, (phase,), {phase: 1.0}, {}, {})

    brackets = [
        EquilibriumBracket(
            "pressure-boundary",
            Experiment(
                "p-low", 3150.0, 650.0, {}, ("lowP",), metadata={"reaction": "pressure"}
            ),
            Experiment(
                "p-high",
                3250.0,
                650.0,
                {},
                ("highP",),
                metadata={"reaction": "pressure"},
            ),
            lower_required_phases=("lowP",),
            upper_required_phases=("highP",),
            lower_forbidden_phases=("highP",),
            upper_forbidden_phases=("lowP",),
        ),
        EquilibriumBracket(
            "temperature-boundary",
            Experiment(
                "t-low",
                3200.0,
                640.0,
                {},
                ("lowT",),
                metadata={"reaction": "temperature"},
            ),
            Experiment(
                "t-high",
                3200.0,
                660.0,
                {},
                ("highT",),
                metadata={"reaction": "temperature"},
            ),
            lower_required_phases=("lowT",),
            upper_required_phases=("highT",),
            lower_forbidden_phases=("highT",),
            upper_forbidden_phases=("lowT",),
        ),
    ]
    problem = PTInversionProblem(
        BoundaryForwardModel(),
        brackets[0].lower,
        ThermodynamicParameter("pressure_bar", "pressure", prior_mean_j_mol=3200.0),
        ThermodynamicParameter("temperature_c", "temperature", prior_mean_j_mol=650.0),
        brackets=brackets,
    )
    from thermoinvert.posterior import MCMCResult

    samples = np.asarray([[3200.0, 650.0], [3200.0, 650.0], [2800.0, 600.0]])
    mcmc_result = MCMCResult(
        ("pressure_bar", "temperature_c"),
        samples,
        np.zeros(3),
        1.0,
        {"pressure_bar": 3200.0, "temperature_c": 650.0},
        {"pressure_bar": 1.0, "temperature_c": 1.0},
        {"pressure_bar": (2800.0, 3200.0), "temperature_c": (600.0, 650.0)},
    )
    result = posterior_predictive_check_pt_brackets(problem, mcmc_result, 3)
    assert result.bracket_ids == ("pressure-boundary", "temperature-boundary")
    assert result.samples_evaluated == 3
    assert result.satisfied_fraction == 2.0 / 3.0
    assert result.mean_misfit > 0.0


def test_numerical_pt_sensitivity_matches_linear_forward_model():
    result = numerical_pt_sensitivity(
        _pt_problem(), [TRUE_PRESSURE_BAR, TRUE_TEMPERATURE_C]
    )
    composition_rows = {
        name: row for name, row in zip(result.observable_names, result.jacobian)
    }
    assert np.allclose(
        composition_rows["natural-sample:composition:Grt:Mg"],
        [1.0 / 50000.0, 0.0],
        atol=1.0e-12,
    )
    assert np.allclose(
        composition_rows["natural-sample:composition:Bt:Fe"],
        [0.0, 1.0 / 5000.0],
        atol=1.0e-12,
    )
    assert result.identifiability.rank == 2


def test_numerical_pt_sensitivity_rejects_discontinuous_bracket_problem():
    bracket = EquilibriumBracket(
        "and-ky",
        Experiment("low", 3150.0, 640.0, {}, ("and",)),
        Experiment("high", 3250.0, 660.0, {}, ("ky",)),
    )
    problem = PTInversionProblem(
        LinearPTForwardModel(),
        bracket.lower,
        ThermodynamicParameter("pressure_bar", "pressure"),
        ThermodynamicParameter("temperature_c", "temperature"),
        bracket=bracket,
    )
    try:
        numerical_pt_sensitivity(problem, [3200.0, 650.0])
    except ValueError:
        pass
    else:
        raise AssertionError("Expected ValueError for a bracket-only P-T problem")

    samples = mcmc_sample(
        _pt_problem(),
        [TRUE_PRESSURE_BAR, TRUE_TEMPERATURE_C],
        n_samples=20,
        burn_in=2,
        proposal_sigma_j_mol=[70.0, 7.0],
        seed=3,
    )
    try:
        posterior_predictive_check_pt(problem, samples)
    except ValueError:
        pass
    else:
        raise AssertionError("Expected ValueError for a bracket-only P-T problem")
