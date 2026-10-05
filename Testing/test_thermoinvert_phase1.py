import json
import sys
from dataclasses import replace
from math import lgamma
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parents[1]))

from thermoinvert.data import (
    DirichletConcentrationParameter,
    EquilibriumBracket,
    Experiment,
    Parameter,
    ThermodynamicParameter,
)
from thermoinvert.forward_model import EquilibriumResult
from thermoinvert.inversion import InversionProblem
from thermoinvert.joint_inversion import JointInversionProblem
from thermoinvert.bias import (
    GroupedObservationBiasReport,
    GroupedObservationBiasSummary,
    ObservationBiasProfileResult,
    apply_group_observation_biases,
    profile_grouped_bias_summary,
    profile_observation_bias_by_group,
    profile_observation_biases_by_group,
    summarize_group_observation_bias,
    summarize_grouped_bias_report,
    summarize_grouped_observation_biases,
)
from thermoinvert.hierarchical import (
    mcmc_sample_combined_grouped_uncertainty,
    mcmc_sample_group_observation_biases,
    mcmc_sample_theoretical_covariance_group_scales,
    mcmc_sample_theoretical_covariance_scale,
    profile_theoretical_covariance_group_scales,
    profile_theoretical_covariance_scale,
)
from thermoinvert.likelihood import (
    bracket_misfit,
    experiment_log_likelihood,
    experiment_misfit,
)
from thermoinvert.posterior import (
    filter_by_acceptance_threshold,
    filter_outliers_mcmc,
    grid_posterior,
    laplace_covariance,
    MCMCResult,
    mcmc_sample,
    mcmc_sample_multi_chain,
    ParallelTemperingResult,
    parallel_tempering_sample,
    analyze_parallel_tempering_modes,
    thermodynamic_integration,
    potential_scale_reduction,
    effective_sample_size,
)
from thermoinvert.validation import (
    database_file_provenance,
    leave_one_group_out_cross_validation,
    leave_one_database_family_out_cross_validation,
    repeated_group_holdout_cross_validation,
    leave_one_out_cross_validation,
    leave_one_phase_assemblage_out_cross_validation,
    leave_one_out_of_domain_pt_cross_validation,
    posterior_predictive_check,
    validate_database_provenance,
)
from thermoinvert.sensitivity import (
    SensitivityResult,
    identify_sensitivity,
    prior_scaled_information,
    read_sensitivity_report,
    write_sensitivity_report,
)
from thermoinvert.thermo import DatabaseOverlay


DB_PATH = (
    Path(__file__).parents[1]
    / "embedded/databases/perplex/DEW17HP622_Zn_2025-reaktoro.json"
)


def test_generic_parameter_uses_native_units_and_name_key():
    parameter = Parameter("pressure", 3000.0, 500.0, 500.0, 8000.0, "bar")

    assert parameter.key == "pressure"
    assert parameter.unit == "bar"
    assert parameter.bounds() == (500.0, 8000.0)


def test_inversion_accepts_unit_neutral_covariance_for_mixed_units():
    class MixedUnitForwardModel:
        def predict(self, experiment, corrections):
            return EquilibriumResult(True, ("A",), {"A": 1.0}, {}, {})

    experiment = Experiment(
        "mixed-units",
        2000.0,
        500.0,
        {"A": 1.0},
        ("A",),
        phase_fractions={"A": 1.0},
        phase_fraction_sigma={"A": 1.0},
    )
    covariance = np.asarray([[100.0, 20.0], [20.0, 250000.0]])
    problem = InversionProblem(
        MixedUnitForwardModel(),
        [experiment],
        [
            ThermodynamicParameter("dH_A", "A", prior_sigma_j_mol=10.0),
            Parameter("pressure", 3000.0, 500.0, 500.0, 8000.0, "bar"),
        ],
        prior_covariance=covariance,
    )

    assert np.allclose(problem.prior_covariance, covariance)
    assert problem.prior_covariance is problem.prior_covariance_j_mol2


def test_mcmc_accepts_unit_neutral_proposals_for_mixed_unit_parameters():
    class MixedUnitForwardModel:
        def predict(self, experiment, corrections):
            value = corrections.get("pressure", 3000.0)
            return EquilibriumResult(True, ("A",), {"A": 1.0}, {}, {"p": value})

    experiment = Experiment(
        "mixed-unit-mcmc",
        3000.0,
        500.0,
        {"A": 1.0},
        ("A",),
        phase_fractions={"A": 1.0},
        phase_fraction_sigma={"A": 1.0},
    )
    parameters = [
        Parameter("pressure", 3000.0, 500.0, 500.0, 8000.0, "bar"),
        Parameter("temperature", 500.0, 50.0, 300.0, 900.0, "degC"),
    ]
    problem = InversionProblem(MixedUnitForwardModel(), [experiment], parameters)
    proposal_covariance = np.asarray([[100.0, 2.0], [2.0, 25.0]])

    result = mcmc_sample(
        problem,
        [3000.0, 500.0],
        n_samples=20,
        burn_in=2,
        proposal_covariance=proposal_covariance,
        seed=13,
    )

    assert result.parameter_units == ("bar", "degC")
    assert np.allclose(result.proposal_covariance, proposal_covariance)


def test_joint_theoretical_scale_mcmc_samples_scale_and_corrections():
    class FakeForwardModel:
        def predict(self, experiment, corrections):
            value = corrections.get("A", 0.0)
            return EquilibriumResult(True, ("A",), {"A": 0.5 + value / 100.0}, {}, {})

    experiment = Experiment(
        "hierarchical",
        1000.0,
        500.0,
        {"A": 1.0},
        ("A",),
        phase_fractions={"A": 0.5},
        phase_fraction_sigma={"A": 0.1},
        theoretical_covariance=[[0.01]],
        theoretical_covariance_keys=("fraction:A",),
    )
    problem = InversionProblem(
        FakeForwardModel(),
        [experiment],
        [ThermodynamicParameter("dH_A", "A", prior_sigma_j_mol=10.0)],
    )

    result = mcmc_sample_theoretical_covariance_scale(
        problem,
        [0.0],
        initial_scale=1.0,
        n_samples=40,
        burn_in=10,
        proposal_sigma_j_mol=[1.0],
        proposal_sigma_scale=0.1,
        seed=4,
    )

    assert result.samples.shape == (30, 2)
    assert result.parameter_names == ("dH_A", "scale_C_T")
    assert result.parameter_units == ("J/mol", "dimensionless")
    assert np.all(result.samples[:, 1] >= 0.0)
    assert np.all(result.samples[:, 1] <= 5.0)

    delegated = problem.sample_theoretical_covariance_scale(
        [0.0], n_samples=12, burn_in=2, seed=5
    )
    assert delegated.parameter_names == result.parameter_names


def test_joint_group_theoretical_scale_mcmc_samples_each_metadata_group():
    class GroupScaleForwardModel:
        def predict(self, experiment, corrections):
            return EquilibriumResult(True, ("A",), {"A": 0.5}, {}, {})

    experiments = [
        Experiment(
            f"group-{group}",
            1000.0,
            500.0,
            {"A": 1.0},
            ("A",),
            phase_fractions={"A": 0.5},
            phase_fraction_sigma={"A": 0.1},
            theoretical_covariance=[[0.01]],
            theoretical_covariance_keys=("fraction:A",),
            metadata={"study": group},
        )
        for group in ("one", "two")
    ]
    problem = InversionProblem(
        GroupScaleForwardModel(),
        experiments,
        [ThermodynamicParameter("dH_A", "A", prior_sigma_j_mol=10.0)],
    )
    result = mcmc_sample_theoretical_covariance_group_scales(
        problem,
        "study",
        [0.0],
        {"one": 1.0, "two": 1.0},
        n_samples=30,
        burn_in=5,
        proposal_sigma_j_mol=[1.0],
        proposal_sigma_scale=0.1,
        seed=6,
    )

    assert result.samples.shape == (25, 3)
    assert result.scale_names == ("one", "two")
    assert result.parameter_names[-2:] == ("scale_C_T[one]", "scale_C_T[two]")


def test_joint_group_observation_bias_mcmc_samples_each_metadata_group():
    class BiasForwardModel:
        def predict(self, experiment, corrections):
            return EquilibriumResult(True, ("A",), {"A": 1.0}, {}, {})

    experiments = [
        Experiment(
            f"bias-{group}",
            1000.0,
            500.0,
            {"A": 1.0},
            ("A",),
            phase_fractions={"A": 0.5},
            phase_fraction_sigma={"A": 0.1},
            metadata={"lab": group},
        )
        for group in ("one", "two")
    ]
    problem = InversionProblem(
        BiasForwardModel(),
        experiments,
        [ThermodynamicParameter("dH_A", "A", prior_sigma_j_mol=10.0)],
    )
    result = mcmc_sample_group_observation_biases(
        problem,
        "lab",
        "fraction:A",
        [0.0],
        {"one": 0.0, "two": 0.0},
        n_samples=30,
        burn_in=5,
        proposal_sigma_j_mol=[1.0],
        proposal_sigma_bias=0.01,
        seed=8,
    )

    assert result.samples.shape == (25, 3)
    assert result.bias_names == ("fraction:A|one", "fraction:A|two")
    assert result.observation_keys == ("fraction:A",)
    assert result.parameter_names[-2:] == (
        "bias[fraction:A|one]",
        "bias[fraction:A|two]",
    )


def test_combined_grouped_uncertainty_mcmc_supports_multiple_observation_keys():
    class MultiKeyCombinedForwardModel:
        def predict(self, experiment, corrections):
            return EquilibriumResult(True, ("A",), {"A": 0.5}, {"A": {"x": 0.5}}, {})

    experiments = [
        Experiment(
            f"multi-combined-{group}",
            1000.0,
            500.0,
            {"A": 1.0},
            ("A",),
            phase_fractions={"A": 0.5},
            phase_fraction_sigma={"A": 0.1},
            phase_compositions={"A": {"x": 0.5}},
            phase_composition_sigma={"A": {"x": 0.1}},
            theoretical_covariance=[[0.01, 0.0], [0.0, 0.01]],
            theoretical_covariance_keys=("fraction:A", "composition:A:x"),
            metadata={"study": group},
        )
        for group in ("one", "two")
    ]
    problem = InversionProblem(
        MultiKeyCombinedForwardModel(),
        experiments,
        [ThermodynamicParameter("dH_A", "A", prior_sigma_j_mol=10.0)],
    )
    result = mcmc_sample_combined_grouped_uncertainty(
        problem,
        "study",
        ("fraction:A", "composition:A:x"),
        [0.0],
        {"one": 1.0, "two": 1.0},
        {
            ("fraction:A", "one"): 0.0,
            ("fraction:A", "two"): 0.0,
            ("composition:A:x", "one"): 0.0,
            ("composition:A:x", "two"): 0.0,
        },
        n_samples=20,
        burn_in=5,
        proposal_sigma_j_mol=[1.0],
        seed=14,
    )

    assert result.samples.shape == (15, 7)
    assert result.observation_keys == ("fraction:A", "composition:A:x")
    assert result.parameter_names[-4:] == (
        "bias[('fraction:A', 'one')]",
        "bias[('fraction:A', 'two')]",
        "bias[('composition:A:x', 'one')]",
        "bias[('composition:A:x', 'two')]",
    )


def test_combined_grouped_uncertainty_mcmc_samples_scales_and_biases_together():
    class CombinedForwardModel:
        def predict(self, experiment, corrections):
            return EquilibriumResult(True, ("A",), {"A": 0.5}, {}, {})

    experiments = [
        Experiment(
            f"combined-{group}",
            1000.0,
            500.0,
            {"A": 1.0},
            ("A",),
            phase_fractions={"A": 0.5},
            phase_fraction_sigma={"A": 0.1},
            theoretical_covariance=[[0.01]],
            theoretical_covariance_keys=("fraction:A",),
            metadata={"study": group},
        )
        for group in ("one", "two")
    ]
    problem = InversionProblem(
        CombinedForwardModel(),
        experiments,
        [ThermodynamicParameter("dH_A", "A", prior_sigma_j_mol=10.0)],
    )
    result = mcmc_sample_combined_grouped_uncertainty(
        problem,
        "study",
        "fraction:A",
        [0.0],
        {"one": 1.0, "two": 1.0},
        {"one": 0.0, "two": 0.0},
        n_samples=25,
        burn_in=5,
        proposal_sigma_j_mol=[1.0],
        seed=12,
    )

    assert result.samples.shape == (20, 5)
    assert result.parameter_names[-4:] == (
        "scale_C_T[one]",
        "scale_C_T[two]",
        "bias[one]",
        "bias[two]",
    )


def test_parallel_tempering_reports_modes_and_swap_diagnostics():
    class BimodalForwardModel:
        def predict(self, experiment, corrections):
            value = corrections.get("A", 0.0)
            return EquilibriumResult(True, ("A",), {"A": value}, {}, {})

    experiment = Experiment(
        "bimodal",
        1000.0,
        500.0,
        {"A": 1.0},
        ("A",),
        phase_fractions={"A": 0.0},
        phase_fraction_sigma={"A": 0.05},
    )
    problem = InversionProblem(
        BimodalForwardModel(),
        [experiment],
        [
            ThermodynamicParameter(
                "p",
                "A",
                prior_sigma_j_mol=10.0,
                lower_bound_j_mol=-20.0,
                upper_bound_j_mol=20.0,
            )
        ],
    )
    result = parallel_tempering_sample(
        problem,
        initial_points=[[-10.0], [10.0]],
        temperatures=[1.0, 4.0],
        n_samples=80,
        burn_in=20,
        proposal_sigma=[2.0],
        mode_separation_fraction=0.2,
        seed=3,
    )

    assert result.samples.shape == (60, 1)
    assert result.parameter_names == ("p",)
    assert len(result.temperatures) == 2
    assert len(result.swap_acceptance_rates) == 1
    assert sum(result.mode_occupancy.values()) == 1.0
    assert result.mode_best_objective
    assert result.mode_transition_count >= 0
    assert 0.0 <= result.mode_transition_rate <= 1.0
    assert set(result.temperature_mode_occupancy) == {1.0, 4.0}
    assert all(
        np.isclose(sum(occupancy.values()), 1.0)
        for occupancy in result.temperature_mode_occupancy.values()
    )
    global_modes = analyze_parallel_tempering_modes(
        result, mode_separation_fraction=0.2
    )
    assert global_modes.temperature_mode_labels.shape == (60, 2)
    assert np.array_equal(
        global_modes.mode_labels, global_modes.temperature_mode_labels[:, 0]
    )
    assert np.isclose(sum(global_modes.mode_probability.values()), 1.0)
    assert len(global_modes.clustering_mode_counts) == 3
    assert set(global_modes.mode_probability_sensitivity) == set(
        global_modes.mode_probability
    )
    assert all(
        low <= global_modes.mode_probability[mode] <= high
        for mode, (low, high) in global_modes.mode_probability_sensitivity.items()
    )
    integration = thermodynamic_integration(result)
    assert np.isfinite(integration.log_normalizer_relative_to_uniform_bounds)
    assert integration.uniform_reference_sample_count == 128


def test_thermodynamic_integration_recovers_constant_objective_normalizer():
    result = ParallelTemperingResult(
        ("x",),
        ("1",),
        (1.0, 2.0, 8.0),
        np.zeros((4, 1)),
        np.full(4, -2.5),
        np.zeros(4, dtype=int),
        {0: 1.0},
        {0: {"x": 0.0}},
        {0: 2.5},
        (0.5, 0.5, 0.5),
        (0.5, 0.5),
        False,
        1,
        0,
        (0.1,),
        all_temperature_objectives=np.full((4, 3), 2.5),
        uniform_reference_objectives=np.full(10, 2.5),
    )

    estimate = thermodynamic_integration(result)

    assert np.isclose(estimate.log_normalizer_relative_to_uniform_bounds, -2.5)
    assert estimate.inverse_temperatures == (0.0, 0.125, 0.5, 1.0)
    assert np.isclose(estimate.monte_carlo_standard_error, 0.0)
    assert np.allclose(estimate.monte_carlo_interval_95, (-2.5, -2.5))
    assert np.isclose(
        estimate.quadratic_ladder_log_normalizer,
        estimate.log_normalizer_relative_to_uniform_bounds,
    )
    assert np.isclose(estimate.quadrature_discrepancy, 0.0)
    assert np.isclose(estimate.maximum_inverse_temperature_step, 0.5)


def test_thermodynamic_integration_reports_reproducible_mc_uncertainty():
    objective_values = np.column_stack(
        (
            1.0 + np.sin(np.arange(80) / 5.0),
            2.0 + np.cos(np.arange(80) / 7.0),
        )
    )
    result = ParallelTemperingResult(
        ("x",),
        ("1",),
        (1.0, 4.0),
        np.zeros((80, 1)),
        -objective_values[:, 0],
        np.zeros(80, dtype=int),
        {0: 1.0},
        {0: {"x": 0.0}},
        {0: 1.0},
        (0.5, 0.5),
        (0.25,),
        False,
        2,
        0,
        (0.1,),
        all_temperature_objectives=objective_values,
        uniform_reference_objectives=np.linspace(0.5, 1.5, 60),
    )

    first = thermodynamic_integration(
        result, bootstrap_replicates=120, block_length=8, seed=52
    )
    repeated = thermodynamic_integration(
        result, bootstrap_replicates=120, block_length=8, seed=52
    )

    assert np.isfinite(first.monte_carlo_standard_error)
    assert first.monte_carlo_standard_error > 0.0
    assert first.monte_carlo_interval_95[0] < first.monte_carlo_interval_95[1]
    assert first.monte_carlo_interval_95 == repeated.monte_carlo_interval_95
    assert first.log_normalizer_relative_to_uniform_bounds == (
        repeated.log_normalizer_relative_to_uniform_bounds
    )
    assert first.quadratic_ladder_log_normalizer is not None
    assert first.quadrature_discrepancy >= 0.0


def test_joint_inversion_recovers_thermo_pressure_and_temperature():
    class JointForwardModel:
        def predict(self, experiment, corrections):
            correction = corrections.get("A", 0.0)
            pressure = experiment.pressure_bar
            temperature = experiment.temperature_c
            return EquilibriumResult(
                True,
                ("A",),
                {"A": 1.0},
                {
                    "A": {
                        "x": 0.5
                        + correction / 1000.0
                        + (pressure - 3000.0) / 10000.0
                        + (temperature - 600.0) / 1000.0
                    }
                },
                {},
            )

    experiment = Experiment(
        "joint",
        2500.0,
        550.0,
        {"A": 1.0},
        ("A",),
        phase_compositions={"A": {"x": 0.5}},
        phase_composition_sigma={"A": {"x": 0.01}},
    )
    problem = JointInversionProblem(
        JointForwardModel(),
        experiment,
        [
            ThermodynamicParameter(
                "dH_A", "A", prior_mean_j_mol=0.0, prior_sigma_j_mol=5000.0
            )
        ],
        Parameter("pressure", 3000.0, 500.0, 500.0, 8000.0, "bar"),
        Parameter("temperature", 600.0, 50.0, 300.0, 900.0, "degC"),
    )
    result = problem.optimize(starts=2, iterations=30, seed=2)

    assert abs(result.corrections_j_mol["pressure_bar"] - 3000.0) <= 250.0
    assert abs(result.corrections_j_mol["temperature_c"] - 600.0) <= 30.0


def test_joint_inversion_works_with_existing_mcmc_sampler():
    class LinearJointForwardModel:
        def predict(self, experiment, corrections):
            return EquilibriumResult(
                True,
                ("A",),
                {"A": 1.0},
                {"A": {"x": 0.5 + corrections.get("A", 0.0) / 1000.0}},
                {},
            )

    experiment = Experiment(
        "joint-mcmc",
        3000.0,
        600.0,
        {"A": 1.0},
        ("A",),
        phase_compositions={"A": {"x": 0.5}},
        phase_composition_sigma={"A": {"x": 0.01}},
    )
    problem = JointInversionProblem(
        LinearJointForwardModel(),
        experiment,
        [ThermodynamicParameter("dH_A", "A", prior_sigma_j_mol=10.0)],
        Parameter("pressure", 3000.0, 500.0, 500.0, 8000.0, "bar"),
        Parameter("temperature", 600.0, 50.0, 300.0, 900.0, "degC"),
    )
    result = mcmc_sample(
        problem,
        [0.0, 3000.0, 600.0],
        n_samples=30,
        burn_in=5,
        proposal_sigma_j_mol=[1.0, 50.0, 5.0],
        seed=1,
    )

    assert result.samples.shape == (25, 3)
    assert result.parameter_names == ("dH_A", "pressure", "temperature")
    assert result.parameter_units == ("J/mol", "bar", "degC")


def test_joint_inversion_correlated_mixed_unit_prior_is_applied():
    class ConstantJointForwardModel:
        def predict(self, experiment, corrections):
            return EquilibriumResult(True, ("A",), {"A": 1.0}, {}, {})

    experiment = Experiment(
        "joint-correlated-prior",
        3000.0,
        600.0,
        {"A": 1.0},
        ("A",),
        phase_fractions={"A": 1.0},
        phase_fraction_sigma={"A": 1.0},
    )
    covariance = np.asarray(
        [[10000.0, 4000.0, 100.0], [4000.0, 250000.0, 500.0], [100.0, 500.0, 2500.0]]
    )
    problem = JointInversionProblem(
        ConstantJointForwardModel(),
        experiment,
        [
            ThermodynamicParameter(
                "dH_A", "A", prior_mean_j_mol=10.0, prior_sigma_j_mol=100.0
            )
        ],
        Parameter("pressure", 3000.0, 500.0, 500.0, 8000.0, "bar"),
        Parameter("temperature", 600.0, 50.0, 300.0, 900.0, "degC"),
        prior_covariance=covariance,
    )
    point = np.asarray([30.0, 3100.0, 620.0])
    deviation = point - np.asarray([10.0, 3000.0, 600.0])
    expected_prior = 0.5 * deviation @ np.linalg.solve(covariance, deviation)

    assert np.allclose(problem.prior_covariance, covariance)
    assert np.isclose(problem.objective_components(point)[1], expected_prior)
    approximation = laplace_covariance(problem, point, step_fraction=1.0e-3)

    assert approximation.parameter_names == ("dH_A", "pressure", "temperature")
    assert approximation.parameter_units == ("J/mol", "bar", "degC")


def test_joint_inversion_shares_thermo_correction_across_multiple_experiments():
    class MultiExperimentForwardModel:
        def predict(self, experiment, corrections):
            correction = corrections.get("A", 0.0)
            value = (
                0.5
                + correction / 1000.0
                + (experiment.pressure_bar - 3000.0) / 10000.0
                + (experiment.temperature_c - 600.0) / 1000.0
            )
            return EquilibriumResult(True, ("A",), {"A": 1.0}, {"A": {"x": value}}, {})

    experiments = [
        Experiment(
            "joint-a",
            3000.0,
            600.0,
            {"A": 1.0},
            ("A",),
            phase_compositions={"A": {"x": 0.5}},
            phase_composition_sigma={"A": {"x": 0.01}},
        ),
        Experiment(
            "joint-b",
            4000.0,
            650.0,
            {"A": 1.0},
            ("A",),
            phase_compositions={"A": {"x": 0.65}},
            phase_composition_sigma={"A": {"x": 0.01}},
        ),
    ]
    problem = JointInversionProblem(
        MultiExperimentForwardModel(),
        experiments,
        [ThermodynamicParameter("dH_A", "A", prior_sigma_j_mol=100.0)],
        [
            Parameter("pressure_a", 3000.0, 500.0, 500.0, 8000.0, "bar"),
            Parameter("pressure_b", 4000.0, 500.0, 500.0, 8000.0, "bar"),
        ],
        [
            Parameter("temperature_a", 600.0, 50.0, 300.0, 900.0, "degC"),
            Parameter("temperature_b", 650.0, 50.0, 300.0, 900.0, "degC"),
        ],
    )
    point = [0.0, 3000.0, 600.0, 4000.0, 650.0]

    data, prior = problem.objective_components(point)

    assert data == 0.0
    assert prior == 0.0
    assert len(problem.parameters) == 5


def test_thermodynamic_parameter_exposes_native_unit_aliases():
    parameter = ThermodynamicParameter(
        "dH_A",
        "A",
        prior_mean_j_mol=25.0,
        prior_sigma_j_mol=50.0,
        lower_bound_j_mol=-100.0,
        upper_bound_j_mol=100.0,
    )

    assert parameter.unit == "J/mol"
    assert isinstance(parameter, Parameter)
    assert parameter.prior_mean_j_mol == parameter.prior_mean
    assert parameter.prior_sigma_j_mol == parameter.prior_sigma
    assert parameter.prior_mean == 25.0
    assert parameter.prior_sigma == 50.0
    assert parameter.bounds() == (-100.0, 100.0)


def test_database_overlay_changes_selected_enthalpy_only_in_copy():
    original = json.loads(DB_PATH.read_text(encoding="utf-8"))
    modified_path = DatabaseOverlay(DB_PATH).create({"Znc": 1250.0})
    try:
        modified = json.loads(modified_path.read_text(encoding="utf-8"))
        assert (
            modified["Species"]["Znc"]["StandardThermoModel"]["Constant"]["H0"]
            == original["Species"]["Znc"]["StandardThermoModel"]["Constant"]["H0"]
            + 1250.0
        )
        assert (
            modified["Species"]["Znc"]["StandardThermoModel"]["Constant"]["G0"]
            == original["Species"]["Znc"]["StandardThermoModel"]["Constant"]["G0"]
        )
        assert json.loads(DB_PATH.read_text(encoding="utf-8")) == original
    finally:
        modified_path.unlink(missing_ok=True)


def test_database_overlay_supports_compound_field_keys():
    original = json.loads(DB_PATH.read_text(encoding="utf-8"))
    modified_path = DatabaseOverlay(DB_PATH).create({"Znc:V0": 0.5, "Znc": 100.0})
    try:
        modified = json.loads(modified_path.read_text(encoding="utf-8"))
        constant = modified["Species"]["Znc"]["StandardThermoModel"]["Constant"]
        original_constant = original["Species"]["Znc"]["StandardThermoModel"][
            "Constant"
        ]
        assert constant["V0"] == original_constant["V0"] + 0.5
        assert constant["H0"] == original_constant["H0"] + 100.0
        assert constant["G0"] == original_constant["G0"]
    finally:
        modified_path.unlink(missing_ok=True)


def test_parameter_key_is_bare_phase_for_default_field_and_compound_otherwise():
    default_field = ThermodynamicParameter("dH_cc", "cc")
    assert default_field.key == "cc"
    volume_field = ThermodynamicParameter("dV_cc", "cc", field="V0")
    assert volume_field.key == "cc:V0"


def test_multi_field_parameters_on_same_phase_are_independently_invertible():
    class FakeForwardModel:
        def predict(self, experiment, corrections):
            enthalpy_error = corrections.get("cc", 0.0) - 500.0
            volume_error = corrections.get("cc:V0", 0.0) - 2.0
            return EquilibriumResult(
                True,
                ("cc",),
                {"cc": 1.0},
                {"cc": {"cc": 0.9 + enthalpy_error / 5000.0 + volume_error / 500.0}},
                {},
            )

    experiment = Experiment(
        id="multi-field",
        pressure_bar=2000.0,
        temperature_c=500.0,
        bulk_composition={"Ca": 1.0},
        observed_phases=("cc",),
        phase_compositions={"cc": {"cc": 0.9}},
        phase_composition_sigma={"cc": {"cc": 0.001}},
    )
    parameters = [
        ThermodynamicParameter(
            "dH_cc",
            "cc",
            field="H0",
            prior_sigma_j_mol=5000.0,
            lower_bound_j_mol=-4000.0,
            upper_bound_j_mol=4000.0,
        ),
        ThermodynamicParameter(
            "dV_cc",
            "cc",
            field="V0",
            prior_sigma_j_mol=5000.0,
            lower_bound_j_mol=-4000.0,
            upper_bound_j_mol=4000.0,
        ),
    ]
    problem = InversionProblem(FakeForwardModel(), [experiment], parameters)
    corrections = {
        parameter.key: parameter.prior_mean_j_mol for parameter in parameters
    }
    assert set(corrections) == {"cc", "cc:V0"}
    # Uniqueness check must accept two parameters on the same phase with
    # different fields (would previously have raised on duplicate phase).
    assert len(problem.parameters) == 2


def test_one_parameter_map_recovers_synthetic_correction():
    class FakeForwardModel:
        def predict(self, experiment, corrections):
            error = corrections.get("Znc", 0.0) - 1500.0
            valid = abs(error) < 2500.0
            return EquilibriumResult(
                valid,
                ("Znc",),
                {"Znc": max(0.0, 1.0 + error / 1000.0)},
                {},
                {"synthetic_error": error},
            )

    experiment = Experiment(
        id="synthetic",
        pressure_bar=2000.0,
        temperature_c=500.0,
        bulk_composition={"Zn": 1.0},
        observed_phases=("Znc",),
        phase_fractions={"Znc": 1.0},
        phase_fraction_sigma={"Znc": 0.01},
    )
    parameter = ThermodynamicParameter(
        name="delta_H_Znc",
        phase="Znc",
        prior_mean_j_mol=0.0,
        prior_sigma_j_mol=5000.0,
        lower_bound_j_mol=-4000.0,
        upper_bound_j_mol=4000.0,
    )
    result = InversionProblem(
        FakeForwardModel(), [experiment], [parameter]
    ).optimize_one_parameter(starts=3, iterations=30)
    assert abs(result.corrections_j_mol["Znc"] - 1500.0) < 100.0


def test_phase_composition_likelihood_penalizes_component_residual():
    experiment = Experiment(
        id="composition",
        pressure_bar=2000.0,
        temperature_c=500.0,
        bulk_composition={"Zn": 1.0},
        observed_phases=("Znc",),
        phase_compositions={"Znc": {"Zn": 0.8}},
        phase_composition_sigma={"Znc": {"Zn": 0.1}},
    )
    prediction = EquilibriumResult(
        True,
        ("Znc",),
        {"Znc": 1.0},
        {"Znc": {"Zn": 0.6}},
        {},
    )
    assert abs(experiment_misfit(experiment, prediction) - 4.0) < 1.0e-12


def test_experiment_misfit_supports_alternative_misfit_forms():
    experiment = Experiment(
        id="composition-forms",
        pressure_bar=2000.0,
        temperature_c=500.0,
        bulk_composition={"Zn": 1.0},
        observed_phases=("Znc",),
        phase_compositions={"Znc": {"Zn": 0.8}},
        phase_composition_sigma={"Znc": {"Zn": 0.1}},
    )
    prediction = EquilibriumResult(
        True, ("Znc",), {"Znc": 1.0}, {"Znc": {"Zn": 0.6}}, {}
    )

    raw_lsq = experiment_misfit(experiment, prediction, misfit_form="raw_lsq")
    assert abs(raw_lsq - (0.6 - 0.8) ** 2) < 1.0e-12

    chi_square = experiment_misfit(experiment, prediction, misfit_form="chi_square")
    assert abs(chi_square - (0.6 - 0.8) ** 2 / 0.6) < 1.0e-12

    try:
        experiment_misfit(experiment, prediction, misfit_form="bogus")
    except ValueError:
        pass
    else:
        raise AssertionError("Expected ValueError for an invalid misfit_form")


def test_predictive_log_likelihood_is_normalized_gaussian_and_phase_score():
    experiment = Experiment(
        "proper-score",
        2000.0,
        500.0,
        {"A": 1.0},
        ("A",),
        absent_phases=("B",),
        phase_fractions={"A": 0.4},
        phase_fraction_sigma={"A": 0.2},
    )
    matching = EquilibriumResult(True, ("A",), {"A": 1.0}, {}, {})
    mismatching = EquilibriumResult(True, ("B",), {"A": 1.0}, {}, {})
    expected_continuous = -0.5 * (
        ((1.0 - 0.4) / 0.2) ** 2 + np.log(2.0 * np.pi * 0.2**2)
    )
    expected_phase = np.log(1.0 - 0.02) * 2.0

    match_score = experiment_log_likelihood(
        experiment, matching, phase_error_probability=0.02
    )
    mismatch_score = experiment_log_likelihood(
        experiment, mismatching, phase_error_probability=0.02
    )

    assert np.isclose(match_score, expected_continuous + expected_phase)
    assert np.isclose(
        mismatch_score,
        expected_continuous + 2.0 * np.log(0.02),
    )


def test_experiment_misfit_supports_logit_likelihood_for_bounded_observations():
    experiment = Experiment(
        id="logit-composition",
        pressure_bar=2000.0,
        temperature_c=500.0,
        bulk_composition={"Zn": 1.0},
        observed_phases=("Znc",),
        phase_compositions={"Znc": {"Zn": 0.8}},
        phase_composition_sigma={"Znc": {"Zn": 0.1}},
    )
    prediction = EquilibriumResult(
        True, ("Znc",), {"Znc": 1.0}, {"Znc": {"Zn": 0.6}}, {}
    )
    logit = lambda value: np.log(value / (1.0 - value))
    expected_sigma = 0.1 / (0.8 * 0.2)
    expected = ((logit(0.6) - logit(0.8)) / expected_sigma) ** 2

    assert np.isclose(
        experiment_misfit(experiment, prediction, misfit_form="logit_lsq"), expected
    )

    invalid = Experiment(
        id="invalid-logit",
        pressure_bar=2000.0,
        temperature_c=500.0,
        bulk_composition={"Zn": 1.0},
        observed_phases=("Znc",),
        phase_compositions={"Znc": {"Zn": 0.0}},
        phase_composition_sigma={"Znc": {"Zn": 0.1}},
    )
    try:
        experiment_misfit(invalid, prediction, misfit_form="logit_lsq")
    except ValueError:
        pass
    else:
        raise AssertionError("Expected ValueError for a boundary logit observation")


def test_experiment_misfit_supports_covariance_aware_alr_likelihood():
    experiment = Experiment(
        id="alr-composition",
        pressure_bar=2000.0,
        temperature_c=500.0,
        bulk_composition={"Zn": 1.0},
        observed_phases=("Znc",),
        phase_compositions={"Znc": {"Zn": 0.8, "Mg": 0.4}},
        observation_covariance=[[0.01, 0.005], [0.005, 0.01]],
        observation_covariance_keys=("composition:Znc:Mg", "composition:Znc:Zn"),
    )
    prediction = EquilibriumResult(
        True, ("Znc",), {"Znc": 1.0}, {"Znc": {"Zn": 0.7, "Mg": 0.3}}, {}
    )
    observed = np.asarray([0.8, 0.4])
    calculated = np.asarray([0.7, 0.3])
    observed /= observed.sum()
    calculated /= calculated.sum()
    residual = np.asarray(
        [np.log(calculated[0] / calculated[1]) - np.log(observed[0] / observed[1])]
    )
    transform = np.asarray([[1.0 / observed[0], -1.0 / observed[1]]])
    covariance = np.asarray([[0.01, 0.005], [0.005, 0.01]])
    transformed_covariance = transform @ covariance @ transform.T

    assert np.isclose(
        experiment_misfit(experiment, prediction, misfit_form="alr_lsq"),
        residual @ np.linalg.solve(transformed_covariance, residual),
    )


def test_alr_likelihood_rejects_nonpositive_or_mixed_observations():
    experiment = Experiment(
        id="invalid-alr",
        pressure_bar=2000.0,
        temperature_c=500.0,
        bulk_composition={"Zn": 1.0},
        observed_phases=("Znc",),
        phase_fractions={"Znc": 0.5},
        phase_fraction_sigma={"Znc": 0.1},
    )
    prediction = EquilibriumResult(True, ("Znc",), {"Znc": 0.5}, {}, {})
    try:
        experiment_misfit(experiment, prediction, misfit_form="alr_lsq")
    except ValueError:
        pass
    else:
        raise AssertionError("Expected ValueError for phase fractions with alr_lsq")


def test_experiment_misfit_supports_logistic_normal_nll():
    experiment = Experiment(
        id="logistic-normal-composition",
        pressure_bar=2000.0,
        temperature_c=500.0,
        bulk_composition={"Zn": 1.0},
        observed_phases=("Znc",),
        phase_compositions={"Znc": {"Zn": 0.8, "Mg": 0.4}},
        observation_covariance=[[0.01, 0.005], [0.005, 0.01]],
        observation_covariance_keys=("composition:Znc:Zn", "composition:Znc:Mg"),
    )
    prediction = EquilibriumResult(
        True, ("Znc",), {"Znc": 1.0}, {"Znc": {"Zn": 0.7, "Mg": 0.3}}, {}
    )
    value = experiment_misfit(experiment, prediction, misfit_form="logistic_normal_nll")

    assert np.isfinite(value)


def test_experiment_misfit_supports_dirichlet_composition_likelihood():
    experiment = Experiment(
        id="dirichlet-composition",
        pressure_bar=2000.0,
        temperature_c=500.0,
        bulk_composition={"Zn": 1.0},
        observed_phases=("Znc",),
        phase_compositions={"Znc": {"Zn": 0.8, "Mg": 0.4}},
        phase_composition_sigma={"Znc": {"Zn": 0.1, "Mg": 0.1}},
    )
    prediction = EquilibriumResult(
        True, ("Znc",), {"Znc": 1.0}, {"Znc": {"Zn": 0.7, "Mg": 0.3}}, {}
    )

    value = experiment_misfit(experiment, prediction, misfit_form="dirichlet_nll")

    assert np.isfinite(value)
    assert value > 0.0

    covariance_experiment = replace(
        experiment,
        observation_covariance=[[0.01, 0.0], [0.0, 0.01]],
        observation_covariance_keys=("composition:Znc:Zn", "composition:Znc:Mg"),
    )
    try:
        experiment_misfit(
            covariance_experiment, prediction, misfit_form="dirichlet_nll"
        )
    except ValueError:
        pass
    else:
        raise AssertionError("Expected ValueError for covariance with dirichlet_nll")


def test_dirichlet_likelihood_accepts_explicit_phase_concentration():
    experiment = Experiment(
        id="explicit-dirichlet-concentration",
        pressure_bar=2000.0,
        temperature_c=500.0,
        bulk_composition={"Zn": 1.0},
        observed_phases=("Znc",),
        phase_compositions={"Znc": {"Zn": 0.8, "Mg": 0.4}},
        phase_composition_concentration={"Znc": 20.0},
    )
    prediction = EquilibriumResult(
        True, ("Znc",), {"Znc": 1.0}, {"Znc": {"Zn": 0.7, "Mg": 0.3}}, {}
    )
    observed = np.asarray([2.0 / 3.0, 1.0 / 3.0])
    alpha = np.asarray([14.0, 6.0])
    expected = lgamma(20.0) - lgamma(14.0) - lgamma(6.0)
    expected -= float(np.sum((alpha - 1.0) * np.log(observed)))

    actual = experiment_misfit(experiment, prediction, misfit_form="dirichlet_nll")

    assert np.isclose(actual, expected)


def test_dirichlet_concentration_can_be_inferred_without_marginal_sigmas():
    class CompositionForwardModel:
        def predict(self, experiment, corrections):
            return EquilibriumResult(
                True,
                ("Znc",),
                {"Znc": 1.0},
                {"Znc": {"Zn": 0.7, "Mg": 0.3}},
                {},
            )

    experiment = Experiment(
        id="inferred-dirichlet-concentration",
        pressure_bar=2000.0,
        temperature_c=500.0,
        bulk_composition={"Zn": 1.0},
        observed_phases=("Znc",),
        phase_compositions={"Znc": {"Zn": 0.8, "Mg": 0.2}},
    )
    problem = InversionProblem(
        CompositionForwardModel(),
        [experiment],
        [
            DirichletConcentrationParameter(
                "Znc",
                prior_mean=20.0,
                prior_sigma=10.0,
                lower_bound=0.1,
                upper_bound=100.0,
            )
        ],
        misfit_form="dirichlet_nll",
    )

    result = mcmc_sample(
        problem,
        [20.0],
        n_samples=30,
        burn_in=5,
        proposal_sigma=[2.0],
        seed=6,
    )
    map_result = problem.optimize_multi_parameter(starts=2, iterations=4, seed=6)

    assert result.parameter_names == ("dirichlet_kappa[Znc]",)
    assert result.samples.shape == (25, 1)
    assert np.all(result.samples[:, 0] > 0.0)
    assert map_result.corrections_j_mol == {}
    assert "dirichlet_kappa[Znc]" in map_result.parameter_values


def test_dirichlet_likelihood_handles_zero_components_as_censored_mass():
    experiment = Experiment(
        id="censored-dirichlet-composition",
        pressure_bar=2000.0,
        temperature_c=500.0,
        bulk_composition={"Zn": 1.0},
        observed_phases=("Znc",),
        phase_compositions={"Znc": {"Zn": 0.7, "Mg": 0.2, "Mn": 0.0}},
        phase_composition_concentration={"Znc": 20.0},
        phase_composition_detection_limit={"Znc": 0.1},
    )
    prediction = EquilibriumResult(
        True,
        ("Znc",),
        {"Znc": 1.0},
        {"Znc": {"Zn": 0.7, "Mg": 0.2, "Mn": 0.1}},
        {},
    )
    alpha_active = np.asarray([14.0, 4.0])
    active_observed = np.asarray([7.0 / 9.0, 2.0 / 9.0])
    conditional_nll = lgamma(18.0) - lgamma(14.0) - lgamma(4.0)
    conditional_nll -= float(np.sum((alpha_active - 1.0) * np.log(active_observed)))
    beta_cdf = 1.0 - 0.9**18 * (1.0 + 18.0 * 0.1)
    expected = conditional_nll - np.log(beta_cdf)

    actual = experiment_misfit(experiment, prediction, misfit_form="dirichlet_nll")

    assert np.isclose(actual, expected)
    try:
        experiment_misfit(
            replace(experiment, phase_composition_detection_limit=None),
            prediction,
            misfit_form="dirichlet_nll",
        )
    except ValueError as error:
        assert "detection limit" in str(error)
    else:
        raise AssertionError("Expected zero components without a limit to be rejected")


def test_experiment_misfit_supports_correlated_observation_covariance():
    experiment = Experiment(
        id="correlated-composition",
        pressure_bar=2000.0,
        temperature_c=500.0,
        bulk_composition={"Zn": 1.0},
        observed_phases=("Znc",),
        phase_compositions={"Znc": {"Zn": 0.8, "Mg": 0.4}},
        observation_covariance=[[0.01, 0.005], [0.005, 0.01]],
        observation_covariance_keys=("composition:Znc:Zn", "composition:Znc:Mg"),
    )
    prediction = EquilibriumResult(
        True, ("Znc",), {"Znc": 1.0}, {"Znc": {"Zn": 0.6, "Mg": 0.3}}, {}
    )
    residual = np.asarray([-0.2, -0.1])
    covariance = np.asarray(experiment.observation_covariance)

    assert np.isclose(
        experiment_misfit(experiment, prediction),
        residual @ np.linalg.solve(covariance, residual),
    )

    try:
        experiment_misfit(experiment, prediction, misfit_form="raw_lsq")
    except ValueError:
        pass
    else:
        raise AssertionError("Expected ValueError for covariance with raw_lsq")


def test_experiment_misfit_combines_keyed_analytical_and_theoretical_covariance():
    experiment = Experiment(
        id="analytical-plus-theoretical",
        pressure_bar=2000.0,
        temperature_c=500.0,
        bulk_composition={"Zn": 1.0},
        observed_phases=("Znc",),
        phase_compositions={"Znc": {"Zn": 0.8, "Mg": 0.4}},
        observation_covariance=[[0.01, 0.002], [0.002, 0.02]],
        observation_covariance_keys=("composition:Znc:Zn", "composition:Znc:Mg"),
        theoretical_covariance=[[0.04, 0.003], [0.003, 0.03]],
        theoretical_covariance_keys=("composition:Znc:Mg", "composition:Znc:Zn"),
    )
    prediction = EquilibriumResult(
        True, ("Znc",), {"Znc": 1.0}, {"Znc": {"Zn": 0.6, "Mg": 0.3}}, {}
    )
    residual = np.asarray([-0.2, -0.1])
    analytical = np.asarray(experiment.observation_covariance)
    theoretical_in_observation_order = np.asarray([[0.03, 0.003], [0.003, 0.04]])

    assert np.isclose(
        experiment_misfit(experiment, prediction),
        residual
        @ np.linalg.solve(analytical + theoretical_in_observation_order, residual),
    )


def test_experiment_misfit_scales_theoretical_covariance():
    experiment = Experiment(
        id="scaled-theoretical",
        pressure_bar=2000.0,
        temperature_c=500.0,
        bulk_composition={"Zn": 1.0},
        observed_phases=("Znc",),
        phase_compositions={"Znc": {"Zn": 0.8}},
        phase_composition_sigma={"Znc": {"Zn": 0.1}},
        theoretical_covariance=[[0.01]],
        theoretical_covariance_keys=("composition:Znc:Zn",),
    )
    prediction = EquilibriumResult(
        True, ("Znc",), {"Znc": 1.0}, {"Znc": {"Zn": 0.6}}, {}
    )

    assert np.isclose(
        experiment_misfit(experiment, prediction, theoretical_covariance_scale=0.0),
        4.0,
    )
    assert np.isclose(
        experiment_misfit(experiment, prediction, theoretical_covariance_scale=1.0),
        2.0,
    )


def test_experiment_misfit_applies_known_observation_bias():
    experiment = Experiment(
        id="biased-observation",
        pressure_bar=2000.0,
        temperature_c=500.0,
        bulk_composition={"Zn": 1.0},
        observed_phases=("Znc",),
        phase_compositions={"Znc": {"Zn": 0.8}},
        phase_composition_sigma={"Znc": {"Zn": 0.1}},
        observation_bias={"composition:Znc:Zn": 0.1},
    )
    prediction = EquilibriumResult(
        True, ("Znc",), {"Znc": 1.0}, {"Znc": {"Zn": 0.7}}, {}
    )

    assert np.isclose(experiment_misfit(experiment, prediction), 0.0)

    try:
        experiment_misfit(
            replace(experiment, observation_bias={"unknown": 0.1}), prediction
        )
    except ValueError:
        pass
    else:
        raise AssertionError("Expected ValueError for an unknown observation bias key")


def test_theoretical_covariance_scale_profile_returns_prior_preferred_scale():
    class MatchingForwardModel:
        def predict(self, experiment, corrections):
            return EquilibriumResult(True, ("A",), {"A": 1.0}, {"A": {"x": 0.8}}, {})

    experiment = Experiment(
        id="profile-scale",
        pressure_bar=2000.0,
        temperature_c=500.0,
        bulk_composition={"A": 1.0},
        observed_phases=("A",),
        phase_compositions={"A": {"x": 0.8}},
        phase_composition_sigma={"A": {"x": 0.1}},
        theoretical_covariance=[[0.01]],
        theoretical_covariance_keys=("composition:A:x",),
    )
    parameter = ThermodynamicParameter("dH_A", "A", prior_sigma_j_mol=1.0)
    problem = InversionProblem(MatchingForwardModel(), [experiment], [parameter])

    result = profile_theoretical_covariance_scale(
        problem,
        scale_bounds=(0.0, 3.0),
        scale_prior_mean=1.0,
        scale_prior_sigma=0.2,
        scale_iterations=8,
        optimize_starts=1,
        optimize_iterations=5,
    )

    assert abs(result.scale - 1.0) < 0.1
    assert result.scale_prior_penalty < 0.2
    assert result.records


def test_inversion_uses_metadata_group_theoretical_covariance_scales():
    class MatchingForwardModel:
        def predict(self, experiment, corrections):
            return EquilibriumResult(True, ("A",), {"A": 1.0}, {"A": {"x": 0.6}}, {})

    experiments = [
        Experiment(
            id="study-a",
            pressure_bar=2000.0,
            temperature_c=500.0,
            bulk_composition={"A": 1.0},
            observed_phases=("A",),
            phase_compositions={"A": {"x": 0.8}},
            phase_composition_sigma={"A": {"x": 0.1}},
            theoretical_covariance=[[0.01]],
            theoretical_covariance_keys=("composition:A:x",),
            metadata={"study": "a"},
        ),
        Experiment(
            id="study-b",
            pressure_bar=2000.0,
            temperature_c=500.0,
            bulk_composition={"A": 1.0},
            observed_phases=("A",),
            phase_compositions={"A": {"x": 0.8}},
            phase_composition_sigma={"A": {"x": 0.1}},
            theoretical_covariance=[[0.01]],
            theoretical_covariance_keys=("composition:A:x",),
            metadata={"study": "b"},
        ),
    ]
    parameter = ThermodynamicParameter("dH_A", "A", prior_sigma_j_mol=1.0)
    problem = InversionProblem(
        MatchingForwardModel(),
        experiments,
        [parameter],
        theoretical_covariance_scale_metadata_key="study",
        theoretical_covariance_scales={"a": 0.0, "b": 2.0},
    )

    data, _ = problem.objective_components([0.0])

    assert np.isclose(data, 4.0 + 0.8)


def test_observation_bias_profile_returns_group_biases():
    class MatchingForwardModel:
        def predict(self, experiment, corrections):
            return EquilibriumResult(True, ("A",), {"A": 1.0}, {"A": {"x": 0.8}}, {})

    experiments = [
        Experiment(
            id=group,
            pressure_bar=2000.0,
            temperature_c=500.0,
            bulk_composition={"A": 1.0},
            observed_phases=("A",),
            phase_compositions={"A": {"x": 0.8}},
            phase_composition_sigma={"A": {"x": 0.1}},
            metadata={"laboratory": group},
        )
        for group in ("lab-a", "lab-b")
    ]
    problem = InversionProblem(
        MatchingForwardModel(),
        experiments,
        [ThermodynamicParameter("dH_A", "A", prior_sigma_j_mol=1.0)],
    )

    result = profile_observation_bias_by_group(
        problem,
        "laboratory",
        "composition:A:x",
        bias_iterations=2,
        optimize_starts=1,
        optimize_iterations=3,
    )

    assert set(result.biases) == {"lab-a", "lab-b"}
    assert result.metadata_key == "laboratory"
    assert result.observation_key == "composition:A:x"
    assert result.records


def test_observation_bias_profile_rejects_missing_observation_key():
    experiment = Experiment(
        id="missing-bias-observation",
        pressure_bar=2000.0,
        temperature_c=500.0,
        bulk_composition={"A": 1.0},
        observed_phases=("A",),
        phase_compositions={"A": {"x": 0.8}},
        phase_composition_sigma={"A": {"x": 0.1}},
        metadata={"laboratory": "lab-a"},
    )
    problem = InversionProblem(
        object(), [experiment], [ThermodynamicParameter("dH_A", "A")]
    )
    try:
        profile_observation_bias_by_group(problem, "laboratory", "fraction:A")
    except ValueError:
        pass
    else:
        raise AssertionError("Expected ValueError for missing bias observation")


def test_profile_observation_biases_by_group_produces_per_key_map():
    class MatchingForwardModel:
        def predict(self, experiment, corrections):
            return EquilibriumResult(
                True, ("A",), {"A": 1.0}, {"A": {"x": 0.8, "y": 0.5}}, {}
            )

    experiments = [
        Experiment(
            id=group,
            pressure_bar=2000.0,
            temperature_c=500.0,
            bulk_composition={"A": 1.0},
            observed_phases=("A",),
            phase_compositions={"A": {"x": 0.8, "y": 0.5}},
            phase_composition_sigma={"A": {"x": 0.1, "y": 0.1}},
            metadata={"laboratory": group},
        )
        for group in ("lab-a", "lab-b")
    ]
    problem = InversionProblem(
        MatchingForwardModel(),
        experiments,
        [ThermodynamicParameter("dH_A", "A", prior_sigma_j_mol=1.0)],
    )

    result = profile_observation_biases_by_group(
        problem,
        "laboratory",
        ("composition:A:x", "composition:A:y"),
        bias_bounds=(-0.2, 0.2),
        bias_prior_mean=0.0,
        bias_prior_sigma=0.1,
        bias_iterations=2,
        optimize_starts=1,
        optimize_iterations=3,
    )

    assert set(result) == {"lab-a", "lab-b"}
    assert set(result["lab-a"]) == {"composition:A:x", "composition:A:y"}
    assert all(
        np.isfinite(value) for values in result.values() for value in values.values()
    )


def test_apply_group_observation_biases_returns_problem_with_biases_applied():
    experiment = Experiment(
        id="lab-a",
        pressure_bar=2000.0,
        temperature_c=500.0,
        bulk_composition={"A": 1.0},
        observed_phases=("A",),
        phase_compositions={"A": {"x": 0.8}},
        phase_composition_sigma={"A": {"x": 0.1}},
        metadata={"laboratory": "lab-a"},
    )
    problem = InversionProblem(
        object(),
        [experiment],
        [ThermodynamicParameter("dH_A", "A", prior_sigma_j_mol=1.0)],
    )

    applied = apply_group_observation_biases(
        problem,
        "laboratory",
        {"lab-a": {"composition:A:x": 0.1}},
    )

    assert applied.experiments[0].observation_bias == {"composition:A:x": 0.1}


def test_summarize_group_observation_bias_preserves_group_estimates_and_prior():
    result = ObservationBiasProfileResult(
        {"lab-a": 0.1, "lab-b": -0.2},
        2.5,
        {"A": 0.0},
        0.25,
        "laboratory",
        "composition:A:x",
        ({"group": "lab-a", "bias": 0.1, "objective": 2.5},),
        (-0.5, 0.5),
        0.0,
        0.2,
    )

    summary = summarize_group_observation_bias(result)
    assert summary.metadata_key == "laboratory"
    assert summary.group_biases == {"lab-a": 0.1, "lab-b": -0.2}
    assert summary.bias_prior_sigma == 0.2
    assert summary.bias_bounds == (-0.5, 0.5)


def test_summarize_grouped_observation_biases_returns_multi_key_report():
    profile_a = ObservationBiasProfileResult(
        {"lab-a": 0.1, "lab-b": -0.2},
        2.5,
        {"A": 0.0},
        0.25,
        "laboratory",
        "composition:A:x",
        ({"group": "lab-a", "bias": 0.1, "objective": 2.5},),
        (-0.5, 0.5),
        0.0,
        0.2,
    )
    profile_b = ObservationBiasProfileResult(
        {"lab-a": 0.05, "lab-b": -0.10},
        1.5,
        {"A": 0.0},
        0.10,
        "laboratory",
        "composition:A:y",
        ({"group": "lab-a", "bias": 0.05, "objective": 1.5},),
        (-0.5, 0.5),
        0.0,
        0.2,
    )

    report = summarize_grouped_observation_biases(
        {
            "composition:A:x": profile_a,
            "composition:A:y": profile_b,
        }
    )

    assert isinstance(report, GroupedObservationBiasReport)
    assert report.observation_keys == ("composition:A:x", "composition:A:y")
    assert report.group_biases["lab-a"]["composition:A:x"] == 0.1
    assert report.group_biases["lab-b"]["composition:A:y"] == -0.10
    assert report.bias_prior_sigma == 0.2


def test_summarize_grouped_bias_report_preserves_metadata():
    report = GroupedObservationBiasReport(
        "laboratory",
        ("composition:A:x", "composition:A:y"),
        {
            "lab-a": {"composition:A:x": 0.10, "composition:A:y": -0.04},
            "lab-b": {"composition:A:x": -0.01, "composition:A:y": 0.03},
        },
        0.0,
        0.2,
        (-0.5, 0.5),
        1.5,
        0.3,
    )

    summary = summarize_grouped_bias_report(report)
    assert isinstance(summary, GroupedObservationBiasSummary)
    assert summary.metadata_key == "laboratory"
    assert summary.observation_keys == ("composition:A:x", "composition:A:y")
    assert summary.group_biases["lab-a"]["composition:A:y"] == -0.04
    assert summary.bias_prior_sigma == 0.2


def test_profile_grouped_bias_summary_builds_directly_from_problem():
    class MatchingForwardModel:
        def predict(self, experiment, corrections):
            return EquilibriumResult(
                True, ("A",), {"A": 1.0}, {"A": {"x": 0.8, "y": 0.5}}, {}
            )

    experiments = [
        Experiment(
            id=group,
            pressure_bar=2000.0,
            temperature_c=500.0,
            bulk_composition={"A": 1.0},
            observed_phases=("A",),
            phase_compositions={"A": {"x": 0.8, "y": 0.5}},
            phase_composition_sigma={"A": {"x": 0.1, "y": 0.1}},
            metadata={"laboratory": group},
        )
        for group in ("lab-a", "lab-b")
    ]
    problem = InversionProblem(
        MatchingForwardModel(),
        experiments,
        [ThermodynamicParameter("dH_A", "A", prior_sigma_j_mol=1.0)],
    )

    summary = profile_grouped_bias_summary(
        problem,
        "laboratory",
        ("composition:A:x", "composition:A:y"),
        bias_bounds=(-0.2, 0.2),
        bias_prior_mean=0.0,
        bias_prior_sigma=0.1,
        bias_iterations=2,
        optimize_starts=1,
        optimize_iterations=3,
    )

    assert isinstance(summary, GroupedObservationBiasSummary)
    assert summary.metadata_key == "laboratory"
    assert set(summary.observation_keys) == {"composition:A:x", "composition:A:y"}
    assert set(summary.group_biases) == {"lab-a", "lab-b"}
    assert summary.bias_prior_sigma == 0.1


def test_group_theoretical_covariance_scale_profile_returns_each_group():
    class MatchingForwardModel:
        def predict(self, experiment, corrections):
            return EquilibriumResult(True, ("A",), {"A": 1.0}, {"A": {"x": 0.8}}, {})

    experiments = [
        Experiment(
            id=group,
            pressure_bar=2000.0,
            temperature_c=500.0,
            bulk_composition={"A": 1.0},
            observed_phases=("A",),
            phase_compositions={"A": {"x": 0.8}},
            phase_composition_sigma={"A": {"x": 0.1}},
            theoretical_covariance=[[0.01]],
            theoretical_covariance_keys=("composition:A:x",),
            metadata={"study": group},
        )
        for group in ("a", "b")
    ]
    problem = InversionProblem(
        MatchingForwardModel(),
        experiments,
        [ThermodynamicParameter("dH_A", "A", prior_sigma_j_mol=1.0)],
    )

    result = profile_theoretical_covariance_group_scales(
        problem,
        "study",
        scale_bounds=(0.0, 2.0),
        scale_prior_mean=1.0,
        scale_prior_sigma=0.2,
        scale_iterations=2,
        optimize_starts=1,
        optimize_iterations=3,
    )

    assert set(result.scales) == {"a", "b"}
    assert result.metadata_key == "study"
    assert result.records


def test_bracket_likelihood_scores_endpoint_phase_violations():
    lower = Experiment("lower", 2000.0, 500.0, {}, ("A",))
    upper = Experiment("upper", 2000.0, 600.0, {}, ("B",))
    bracket = EquilibriumBracket(
        "A-to-B",
        lower,
        upper,
        lower_required_phases=("A",),
        upper_required_phases=("B",),
        lower_forbidden_phases=("B",),
        upper_forbidden_phases=("A",),
    )
    valid_lower = EquilibriumResult(True, ("A",), {}, {}, {})
    valid_upper = EquilibriumResult(True, ("B",), {}, {}, {})
    invalid_upper = EquilibriumResult(True, ("A",), {}, {}, {})
    assert bracket_misfit(bracket, valid_lower, valid_upper) == 0.0
    assert bracket_misfit(bracket, valid_lower, invalid_upper) == 200.0


def test_condition_uncertainty_uses_deterministic_quadrature_nodes():
    experiment = Experiment(
        "uncertain",
        2000.0,
        500.0,
        {},
        ("Znc",),
        pressure_sigma_bar=50.0,
        temperature_sigma_c=10.0,
    )

    class FakeForwardModel:
        def predict(self, experiment, corrections):
            return EquilibriumResult(True, ("Znc",), {"Znc": 1.0}, {}, {})

    problem = InversionProblem(
        FakeForwardModel(),
        [experiment],
        [ThermodynamicParameter("dH_Znc", "Znc")],
        marginalize_condition_uncertainty=True,
        condition_quadrature_order=3,
    )
    nodes = problem._condition_nodes(experiment)
    assert len(nodes) == 9
    assert abs(sum(weight for _, weight in nodes) - 1.0) < 1.0e-12
    assert any(
        node.pressure_bar != experiment.pressure_bar
        or node.temperature_c != experiment.temperature_c
        for node, _ in nodes
    )


def test_sensitivity_identifies_full_rank_parameters():
    result = identify_sensitivity([[1.0, 0.0], [0.0, 2.0]], ["dH_A", "dH_B"])
    assert result.rank == 2
    assert result.condition_number == 2.0
    assert result.weak_directions == ()


def test_sensitivity_reports_unidentifiable_parameter_direction():
    result = identify_sensitivity([[1.0, 1.0], [2.0, 2.0]], ["dH_A", "dH_B"])
    assert result.rank == 1
    assert len(result.weak_directions) == 1
    direction = result.weak_directions[0]
    assert abs(abs(direction["dH_A"]) - abs(direction["dH_B"])) < 1.0e-12


def test_prior_scaled_information_preserves_correlated_prior_geometry():
    jacobian = np.asarray([[1.0, 0.0]])
    parameters = [
        ThermodynamicParameter("dH_A", "A", prior_sigma_j_mol=1.0),
        ThermodynamicParameter("dH_B", "B", prior_sigma_j_mol=1.0),
    ]
    covariance = np.asarray([[4.0, 3.0], [3.0, 9.0]])

    information = prior_scaled_information(jacobian, parameters, covariance)

    assert np.allclose(information, [[4.0, 0.0], [0.0, 0.0]])


def test_sensitivity_report_writes_machine_readable_diagnostics(tmp_path):
    parameters = [
        ThermodynamicParameter("dH_A", "A", prior_sigma_j_mol=2.0),
        ThermodynamicParameter("dH_B", "B", prior_sigma_j_mol=3.0),
    ]
    result = SensitivityResult(
        ("obs_a", "obs_b"),
        ("dH_A", "dH_B"),
        np.zeros(2),
        np.asarray([[1.0, 0.0], [0.0, 2.0]]),
        identify_sensitivity([[1.0, 0.0], [0.0, 2.0]], ["dH_A", "dH_B"]),
    )
    write_sensitivity_report(result, parameters, tmp_path)
    assert (tmp_path / "sensitivity_matrix.csv").is_file()
    assert (tmp_path / "identifiability.json").is_file()
    restored = read_sensitivity_report(tmp_path)
    assert restored.observable_names == result.observable_names
    assert np.allclose(restored.values, result.values)
    assert np.allclose(restored.jacobian, result.jacobian)
    correlated_information = prior_scaled_information(
        result.jacobian,
        parameters,
        np.asarray([[4.0, 3.0], [3.0, 9.0]]),
    )
    write_sensitivity_report(
        result,
        parameters,
        tmp_path / "correlated",
        prior_covariance=np.asarray([[4.0, 3.0], [3.0, 9.0]]),
    )
    assert np.allclose(
        np.loadtxt(
            tmp_path / "correlated" / "prior_scaled_information.csv", delimiter=","
        ),
        correlated_information,
    )
    report = json.loads(
        (tmp_path / "correlated" / "identifiability.json").read_text(encoding="utf-8")
    )
    assert report["prior_covariance_j_mol2"] == [[4.0, 3.0], [3.0, 9.0]]
    assert np.allclose(
        prior_scaled_information(result.jacobian, parameters), [[4.0, 0.0], [0.0, 36.0]]
    )


def test_multi_parameter_map_recovers_two_synthetic_corrections():
    class FakeForwardModel:
        def predict(self, experiment, corrections):
            value_a = corrections.get("A", 0.0)
            value_b = corrections.get("B", 0.0)
            return EquilibriumResult(
                True,
                ("A", "B"),
                {"A": 0.5, "B": 0.5},
                {},
                {"synthetic": (value_a, value_b)},
            )

    experiment = Experiment(
        id="synthetic-two",
        pressure_bar=2000.0,
        temperature_c=500.0,
        bulk_composition={"A": 1.0},
        observed_phases=("A", "B"),
        phase_compositions={"A": {"A": 0.65}, "B": {"B": 0.65}},
        phase_composition_sigma={"A": {"A": 0.001}, "B": {"B": 0.001}},
    )

    class CompositionForwardModel(FakeForwardModel):
        def predict(self, experiment, corrections):
            value_a = corrections.get("A", 0.0)
            value_b = corrections.get("B", 0.0)
            return EquilibriumResult(
                True,
                ("A", "B"),
                {"A": 0.5, "B": 0.5},
                {
                    "A": {"A": 0.65 + (value_a - 1200.0) / 10000.0},
                    "B": {"B": 0.65 + (value_b + 800.0) / 10000.0},
                },
                {},
            )

    parameters = [
        ThermodynamicParameter(
            "dH_A",
            "A",
            prior_sigma_j_mol=5000.0,
            lower_bound_j_mol=-4000.0,
            upper_bound_j_mol=4000.0,
        ),
        ThermodynamicParameter(
            "dH_B",
            "B",
            prior_sigma_j_mol=5000.0,
            lower_bound_j_mol=-4000.0,
            upper_bound_j_mol=4000.0,
        ),
    ]
    result = InversionProblem(
        CompositionForwardModel(), [experiment], parameters
    ).optimize_multi_parameter(starts=4, iterations=80)
    assert abs(result.corrections_j_mol["A"] - 1200.0) < 100.0
    assert abs(result.corrections_j_mol["B"] + 800.0) < 100.0


def test_hot_start_multi_parameter_map_converges_near_central_point():
    class CompositionForwardModel:
        def predict(self, experiment, corrections):
            value_a = corrections.get("A", 0.0)
            value_b = corrections.get("B", 0.0)
            return EquilibriumResult(
                True,
                ("A", "B"),
                {"A": 0.5, "B": 0.5},
                {
                    "A": {"A": 0.65 + (value_a - 1200.0) / 10000.0},
                    "B": {"B": 0.65 + (value_b + 800.0) / 10000.0},
                },
                {},
            )

    experiment = Experiment(
        id="hot-start-two",
        pressure_bar=2000.0,
        temperature_c=500.0,
        bulk_composition={"A": 1.0},
        observed_phases=("A", "B"),
        phase_compositions={"A": {"A": 0.65}, "B": {"B": 0.65}},
        phase_composition_sigma={"A": {"A": 0.001}, "B": {"B": 0.001}},
    )
    parameters = [
        ThermodynamicParameter(
            "dH_A",
            "A",
            prior_sigma_j_mol=5000.0,
            lower_bound_j_mol=-4000.0,
            upper_bound_j_mol=4000.0,
        ),
        ThermodynamicParameter(
            "dH_B",
            "B",
            prior_sigma_j_mol=5000.0,
            lower_bound_j_mol=-4000.0,
            upper_bound_j_mol=4000.0,
        ),
    ]
    problem = InversionProblem(CompositionForwardModel(), [experiment], parameters)
    result = problem.optimize_multi_parameter(
        starts=4,
        iterations=80,
        central_point=[1150.0, -750.0],
        hot_start_spread_fraction=0.05,
    )
    assert abs(result.corrections_j_mol["A"] - 1200.0) < 100.0
    assert abs(result.corrections_j_mol["B"] + 800.0) < 100.0

    try:
        problem.optimize_multi_parameter(
            central_point=[0.0, 0.0], hot_start_spread_fraction=0.0
        )
    except ValueError:
        pass
    else:
        raise AssertionError("Expected ValueError for a non-positive spread fraction")

    try:
        problem.optimize_multi_parameter(central_point=[0.0])
    except ValueError:
        pass
    else:
        raise AssertionError(
            "Expected ValueError for a mismatched central_point length"
        )


def test_bayes_objective_criterion_adds_fixed_complexity_penalty():
    class FakeForwardModel:
        def predict(self, experiment, corrections):
            return EquilibriumResult(True, ("A",), {"A": 1.0}, {}, {})

    experiments = [
        Experiment(
            id=f"exp-{i}",
            pressure_bar=2000.0,
            temperature_c=500.0,
            bulk_composition={"A": 1.0},
            observed_phases=("A",),
            phase_fractions={"A": 1.0},
            phase_fraction_sigma={"A": 0.1},
        )
        for i in range(3)
    ]
    parameters = [ThermodynamicParameter("dH_A", "A", prior_sigma_j_mol=1000.0)]

    map_problem = InversionProblem(FakeForwardModel(), experiments, parameters)
    bayes_problem = InversionProblem(
        FakeForwardModel(), experiments, parameters, objective_criterion="bayes"
    )
    point = [50.0]
    expected_penalty = 0.5 * len(parameters) * np.log(len(experiments))
    assert abs(bayes_problem.bayes_penalty() - expected_penalty) < 1.0e-12
    assert (
        abs(
            bayes_problem.objective(point)
            - map_problem.objective(point)
            - expected_penalty
        )
        < 1.0e-9
    )

    try:
        InversionProblem(
            FakeForwardModel(), experiments, parameters, objective_criterion="bogus"
        )
    except ValueError:
        pass
    else:
        raise AssertionError("Expected ValueError for an invalid objective_criterion")


def test_correlated_prior_covariance_changes_joint_prior_penalty():
    class FakeForwardModel:
        def predict(self, experiment, corrections):
            return EquilibriumResult(True, ("A", "B"), {"A": 0.5, "B": 0.5}, {}, {})

    experiment = Experiment(
        "correlated-prior",
        2000.0,
        500.0,
        {"A": 1.0},
        ("A", "B"),
    )
    parameters = [
        ThermodynamicParameter("dH_A", "A", prior_sigma_j_mol=2.0),
        ThermodynamicParameter("dH_B", "B", prior_sigma_j_mol=2.0),
    ]
    covariance = np.asarray([[4.0, 3.0], [3.0, 4.0]])
    problem = InversionProblem(
        FakeForwardModel(),
        [experiment],
        parameters,
        prior_covariance_j_mol2=covariance,
    )
    _, prior = problem.objective_components([1.0, 1.0])
    expected = (
        0.5
        * np.asarray([1.0, 1.0])
        @ np.linalg.inv(covariance)
        @ np.asarray([1.0, 1.0])
    )
    assert abs(prior - expected) < 1.0e-12

    try:
        InversionProblem(
            FakeForwardModel(),
            [experiment],
            parameters,
            prior_covariance_j_mol2=[[1.0, 2.0], [2.0, 1.0]],
        )
    except ValueError:
        pass
    else:
        raise AssertionError(
            "Expected ValueError for a non-positive-definite prior covariance"
        )


def test_laplace_covariance_matches_analytic_quadratic_objective():
    # Constructed so experiment_misfit == ((corrections - target) / sigma) ** 2
    # for each independent phase, giving an analytically known diagonal Hessian
    # of 2 / sigma_data**2 (plus the Gaussian prior's own precision).
    sigma_data = 200.0
    target_a, target_b = 900.0, -400.0

    class QuadraticForwardModel:
        def predict(self, experiment, corrections):
            value_a = corrections.get("A", 0.0)
            value_b = corrections.get("B", 0.0)
            return EquilibriumResult(
                True,
                ("A", "B"),
                {
                    "A": 1.0 + (value_a - target_a) / sigma_data,
                    "B": 1.0 + (value_b - target_b) / sigma_data,
                },
                {},
                {},
            )

    experiment = Experiment(
        id="quadratic",
        pressure_bar=2000.0,
        temperature_c=500.0,
        bulk_composition={"A": 1.0},
        observed_phases=("A", "B"),
        phase_fractions={"A": 1.0, "B": 1.0},
        phase_fraction_sigma={"A": 1.0, "B": 1.0},
    )
    prior_sigma = 5000.0
    parameters = [
        ThermodynamicParameter("dH_A", "A", prior_sigma_j_mol=prior_sigma),
        ThermodynamicParameter("dH_B", "B", prior_sigma_j_mol=prior_sigma),
    ]
    problem = InversionProblem(QuadraticForwardModel(), [experiment], parameters)
    approximation = laplace_covariance(problem, [target_a, target_b])

    assert approximation.is_positive_definite
    expected_precision = 2.0 / sigma_data**2 + 1.0 / prior_sigma**2
    assert abs(approximation.hessian[0, 0] - expected_precision) < 1.0e-6
    assert abs(approximation.hessian[1, 1] - expected_precision) < 1.0e-6
    assert abs(approximation.hessian[0, 1]) < 1.0e-6
    assert abs(approximation.covariance[0, 0] - 1.0 / expected_precision) < 1.0e-6
    assert abs(approximation.correlation[0, 1]) < 1.0e-6


def test_laplace_covariance_uncertainty_source_isolates_data_and_prior():
    sigma_data = 200.0
    target_a, target_b = 900.0, -400.0

    class QuadraticForwardModel:
        def predict(self, experiment, corrections):
            value_a = corrections.get("A", 0.0)
            value_b = corrections.get("B", 0.0)
            return EquilibriumResult(
                True,
                ("A", "B"),
                {
                    "A": 1.0 + (value_a - target_a) / sigma_data,
                    "B": 1.0 + (value_b - target_b) / sigma_data,
                },
                {},
                {},
            )

    experiment = Experiment(
        id="quadratic-source",
        pressure_bar=2000.0,
        temperature_c=500.0,
        bulk_composition={"A": 1.0},
        observed_phases=("A", "B"),
        phase_fractions={"A": 1.0, "B": 1.0},
        phase_fraction_sigma={"A": 1.0, "B": 1.0},
    )
    prior_sigma = 5000.0
    parameters = [
        ThermodynamicParameter("dH_A", "A", prior_sigma_j_mol=prior_sigma),
        ThermodynamicParameter("dH_B", "B", prior_sigma_j_mol=prior_sigma),
    ]
    problem = InversionProblem(QuadraticForwardModel(), [experiment], parameters)

    data_only = laplace_covariance(
        problem, [target_a, target_b], uncertainty_source="data"
    )
    assert data_only.uncertainty_source == "data"
    expected_data_precision = 2.0 / sigma_data**2
    assert abs(data_only.hessian[0, 0] - expected_data_precision) < 1.0e-6

    prior_only = laplace_covariance(problem, [0.0, 0.0], uncertainty_source="prior")
    assert prior_only.uncertainty_source == "prior"
    expected_prior_precision = 1.0 / prior_sigma**2
    assert abs(prior_only.hessian[0, 0] - expected_prior_precision) < 1.0e-9

    try:
        laplace_covariance(problem, [0.0, 0.0], uncertainty_source="bogus")
    except ValueError:
        pass
    else:
        raise AssertionError("Expected ValueError for an invalid uncertainty_source")


def test_laplace_covariance_reports_hard_bound_limitation():
    parameter = ThermodynamicParameter(
        "dH_A",
        "A",
        prior_sigma_j_mol=2.0,
        lower_bound_j_mol=-1.0,
        upper_bound_j_mol=1.0,
    )
    problem = InversionProblem(object(), [], [parameter])

    boundary_result = laplace_covariance(problem, [1.0])
    interior_result = laplace_covariance(problem, [0.99])

    assert boundary_result.boundary_limited
    assert boundary_result.boundary_parameters == ("dH_A",)
    assert boundary_result.covariance is None
    assert not interior_result.boundary_limited
    assert np.isclose(interior_result.hessian[0, 0], 0.25)


def test_grid_posterior_matches_analytic_gaussian_moments():
    sigma_data = 200.0
    target = 900.0
    prior_sigma = 5000.0

    class QuadraticForwardModel:
        def predict(self, experiment, corrections):
            value = corrections.get("A", 0.0)
            return EquilibriumResult(
                True, ("A",), {"A": 1.0 + (value - target) / sigma_data}, {}, {}
            )

    experiment = Experiment(
        id="quadratic-1d",
        pressure_bar=2000.0,
        temperature_c=500.0,
        bulk_composition={"A": 1.0},
        observed_phases=("A",),
        phase_fractions={"A": 1.0},
        phase_fraction_sigma={"A": 1.0},
    )
    parameters = [ThermodynamicParameter("dH_A", "A", prior_sigma_j_mol=prior_sigma)]
    problem = InversionProblem(QuadraticForwardModel(), [experiment], parameters)

    expected_precision = 2.0 / sigma_data**2 + 1.0 / prior_sigma**2
    expected_std = 1.0 / np.sqrt(expected_precision)
    expected_mean = (2.0 * target / sigma_data**2) / expected_precision

    bounds = [(expected_mean - 6.0 * expected_std, expected_mean + 6.0 * expected_std)]
    result = grid_posterior(problem, bounds, resolution=401)

    assert abs(result.marginal_mean["dH_A"] - expected_mean) < 1.0
    assert abs(result.marginal_std["dH_A"] - expected_std) < 1.0
    assert np.isclose(
        np.sum(result.posterior_density) * (result.axes[0][1] - result.axes[0][0]),
        1.0,
        atol=1.0e-6,
    )


def test_grid_posterior_rejects_out_of_support_bounds():
    parameter = ThermodynamicParameter(
        "dH_A", "A", lower_bound_j_mol=-1.0, upper_bound_j_mol=1.0
    )
    problem = InversionProblem(object(), [], [parameter])
    try:
        grid_posterior(problem, [(-2.0, 0.5)], resolution=5)
    except ValueError:
        pass
    else:
        raise AssertionError("Expected ValueError for out-of-support grid bounds")


def test_mcmc_sample_matches_analytic_gaussian_moments():
    sigma_data = 200.0
    target = 900.0
    prior_sigma = 5000.0

    class QuadraticForwardModel:
        def predict(self, experiment, corrections):
            value = corrections.get("A", 0.0)
            return EquilibriumResult(
                True, ("A",), {"A": 1.0 + (value - target) / sigma_data}, {}, {}
            )

    experiment = Experiment(
        id="quadratic-mcmc",
        pressure_bar=2000.0,
        temperature_c=500.0,
        bulk_composition={"A": 1.0},
        observed_phases=("A",),
        phase_fractions={"A": 1.0},
        phase_fraction_sigma={"A": 1.0},
    )
    parameters = [ThermodynamicParameter("dH_A", "A", prior_sigma_j_mol=prior_sigma)]
    problem = InversionProblem(QuadraticForwardModel(), [experiment], parameters)

    expected_precision = 2.0 / sigma_data**2 + 1.0 / prior_sigma**2
    expected_std = 1.0 / np.sqrt(expected_precision)
    expected_mean = (2.0 * target / sigma_data**2) / expected_precision

    result = mcmc_sample(
        problem,
        [expected_mean],
        n_samples=20000,
        burn_in=2000,
        proposal_sigma_j_mol=[expected_std * 1.5],
        seed=0,
    )

    assert 0.1 < result.acceptance_rate < 0.9
    assert abs(result.posterior_mean["dH_A"] - expected_mean) < 10.0
    assert abs(result.posterior_std["dH_A"] - expected_std) < 15.0
    low, high = result.credible_interval_95["dH_A"]
    assert low < expected_mean < high


def test_mcmc_sample_prior_only_recovers_the_prior_distribution():
    prior_sigma = 5000.0

    class QuadraticForwardModel:
        def predict(self, experiment, corrections):
            value = corrections.get("A", 0.0)
            # Data term intentionally still varies with the correction, but
            # uncertainty_source="prior" must ignore it entirely.
            return EquilibriumResult(True, ("A",), {"A": 1.0 + value / 10.0}, {}, {})

    experiment = Experiment(
        id="prior-only",
        pressure_bar=2000.0,
        temperature_c=500.0,
        bulk_composition={"A": 1.0},
        observed_phases=("A",),
        phase_fractions={"A": 1.0},
        phase_fraction_sigma={"A": 1.0},
    )
    parameters = [ThermodynamicParameter("dH_A", "A", prior_sigma_j_mol=prior_sigma)]
    problem = InversionProblem(QuadraticForwardModel(), [experiment], parameters)

    result = mcmc_sample(
        problem,
        [0.0],
        n_samples=20000,
        burn_in=2000,
        proposal_sigma_j_mol=[prior_sigma * 0.8],
        seed=4,
        uncertainty_source="prior",
    )

    assert result.uncertainty_source == "prior"
    assert abs(result.posterior_mean["dH_A"] - 0.0) < 200.0
    assert abs(result.posterior_std["dH_A"] - prior_sigma) < 300.0


def test_mcmc_sample_supports_a_correlated_proposal_covariance():
    parameters = [
        ThermodynamicParameter("dH_A", "A", prior_sigma_j_mol=1.0),
        ThermodynamicParameter("dH_B", "B", prior_sigma_j_mol=1.0),
    ]
    covariance = np.asarray([[1.0, 0.8], [0.8, 1.0]])
    problem = InversionProblem(
        object(), [], parameters, prior_covariance_j_mol2=covariance
    )

    result = mcmc_sample(
        problem,
        [0.0, 0.0],
        n_samples=10000,
        burn_in=1000,
        proposal_covariance_j_mol2=covariance * 0.8,
        seed=6,
    )

    assert 0.1 < result.acceptance_rate < 0.9
    assert np.corrcoef(result.samples.T)[0, 1] > 0.65
    assert np.allclose(result.proposal_covariance_j_mol2, covariance * 0.8)
    assert result.seed == 6
    assert result.burn_in == 1000

    try:
        mcmc_sample(
            problem,
            [0.0, 0.0],
            n_samples=10,
            burn_in=1,
            proposal_sigma_j_mol=[1.0, 1.0],
            proposal_covariance_j_mol2=covariance,
        )
    except ValueError:
        pass
    else:
        raise AssertionError("Expected ValueError for two proposal specifications")


def test_objective_and_mcmc_enforce_declared_parameter_bounds():
    parameter = ThermodynamicParameter(
        "dH_A",
        "A",
        prior_sigma_j_mol=10.0,
        lower_bound_j_mol=-1.0,
        upper_bound_j_mol=1.0,
    )
    problem = InversionProblem(object(), [], [parameter])

    assert np.isinf(problem.objective([1.1]))
    result = mcmc_sample(
        problem,
        [0.0],
        n_samples=200,
        burn_in=10,
        proposal_sigma_j_mol=[5.0],
        seed=5,
    )
    assert np.all(result.samples >= -1.0)
    assert np.all(result.samples <= 1.0)

    try:
        mcmc_sample(problem, [1.1], n_samples=10, burn_in=1)
    except ValueError:
        pass
    else:
        raise AssertionError("Expected ValueError for an out-of-bounds MCMC start")


def test_covariance_paths_reject_nonfinite_matrices():
    parameter = ThermodynamicParameter("dH_A", "A", prior_sigma_j_mol=1.0)
    try:
        InversionProblem(object(), [], [parameter], prior_covariance_j_mol2=[[np.nan]])
    except ValueError:
        pass
    else:
        raise AssertionError("Expected ValueError for non-finite prior covariance")

    problem = InversionProblem(object(), [], [parameter])
    try:
        mcmc_sample(
            problem,
            [0.0],
            n_samples=10,
            burn_in=1,
            proposal_covariance_j_mol2=[[np.nan]],
        )
    except ValueError:
        pass
    else:
        raise AssertionError("Expected ValueError for non-finite proposal covariance")


def test_mcmc_rejects_nonfinite_initial_points_and_proposal_sigmas():
    parameter = ThermodynamicParameter("dH_A", "A", prior_sigma_j_mol=1.0)
    problem = InversionProblem(object(), [], [parameter])
    for initial_point, proposal_sigma in [([np.nan], [1.0]), ([0.0], [np.inf])]:
        try:
            mcmc_sample(
                problem,
                initial_point,
                n_samples=10,
                burn_in=1,
                proposal_sigma_j_mol=proposal_sigma,
            )
        except ValueError:
            pass
        else:
            raise AssertionError("Expected ValueError for non-finite MCMC input")


def test_leave_one_out_cross_validation_recovers_consistent_correction():
    class FakeForwardModel:
        def predict(self, experiment, corrections):
            error = corrections.get("Znc", 0.0) - 1500.0
            valid = abs(error) < 2500.0
            return EquilibriumResult(
                valid, ("Znc",), {"Znc": max(0.0, 1.0 + error / 1000.0)}, {}, {}
            )

    experiments = [
        Experiment(
            id=f"loocv-{i}",
            pressure_bar=2000.0,
            temperature_c=500.0,
            bulk_composition={"Zn": 1.0},
            observed_phases=("Znc",),
            phase_fractions={"Znc": 1.0},
            phase_fraction_sigma={"Znc": 0.01},
        )
        for i in range(3)
    ]
    parameter = ThermodynamicParameter(
        name="delta_H_Znc",
        phase="Znc",
        prior_mean_j_mol=0.0,
        prior_sigma_j_mol=5000.0,
        lower_bound_j_mol=-4000.0,
        upper_bound_j_mol=4000.0,
    )
    problem = InversionProblem(FakeForwardModel(), experiments, [parameter])
    result = leave_one_out_cross_validation(
        problem, optimize_starts=3, optimize_iterations=30
    )

    assert len(result.folds) == 3
    assert result.mean_misfit < 1.0
    for fold in result.folds:
        assert abs(fold.corrections_j_mol["Znc"] - 1500.0) < 100.0


def test_leave_one_group_out_cross_validation_holds_out_studies():
    class StudyForwardModel:
        def predict(self, experiment, corrections):
            return EquilibriumResult(
                True,
                ("A",),
                {"A": corrections.get("A", 0.0) / 1000.0},
                {},
                {},
            )

    experiments = [
        Experiment(
            id=f"study-{study}-{index}",
            pressure_bar=2000.0,
            temperature_c=500.0,
            bulk_composition={"A": 1.0},
            observed_phases=("A",),
            phase_fractions={"A": 1.0},
            phase_fraction_sigma={"A": 1.0},
            metadata={"study": study},
        )
        for study in ("one", "two")
        for index in range(2)
    ]
    parameter = ThermodynamicParameter(
        "dH_A",
        "A",
        prior_sigma_j_mol=5000.0,
        lower_bound_j_mol=-4000.0,
        upper_bound_j_mol=4000.0,
    )
    problem = InversionProblem(
        StudyForwardModel(), experiments, [parameter], misfit_form="raw_lsq"
    )

    result = leave_one_group_out_cross_validation(
        problem, "study", optimize_starts=2, optimize_iterations=20
    )

    assert {fold.held_out_experiment_id for fold in result.folds} == {
        "study=one",
        "study=two",
    }
    assert len(result.folds) == 2
    assert all(len(fold.held_out_experiment_ids) == 2 for fold in result.folds)
    assert all(len(fold.held_out_misfits) == 2 for fold in result.folds)
    assert result.mean_misfit < 1.0e-3


def test_leave_one_phase_assemblage_out_cross_validation_groups_observed_phases():
    class AssemblageForwardModel:
        def predict(self, experiment, corrections):
            return EquilibriumResult(
                True,
                experiment.observed_phases,
                {phase: 1.0 for phase in experiment.observed_phases},
                {},
                {},
            )

    experiments = [
        Experiment(
            f"assemblage-{index}",
            2000.0,
            500.0,
            {"A": 1.0},
            phases,
            phase_fractions={phases[0]: 1.0},
            phase_fraction_sigma={phases[0]: 1.0},
        )
        for index, phases in enumerate((("A",), ("A",), ("B",), ("B",)))
    ]
    problem = InversionProblem(
        AssemblageForwardModel(),
        experiments,
        [ThermodynamicParameter("dH_A", "A", prior_sigma_j_mol=10.0)],
        misfit_form="raw_lsq",
    )

    result = leave_one_phase_assemblage_out_cross_validation(
        problem, optimize_starts=1, optimize_iterations=5
    )

    assert len(result.folds) == 2
    assert result.group_metadata_key == "__thermoinvert_phase_assemblage__"
    assert all(len(fold.held_out_experiment_ids) == 2 for fold in result.folds)


def test_out_of_domain_pt_cross_validation_holds_out_external_rectangle():
    class DomainForwardModel:
        def predict(self, experiment, corrections):
            return EquilibriumResult(
                True,
                ("A",),
                {"A": 1.0},
                {},
                {},
            )

    experiments = [
        Experiment(
            "inside",
            2000.0,
            500.0,
            {"A": 1.0},
            ("A",),
            phase_fractions={"A": 1.0},
            phase_fraction_sigma={"A": 1.0},
        ),
        Experiment(
            "outside",
            5000.0,
            800.0,
            {"A": 1.0},
            ("A",),
            phase_fractions={"A": 1.0},
            phase_fraction_sigma={"A": 1.0},
        ),
    ]
    problem = InversionProblem(
        DomainForwardModel(),
        experiments,
        [ThermodynamicParameter("dH_A", "A", prior_sigma_j_mol=10.0)],
    )

    result = leave_one_out_of_domain_pt_cross_validation(
        problem,
        (1000.0, 3000.0),
        (400.0, 600.0),
        optimize_starts=1,
        optimize_iterations=5,
    )

    assert len(result.folds) == 1
    assert result.folds[0].held_out_experiment_ids == ("outside",)
    assert result.group_metadata_key == "out_of_domain_pt"


def test_leave_one_database_family_out_cross_validation_holds_out_families():
    class DatabaseFamilyForwardModel:
        def predict(self, experiment, corrections):
            return EquilibriumResult(True, ("A",), {"A": 1.0}, {}, {})

    experiments = [
        Experiment(
            f"family-{family}",
            2000.0,
            500.0,
            {"A": 1.0},
            ("A",),
            phase_fractions={"A": 1.0},
            phase_fraction_sigma={"A": 1.0},
            metadata={"database_family": family},
        )
        for family in ("legacy", "new")
    ]
    problem = InversionProblem(
        DatabaseFamilyForwardModel(),
        experiments,
        [ThermodynamicParameter("dH_A", "A", prior_sigma_j_mol=10.0)],
    )

    result = leave_one_database_family_out_cross_validation(
        problem, optimize_starts=1, optimize_iterations=5
    )

    assert result.group_metadata_key == "database_family"
    assert {fold.held_out_experiment_id for fold in result.folds} == {
        "database_family=legacy",
        "database_family=new",
    }
    repeated = repeated_group_holdout_cross_validation(
        problem,
        "database_family",
        repeats=3,
        test_fraction=0.5,
        seed=4,
        optimize_starts=1,
        optimize_iterations=5,
    )
    assert len(repeated.folds) == 3
    assert all(fold.held_out_experiment_ids for fold in repeated.folds)


def test_database_provenance_derives_family_from_file_hash(tmp_path):
    first = tmp_path / "first.json"
    second = tmp_path / "second.json"
    first.write_text('{"database": "one"}', encoding="utf-8")
    second.write_text('{"database": "two"}', encoding="utf-8")
    first_provenance = database_file_provenance(first)
    experiments = [
        Experiment(
            f"db-{index}",
            2000.0,
            500.0,
            {"A": 1.0},
            ("A",),
            phase_fractions={"A": 1.0},
            phase_fraction_sigma={"A": 1.0},
            metadata={"database_path": str(path)},
        )
        for index, path in enumerate((first, second))
    ]

    records = validate_database_provenance(experiments)

    assert records[0]["database_sha256"] == first_provenance["database_sha256"]
    assert records[0]["database_sha256"] != records[1]["database_sha256"]

    class FileDatabaseForwardModel:
        def predict(self, experiment, corrections):
            return EquilibriumResult(True, ("A",), {"A": 1.0}, {}, {})

    problem = InversionProblem(
        FileDatabaseForwardModel(),
        experiments,
        [ThermodynamicParameter("dH_A", "A", prior_sigma_j_mol=10.0)],
    )
    result = leave_one_database_family_out_cross_validation(
        problem, optimize_starts=1, optimize_iterations=3
    )
    assert len(result.folds) == 2
    assert all(
        fold.held_out_experiment_id.startswith("database_family=")
        for fold in result.folds
    )


def test_loocv_preserves_parent_misfit_form():
    class LinearForwardModel:
        def predict(self, experiment, corrections):
            return EquilibriumResult(
                True,
                ("A",),
                {"A": 1.0},
                {},
                {"A": corrections.get("A", 0.0) / 1000.0},
            )

    experiments = [
        Experiment(
            id="raw-lsq-1",
            pressure_bar=2000.0,
            temperature_c=500.0,
            bulk_composition={"A": 1.0},
            observed_phases=("A",),
            phase_fractions={"A": 0.1},
            phase_fraction_sigma={"A": 0.1},
        ),
        Experiment(
            id="raw-lsq-2",
            pressure_bar=2000.0,
            temperature_c=500.0,
            bulk_composition={"A": 1.0},
            observed_phases=("A",),
            phase_fractions={"A": 0.9},
            phase_fraction_sigma={"A": 0.1},
        ),
    ]
    parameter = ThermodynamicParameter(
        "dH_A",
        "A",
        prior_sigma_j_mol=1.0e6,
        lower_bound_j_mol=-1.0e6,
        upper_bound_j_mol=1.0e6,
    )
    problem = InversionProblem(
        LinearForwardModel(),
        experiments,
        [parameter],
        objective_criterion="bayes",
        misfit_form="raw_lsq",
    )

    result = leave_one_out_cross_validation(
        problem, optimize_starts=2, optimize_iterations=40
    )

    assert 0.35 < result.mean_misfit < 0.45


def test_posterior_predictive_check_covers_observed_value():
    sigma_data = 200.0
    target = 900.0
    prior_sigma = 5000.0

    class QuadraticForwardModel:
        def predict(self, experiment, corrections):
            value = corrections.get("A", 0.0)
            return EquilibriumResult(
                True, ("A",), {"A": 1.0 + (value - target) / sigma_data}, {}, {}
            )

    experiment = Experiment(
        id="quadratic-ppc",
        pressure_bar=2000.0,
        temperature_c=500.0,
        bulk_composition={"A": 1.0},
        observed_phases=("A",),
        phase_fractions={"A": 1.0},
        phase_fraction_sigma={"A": 1.0},
    )
    parameters = [ThermodynamicParameter("dH_A", "A", prior_sigma_j_mol=prior_sigma)]
    problem = InversionProblem(QuadraticForwardModel(), [experiment], parameters)

    expected_precision = 2.0 / sigma_data**2 + 1.0 / prior_sigma**2
    expected_std = 1.0 / np.sqrt(expected_precision)
    expected_mean = (2.0 * target / sigma_data**2) / expected_precision

    mcmc_result = mcmc_sample(
        problem,
        [expected_mean],
        n_samples=6000,
        burn_in=1000,
        proposal_sigma_j_mol=[expected_std * 1.5],
        seed=1,
    )
    ppc = posterior_predictive_check(
        problem, experiment, mcmc_result, n_predictive_samples=100
    )

    assert ppc.observable_names == ("phase_present:A", "phase_fraction:A")
    fraction_index = ppc.observable_names.index("phase_fraction:A")
    assert abs(ppc.predictive_mean[fraction_index] - 1.0) < 0.2
    assert ppc.within_95["phase_fraction:A"]
    assert ppc.within_95["phase_present:A"]


def test_posterior_predictive_check_can_include_observation_uncertainty():
    class ConstantForwardModel:
        def predict(self, experiment, corrections):
            return EquilibriumResult(True, ("A",), {"A": 1.0}, {}, {"A": 0.5})

    experiment = Experiment(
        id="predictive-noise",
        pressure_bar=2000.0,
        temperature_c=500.0,
        bulk_composition={"A": 1.0},
        observed_phases=("A",),
        phase_fractions={"A": 0.5},
        phase_fraction_sigma={"A": 0.2},
    )
    parameter = ThermodynamicParameter("dH_A", "A", prior_sigma_j_mol=1.0)
    problem = InversionProblem(ConstantForwardModel(), [experiment], [parameter])
    mcmc_result = MCMCResult(
        ("dH_A",),
        np.zeros((100, 1)),
        np.zeros(100),
        1.0,
        {"dH_A": 0.0},
        {"dH_A": 0.0},
        {"dH_A": (0.0, 0.0)},
    )

    without_noise = posterior_predictive_check(problem, experiment, mcmc_result, 100)
    with_noise = posterior_predictive_check(
        problem,
        experiment,
        mcmc_result,
        100,
        include_observation_uncertainty=True,
        seed=12,
    )
    repeated = posterior_predictive_check(
        problem,
        experiment,
        mcmc_result,
        100,
        include_observation_uncertainty=True,
        seed=12,
    )

    assert without_noise.predictive_std[1] == 0.0
    assert with_noise.predictive_std[1] > 0.1
    assert np.allclose(with_noise.predictive_mean, repeated.predictive_mean)
    assert with_noise.include_observation_uncertainty
    assert with_noise.seed == 12


def test_posterior_predictive_check_supports_seeded_random_sample_selection():
    class SampleForwardModel:
        def predict(self, experiment, corrections):
            value = corrections.get("A", 0.0)
            return EquilibriumResult(True, ("A",), {"A": 1.0}, {}, {"A": value})

    experiment = Experiment(
        id="random-selection",
        pressure_bar=2000.0,
        temperature_c=500.0,
        bulk_composition={"A": 1.0},
        observed_phases=("A",),
        phase_fractions={"A": 1.0},
        phase_fraction_sigma={"A": 1.0},
    )
    parameter = ThermodynamicParameter("dH_A", "A", prior_sigma_j_mol=1.0)
    problem = InversionProblem(SampleForwardModel(), [experiment], [parameter])
    samples = np.arange(100.0).reshape(-1, 1)
    mcmc_result = MCMCResult(("dH_A",), samples, np.zeros(100), 1.0, {}, {}, {})

    first = posterior_predictive_check(
        problem, experiment, mcmc_result, 20, sample_selection="random", seed=9
    )
    second = posterior_predictive_check(
        problem, experiment, mcmc_result, 20, sample_selection="random", seed=9
    )

    assert first.sample_selection == "random"
    assert first.seed == 9
    assert np.allclose(first.predictive_mean, second.predictive_mean)


def test_posterior_predictive_noise_rejects_incomplete_covariance_keys():
    class ConstantForwardModel:
        def predict(self, experiment, corrections):
            return EquilibriumResult(
                True, ("A",), {"A": 1.0}, {"A": {"x": 0.2}}, {"A": 0.5}
            )

    experiment = Experiment(
        id="invalid-predictive-covariance",
        pressure_bar=2000.0,
        temperature_c=500.0,
        bulk_composition={"A": 1.0},
        observed_phases=("A",),
        phase_fractions={"A": 0.5},
        phase_compositions={"A": {"x": 0.2}},
        observation_covariance=[[0.01]],
        observation_covariance_keys=("fraction:A",),
    )
    parameter = ThermodynamicParameter("dH_A", "A", prior_sigma_j_mol=1.0)
    problem = InversionProblem(ConstantForwardModel(), [experiment], [parameter])
    mcmc_result = MCMCResult(("dH_A",), np.zeros((4, 1)), np.zeros(4), 1.0, {}, {}, {})
    try:
        posterior_predictive_check(
            problem,
            experiment,
            mcmc_result,
            include_observation_uncertainty=True,
        )
    except ValueError:
        pass
    else:
        raise AssertionError(
            "Expected ValueError for incomplete predictive covariance keys"
        )


def test_posterior_predictive_noise_accepts_rank_deficient_covariance():
    class ConstantForwardModel:
        def predict(self, experiment, corrections):
            return EquilibriumResult(
                True,
                ("A",),
                {"A": 1.0},
                {"A": {"x": 0.2, "y": 0.4}},
                {"A": 0.5},
            )

    experiment = Experiment(
        id="rank-deficient-predictive-covariance",
        pressure_bar=2000.0,
        temperature_c=500.0,
        bulk_composition={"A": 1.0},
        observed_phases=("A",),
        phase_compositions={"A": {"x": 0.2, "y": 0.4}},
        observation_covariance=[[0.01, 0.01], [0.01, 0.01]],
        observation_covariance_keys=("composition:A:x", "composition:A:y"),
    )
    parameter = ThermodynamicParameter("dH_A", "A", prior_sigma_j_mol=1.0)
    problem = InversionProblem(ConstantForwardModel(), [experiment], [parameter])
    mcmc_result = MCMCResult(
        ("dH_A",), np.zeros((20, 1)), np.zeros(20), 1.0, {}, {}, {}
    )

    result = posterior_predictive_check(
        problem,
        experiment,
        mcmc_result,
        include_observation_uncertainty=True,
        seed=4,
    )

    assert result.predictive_std[1] > 0.0


def _fake_mcmc_result_with_outliers():
    from thermoinvert.posterior import MCMCResult

    rng = np.random.default_rng(3)
    inliers = rng.normal(loc=0.0, scale=1.0, size=(200, 1))
    outliers = np.array([[50.0], [60.0], [70.0]])
    samples = np.vstack([inliers, outliers])
    objective = np.sum(samples**2, axis=1)
    return MCMCResult(
        parameter_names=("p0",),
        samples=samples,
        log_posterior=-objective,
        acceptance_rate=0.4,
        posterior_mean={"p0": float(np.mean(samples))},
        posterior_std={"p0": float(np.std(samples))},
        credible_interval_95={"p0": (-2.0, 2.0)},
    )


def test_filter_outliers_mcmc_drops_high_objective_tail():
    result = _fake_mcmc_result_with_outliers()
    filtered = filter_outliers_mcmc(result, z_threshold=3.0)
    assert filtered.samples.shape[0] < result.samples.shape[0]
    assert np.max(filtered.samples) < 10.0


def test_filter_by_acceptance_threshold_keeps_only_low_objective_samples():
    result = _fake_mcmc_result_with_outliers()
    max_objective = 10.0
    filtered = filter_by_acceptance_threshold(result, max_objective)
    assert np.all(-filtered.log_posterior <= max_objective)
    assert filtered.samples.shape[0] < result.samples.shape[0]


def test_filter_by_acceptance_threshold_raises_when_nothing_passes():
    result = _fake_mcmc_result_with_outliers()
    try:
        filter_by_acceptance_threshold(result, max_objective=-1.0)
    except ValueError:
        return
    raise AssertionError("Expected ValueError when no samples satisfy the threshold")


def test_effective_sample_size_is_near_n_for_iid_samples():
    rng = np.random.default_rng(5)
    iid_samples = rng.normal(size=5000)
    ess = effective_sample_size(iid_samples)
    assert ess > 0.7 * len(iid_samples)


def test_effective_sample_size_is_much_smaller_for_a_random_walk():
    rng = np.random.default_rng(5)
    random_walk = np.cumsum(rng.normal(scale=0.1, size=5000))
    ess = effective_sample_size(random_walk)
    assert ess < 0.1 * len(random_walk)


def test_potential_scale_reduction_near_one_for_well_mixed_chains():
    rng = np.random.default_rng(6)
    chains = rng.normal(loc=0.0, scale=1.0, size=(4, 2000))
    r_hat = potential_scale_reduction(chains)
    assert abs(r_hat - 1.0) < 0.05


def test_potential_scale_reduction_is_large_for_poorly_mixed_chains():
    rng = np.random.default_rng(6)
    chains = np.stack(
        [
            rng.normal(loc=offset, scale=0.1, size=2000)
            for offset in (-5.0, -2.0, 2.0, 5.0)
        ]
    )
    r_hat = potential_scale_reduction(chains)
    assert r_hat > 1.5


def test_mcmc_sample_multi_chain_reports_good_convergence_on_quadratic_objective():
    sigma_data = 200.0
    target = 900.0
    prior_sigma = 5000.0

    class QuadraticForwardModel:
        def predict(self, experiment, corrections):
            value = corrections.get("A", 0.0)
            return EquilibriumResult(
                True, ("A",), {"A": 1.0 + (value - target) / sigma_data}, {}, {}
            )

    experiment = Experiment(
        id="quadratic-multichain",
        pressure_bar=2000.0,
        temperature_c=500.0,
        bulk_composition={"A": 1.0},
        observed_phases=("A",),
        phase_fractions={"A": 1.0},
        phase_fraction_sigma={"A": 1.0},
    )
    parameters = [ThermodynamicParameter("dH_A", "A", prior_sigma_j_mol=prior_sigma)]
    problem = InversionProblem(QuadraticForwardModel(), [experiment], parameters)

    expected_precision = 2.0 / sigma_data**2 + 1.0 / prior_sigma**2
    expected_std = 1.0 / np.sqrt(expected_precision)
    expected_mean = (2.0 * target / sigma_data**2) / expected_precision

    result = mcmc_sample_multi_chain(
        problem,
        initial_points=[
            [expected_mean - 2 * expected_std],
            [expected_mean + 2 * expected_std],
        ],
        n_samples=6000,
        burn_in=1000,
        proposal_sigma_j_mol=[expected_std * 1.5],
        seeds=[10, 11],
    )
    assert len(result.chains) == 2
    assert result.initial_points == (
        (expected_mean - 2 * expected_std,),
        (expected_mean + 2 * expected_std,),
    )
    assert result.seeds == (10, 11)
    assert result.parameter_units == ("J/mol",)
    assert result.r_hat["dH_A"] < 1.1
    assert result.effective_sample_size["dH_A"] > 100.0
