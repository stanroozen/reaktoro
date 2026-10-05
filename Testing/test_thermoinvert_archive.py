"""Round-trip tests for persistent CSV/JSON archival of MCMC/LOOCV/grid results."""

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parents[1]))

from thermoinvert.archive import (
    read_grid_posterior,
    read_hierarchical_mcmc_result,
    read_grouped_hierarchical_mcmc_result,
    read_grouped_observation_bias_mcmc_result,
    read_combined_grouped_uncertainty_mcmc_result,
    read_bracket_posterior_predictive_result,
    read_grouped_observation_bias_report,
    read_grouped_observation_bias_summary,
    read_inversion_configuration,
    read_joint_inversion_configuration,
    read_pt_inversion_configuration,
    read_loocv_result,
    read_mcmc_result,
    read_multi_chain_mcmc_result,
    read_parallel_tempering_result,
    read_observation_bias_profile_result,
    read_posterior_predictive_result,
    read_posterior_resolution_result,
    write_grouped_observation_bias_report,
    write_grouped_observation_bias_summary,
    write_grid_posterior,
    write_hierarchical_mcmc_result,
    write_grouped_hierarchical_mcmc_result,
    write_grouped_observation_bias_mcmc_result,
    write_combined_grouped_uncertainty_mcmc_result,
    write_bracket_posterior_predictive_result,
    write_inversion_configuration,
    write_joint_inversion_configuration,
    write_loocv_result,
    write_mcmc_result,
    write_multi_chain_mcmc_result,
    write_parallel_tempering_result,
    write_observation_bias_profile_result,
    write_posterior_predictive_result,
    write_posterior_resolution_result,
)
from thermoinvert.data import (
    DirichletConcentrationParameter,
    Experiment,
    Parameter,
    ThermodynamicParameter,
)
from thermoinvert.hierarchical import (
    CombinedGroupedUncertaintyMCMCResult,
    GroupedHierarchicalMCMCResult,
    GroupedObservationBiasMCMCResult,
    HierarchicalMCMCResult,
)
from thermoinvert.joint_inversion import JointInversionProblem
from thermoinvert.posterior import ParallelTemperingResult
from thermoinvert.bias import (
    GroupedObservationBiasReport,
    GroupedObservationBiasSummary,
    ObservationBiasProfileResult,
)
from thermoinvert.inversion import InversionProblem
from thermoinvert.posterior import GridPosteriorResult, MCMCResult, MultiChainMCMCResult
from thermoinvert.sensitivity import PosteriorResolutionResult
from thermoinvert.validation import (
    BracketPosteriorPredictiveResult,
    LOOCVFold,
    LOOCVResult,
    PosteriorPredictiveResult,
)


def test_mcmc_result_round_trips_through_disk(tmp_path):
    rng = np.random.default_rng(0)
    samples = rng.normal(size=(50, 2))
    original = MCMCResult(
        parameter_names=("dH_A", "dH_B"),
        samples=samples,
        log_posterior=-np.sum(samples**2, axis=1),
        acceptance_rate=0.42,
        posterior_mean={
            "dH_A": float(np.mean(samples[:, 0])),
            "dH_B": float(np.mean(samples[:, 1])),
        },
        posterior_std={
            "dH_A": float(np.std(samples[:, 0])),
            "dH_B": float(np.std(samples[:, 1])),
        },
        credible_interval_95={"dH_A": (-2.0, 2.0), "dH_B": (-3.0, 3.0)},
        uncertainty_source="data",
        proposal_covariance_j_mol2=np.asarray([[0.25, 0.1], [0.1, 0.5]]),
        seed=17,
        burn_in=10,
        parameter_units=("J/mol", "J/mol"),
    )
    write_mcmc_result(original, tmp_path)
    restored = read_mcmc_result(tmp_path)

    assert restored.parameter_names == original.parameter_names
    assert np.allclose(restored.samples, original.samples, atol=1.0e-6)
    assert np.allclose(restored.log_posterior, original.log_posterior, atol=1.0e-6)
    assert restored.acceptance_rate == original.acceptance_rate
    assert restored.posterior_mean == original.posterior_mean
    assert restored.uncertainty_source == "data"
    assert restored.credible_interval_95["dH_A"] == (-2.0, 2.0)
    assert restored.proposal_covariance_j_mol2 is not None
    assert original.proposal_covariance_j_mol2 is not None
    assert np.allclose(
        restored.proposal_covariance_j_mol2, original.proposal_covariance_j_mol2
    )
    assert restored.proposal_covariance is not None
    assert original.proposal_covariance is not None
    assert np.allclose(restored.proposal_covariance, original.proposal_covariance)
    assert restored.seed == 17
    assert restored.burn_in == 10
    assert restored.parameter_units == ("J/mol", "J/mol")

    legacy_summary = json.loads((tmp_path / "summary.json").read_text(encoding="utf-8"))
    assert (
        legacy_summary["proposal_covariance"]
        == legacy_summary["proposal_covariance_j_mol2"]
    )
    legacy_dir = tmp_path / "legacy"
    legacy_dir.mkdir()
    (legacy_dir / "samples.csv").write_bytes((tmp_path / "samples.csv").read_bytes())
    (legacy_dir / "log_posterior.csv").write_bytes(
        (tmp_path / "log_posterior.csv").read_bytes()
    )
    legacy_summary.pop("proposal_covariance")
    (legacy_dir / "summary.json").write_text(
        json.dumps(legacy_summary), encoding="utf-8"
    )
    legacy_restored = read_mcmc_result(legacy_dir)
    assert legacy_restored.proposal_covariance is not None
    assert original.proposal_covariance is not None
    assert np.allclose(
        legacy_restored.proposal_covariance, original.proposal_covariance
    )


def test_hierarchical_mcmc_result_round_trips_through_disk(tmp_path):
    samples = np.asarray([[1.0, 0.8], [2.0, 1.2], [1.5, 1.0]])
    original = HierarchicalMCMCResult(
        ("dH_A", "scale_C_T"),
        ("J/mol", "dimensionless"),
        samples,
        -np.sum(samples**2, axis=1),
        0.4,
        {"dH_A": 1.5, "scale_C_T": 1.0},
        {"dH_A": 0.5, "scale_C_T": 0.2},
        {"dH_A": (1.0, 2.0), "scale_C_T": (0.7, 1.3)},
        "scale_C_T",
        1.0,
        0.5,
        (0.0, 5.0),
        7,
        1,
    )
    write_hierarchical_mcmc_result(original, tmp_path)
    restored = read_hierarchical_mcmc_result(tmp_path)

    assert restored.parameter_names == original.parameter_names
    assert restored.parameter_units == original.parameter_units
    assert np.allclose(restored.samples, original.samples)
    assert np.allclose(restored.log_posterior, original.log_posterior)
    assert restored.scale_name == "scale_C_T"
    assert restored.scale_bounds == (0.0, 5.0)


def test_grouped_hierarchical_mcmc_result_round_trips_through_disk(tmp_path):
    samples = np.asarray([[1.0, 0.8, 1.2], [2.0, 1.1, 0.9]])
    original = GroupedHierarchicalMCMCResult(
        ("dH_A", "scale_C_T[one]", "scale_C_T[two]"),
        ("J/mol", "dimensionless", "dimensionless"),
        samples,
        -np.sum(samples**2, axis=1),
        0.4,
        {"dH_A": 1.5, "scale_C_T[one]": 0.95, "scale_C_T[two]": 1.05},
        {"dH_A": 0.5, "scale_C_T[one]": 0.2, "scale_C_T[two]": 0.2},
        {name: (0.1, 2.0) for name in ("dH_A", "scale_C_T[one]", "scale_C_T[two]")},
        "study",
        ("one", "two"),
        {"one": 1.0, "two": 1.0},
        {"one": 0.5, "two": 0.5},
        (0.0, 5.0),
        9,
        1,
    )
    write_grouped_hierarchical_mcmc_result(original, tmp_path)
    restored = read_grouped_hierarchical_mcmc_result(tmp_path)

    assert restored.metadata_key == "study"
    assert restored.scale_names == ("one", "two")
    assert np.allclose(restored.samples, original.samples)


def test_grouped_observation_bias_mcmc_result_round_trips_through_disk(tmp_path):
    samples = np.asarray([[1.0, 0.01, -0.02], [2.0, 0.02, 0.01]])
    original = GroupedObservationBiasMCMCResult(
        ("dH_A", "bias[one]", "bias[two]"),
        ("J/mol", "native observation unit", "native observation unit"),
        samples,
        -np.sum(samples**2, axis=1),
        0.4,
        {"dH_A": 1.5, "bias[one]": 0.015, "bias[two]": -0.005},
        {"dH_A": 0.5, "bias[one]": 0.005, "bias[two]": 0.015},
        {name: (-0.1, 0.1) for name in ("dH_A", "bias[one]", "bias[two]")},
        "laboratory",
        "fraction:A",
        ("one", "two"),
        0.0,
        0.1,
        (-0.1, 0.1),
        11,
        1,
    )
    write_grouped_observation_bias_mcmc_result(original, tmp_path)
    restored = read_grouped_observation_bias_mcmc_result(tmp_path)

    assert restored.metadata_key == "laboratory"
    assert restored.observation_key == "fraction:A"
    assert restored.bias_names == ("one", "two")
    assert restored.observation_keys == ("fraction:A",)
    assert np.allclose(restored.samples, original.samples)


def test_combined_grouped_uncertainty_mcmc_result_round_trips_through_disk(tmp_path):
    samples = np.asarray([[1.0, 0.9, 1.1, 0.01, -0.02]])
    original = CombinedGroupedUncertaintyMCMCResult(
        ("dH_A", "scale_C_T[one]", "scale_C_T[two]", "bias[one]", "bias[two]"),
        ("J/mol", "dimensionless", "dimensionless", "native", "native"),
        samples,
        np.asarray([-1.0]),
        0.5,
        {
            name: float(samples[0, index])
            for index, name in enumerate(
                ("dH_A", "scale_C_T[one]", "scale_C_T[two]", "bias[one]", "bias[two]")
            )
        },
        {
            name: 0.0
            for name in (
                "dH_A",
                "scale_C_T[one]",
                "scale_C_T[two]",
                "bias[one]",
                "bias[two]",
            )
        },
        {
            name: (-1.0, 1.0)
            for name in (
                "dH_A",
                "scale_C_T[one]",
                "scale_C_T[two]",
                "bias[one]",
                "bias[two]",
            )
        },
        "study",
        "fraction:A",
        ("one", "two"),
        (0.0, 5.0),
        1.0,
        1.0,
        (-0.1, 0.1),
        0.0,
        0.1,
        4,
        0,
    )
    write_combined_grouped_uncertainty_mcmc_result(original, tmp_path)
    restored = read_combined_grouped_uncertainty_mcmc_result(tmp_path)

    assert restored.metadata_key == "study"
    assert restored.observation_key == "fraction:A"
    assert restored.groups == ("one", "two")
    assert restored.observation_keys == ("fraction:A",)
    assert restored.bias_labels == ()
    assert np.allclose(restored.samples, original.samples)


def test_inversion_archive_preserves_generic_parameter_metadata(tmp_path):
    experiment = Experiment(
        "generic-parameter",
        1000.0,
        500.0,
        {"A": 1.0},
        ("A",),
        phase_fractions={"A": 0.5},
    )
    problem = InversionProblem(
        object(),
        [experiment],
        [Parameter("pressure", 3000.0, 500.0, 500.0, 8000.0, "bar")],
    )

    write_inversion_configuration(problem, tmp_path)
    restored = read_inversion_configuration(tmp_path)

    record = restored["configuration"]["parameters"][0]
    assert record["parameter_kind"] == "scalar"
    assert record["name"] == "pressure"
    assert record["unit"] == "bar"
    assert record["key"] == "pressure"


def test_parallel_tempering_result_round_trips_through_disk(tmp_path):
    original = ParallelTemperingResult(
        ("p",),
        ("J/mol",),
        (1.0, 4.0),
        np.asarray([[-1.0], [1.0], [-0.9]]),
        np.asarray([-2.0, -3.0, -2.1]),
        np.asarray([0, 1, 0]),
        {0: 2.0 / 3.0, 1: 1.0 / 3.0},
        {0: {"p": -0.95}, 1: {"p": 1.0}},
        {0: 2.0, 1: 3.0},
        (0.5, 0.7),
        (0.25,),
        True,
        3,
        1,
        (0.5,),
        {},
        0,
        0.0,
        np.asarray([[[-1.0], [1.0]], [[-0.9], [0.9]]]),
        np.asarray([[0, 1], [0, 1]]),
        np.asarray([[2.0, 3.0], [2.1, 3.1]]),
        np.asarray([1.8, 2.2, 2.7]),
    )
    write_parallel_tempering_result(original, tmp_path)
    restored = read_parallel_tempering_result(tmp_path)

    assert restored.parameter_names == original.parameter_names
    assert restored.temperatures == original.temperatures
    assert np.allclose(restored.samples, original.samples)
    assert np.array_equal(restored.mode_labels, original.mode_labels)
    assert restored.mode_occupancy == original.mode_occupancy
    assert restored.crossed_modes is True
    assert restored.temperature_mode_occupancy == (
        original.temperature_mode_occupancy or {}
    )
    assert restored.mode_transition_count == original.mode_transition_count
    assert restored.mode_transition_rate == original.mode_transition_rate
    assert np.array_equal(
        restored.all_temperature_samples, original.all_temperature_samples
    )
    assert np.array_equal(
        restored.temperature_mode_labels, original.temperature_mode_labels
    )
    assert np.array_equal(
        restored.all_temperature_objectives, original.all_temperature_objectives
    )
    assert np.array_equal(
        restored.uniform_reference_objectives,
        original.uniform_reference_objectives,
    )


def test_joint_inversion_configuration_round_trips_metadata(tmp_path):
    experiment = Experiment(
        "joint-archive",
        2500.0,
        550.0,
        {"A": 1.0},
        ("A",),
        phase_compositions={"A": {"x": 0.5}},
    )
    problem = JointInversionProblem(
        object(),
        experiment,
        [ThermodynamicParameter("dH_A", "A")],
        Parameter("pressure", 3000.0, 500.0, 500.0, 8000.0, "bar"),
        Parameter("temperature", 600.0, 50.0, 300.0, 900.0, "degC"),
    )
    write_joint_inversion_configuration(problem, tmp_path)
    restored = read_joint_inversion_configuration(tmp_path)

    assert restored["configuration"]["parameter_order"] == [
        "dH_A",
        "pressure",
        "temperature",
    ]
    assert restored["configuration"]["parameters"][1]["unit"] == "bar"


def test_multi_experiment_joint_configuration_round_trips_all_experiments(tmp_path):
    experiments = (
        Experiment(
            "joint-a",
            3000.0,
            600.0,
            {"A": 1.0},
            ("A",),
            phase_fractions={"A": 1.0},
        ),
        Experiment(
            "joint-b",
            4000.0,
            650.0,
            {"A": 1.0},
            ("A",),
            phase_fractions={"A": 1.0},
        ),
    )
    problem = JointInversionProblem(
        object(),
        experiments,
        [ThermodynamicParameter("dH_A", "A")],
        [
            Parameter("pressure_a", 3000.0, 500.0, 500.0, 8000.0, "bar"),
            Parameter("pressure_b", 4000.0, 500.0, 500.0, 8000.0, "bar"),
        ],
        [
            Parameter("temperature_a", 600.0, 50.0, 300.0, 900.0, "degC"),
            Parameter("temperature_b", 650.0, 50.0, 300.0, 900.0, "degC"),
        ],
    )
    write_joint_inversion_configuration(problem, tmp_path)
    restored = read_joint_inversion_configuration(tmp_path)

    assert restored["configuration"]["experiment_count"] == 2
    assert restored["configuration"]["experiment_ids"] == ["joint-a", "joint-b"]
    assert [item["id"] for item in restored["experiments"]] == ["joint-a", "joint-b"]


def test_mcmc_result_round_trips_for_a_single_parameter(tmp_path):
    rng = np.random.default_rng(1)
    samples = rng.normal(size=(30, 1))
    original = MCMCResult(
        parameter_names=("dH_A",),
        samples=samples,
        log_posterior=-np.sum(samples**2, axis=1),
        acceptance_rate=0.3,
        posterior_mean={"dH_A": float(np.mean(samples))},
        posterior_std={"dH_A": float(np.std(samples))},
        credible_interval_95={"dH_A": (-1.0, 1.0)},
    )
    write_mcmc_result(original, tmp_path)
    restored = read_mcmc_result(tmp_path)
    assert restored.samples.shape == (30, 1)
    assert np.allclose(restored.samples, original.samples, atol=1.0e-6)


def test_mcmc_archive_rejects_invalid_parameter_units(tmp_path):
    (tmp_path / "summary.json").write_text(
        json.dumps({"parameter_names": ["dH_A"], "parameter_units": [""]}),
        encoding="utf-8",
    )
    np.savetxt(
        tmp_path / "samples.csv", [[0.0]], delimiter=",", header="dH_A", comments=""
    )
    np.savetxt(tmp_path / "log_posterior.csv", [0.0], delimiter=",")
    try:
        read_mcmc_result(tmp_path)
    except ValueError:
        pass
    else:
        raise AssertionError("Expected ValueError for invalid MCMC units")


def test_multi_chain_mcmc_result_round_trips_through_disk(tmp_path):
    chains = tuple(
        MCMCResult(
            ("dH_A",),
            np.asarray([[float(seed)], [float(seed + 1)]]),
            np.asarray([-1.0, -2.0]),
            0.4,
            {"dH_A": float(seed) + 0.5},
            {"dH_A": 0.5},
            {"dH_A": (float(seed), float(seed + 1))},
            proposal_covariance_j_mol2=np.asarray([[0.25]]),
            seed=seed,
            burn_in=10,
        )
        for seed in (5, 6)
    )
    original = MultiChainMCMCResult(
        ("dH_A",),
        chains,
        {"dH_A": 1.02},
        {"dH_A": 3.5},
        ((-2.0,), (2.0,)),
        (5, 6),
        ("J/mol",),
    )

    write_multi_chain_mcmc_result(original, tmp_path)
    restored = read_multi_chain_mcmc_result(tmp_path)

    assert restored.r_hat == original.r_hat
    assert restored.effective_sample_size == original.effective_sample_size
    assert restored.initial_points == original.initial_points
    assert restored.seeds == original.seeds
    assert restored.parameter_units == original.parameter_units
    assert len(restored.chains) == 2
    assert restored.chains[1].seed == 6
    assert np.allclose(restored.chains[0].proposal_covariance_j_mol2, [[0.25]])


def test_posterior_predictive_results_round_trip_through_disk(tmp_path):
    continuous = PosteriorPredictiveResult(
        ("phase_fraction:A",),
        np.asarray([0.5]),
        np.asarray([0.51]),
        np.asarray([0.1]),
        {"phase_fraction:A": (0.3, 0.7)},
        {"phase_fraction:A": True},
        True,
        12,
        "random",
    )
    write_posterior_predictive_result(continuous, tmp_path / "continuous")
    restored_continuous = read_posterior_predictive_result(tmp_path / "continuous")
    assert np.allclose(restored_continuous.predictive_mean, continuous.predictive_mean)
    assert restored_continuous.within_95 == continuous.within_95
    assert restored_continuous.include_observation_uncertainty
    assert restored_continuous.sample_selection == "random"
    assert restored_continuous.seed == 12

    bracket = BracketPosteriorPredictiveResult(
        ("bracket-1",), 4, 0.75, 0.25, (0.0, 1.0), (0.0, 1.0, 0.0, 0.0)
    )
    write_bracket_posterior_predictive_result(bracket, tmp_path / "bracket")
    restored_bracket = read_bracket_posterior_predictive_result(tmp_path / "bracket")
    assert restored_bracket.bracket_ids == bracket.bracket_ids
    assert restored_bracket.misfits == bracket.misfits


def test_observation_bias_profile_round_trips_through_disk(tmp_path):
    original = ObservationBiasProfileResult(
        {"lab-a": 0.05, "lab-b": -0.02},
        1.25,
        {"A": 100.0},
        0.5,
        "laboratory",
        "composition:A:x",
        ({"group": "lab-a", "bias": 0.05, "objective": 1.25},),
        (-0.5, 0.5),
        0.0,
        0.2,
    )

    write_observation_bias_profile_result(original, tmp_path)
    restored = read_observation_bias_profile_result(tmp_path)

    assert restored.biases == original.biases
    assert restored.metadata_key == "laboratory"
    assert restored.observation_key == "composition:A:x"
    assert restored.records == original.records
    assert restored.bias_bounds == (-0.5, 0.5)
    assert restored.bias_prior_sigma == 0.2


def test_grouped_observation_bias_report_round_trips_through_disk(tmp_path):
    original = GroupedObservationBiasReport(
        metadata_key="laboratory",
        observation_keys=("composition:A:x", "composition:A:y"),
        group_biases={
            "lab-a": {"composition:A:x": 0.05, "composition:A:y": -0.02},
            "lab-b": {"composition:A:x": -0.01, "composition:A:y": 0.03},
        },
        bias_prior_mean=0.0,
        bias_prior_sigma=0.2,
        bias_bounds=(-0.5, 0.5),
        objective=1.3,
        bias_prior_penalty=0.4,
    )

    write_grouped_observation_bias_report(original, tmp_path)
    restored = read_grouped_observation_bias_report(tmp_path)

    assert restored.metadata_key == original.metadata_key
    assert restored.observation_keys == original.observation_keys
    assert restored.group_biases == original.group_biases
    assert restored.bias_prior_sigma == 0.2
    assert restored.bias_bounds == (-0.5, 0.5)
    assert restored.objective == 1.3


def test_grouped_observation_bias_summary_round_trips_through_disk(tmp_path):
    original = GroupedObservationBiasSummary(
        metadata_key="laboratory",
        observation_keys=("composition:A:x", "composition:A:y"),
        group_biases={
            "lab-a": {"composition:A:x": 0.05, "composition:A:y": -0.02},
            "lab-b": {"composition:A:x": -0.01, "composition:A:y": 0.03},
        },
        bias_prior_mean=0.0,
        bias_prior_sigma=0.2,
        bias_bounds=(-0.5, 0.5),
        objective=1.3,
        bias_prior_penalty=0.4,
    )

    write_grouped_observation_bias_summary(original, tmp_path)
    restored = read_grouped_observation_bias_summary(tmp_path)

    assert restored.metadata_key == original.metadata_key
    assert restored.observation_keys == original.observation_keys
    assert restored.group_biases == original.group_biases
    assert restored.bias_prior_sigma == original.bias_prior_sigma
    assert restored.bias_bounds == original.bias_bounds


def test_loocv_result_round_trips_through_disk(tmp_path):
    original = LOOCVResult(
        folds=(
            LOOCVFold(
                "exp-1",
                {"A": 100.0, "B": -50.0},
                0.5,
                ("held-1",),
                (0.5,),
                (-2.3,),
            ),
            LOOCVFold("exp-2", {"A": 110.0, "B": -40.0}, 0.7),
        ),
        mean_misfit=0.6,
        group_metadata_key="study",
        mean_log_predictive_density=-2.3,
    )
    write_loocv_result(original, tmp_path)
    restored = read_loocv_result(tmp_path)

    assert restored.mean_misfit == original.mean_misfit
    assert len(restored.folds) == 2
    assert restored.folds[0].held_out_experiment_id == "exp-1"
    assert restored.folds[0].corrections_j_mol == {"A": 100.0, "B": -50.0}
    assert restored.folds[1].misfit == 0.7
    assert restored.folds[0].held_out_experiment_ids == ("held-1",)
    assert restored.folds[0].held_out_misfits == (0.5,)
    assert restored.folds[0].held_out_log_predictive_densities == (-2.3,)
    assert restored.group_metadata_key == "study"
    assert restored.mean_log_predictive_density == -2.3


def test_loocv_archive_rejects_mismatched_group_detail_lengths(tmp_path):
    (tmp_path / "folds.csv").write_text(
        "held_out_experiment_id,held_out_experiment_ids,held_out_misfits,misfit,corrections_json\n"
        'study=one,"["one-1", "one-2"]","[1.0]",1.0,"{}"\n',
        encoding="utf-8",
    )
    (tmp_path / "summary.json").write_text(
        json.dumps({"mean_misfit": 1.0, "group_metadata_key": "study"}),
        encoding="utf-8",
    )
    try:
        read_loocv_result(tmp_path)
    except ValueError:
        pass
    else:
        raise AssertionError("Expected ValueError for mismatched LOOCV detail lengths")


def test_phase_assemblage_holdout_archive_preserves_provenance(tmp_path):
    original = LOOCVResult(
        folds=(
            LOOCVFold(
                "__thermoinvert_phase_assemblage__=A",
                {"A": 10.0},
                0.4,
                ("exp-a1", "exp-a2"),
                (0.3, 0.5),
            ),
            LOOCVFold(
                "__thermoinvert_phase_assemblage__=B",
                {"A": 12.0},
                0.8,
                ("exp-b1",),
                (0.8,),
            ),
        ),
        mean_misfit=0.6,
        group_metadata_key="__thermoinvert_phase_assemblage__",
    )
    write_loocv_result(original, tmp_path)
    restored = read_loocv_result(tmp_path)

    assert restored.group_metadata_key == original.group_metadata_key
    assert restored.folds[0].held_out_experiment_ids == ("exp-a1", "exp-a2")
    assert restored.folds[0].held_out_misfits == (0.3, 0.5)


def test_out_of_domain_holdout_archive_preserves_provenance(tmp_path):
    original = LOOCVResult(
        folds=(
            LOOCVFold(
                "out_of_domain_pt",
                {"A": 5.0},
                1.2,
                ("outside-pt",),
                (1.2,),
            ),
        ),
        mean_misfit=1.2,
        group_metadata_key="out_of_domain_pt",
    )
    write_loocv_result(original, tmp_path)
    restored = read_loocv_result(tmp_path)

    assert restored.group_metadata_key == "out_of_domain_pt"
    assert restored.folds[0].held_out_experiment_ids == ("outside-pt",)
    assert restored.folds[0].held_out_misfits == (1.2,)


def test_grid_posterior_round_trips_for_one_and_two_dimensions(tmp_path):
    axis = np.linspace(-1.0, 1.0, 11)
    density_1d = np.exp(-(axis**2))
    density_1d /= np.sum(density_1d) * (axis[1] - axis[0])
    original_1d = GridPosteriorResult(
        parameter_names=("p0",),
        axes=(axis,),
        log_posterior=-(axis**2),
        posterior_density=density_1d,
        marginal_mean={"p0": 0.0},
        marginal_std={"p0": 0.5},
    )
    write_grid_posterior(original_1d, tmp_path / "grid1d")
    restored_1d = read_grid_posterior(tmp_path / "grid1d")
    assert restored_1d.posterior_density.shape == (11,)
    assert np.allclose(restored_1d.posterior_density, density_1d, atol=1.0e-9)
    assert restored_1d.marginal_mean == {"p0": 0.0}

    density_2d = np.outer(density_1d, density_1d)
    original_2d = GridPosteriorResult(
        parameter_names=("p0", "p1"),
        axes=(axis, axis),
        log_posterior=np.zeros((11, 11)),
        posterior_density=density_2d,
        marginal_mean={"p0": 0.0, "p1": 0.0},
        marginal_std={"p0": 0.5, "p1": 0.5},
    )
    write_grid_posterior(original_2d, tmp_path / "grid2d")
    restored_2d = read_grid_posterior(tmp_path / "grid2d")
    assert restored_2d.posterior_density.shape == (11, 11)
    assert np.allclose(restored_2d.posterior_density, density_2d, atol=1.0e-9)


def test_posterior_resolution_result_round_trips_through_disk(tmp_path):
    original = PosteriorResolutionResult(
        ("p0", "p1"),
        3,
        np.asarray([[0.7, 0.1], [0.2, 0.6]]),
        np.asarray([[0.4, 0.0], [0.1, 0.3]]),
        np.asarray([[0.9, 0.3], [0.4, 0.8]]),
        np.asarray([0.7, 0.6]),
        np.asarray([0.4, 0.3]),
        np.asarray([0.9, 0.8]),
        (0, 1, 1),
        (0, 4, 9),
        "evenly_spaced_balanced_chains",
        ("exp:fraction:A", "exp:fraction:B", "exp:composition:A:X"),
        (0.1, 0.9),
        0.002,
        "explicit",
        "problem_assembled",
        np.asarray([[0.04, 0.01, 0.0], [0.01, 0.09, 0.0], [0.0, 0.0, 0.01]]),
        np.asarray([[0.08, 0.02, 0.0], [0.02, 0.18, 0.0], [0.0, 0.0, 0.04]]),
        np.asarray([[1.0, 0.2], [0.2, 2.0]]),
    )
    write_posterior_resolution_result(original, tmp_path)
    restored = read_posterior_resolution_result(tmp_path)

    assert restored.parameter_names == original.parameter_names
    assert restored.sample_count == original.sample_count
    assert restored.selected_chain_indices == original.selected_chain_indices
    assert restored.selected_draw_indices == original.selected_draw_indices
    assert restored.selection_strategy == original.selection_strategy
    assert restored.observable_names == original.observable_names
    assert restored.credible_interval == original.credible_interval
    assert restored.step_fraction == original.step_fraction
    assert (
        restored.observation_covariance_source == original.observation_covariance_source
    )
    assert (
        restored.theoretical_covariance_source == original.theoretical_covariance_source
    )
    assert restored.observation_covariance is not None
    assert original.observation_covariance is not None
    assert restored.theoretical_covariance is not None
    assert original.theoretical_covariance is not None
    assert np.allclose(restored.observation_covariance, original.observation_covariance)
    assert np.allclose(restored.theoretical_covariance, original.theoretical_covariance)
    assert restored.prior_covariance is not None
    assert original.prior_covariance is not None
    assert np.allclose(restored.prior_covariance, original.prior_covariance)
    for field in (
        "mean_matrix",
        "lower_matrix",
        "upper_matrix",
        "mean_diagonal",
        "lower_diagonal",
        "upper_diagonal",
    ):
        assert np.allclose(getattr(restored, field), getattr(original, field))


def test_posterior_resolution_archive_rejects_nonfinite_matrix(tmp_path):
    archive_path = tmp_path / "posterior_resolution.json"
    archive_path.write_text(
        json.dumps(
            {
                "parameter_names": ["p0"],
                "sample_count": 2,
                "mean_matrix": [[float("nan")]],
                "lower_matrix": [[0.0]],
                "upper_matrix": [[1.0]],
                "mean_diagonal": [0.5],
                "lower_diagonal": [0.0],
                "upper_diagonal": [1.0],
            }
        ),
        encoding="utf-8",
    )

    try:
        read_posterior_resolution_result(tmp_path)
    except ValueError as error:
        assert "mean_matrix" in str(error)
    else:
        raise AssertionError("Expected ValueError for nonfinite resolution matrix")


def test_posterior_resolution_archive_rejects_mismatched_selection_provenance(
    tmp_path,
):
    archive_path = tmp_path / "posterior_resolution.json"
    archive_path.write_text(
        json.dumps(
            {
                "parameter_names": ["p0"],
                "sample_count": 2,
                "mean_matrix": [[0.5]],
                "lower_matrix": [[0.2]],
                "upper_matrix": [[0.8]],
                "mean_diagonal": [0.5],
                "lower_diagonal": [0.2],
                "upper_diagonal": [0.8],
                "selected_chain_indices": [0],
                "selected_draw_indices": [1],
            }
        ),
        encoding="utf-8",
    )

    try:
        read_posterior_resolution_result(tmp_path)
    except ValueError as error:
        assert "sample_count" in str(error)
    else:
        raise AssertionError("Expected ValueError for incomplete selection provenance")


def test_posterior_resolution_archive_rejects_invalid_credible_interval(tmp_path):
    archive_path = tmp_path / "posterior_resolution.json"
    archive_path.write_text(
        json.dumps(
            {
                "parameter_names": ["p0"],
                "sample_count": 1,
                "mean_matrix": [[0.5]],
                "lower_matrix": [[0.2]],
                "upper_matrix": [[0.8]],
                "mean_diagonal": [0.5],
                "lower_diagonal": [0.2],
                "upper_diagonal": [0.8],
                "credible_interval": [0.95, 0.05],
            }
        ),
        encoding="utf-8",
    )

    try:
        read_posterior_resolution_result(tmp_path)
    except ValueError as error:
        assert "credible_interval" in str(error)
    else:
        raise AssertionError("Expected ValueError for reversed credible interval")


def test_posterior_resolution_archive_rejects_wrong_prior_covariance_dimension(
    tmp_path,
):
    archive_path = tmp_path / "posterior_resolution.json"
    archive_path.write_text(
        json.dumps(
            {
                "parameter_names": ["p0"],
                "sample_count": 1,
                "mean_matrix": [[0.5]],
                "lower_matrix": [[0.2]],
                "upper_matrix": [[0.8]],
                "mean_diagonal": [0.5],
                "lower_diagonal": [0.2],
                "upper_diagonal": [0.8],
                "prior_covariance": [[1.0, 0.0], [0.0, 1.0]],
            }
        ),
        encoding="utf-8",
    )

    try:
        read_posterior_resolution_result(tmp_path)
    except ValueError as error:
        assert "parameter_names length" in str(error)
    else:
        raise AssertionError("Expected ValueError for wrong prior covariance dimension")


def test_inversion_configuration_archives_correlated_prior_and_observations(tmp_path):
    experiment = Experiment(
        "exp",
        2000.0,
        500.0,
        {"A": 1.0},
        ("A",),
        phase_fractions={"A": 1.0},
        observation_covariance=[[0.01]],
        observation_covariance_keys=("fraction:A",),
        theoretical_covariance=[[0.02]],
        theoretical_covariance_keys=("fraction:A",),
    )
    parameter = ThermodynamicParameter("dH_A", "A", prior_sigma_j_mol=2.0)
    problem = InversionProblem(
        object(),
        [experiment],
        [parameter],
        prior_covariance_j_mol2=[[4.0]],
        theoretical_covariance_scale=2.0,
    )

    write_inversion_configuration(problem, tmp_path)

    configuration = json.loads(
        (tmp_path / "inversion_configuration.json").read_text(encoding="utf-8")
    )
    assert configuration["prior_covariance"] == [[4.0]]
    experiments = json.loads(
        (tmp_path / "experiments.json").read_text(encoding="utf-8")
    )
    assert configuration["parameter_order"] == ["A"]
    assert configuration["prior_covariance_j_mol2"] == [[4.0]]
    parameter_record = configuration["parameters"][0]
    assert parameter_record["prior_mean"] == parameter.prior_mean
    assert parameter_record["prior_sigma"] == parameter.prior_sigma
    assert parameter_record["prior_sigma_j_mol"] == parameter.prior_sigma
    assert experiments[0]["theoretical_covariance"] == [[0.02]]
    assert configuration["theoretical_covariance_scale"] == 2.0
    restored = read_inversion_configuration(tmp_path)
    assert restored["configuration"]["parameter_order"] == ["A"]
    assert restored["experiments"][0]["id"] == "exp"


def test_inversion_configuration_archives_dirichlet_concentration_parameter(tmp_path):
    experiment = Experiment(
        "composition",
        2000.0,
        500.0,
        {"Zn": 1.0},
        ("Znc",),
        phase_compositions={"Znc": {"Zn": 0.8, "Mg": 0.2, "Mn": 0.0}},
        phase_composition_concentration={"Znc": 25.0},
        phase_composition_detection_limit={"Znc": 0.02},
    )
    parameter = DirichletConcentrationParameter(
        "Znc", prior_mean=30.0, prior_sigma=8.0, lower_bound=0.5, upper_bound=120.0
    )
    problem = InversionProblem(
        object(), [experiment], [parameter], misfit_form="dirichlet_nll"
    )

    write_inversion_configuration(problem, tmp_path)
    restored = read_inversion_configuration(tmp_path)
    record = restored["configuration"]["parameters"][0]

    assert record["parameter_kind"] == "dirichlet_concentration"
    assert record["phase"] == "Znc"
    assert record["key"] == "dirichlet_concentration:Znc"
    assert record["prior_mean"] == 30.0


def test_inversion_configuration_reader_rejects_invalid_prior(tmp_path):
    (tmp_path / "inversion_configuration.json").write_text(
        json.dumps(
            {
                "parameter_order": ["A", "B"],
                "parameters": [{"key": "A"}, {"key": "B"}],
                "prior_covariance_j_mol2": [[1.0, 2.0], [2.0, 1.0]],
            }
        ),
        encoding="utf-8",
    )
    (tmp_path / "experiments.json").write_text("[]", encoding="utf-8")
    (tmp_path / "brackets.json").write_text("[]", encoding="utf-8")
    try:
        read_inversion_configuration(tmp_path)
    except ValueError:
        pass
    else:
        raise AssertionError("Expected ValueError for indefinite archived prior")


def test_inversion_configuration_reader_rejects_duplicate_parameter_keys(tmp_path):
    (tmp_path / "inversion_configuration.json").write_text(
        json.dumps(
            {
                "parameter_order": ["A", "A"],
                "parameters": [{"key": "A"}, {"key": "A"}],
                "prior_covariance_j_mol2": [[1.0, 0.0], [0.0, 1.0]],
            }
        ),
        encoding="utf-8",
    )
    (tmp_path / "experiments.json").write_text("[]", encoding="utf-8")
    (tmp_path / "brackets.json").write_text("[]", encoding="utf-8")
    try:
        read_inversion_configuration(tmp_path)
    except ValueError:
        pass
    else:
        raise AssertionError(
            "Expected ValueError for duplicate archived parameter keys"
        )


def test_inversion_configuration_reader_rejects_missing_required_covariance(tmp_path):
    (tmp_path / "inversion_configuration.json").write_text(
        json.dumps(
            {
                "parameter_order": ["A"],
                "parameters": [{"key": "A"}],
            }
        ),
        encoding="utf-8",
    )
    (tmp_path / "experiments.json").write_text("[]", encoding="utf-8")
    (tmp_path / "brackets.json").write_text("[]", encoding="utf-8")
    try:
        read_inversion_configuration(tmp_path)
    except ValueError:
        pass
    else:
        raise AssertionError("Expected ValueError for missing archived covariance")


def test_configuration_readers_reject_wrong_json_container_types(tmp_path):
    (tmp_path / "inversion_configuration.json").write_text("[]", encoding="utf-8")
    try:
        read_inversion_configuration(tmp_path)
    except ValueError:
        pass
    else:
        raise AssertionError(
            "Expected ValueError for non-object inversion configuration"
        )

    (tmp_path / "pt_inversion_configuration.json").write_text("[]", encoding="utf-8")
    try:
        read_pt_inversion_configuration(tmp_path)
    except ValueError:
        pass
    else:
        raise AssertionError("Expected ValueError for non-object P-T configuration")


def test_pt_configuration_reader_rejects_incomplete_parameter_metadata(tmp_path):
    (tmp_path / "pt_inversion_configuration.json").write_text(
        json.dumps(
            {
                "parameter_order": ["pressure_bar", "temperature_c"],
                "parameters": [{"name": "pressure_bar"}, {"name": "wrong"}],
                "prior_covariance_j_mol2": [[1.0, 0.0], [0.0, 1.0]],
                "problem_mode": "composition",
                "bracket_ids": [],
            }
        ),
        encoding="utf-8",
    )
    (tmp_path / "experiment.json").write_text('{"id": "exp"}', encoding="utf-8")
    (tmp_path / "brackets.json").write_text("[]", encoding="utf-8")
    try:
        read_pt_inversion_configuration(tmp_path)
    except ValueError:
        pass
    else:
        raise AssertionError(
            "Expected ValueError for incomplete P-T parameter metadata"
        )


def test_pt_configuration_reader_rejects_invalid_corrections(tmp_path):
    (tmp_path / "pt_inversion_configuration.json").write_text(
        json.dumps(
            {
                "parameter_order": ["pressure_bar", "temperature_c"],
                "parameters": [
                    {
                        "name": "pressure_bar",
                        "prior_mean_j_mol": 1.0,
                        "prior_sigma_j_mol": 1.0,
                    },
                    {
                        "name": "temperature_c",
                        "prior_mean_j_mol": 1.0,
                        "prior_sigma_j_mol": 1.0,
                    },
                ],
                "prior_covariance_j_mol2": [[1.0, 0.0], [0.0, 1.0]],
                "missing_phase_penalty": 1.0,
                "extra_phase_penalty": 1.0,
                "misfit_form": "weighted_lsq",
                "corrections_j_mol": {"A": "bad"},
                "problem_mode": "composition",
                "bracket_ids": [],
            }
        ),
        encoding="utf-8",
    )
    (tmp_path / "experiment.json").write_text('{"id": "exp"}', encoding="utf-8")
    (tmp_path / "brackets.json").write_text("[]", encoding="utf-8")
    try:
        read_pt_inversion_configuration(tmp_path)
    except ValueError:
        pass
    else:
        raise AssertionError("Expected ValueError for invalid P-T corrections")


def test_pt_configuration_reader_rejects_bracket_id_mismatch(tmp_path):
    (tmp_path / "pt_inversion_configuration.json").write_text(
        json.dumps(
            {
                "parameter_order": ["pressure_bar", "temperature_c"],
                "parameters": [
                    {
                        "name": "pressure_bar",
                        "prior_mean_j_mol": 1.0,
                        "prior_sigma_j_mol": 1.0,
                    },
                    {
                        "name": "temperature_c",
                        "prior_mean_j_mol": 1.0,
                        "prior_sigma_j_mol": 1.0,
                    },
                ],
                "prior_covariance_j_mol2": [[1.0, 0.0], [0.0, 1.0]],
                "missing_phase_penalty": 1.0,
                "extra_phase_penalty": 1.0,
                "misfit_form": "weighted_lsq",
                "corrections_j_mol": {},
                "problem_mode": "bracket",
                "bracket_ids": ["declared"],
            }
        ),
        encoding="utf-8",
    )
    (tmp_path / "experiment.json").write_text('{"id": "exp"}', encoding="utf-8")
    (tmp_path / "brackets.json").write_text(
        json.dumps([{"id": "stored"}]), encoding="utf-8"
    )
    try:
        read_pt_inversion_configuration(tmp_path)
    except ValueError:
        pass
    else:
        raise AssertionError("Expected ValueError for bracket ID mismatch")


def test_inversion_configuration_reader_rejects_duplicate_record_ids(tmp_path):
    (tmp_path / "inversion_configuration.json").write_text(
        json.dumps(
            {
                "parameter_order": ["A"],
                "parameters": [{"key": "A"}],
                "prior_covariance_j_mol2": [[1.0]],
            }
        ),
        encoding="utf-8",
    )
    (tmp_path / "experiments.json").write_text(
        json.dumps([{"id": "same"}, {"id": "same"}]), encoding="utf-8"
    )
    (tmp_path / "brackets.json").write_text("[]", encoding="utf-8")
    try:
        read_inversion_configuration(tmp_path)
    except ValueError:
        pass
    else:
        raise AssertionError("Expected ValueError for duplicate experiment ids")


def test_inversion_configuration_reader_rejects_malformed_experiment_schema(tmp_path):
    (tmp_path / "inversion_configuration.json").write_text(
        json.dumps(
            {
                "parameter_order": ["A"],
                "parameters": [{"key": "A"}],
                "prior_covariance_j_mol2": [[1.0]],
                "missing_phase_penalty": 1.0,
                "extra_phase_penalty": 1.0,
                "condition_quadrature_order": 3,
                "objective_criterion": "map",
                "misfit_form": "weighted_lsq",
            }
        ),
        encoding="utf-8",
    )
    (tmp_path / "experiments.json").write_text(
        json.dumps([{"id": "exp", "pressure_bar": 1.0}]), encoding="utf-8"
    )
    (tmp_path / "brackets.json").write_text("[]", encoding="utf-8")
    try:
        read_inversion_configuration(tmp_path)
    except ValueError:
        pass
    else:
        raise AssertionError("Expected ValueError for malformed experiment schema")


def test_inversion_configuration_reader_rejects_nonstring_parameter_keys(tmp_path):
    (tmp_path / "inversion_configuration.json").write_text(
        json.dumps(
            {
                "parameter_order": [{"key": "A"}],
                "parameters": [{"key": {"key": "A"}}],
                "prior_covariance_j_mol2": [[1.0]],
            }
        ),
        encoding="utf-8",
    )
    (tmp_path / "experiments.json").write_text("[]", encoding="utf-8")
    (tmp_path / "brackets.json").write_text("[]", encoding="utf-8")
    try:
        read_inversion_configuration(tmp_path)
    except ValueError:
        pass
    else:
        raise AssertionError("Expected ValueError for non-string parameter key")


def test_inversion_configuration_reader_rejects_invalid_objective_settings(tmp_path):
    (tmp_path / "inversion_configuration.json").write_text(
        json.dumps(
            {
                "parameter_order": ["A"],
                "parameters": [{"key": "A"}],
                "prior_covariance_j_mol2": [[1.0]],
                "missing_phase_penalty": -1.0,
                "extra_phase_penalty": 1.0,
                "condition_quadrature_order": 3,
                "objective_criterion": "map",
                "misfit_form": "weighted_lsq",
            }
        ),
        encoding="utf-8",
    )
    (tmp_path / "experiments.json").write_text("[]", encoding="utf-8")
    (tmp_path / "brackets.json").write_text("[]", encoding="utf-8")
    try:
        read_inversion_configuration(tmp_path)
    except ValueError:
        pass
    else:
        raise AssertionError("Expected ValueError for invalid objective settings")


def test_inversion_configuration_reader_rejects_invalid_parameter_metadata(tmp_path):
    (tmp_path / "inversion_configuration.json").write_text(
        json.dumps(
            {
                "parameter_order": ["A"],
                "parameters": [
                    {
                        "key": "A",
                        "name": "dH_A",
                        "phase": "A",
                        "field": "H0",
                        "parameter_type": "dHf",
                        "prior_mean_j_mol": 0.0,
                        "prior_sigma_j_mol": -1.0,
                        "lower_bound_j_mol": -1.0,
                        "upper_bound_j_mol": 1.0,
                    }
                ],
                "prior_covariance_j_mol2": [[1.0]],
                "missing_phase_penalty": 1.0,
                "extra_phase_penalty": 1.0,
                "condition_quadrature_order": 3,
                "objective_criterion": "map",
                "misfit_form": "weighted_lsq",
            }
        ),
        encoding="utf-8",
    )
    (tmp_path / "experiments.json").write_text("[]", encoding="utf-8")
    (tmp_path / "brackets.json").write_text("[]", encoding="utf-8")
    try:
        read_inversion_configuration(tmp_path)
    except ValueError:
        pass
    else:
        raise AssertionError("Expected ValueError for invalid parameter metadata")


def test_inversion_configuration_reader_rejects_inconsistent_parameter_key(tmp_path):
    (tmp_path / "inversion_configuration.json").write_text(
        json.dumps(
            {
                "parameter_order": ["wrong"],
                "parameters": [
                    {
                        "key": "wrong",
                        "name": "dH_A",
                        "phase": "A",
                        "field": "H0",
                        "parameter_type": "dHf",
                        "prior_mean_j_mol": 0.0,
                        "prior_sigma_j_mol": 1.0,
                    }
                ],
                "prior_covariance_j_mol2": [[1.0]],
                "missing_phase_penalty": 1.0,
                "extra_phase_penalty": 1.0,
                "condition_quadrature_order": 3,
                "objective_criterion": "map",
                "misfit_form": "weighted_lsq",
            }
        ),
        encoding="utf-8",
    )
    (tmp_path / "experiments.json").write_text("[]", encoding="utf-8")
    (tmp_path / "brackets.json").write_text("[]", encoding="utf-8")
    try:
        read_inversion_configuration(tmp_path)
    except ValueError:
        pass
    else:
        raise AssertionError("Expected ValueError for inconsistent parameter key")


def test_inversion_configuration_reader_rejects_unsupported_field(tmp_path):
    (tmp_path / "inversion_configuration.json").write_text(
        json.dumps(
            {
                "parameter_order": ["A:bad"],
                "parameters": [
                    {
                        "key": "A:bad",
                        "name": "dH_A",
                        "phase": "A",
                        "field": "bad",
                        "parameter_type": "dHf",
                        "prior_mean_j_mol": 0.0,
                        "prior_sigma_j_mol": 1.0,
                    }
                ],
                "prior_covariance_j_mol2": [[1.0]],
                "missing_phase_penalty": 1.0,
                "extra_phase_penalty": 1.0,
                "condition_quadrature_order": 3,
                "objective_criterion": "map",
                "misfit_form": "weighted_lsq",
            }
        ),
        encoding="utf-8",
    )
    (tmp_path / "experiments.json").write_text("[]", encoding="utf-8")
    (tmp_path / "brackets.json").write_text("[]", encoding="utf-8")
    try:
        read_inversion_configuration(tmp_path)
    except ValueError:
        pass
    else:
        raise AssertionError("Expected ValueError for unsupported archived field")
