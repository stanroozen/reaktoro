"""Tests for the shrinkage-based + structural identifiability diagnostic."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parents[1]))

from thermoinvert.data import Experiment, ThermodynamicParameter
from thermoinvert.forward_model import EquilibriumResult
from thermoinvert.identifiability import (
    diagnose_grid_modes,
    diagnose_map_globality,
    diagnose_posterior_shape,
    gaussian_kl_divergence,
    gaussian_multivariate_kl_divergence,
    identifiability_summary,
    information_gain_summary,
    grid_information_gain,
    multivariate_information_gain,
    posterior_correlation_matrix,
)
from thermoinvert.inversion import InversionProblem, MAPResult
from thermoinvert.posterior import (
    GridPosteriorResult,
    MCMCResult,
    MultiChainMCMCResult,
    LaplaceApproximation,
    laplace_covariance,
)
from thermoinvert.sensitivity import (
    identify_sensitivity,
    local_resolution_matrix,
    posterior_averaged_resolution,
    posterior_resolution_from_mcmc,
)


def test_identifiability_summary_flags_data_constrained_parameter():
    # Data term is sharp (small sigma_data) relative to a wide prior, so the
    # data should dominate: shrinkage should be close to 1.
    sigma_data = 50.0
    target = 900.0
    prior_sigma = 5000.0

    class QuadraticForwardModel:
        def predict(self, experiment, corrections):
            value = corrections.get("A", 0.0)
            return EquilibriumResult(
                True, ("A",), {"A": 1.0 + (value - target) / sigma_data}, {}, {}
            )

    experiment = Experiment(
        id="constrained",
        pressure_bar=2000.0,
        temperature_c=500.0,
        bulk_composition={"A": 1.0},
        observed_phases=("A",),
        phase_fractions={"A": 1.0},
        phase_fraction_sigma={"A": 1.0},
    )
    parameters = [ThermodynamicParameter("dH_A", "A", prior_sigma_j_mol=prior_sigma)]
    problem = InversionProblem(QuadraticForwardModel(), [experiment], parameters)

    all_result = laplace_covariance(problem, [target], uncertainty_source="all")
    prior_result = laplace_covariance(problem, [0.0], uncertainty_source="prior")

    summary = identifiability_summary(problem, all_result, prior_result)
    entry = summary["dH_A"]
    assert entry.status == "data_constrained"
    assert entry.shrinkage > 0.9


def test_identifiability_summary_flags_prior_dominated_parameter():
    # Data term is very weak (huge sigma_data) relative to the prior, so the
    # data should barely shrink the prior: shrinkage should be close to 0.
    sigma_data = 1.0e8
    target = 900.0
    prior_sigma = 5000.0

    class WeakDataForwardModel:
        def predict(self, experiment, corrections):
            value = corrections.get("A", 0.0)
            return EquilibriumResult(
                True, ("A",), {"A": 1.0 + (value - target) / sigma_data}, {}, {}
            )

    experiment = Experiment(
        id="unconstrained",
        pressure_bar=2000.0,
        temperature_c=500.0,
        bulk_composition={"A": 1.0},
        observed_phases=("A",),
        phase_fractions={"A": 1.0},
        phase_fraction_sigma={"A": 1.0},
    )
    parameters = [ThermodynamicParameter("dH_A", "A", prior_sigma_j_mol=prior_sigma)]
    problem = InversionProblem(WeakDataForwardModel(), [experiment], parameters)

    all_result = laplace_covariance(problem, [0.0], uncertainty_source="all")
    prior_result = laplace_covariance(problem, [0.0], uncertainty_source="prior")

    summary = identifiability_summary(problem, all_result, prior_result)
    entry = summary["dH_A"]
    assert entry.status == "prior_dominated"
    assert entry.shrinkage < 0.1


def test_identifiability_summary_cross_references_structural_weak_directions():
    class FakeForwardModel:
        def predict(self, experiment, corrections):
            value_a = corrections.get("A", 0.0)
            value_b = corrections.get("B", 0.0)
            return EquilibriumResult(
                True, ("A",), {"A": 1.0 + (value_a + value_b) / 1000.0}, {}, {}
            )

    experiment = Experiment(
        id="collinear",
        pressure_bar=2000.0,
        temperature_c=500.0,
        bulk_composition={"A": 1.0},
        observed_phases=("A",),
        phase_fractions={"A": 1.0},
        phase_fraction_sigma={"A": 1.0},
    )
    parameters = [
        ThermodynamicParameter("dH_A", "A", prior_sigma_j_mol=5000.0),
        ThermodynamicParameter("dH_B", "B", prior_sigma_j_mol=5000.0),
    ]
    problem = InversionProblem(FakeForwardModel(), [experiment], parameters)

    # Jacobian columns are identical (A and B enter only as their sum), so
    # the pair is collinear -- one weak direction should show up.
    sensitivity = identify_sensitivity([[1.0, 1.0]], ["dH_A", "dH_B"])
    assert sensitivity.rank == 1

    all_result = laplace_covariance(problem, [0.0, 0.0], uncertainty_source="all")
    prior_result = laplace_covariance(problem, [0.0, 0.0], uncertainty_source="prior")
    summary = identifiability_summary(
        problem, all_result, prior_result, sensitivity=sensitivity
    )
    assert summary["dH_A"].structurally_weak


def test_gaussian_kl_divergence_is_zero_when_posterior_equals_prior():
    assert abs(gaussian_kl_divergence(0.0, 100.0, 0.0, 100.0)) < 1.0e-12

    try:
        gaussian_kl_divergence(0.0, -1.0, 0.0, 100.0)
    except ValueError:
        pass
    else:
        raise AssertionError("Expected ValueError for a non-positive variance")


def test_gaussian_multivariate_kl_rejects_nonfinite_inputs():
    try:
        gaussian_multivariate_kl_divergence([0.0], [[1.0]], [0.0], [[np.nan]])
    except ValueError:
        pass
    else:
        raise AssertionError("Expected ValueError for non-finite covariance")


def test_information_gain_summary_is_large_for_data_constrained_parameter():
    sigma_data = 50.0
    target = 900.0
    prior_sigma = 5000.0

    class QuadraticForwardModel:
        def predict(self, experiment, corrections):
            value = corrections.get("A", 0.0)
            return EquilibriumResult(
                True, ("A",), {"A": 1.0 + (value - target) / sigma_data}, {}, {}
            )

    experiment = Experiment(
        id="info-gain-constrained",
        pressure_bar=2000.0,
        temperature_c=500.0,
        bulk_composition={"A": 1.0},
        observed_phases=("A",),
        phase_fractions={"A": 1.0},
        phase_fraction_sigma={"A": 1.0},
    )
    parameters = [ThermodynamicParameter("dH_A", "A", prior_sigma_j_mol=prior_sigma)]
    problem = InversionProblem(QuadraticForwardModel(), [experiment], parameters)
    all_result = laplace_covariance(problem, [target], uncertainty_source="all")

    gains = information_gain_summary(problem, all_result)
    assert gains["dH_A"] > 2.0


def test_information_gain_summary_is_near_zero_for_prior_dominated_parameter():
    sigma_data = 1.0e8
    target = 900.0
    prior_sigma = 5000.0

    class WeakDataForwardModel:
        def predict(self, experiment, corrections):
            value = corrections.get("A", 0.0)
            return EquilibriumResult(
                True, ("A",), {"A": 1.0 + (value - target) / sigma_data}, {}, {}
            )

    experiment = Experiment(
        id="info-gain-unconstrained",
        pressure_bar=2000.0,
        temperature_c=500.0,
        bulk_composition={"A": 1.0},
        observed_phases=("A",),
        phase_fractions={"A": 1.0},
        phase_fraction_sigma={"A": 1.0},
    )
    parameters = [ThermodynamicParameter("dH_A", "A", prior_sigma_j_mol=prior_sigma)]
    problem = InversionProblem(WeakDataForwardModel(), [experiment], parameters)
    all_result = laplace_covariance(problem, [0.0], uncertainty_source="all")

    gains = information_gain_summary(problem, all_result)
    assert gains["dH_A"] < 0.01


def test_multivariate_kl_and_posterior_correlation_include_parameter_tradeoffs():
    prior_mean = np.zeros(2)
    prior_covariance = np.eye(2)
    posterior_covariance = np.asarray([[0.25, 0.2], [0.2, 0.25]])
    posterior_mean = np.asarray([0.5, -0.25])
    information = gaussian_multivariate_kl_divergence(
        posterior_mean, posterior_covariance, prior_mean, prior_covariance
    )
    independent_marginal_information = sum(
        gaussian_kl_divergence(
            posterior_mean[index],
            posterior_covariance[index, index],
            prior_mean[index],
            prior_covariance[index, index],
        )
        for index in range(2)
    )
    assert information > independent_marginal_information
    correlation = posterior_correlation_matrix(
        LaplaceApproximation(
            ("a", "b"),
            posterior_mean,
            np.eye(2),
            posterior_covariance,
            posterior_covariance / 0.25,
            True,
        )
    )
    assert np.isclose(correlation[0, 1], 0.8)


def test_multivariate_information_gain_uses_correlated_problem_prior():
    parameters = [
        ThermodynamicParameter("dH_A", "A", prior_sigma_j_mol=2.0),
        ThermodynamicParameter("dH_B", "B", prior_sigma_j_mol=3.0),
    ]
    prior_covariance = np.asarray([[4.0, 3.0], [3.0, 9.0]])
    problem = InversionProblem(
        object(), [], parameters, prior_covariance_j_mol2=prior_covariance
    )
    posterior_mean = np.asarray([0.5, -0.25])
    posterior_covariance = np.asarray([[1.0, 0.2], [0.2, 2.0]])
    result = LaplaceApproximation(
        ("dH_A", "dH_B"),
        posterior_mean,
        np.eye(2),
        posterior_covariance,
        np.linalg.inv(posterior_covariance),
        True,
    )

    assert np.isclose(
        multivariate_information_gain(problem, result),
        gaussian_multivariate_kl_divergence(
            posterior_mean, posterior_covariance, np.zeros(2), prior_covariance
        ),
    )


def test_map_globality_diagnostic_clusters_agreeing_cold_starts():
    parameters = [
        ThermodynamicParameter(
            "dH_A", "A", lower_bound_j_mol=-100.0, upper_bound_j_mol=100.0
        )
    ]
    problem = InversionProblem(object(), [], parameters)
    result = MAPResult(
        {"A": 10.0},
        2.0,
        [
            {"start_j_mol": [-50.0], "map_j_mol": [10.0], "objective": 2.0},
            {"start_j_mol": [50.0], "map_j_mol": [10.5], "objective": 2.0001},
        ],
    )

    diagnostic = diagnose_map_globality(problem, result)

    assert len(diagnostic.basins) == 1
    assert diagnostic.near_optimal_basin_count == 1
    assert diagnostic.consistent_cold_starts
    assert diagnostic.basins[0].start_count == 2


def test_map_globality_diagnostic_flags_competing_near_optimal_basins():
    parameters = [
        ThermodynamicParameter(
            "dH_A", "A", lower_bound_j_mol=-100.0, upper_bound_j_mol=100.0
        )
    ]
    problem = InversionProblem(object(), [], parameters)
    result = MAPResult(
        {"A": -60.0},
        2.0,
        [
            {"start_j_mol": [-80.0], "map_j_mol": [-60.0], "objective": 2.0},
            {"start_j_mol": [80.0], "map_j_mol": [60.0], "objective": 2.0005},
            {"start_j_mol": [20.0], "map_j_mol": [10.0], "objective": 4.0},
        ],
    )

    diagnostic = diagnose_map_globality(problem, result)

    assert len(diagnostic.basins) == 3
    assert diagnostic.near_optimal_basin_count == 2
    assert not diagnostic.consistent_cold_starts


def test_local_resolution_matrix_quantifies_data_vs_prior_resolution():
    result = local_resolution_matrix(
        np.asarray([[2.0, 0.0]]),
        ("data_resolved", "prior_only"),
        np.eye(2),
        np.asarray([[4.0]]),
    )

    assert np.allclose(result.matrix, np.asarray([[0.5, 0.0], [0.0, 0.0]]))
    assert np.allclose(result.diagonal, np.asarray([0.5, 0.0]))
    assert np.allclose(
        result.posterior_covariance, np.asarray([[0.5, 0.0], [0.0, 1.0]])
    )


def test_local_resolution_matrix_combines_analytical_and_theoretical_covariance():
    result = local_resolution_matrix(
        np.asarray([[1.0]]),
        ("dH_A",),
        np.asarray([[1.0]]),
        observation_covariance=np.asarray([[1.0]]),
        theoretical_covariance=np.asarray([[1.0]]),
    )

    assert np.allclose(result.matrix, [[1.0 / 3.0]])
    assert np.allclose(result.posterior_covariance, [[2.0 / 3.0]])


def test_posterior_averaged_resolution_summarizes_sampled_jacobians():
    result = posterior_averaged_resolution(
        np.asarray([[[1.0]], [[2.0]], [[3.0]]]),
        ("dH_A",),
        np.asarray([[1.0]]),
        credible_interval=(0.0, 1.0),
    )
    expected = np.asarray([0.5, 0.8, 0.9])

    assert result.sample_count == 3
    assert result.parameter_names == ("dH_A",)
    assert np.allclose(result.mean_matrix, [[np.mean(expected)]])
    assert np.allclose(result.lower_matrix, [[expected.min()]])
    assert np.allclose(result.upper_matrix, [[expected.max()]])
    assert np.allclose(result.mean_diagonal, [np.mean(expected)])
    assert np.allclose(result.lower_diagonal, [expected.min()])
    assert np.allclose(result.upper_diagonal, [expected.max()])


def test_posterior_averaged_resolution_rejects_empty_jacobian_stack():
    try:
        posterior_averaged_resolution(
            np.empty((0, 1, 1)), ("dH_A",), np.asarray([[1.0]])
        )
    except ValueError as error:
        assert "jacobians" in str(error)
    else:
        raise AssertionError("Expected ValueError for empty Jacobian samples")


def test_posterior_resolution_from_mcmc_evaluates_each_draw():
    class NonlinearForwardModel:
        def predict(self, experiment, corrections):
            value = corrections.get("A", 0.0)
            return EquilibriumResult(
                True, ("A",), {"A": 0.5 + value**2 / 100.0}, {}, {}
            )

    experiment = Experiment(
        id="nonlinear",
        pressure_bar=1000.0,
        temperature_c=500.0,
        bulk_composition={"A": 1.0},
        observed_phases=("A",),
        phase_fractions={"A": 0.5},
        phase_fraction_sigma={"A": 0.1},
    )
    problem = InversionProblem(
        NonlinearForwardModel(),
        [experiment],
        [
            ThermodynamicParameter(
                "p",
                "A",
                prior_sigma_j_mol=1.0,
                lower_bound_j_mol=-3.0,
                upper_bound_j_mol=3.0,
            )
        ],
    )
    posterior = MCMCResult(
        ("p",), np.asarray([[1.0], [2.0]]), np.zeros(2), 0.5, {}, {}, {}
    )

    result = posterior_resolution_from_mcmc(
        problem, posterior, step_fraction=0.002, credible_interval=(0.1, 0.9)
    )

    assert result.sample_count == 2
    assert result.mean_diagonal[0] > 0.0
    assert result.lower_diagonal[0] < result.upper_diagonal[0]
    assert result.selected_chain_indices == (0, 0)
    assert result.selected_draw_indices == (0, 1)
    assert result.selection_strategy == "evenly_spaced_single_chain"
    assert result.observable_names == ("nonlinear:fraction:A",)
    assert result.credible_interval == (0.1, 0.9)
    assert result.step_fraction == 0.002
    assert result.observation_covariance_source == "problem_assembled"
    assert result.theoretical_covariance_source == "problem_assembled"
    assert np.allclose(result.observation_covariance, [[0.01]])
    assert np.allclose(result.theoretical_covariance, [[0.0]])
    assert np.allclose(result.prior_covariance, [[1.0]])

    explicit = posterior_resolution_from_mcmc(
        problem,
        posterior,
        observation_covariance=np.asarray([[1.0]]),
        theoretical_covariance=np.asarray([[0.5]]),
    )
    assert explicit.observation_covariance_source == "explicit"
    assert explicit.theoretical_covariance_source == "explicit"
    assert np.allclose(explicit.observation_covariance, [[1.0]])
    assert np.allclose(explicit.theoretical_covariance, [[0.5]])
    assert np.allclose(explicit.prior_covariance, [[1.0]])


def test_posterior_resolution_from_mcmc_uses_one_sided_difference_at_bound():
    class LinearForwardModel:
        def predict(self, experiment, corrections):
            value = corrections.get("A", 0.0)
            return EquilibriumResult(True, ("A",), {"A": value / 10.0}, {}, {})

    experiment = Experiment(
        "bounded",
        1000.0,
        500.0,
        {"A": 1.0},
        ("A",),
        phase_fractions={"A": 0.5},
        phase_fraction_sigma={"A": 0.1},
    )
    problem = InversionProblem(
        LinearForwardModel(),
        [experiment],
        [
            ThermodynamicParameter(
                "p",
                "A",
                prior_sigma_j_mol=1.0,
                lower_bound_j_mol=-1.0,
                upper_bound_j_mol=1.0,
            )
        ],
    )
    posterior = MCMCResult(("p",), np.asarray([[1.0]]), np.zeros(1), 0.5, {}, {}, {})

    result = posterior_resolution_from_mcmc(problem, posterior)

    assert result.sample_count == 1
    assert result.mean_diagonal[0] > 0.0


def test_posterior_resolution_reorders_experiment_covariance_keys():
    class TwoFractionForwardModel:
        def predict(self, experiment, corrections):
            value = corrections.get("A", 0.0)
            return EquilibriumResult(
                True,
                ("A", "B"),
                {"A": 0.5 + value / 10.0, "B": 0.4 + value / 20.0},
                {},
                {},
            )

    experiment = Experiment(
        "covariance-order",
        1000.0,
        500.0,
        {"A": 0.5, "B": 0.5},
        ("A", "B"),
        phase_fractions={"B": 0.4, "A": 0.5},
        observation_covariance=[[0.09, 0.01], [0.01, 0.04]],
        observation_covariance_keys=("fraction:B", "fraction:A"),
        theoretical_covariance=[[0.18, 0.02], [0.02, 0.08]],
        theoretical_covariance_keys=("fraction:B", "fraction:A"),
    )
    problem = InversionProblem(
        TwoFractionForwardModel(),
        [experiment],
        [
            ThermodynamicParameter(
                "p",
                "A",
                prior_sigma_j_mol=1.0,
                lower_bound_j_mol=-2.0,
                upper_bound_j_mol=2.0,
            )
        ],
        theoretical_covariance_scale=2.0,
    )
    posterior = MCMCResult(("p",), np.asarray([[0.0]]), np.zeros(1), 0.5, {}, {}, {})

    result = posterior_resolution_from_mcmc(problem, posterior)
    expected = local_resolution_matrix(
        np.asarray([[0.1], [0.05]]),
        ("p",),
        problem.prior_covariance_j_mol2,
        np.asarray([[0.04, 0.01], [0.01, 0.09]])
        + 4.0 * np.asarray([[0.08, 0.02], [0.02, 0.18]]),
    )

    assert result.mean_matrix.shape == (1, 1)
    assert np.allclose(result.mean_matrix, expected.matrix, rtol=1.0e-5)


def test_posterior_resolution_from_mcmc_rejects_parameter_order_mismatch():
    problem = InversionProblem(
        object(),
        [
            Experiment(
                "order-check",
                1000.0,
                500.0,
                {"A": 1.0},
                ("A",),
                phase_fractions={"A": 0.5},
            )
        ],
        [ThermodynamicParameter("p", "A", prior_sigma_j_mol=1.0)],
    )
    posterior = MCMCResult(
        ("different",), np.asarray([[0.0]]), np.zeros(1), 0.5, {}, {}, {}
    )

    try:
        posterior_resolution_from_mcmc(problem, posterior)
    except ValueError as error:
        assert "parameter names" in str(error)
    else:
        raise AssertionError("Expected ValueError for parameter order mismatch")


def test_posterior_resolution_from_mcmc_rejects_invalid_forward_solve():
    class FailedForwardModel:
        def predict(self, experiment, corrections):
            return EquilibriumResult(False, (), {}, {}, {"solver": "failed"})

    experiment = Experiment(
        "failed-solve",
        1000.0,
        500.0,
        {"A": 1.0},
        ("A",),
        phase_fractions={"A": 0.5},
        phase_fraction_sigma={"A": 0.1},
    )
    problem = InversionProblem(
        FailedForwardModel(),
        [experiment],
        [ThermodynamicParameter("p", "A", prior_sigma_j_mol=1.0)],
    )
    posterior = MCMCResult(("p",), np.asarray([[0.0]]), np.zeros(1), 0.5, {}, {}, {})

    try:
        posterior_resolution_from_mcmc(problem, posterior)
    except ValueError as error:
        assert "forward model failed" in str(error)
    else:
        raise AssertionError("Expected ValueError for invalid forward solve")


def test_posterior_resolution_from_mcmc_rejects_phase_assemblage_switch():
    class PhaseSwitchForwardModel:
        def predict(self, experiment, corrections):
            value = corrections.get("A", 0.0)
            phases = ("A",) if value <= 0.0 else ("B",)
            return EquilibriumResult(True, phases, {"A": 0.5 + value / 10.0}, {}, {})

    experiment = Experiment(
        "phase-switch",
        1000.0,
        500.0,
        {"A": 1.0},
        ("A",),
        phase_fractions={"A": 0.5},
        phase_fraction_sigma={"A": 0.1},
    )
    problem = InversionProblem(
        PhaseSwitchForwardModel(),
        [experiment],
        [ThermodynamicParameter("p", "A", prior_sigma_j_mol=1.0)],
    )
    posterior = MCMCResult(("p",), np.asarray([[0.0]]), np.zeros(1), 0.5, {}, {}, {})

    try:
        posterior_resolution_from_mcmc(problem, posterior)
    except ValueError as error:
        assert "assemblage changed" in str(error)
    else:
        raise AssertionError("Expected ValueError for phase-assemblage switch")


def test_posterior_resolution_from_mcmc_balances_draws_across_chains():
    class NonlinearForwardModel:
        def predict(self, experiment, corrections):
            value = corrections.get("A", 0.0)
            return EquilibriumResult(
                True, ("A",), {"A": 0.5 + value**2 / 100.0}, {}, {}
            )

    experiment = Experiment(
        "multi-chain",
        1000.0,
        500.0,
        {"A": 1.0},
        ("A",),
        phase_fractions={"A": 0.5},
        phase_fraction_sigma={"A": 0.1},
    )
    problem = InversionProblem(
        NonlinearForwardModel(),
        [experiment],
        [ThermodynamicParameter("p", "A", prior_sigma_j_mol=1.0)],
    )
    first_chain = MCMCResult(
        ("p",), np.asarray([[0.5], [1.0]]), np.zeros(2), 0.5, {}, {}, {}
    )
    second_chain = MCMCResult(
        ("p",), np.asarray([[2.0], [2.5]]), np.zeros(2), 0.5, {}, {}, {}
    )
    posterior = MultiChainMCMCResult(
        ("p",),
        (first_chain, second_chain),
        {"p": 1.0},
        {"p": 4.0},
    )

    result = posterior_resolution_from_mcmc(problem, posterior, max_samples=2)

    assert result.sample_count == 2
    assert result.lower_diagonal[0] < result.upper_diagonal[0]
    assert result.selected_chain_indices == (0, 1)
    assert result.selected_draw_indices == (0, 0)
    assert result.selection_strategy == "evenly_spaced_balanced_chains"


def test_posterior_resolution_from_mcmc_requires_budget_for_each_chain():
    problem = InversionProblem(
        object(),
        [
            Experiment(
                "budget",
                1000.0,
                500.0,
                {"A": 1.0},
                ("A",),
                phase_fractions={"A": 0.5},
                phase_fraction_sigma={"A": 0.1},
            )
        ],
        [ThermodynamicParameter("p", "A", prior_sigma_j_mol=1.0)],
    )
    chain = MCMCResult(("p",), np.asarray([[0.0]]), np.zeros(1), 0.5, {}, {}, {})
    posterior = MultiChainMCMCResult(("p",), (chain, chain), {"p": 1.0}, {"p": 2.0})

    try:
        posterior_resolution_from_mcmc(problem, posterior, max_samples=1)
    except ValueError as error:
        assert "number of chains" in str(error)
    else:
        raise AssertionError("Expected ValueError when sample budget excludes a chain")


def test_posterior_resolution_from_mcmc_rejects_mismatched_chain_units():
    experiment = Experiment(
        "unit-check",
        1000.0,
        500.0,
        {"A": 1.0},
        ("A",),
        phase_fractions={"A": 0.5},
        phase_fraction_sigma={"A": 0.1},
    )
    problem = InversionProblem(
        object(),
        [experiment],
        [ThermodynamicParameter("p", "A", prior_sigma_j_mol=1.0)],
    )
    first_chain = MCMCResult(
        ("p",),
        np.asarray([[0.0]]),
        np.zeros(1),
        0.5,
        {},
        {},
        {},
        parameter_units=("J/mol",),
    )
    second_chain = MCMCResult(
        ("p",),
        np.asarray([[0.1]]),
        np.zeros(1),
        0.5,
        {},
        {},
        {},
        parameter_units=("kJ/mol",),
    )
    posterior = MultiChainMCMCResult(
        ("p",),
        (first_chain, second_chain),
        {"p": 1.0},
        {"p": 2.0},
        parameter_units=("J/mol",),
    )

    try:
        posterior_resolution_from_mcmc(problem, posterior)
    except ValueError as error:
        assert "chain parameter units" in str(error)
    else:
        raise AssertionError("Expected ValueError for mismatched chain units")


def test_posterior_resolution_from_mcmc_rejects_mismatched_aggregate_units():
    experiment = Experiment(
        "aggregate-unit-check",
        1000.0,
        500.0,
        {"A": 1.0},
        ("A",),
        phase_fractions={"A": 0.5},
        phase_fraction_sigma={"A": 0.1},
    )
    problem = InversionProblem(
        object(),
        [experiment],
        [ThermodynamicParameter("p", "A", prior_sigma_j_mol=1.0)],
    )
    chain = MCMCResult(
        ("p",),
        np.asarray([[0.0]]),
        np.zeros(1),
        0.5,
        {},
        {},
        {},
        parameter_units=("J/mol",),
    )
    posterior = MultiChainMCMCResult(
        ("p",),
        (chain, chain),
        {"p": 1.0},
        {"p": 2.0},
        parameter_units=("kJ/mol",),
    )

    try:
        posterior_resolution_from_mcmc(problem, posterior)
    except ValueError as error:
        assert "posterior parameter units" in str(error)
    else:
        raise AssertionError("Expected ValueError for mismatched aggregate units")


def test_posterior_shape_diagnostic_distinguishes_normal_and_skewed_samples():
    rng = np.random.default_rng(42)
    normal = rng.normal(size=(20000, 1))
    skewed = rng.exponential(size=(20000, 1))

    def result_for(samples):
        return MCMCResult(("dH_A",), samples, np.zeros(len(samples)), 0.25, {}, {}, {})

    normal_shape = diagnose_posterior_shape(result_for(normal))
    skewed_shape = diagnose_posterior_shape(result_for(skewed))

    assert not normal_shape.non_gaussian_indicated
    assert skewed_shape.non_gaussian_indicated
    assert skewed_shape.parameters[0].skewness > 1.5


def test_grid_mode_diagnostic_flags_two_separated_posterior_modes():
    axis = np.linspace(-2.0, 2.0, 9)
    density = np.asarray([0.01, 0.1, 0.8, 0.1, 0.01, 0.1, 1.0, 0.1, 0.01])
    grid = GridPosteriorResult(
        ("dH_A",), (axis,), np.log(density), density, {"dH_A": 0.0}, {"dH_A": 1.0}
    )

    diagnostic = diagnose_grid_modes(grid, minimum_relative_density=0.5)

    assert diagnostic.multimodal
    assert len(diagnostic.modes) == 2
    assert {mode.coordinates["dH_A"] for mode in diagnostic.modes} == {-1.0, 1.0}


def test_grid_information_gain_uses_discrete_correlated_problem_prior():
    parameter = ThermodynamicParameter(
        "dH_A",
        "A",
        prior_sigma_j_mol=1.0,
        lower_bound_j_mol=-1.0,
        upper_bound_j_mol=1.0,
    )
    problem = InversionProblem(object(), [], [parameter])
    axis = np.linspace(-1.0, 1.0, 101)
    prior_density = np.exp(-0.5 * axis**2)
    prior_density /= np.sum(prior_density) * (axis[1] - axis[0])
    grid = GridPosteriorResult(
        ("dH_A",), (axis,), np.log(prior_density), prior_density, {}, {}
    )

    information = grid_information_gain(problem, grid)

    assert abs(information.nats) < 1.0e-12
    assert information.grid_covers_parameter_bounds
