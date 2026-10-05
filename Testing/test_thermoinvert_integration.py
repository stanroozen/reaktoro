"""End-to-end integration test exercising the full thermoinvert pipeline together:
multi-parameter MAP with brackets and condition-uncertainty marginalization,
leave-one-out cross-validation, MCMC sampling, Laplace covariance, and a
posterior predictive check -- all on one synthetic two-parameter scenario.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parents[1]))

from thermoinvert.data import EquilibriumBracket, Experiment, ThermodynamicParameter
from thermoinvert.forward_model import EquilibriumResult
from thermoinvert.inversion import InversionProblem
from thermoinvert.posterior import laplace_covariance, mcmc_sample
from thermoinvert.validation import (
    leave_one_out_cross_validation,
    posterior_predictive_check,
)

TRUE_A = 1100.0
TRUE_B = -600.0
THRESHOLD_C = 550.0


class TwoPhaseForwardModel:
    """Phase stability depends only on temperature; corrections only shift
    the predicted composition of whichever phase is stable, keeping the
    objective smooth enough for gradient-sensitive diagnostics."""

    def predict(self, experiment, corrections):
        value_a = corrections.get("A", 0.0)
        value_b = corrections.get("B", 0.0)
        if experiment.temperature_c < THRESHOLD_C:
            stable, fractions = ("A",), {"A": 1.0}
        else:
            stable, fractions = ("B",), {"B": 1.0}
        compositions = {
            "A": {"A": 0.7 + (value_a - TRUE_A) / 10000.0},
            "B": {"B": 0.7 + (value_b - TRUE_B) / 10000.0},
        }
        return EquilibriumResult(True, stable, fractions, compositions, {})


def _build_problem():
    experiments = [
        Experiment(
            id=f"A-{temperature}",
            pressure_bar=2000.0,
            temperature_c=temperature,
            bulk_composition={"A": 1.0},
            observed_phases=("A",),
            phase_compositions={"A": {"A": 0.7}},
            phase_composition_sigma={"A": {"A": 0.01}},
            temperature_sigma_c=2.0 if temperature == 520.0 else None,
            pressure_sigma_bar=20.0 if temperature == 520.0 else None,
        )
        for temperature in (500.0, 510.0, 520.0)
    ] + [
        Experiment(
            id=f"B-{temperature}",
            pressure_bar=2000.0,
            temperature_c=temperature,
            bulk_composition={"B": 1.0},
            observed_phases=("B",),
            phase_compositions={"B": {"B": 0.7}},
            phase_composition_sigma={"B": {"B": 0.01}},
            temperature_sigma_c=2.0 if temperature == 580.0 else None,
            pressure_sigma_bar=20.0 if temperature == 580.0 else None,
        )
        for temperature in (580.0, 590.0, 600.0)
    ]
    bracket = EquilibriumBracket(
        "A-to-B",
        lower=Experiment("bracket-lower", 2000.0, 540.0, {}, ("A",)),
        upper=Experiment("bracket-upper", 2000.0, 560.0, {}, ("B",)),
        lower_required_phases=("A",),
        upper_required_phases=("B",),
        lower_forbidden_phases=("B",),
        upper_forbidden_phases=("A",),
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
    return InversionProblem(
        TwoPhaseForwardModel(),
        experiments,
        parameters,
        brackets=[bracket],
        marginalize_condition_uncertainty=True,
        condition_quadrature_order=3,
    )


def test_full_pipeline_runs_consistently_end_to_end():
    problem = _build_problem()

    map_result = problem.optimize_multi_parameter(starts=6, iterations=100)
    assert abs(map_result.corrections_j_mol["A"] - TRUE_A) < 150.0
    assert abs(map_result.corrections_j_mol["B"] - TRUE_B) < 150.0

    loocv = leave_one_out_cross_validation(
        problem, optimize_starts=3, optimize_iterations=40
    )
    assert len(loocv.folds) == len(problem.experiments)
    assert loocv.mean_misfit < 5.0

    map_point = [map_result.corrections_j_mol["A"], map_result.corrections_j_mol["B"]]
    approximation = laplace_covariance(problem, map_point)
    assert approximation.hessian.shape == (2, 2)

    mcmc_result = mcmc_sample(
        problem,
        map_point,
        n_samples=4000,
        burn_in=500,
        proposal_sigma_j_mol=[80.0, 80.0],
        seed=2,
    )
    assert 0.05 < mcmc_result.acceptance_rate < 0.95
    assert abs(mcmc_result.posterior_mean["dH_A"] - TRUE_A) < 400.0
    assert abs(mcmc_result.posterior_mean["dH_B"] - TRUE_B) < 400.0

    ppc = posterior_predictive_check(
        problem, problem.experiments[0], mcmc_result, n_predictive_samples=50
    )
    assert ppc.within_95["phase_present:A"]
    assert ppc.within_95["composition:A:A"]
