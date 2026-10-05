"""Tests for input-space analytical/thermodynamic perturbation propagation."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parents[1]))

from thermoinvert.data import Experiment, ThermodynamicParameter
from thermoinvert.forward_model import EquilibriumResult
from thermoinvert.input_perturbation import (
    InputPerturbationConfig,
    _perturb_experiment,
    _sample_covariance,
    covariance_diagnostics,
    load_packed_covariance_json,
    run_input_perturbation,
)
from thermoinvert.archive import (
    read_input_perturbation_result,
    write_input_perturbation_result,
)
from thermoinvert.inversion import InversionProblem


class SyntheticForwardModel:
    def predict(self, experiment, corrections):
        thermo_offset = corrections.get("A", 0.0)
        target = experiment.phase_compositions["A"]["A"]
        value = 0.5 + thermo_offset / 10000.0
        return EquilibriumResult(
            True,
            ("A",),
            {"A": 1.0},
            {"A": {"A": value}},
            {"target": target},
        )


def _problem() -> InversionProblem:
    experiment = Experiment(
        id="synthetic-input-uncertainty",
        pressure_bar=2000.0,
        temperature_c=500.0,
        bulk_composition={"A": 1.0},
        observed_phases=("A",),
        phase_compositions={"A": {"A": 0.62}},
        phase_composition_sigma={"A": {"A": 0.01}},
    )
    parameter = ThermodynamicParameter(
        "dH_A",
        "A",
        prior_sigma_j_mol=5000.0,
        lower_bound_j_mol=-4000.0,
        upper_bound_j_mol=4000.0,
    )
    return InversionProblem(SyntheticForwardModel(), [experiment], [parameter])


def test_analytical_input_perturbations_produce_map_distribution():
    result = run_input_perturbation(
        _problem(),
        InputPerturbationConfig(
            n_iterations=8,
            source="analytical",
            seed=2,
            optimize_starts=2,
            optimize_iterations=25,
        ),
    )
    assert result.source == "analytical"
    assert result.success_count == 8
    assert result.corrections_samples.shape == (8, 1)
    assert np.isfinite(result.correction_mean["dH_A"])
    assert result.correction_std["dH_A"] > 0.0


def test_thermodynamic_input_perturbations_use_nuisance_offsets():
    result = run_input_perturbation(
        _problem(),
        InputPerturbationConfig(
            n_iterations=8,
            source="thermodynamic",
            thermodynamic_sigma_j_mol={"A": 500.0},
            seed=3,
            optimize_starts=2,
            optimize_iterations=25,
        ),
    )
    assert result.source == "thermodynamic"
    assert result.success_count == 8
    assert result.correction_std["dH_A"] > 0.0


def test_all_input_perturbation_requires_thermodynamic_sigma():
    try:
        run_input_perturbation(
            _problem(), InputPerturbationConfig(n_iterations=1, source="all")
        )
    except ValueError:
        pass
    else:
        raise AssertionError("Expected ValueError without thermodynamic sigma")


def test_correlated_covariance_sampling_preserves_off_diagonal_correlation():
    covariance = np.asarray([[1.0, 0.8], [0.8, 1.0]])
    samples = np.asarray(
        [
            _sample_covariance(covariance, np.random.default_rng(seed))
            for seed in range(5000)
        ]
    )
    assert np.corrcoef(samples.T)[0, 1] > 0.75
    diagnostics = covariance_diagnostics(covariance)
    assert diagnostics["is_psd"]
    assert diagnostics["effective_rank"] == 2


def test_covariance_sampling_accepts_semidefinite_and_rejects_nonfinite():
    covariance = np.asarray([[1.0, 1.0], [1.0, 1.0]])
    sample = _sample_covariance(covariance, np.random.default_rng(4))
    assert np.isclose(sample[0], sample[1])
    assert covariance_diagnostics(covariance)["effective_rank"] == 1

    try:
        covariance_diagnostics(np.asarray([[1.0, np.nan], [np.nan, 1.0]]))
    except ValueError:
        pass
    else:
        raise AssertionError("Expected ValueError for non-finite covariance")


def test_analytical_perturbations_use_experiment_observation_covariance():
    experiment = Experiment(
        id="correlated-observations",
        pressure_bar=2000.0,
        temperature_c=500.0,
        bulk_composition={"A": 1.0},
        observed_phases=("A",),
        phase_compositions={"A": {"x": 0.2, "y": 0.4}},
        observation_covariance=[[0.01, 0.008], [0.008, 0.01]],
        observation_covariance_keys=("composition:A:x", "composition:A:y"),
    )
    rng = np.random.default_rng(3)
    perturbations = []
    for _ in range(2000):
        perturbed = _perturb_experiment(experiment, rng, 1.0)
        perturbations.append(
            [
                perturbed.phase_compositions["A"][key]
                - experiment.phase_compositions["A"][key]
                for key in ("x", "y")
            ]
        )
    perturbations = np.asarray(perturbations)

    assert np.corrcoef(perturbations.T)[0, 1] > 0.7


def test_packed_holland_powell_covariance_json_loader(tmp_path):
    path = tmp_path / "covariance.json"
    path.write_text(
        '{"Entities": ["A", "B"], "PackedUpperTriangle": [1.0, 0.5, 4.0]}',
        encoding="utf-8",
    )
    keys, covariance = load_packed_covariance_json(path, scale=1000.0)
    assert keys == ("A", "B")
    assert np.allclose(covariance, [[1000.0, 500.0], [500.0, 4000.0]])


def test_correlated_perturbation_result_retains_covariance_provenance():
    covariance = np.asarray([[100.0, 80.0], [80.0, 100.0]])
    result = run_input_perturbation(
        _problem(),
        InputPerturbationConfig(
            n_iterations=2,
            source="thermodynamic",
            thermodynamic_covariance_j_mol2=covariance,
            thermodynamic_covariance_keys=("A", "unused"),
            seed=8,
            optimize_starts=1,
            optimize_iterations=10,
        ),
    )
    assert result.thermodynamic_covariance_keys == ("A", "unused")
    assert result.thermodynamic_covariance_diagnostics["is_psd"]
    assert result.thermodynamic_covariance_diagnostics["effective_rank"] == 2


def test_input_perturbation_result_round_trips_with_covariance_provenance(tmp_path):
    covariance = np.asarray([[100.0, 80.0], [80.0, 100.0]])
    original = run_input_perturbation(
        _problem(),
        InputPerturbationConfig(
            n_iterations=2,
            source="thermodynamic",
            thermodynamic_covariance_j_mol2=covariance,
            thermodynamic_covariance_keys=("A", "unused"),
            seed=9,
            optimize_starts=1,
            optimize_iterations=10,
        ),
    )
    write_input_perturbation_result(original, tmp_path)
    restored = read_input_perturbation_result(tmp_path)
    assert restored.source == original.source
    assert np.allclose(restored.corrections_samples, original.corrections_samples)
    assert restored.thermodynamic_covariance_keys == ("A", "unused")
    assert restored.thermodynamic_covariance_diagnostics["is_psd"]
