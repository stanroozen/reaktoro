"""Round-trip tests for JSON serialization of Experiment/EquilibriumBracket datasets."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parents[1]))

from thermoinvert.data import EquilibriumBracket, Experiment
from thermoinvert.io import (
    load_brackets,
    load_experiments,
    save_brackets,
    save_experiments,
)


def test_experiments_round_trip_through_json(tmp_path):
    experiments = [
        Experiment(
            id="exp-1",
            pressure_bar=2000.0,
            temperature_c=500.0,
            bulk_composition={"Zn": 1.0, "Ca": 0.5},
            observed_phases=("Znc", "cc"),
            pressure_sigma_bar=20.0,
            temperature_sigma_c=5.0,
            absent_phases=("Wlm",),
            phase_compositions={"Znc": {"Zn": 0.9, "Mg": 0.0}},
            phase_composition_sigma={"Znc": {"Zn": 0.05}},
            phase_composition_concentration={"Znc": 24.0},
            phase_composition_detection_limit={"Znc": 0.02},
            phase_fractions={"Znc": 0.6, "cc": 0.4},
            phase_fraction_sigma={"Znc": 0.02, "cc": 0.02},
            metadata={"source": "unit-test"},
        ),
        Experiment(
            id="exp-2",
            pressure_bar=1500.0,
            temperature_c=400.0,
            bulk_composition={"Mg": 1.0},
            observed_phases=("mag",),
        ),
    ]
    path = tmp_path / "experiments.json"
    save_experiments(experiments, path)
    restored = load_experiments(path)

    assert len(restored) == 2
    assert restored[0] == experiments[0]
    assert restored[1] == experiments[1]


def test_brackets_round_trip_through_json(tmp_path):
    lower = Experiment("lower", 2000.0, 540.0, {}, ("A",))
    upper = Experiment("upper", 2000.0, 560.0, {}, ("B",))
    brackets = [
        EquilibriumBracket(
            "A-to-B",
            lower,
            upper,
            lower_required_phases=("A",),
            upper_required_phases=("B",),
            lower_forbidden_phases=("B",),
            upper_forbidden_phases=("A",),
        )
    ]
    path = tmp_path / "brackets.json"
    save_brackets(brackets, path)
    restored = load_brackets(path)

    assert len(restored) == 1
    assert restored[0] == brackets[0]
