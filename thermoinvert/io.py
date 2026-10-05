"""JSON serialization for Experiment/EquilibriumBracket datasets.

mc_fit reads observational data from `.imc` text files; thermoinvert has no
equivalent on-disk dataset format yet, so literature/experimental data can
only be constructed inline in Python. This module fills that gap with a
plain JSON representation so a dataset can be curated once and reloaded by
any script, independent of forward-model or optimizer code.
"""

import json
from pathlib import Path
from typing import Any

from .data import EquilibriumBracket, Experiment


def experiment_to_dict(experiment: Experiment) -> dict[str, Any]:
    return {
        "id": experiment.id,
        "pressure_bar": experiment.pressure_bar,
        "temperature_c": experiment.temperature_c,
        "bulk_composition": experiment.bulk_composition,
        "observed_phases": list(experiment.observed_phases),
        "pressure_sigma_bar": experiment.pressure_sigma_bar,
        "temperature_sigma_c": experiment.temperature_sigma_c,
        "absent_phases": list(experiment.absent_phases),
        "phase_compositions": experiment.phase_compositions,
        "phase_composition_sigma": experiment.phase_composition_sigma,
        "phase_composition_concentration": experiment.phase_composition_concentration,
        "phase_composition_detection_limit": experiment.phase_composition_detection_limit,
        "phase_fractions": experiment.phase_fractions,
        "phase_fraction_sigma": experiment.phase_fraction_sigma,
        "observation_covariance": experiment.observation_covariance,
        "observation_covariance_keys": list(experiment.observation_covariance_keys),
        "theoretical_covariance": experiment.theoretical_covariance,
        "theoretical_covariance_keys": list(experiment.theoretical_covariance_keys),
        "observation_bias": experiment.observation_bias,
        "metadata": experiment.metadata,
    }


def experiment_from_dict(data: dict[str, Any]) -> Experiment:
    return Experiment(
        id=data["id"],
        pressure_bar=data["pressure_bar"],
        temperature_c=data["temperature_c"],
        bulk_composition=data["bulk_composition"],
        observed_phases=tuple(data["observed_phases"]),
        pressure_sigma_bar=data.get("pressure_sigma_bar"),
        temperature_sigma_c=data.get("temperature_sigma_c"),
        absent_phases=tuple(data.get("absent_phases", ())),
        phase_compositions=data.get("phase_compositions"),
        phase_composition_sigma=data.get("phase_composition_sigma"),
        phase_composition_concentration=data.get("phase_composition_concentration"),
        phase_composition_detection_limit=data.get("phase_composition_detection_limit"),
        phase_fractions=data.get("phase_fractions"),
        phase_fraction_sigma=data.get("phase_fraction_sigma"),
        observation_covariance=data.get("observation_covariance"),
        observation_covariance_keys=tuple(data.get("observation_covariance_keys", ())),
        theoretical_covariance=data.get("theoretical_covariance"),
        theoretical_covariance_keys=tuple(data.get("theoretical_covariance_keys", ())),
        observation_bias=data.get("observation_bias"),
        metadata=data.get("metadata", {}),
    )


def bracket_to_dict(bracket: EquilibriumBracket) -> dict[str, Any]:
    return {
        "id": bracket.id,
        "lower": experiment_to_dict(bracket.lower),
        "upper": experiment_to_dict(bracket.upper),
        "lower_required_phases": list(bracket.lower_required_phases),
        "upper_required_phases": list(bracket.upper_required_phases),
        "lower_forbidden_phases": list(bracket.lower_forbidden_phases),
        "upper_forbidden_phases": list(bracket.upper_forbidden_phases),
    }


def bracket_from_dict(data: dict[str, Any]) -> EquilibriumBracket:
    return EquilibriumBracket(
        id=data["id"],
        lower=experiment_from_dict(data["lower"]),
        upper=experiment_from_dict(data["upper"]),
        lower_required_phases=tuple(data.get("lower_required_phases", ())),
        upper_required_phases=tuple(data.get("upper_required_phases", ())),
        lower_forbidden_phases=tuple(data.get("lower_forbidden_phases", ())),
        upper_forbidden_phases=tuple(data.get("upper_forbidden_phases", ())),
    )


def save_experiments(experiments: list[Experiment], path: str | Path) -> None:
    """Write a list of experiments to one JSON file."""
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(
            [experiment_to_dict(experiment) for experiment in experiments], indent=2
        ),
        encoding="utf-8",
    )


def load_experiments(path: str | Path) -> list[Experiment]:
    """Read a list of experiments previously written by `save_experiments`."""
    records = json.loads(Path(path).read_text(encoding="utf-8"))
    return [experiment_from_dict(record) for record in records]


def save_brackets(brackets: list[EquilibriumBracket], path: str | Path) -> None:
    """Write a list of equilibrium brackets to one JSON file."""
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps([bracket_to_dict(bracket) for bracket in brackets], indent=2),
        encoding="utf-8",
    )


def load_brackets(path: str | Path) -> list[EquilibriumBracket]:
    """Read a list of equilibrium brackets previously written by `save_brackets`."""
    records = json.loads(Path(path).read_text(encoding="utf-8"))
    return [bracket_from_dict(record) for record in records]
