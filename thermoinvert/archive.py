"""Persistent CSV/JSON archival for MCMC, LOOCV, and grid-posterior results.

Extends the `sensitivity.write_sensitivity_report` pattern so results can be
reloaded and re-plotted later without rerunning the (potentially expensive)
forward model.
"""

import csv
import json
from pathlib import Path

import numpy as np

from .data import (
    DirichletConcentrationParameter,
    Parameter,
    ThermodynamicParameter,
)
from .forward_model import ReaktoroForwardModel
from .bias import (
    GroupedObservationBiasReport,
    GroupedObservationBiasSummary,
    ObservationBiasProfileResult,
)
from .hierarchical import (
    CombinedGroupedUncertaintyMCMCResult,
    GroupedHierarchicalMCMCResult,
    GroupedObservationBiasMCMCResult,
    HierarchicalMCMCResult,
)
from .input_perturbation import InputPerturbationResult
from .inversion import InversionProblem
from .joint_inversion import JointInversionProblem
from .io import bracket_to_dict, experiment_to_dict
from .forward_model import ReaktoroForwardModel
from .posterior import (
    GridPosteriorResult,
    MCMCResult,
    MultiChainMCMCResult,
    ParallelTemperingResult,
)
from .pt_inversion import PTCalibrationBiasProblem, PTInversionProblem
from .pt_calibration import GroupedPTCalibrationProblem
from .sensitivity import PosteriorResolutionResult
from .validation import (
    BracketPosteriorPredictiveResult,
    LOOCVFold,
    LOOCVResult,
    PosteriorPredictiveResult,
)


def _forward_model_provenance(forward_model: object) -> dict[str, object]:
    implementation = type(forward_model)
    record: dict[str, object] = {
        "implementation": f"{implementation.__module__}.{implementation.__qualname__}",
        "settings": None,
    }
    if isinstance(forward_model, ReaktoroForwardModel):
        record["settings"] = forward_model.provenance()
    return record


def _required_configuration_value(configuration: dict[str, object], key: str) -> object:
    if key not in configuration:
        raise ValueError(f"configuration archive is missing required field: {key}")
    return configuration[key]


def _validate_records(records: object, name: str) -> list[dict[str, object]]:
    if not isinstance(records, list) or any(
        not isinstance(record, dict) for record in records
    ):
        raise ValueError(f"{name} must be an array of JSON objects")
    ids = [record.get("id") for record in records]
    if any(not isinstance(identifier, str) or not identifier for identifier in ids):
        raise ValueError(f"{name} records must have non-empty string ids")
    if len(set(ids)) != len(ids):
        raise ValueError(f"{name} record ids must be unique")
    return records


def _validate_experiment_records(records: list[dict[str, object]], name: str) -> None:
    required = {
        "id",
        "pressure_bar",
        "temperature_c",
        "bulk_composition",
        "observed_phases",
    }
    for record in records:
        if not required.issubset(record):
            raise ValueError(f"{name} records are missing required experiment fields")
        for field in ("pressure_bar", "temperature_c"):
            value = record[field]
            if not isinstance(value, (int, float)) or not np.isfinite(value):
                raise ValueError(f"{name} field {field} must be finite")
        if not isinstance(record["bulk_composition"], dict) or not isinstance(
            record["observed_phases"], list
        ):
            raise ValueError(f"{name} composition and phase fields have invalid types")
        concentrations = record.get("phase_composition_concentration")
        if concentrations is not None:
            compositions = record.get("phase_compositions") or {}
            if not isinstance(concentrations, dict) or not isinstance(
                compositions, dict
            ):
                raise ValueError(
                    f"{name} phase composition concentrations must be an object"
                )
            if any(
                phase not in compositions
                or not isinstance(value, (int, float))
                or not np.isfinite(value)
                or value <= 0.0
                for phase, value in concentrations.items()
            ):
                raise ValueError(
                    f"{name} phase composition concentrations must be positive "
                    "and refer to recorded phase compositions"
                )
        detection_limits = record.get("phase_composition_detection_limit")
        if detection_limits is not None:
            compositions = record.get("phase_compositions") or {}
            if not isinstance(detection_limits, dict) or not isinstance(
                compositions, dict
            ):
                raise ValueError(f"{name} phase detection limits must be an object")
            if any(
                phase not in compositions
                or not isinstance(compositions[phase], dict)
                or not any(value == 0.0 for value in compositions[phase].values())
                or not isinstance(limit, (int, float))
                or not np.isfinite(limit)
                or not 0.0 < limit < 1.0
                for phase, limit in detection_limits.items()
            ):
                raise ValueError(
                    f"{name} phase detection limits must lie in (0, 1) and "
                    "refer to zero-marked compositions"
                )


def _validate_bracket_records(records: list[dict[str, object]], name: str) -> None:
    for record in records:
        if not isinstance(record.get("lower"), dict) or not isinstance(
            record.get("upper"), dict
        ):
            raise ValueError(f"{name} records must contain lower and upper experiments")


def _validate_inversion_settings(configuration: dict[str, object]) -> None:
    for field in ("missing_phase_penalty", "extra_phase_penalty"):
        value = configuration.get(field)
        if not isinstance(value, (int, float)) or not np.isfinite(value) or value < 0.0:
            raise ValueError(
                f"inversion setting {field} must be finite and non-negative"
            )
    criterion = configuration.get("objective_criterion")
    if criterion not in {"map", "bayes"}:
        raise ValueError("invalid inversion objective criterion")
    if configuration.get("misfit_form") not in {
        "weighted_lsq",
        "logit_lsq",
        "alr_lsq",
        "dirichlet_nll",
        "chi_square",
        "raw_lsq",
    }:
        raise ValueError("invalid inversion misfit form")
    order = configuration.get("condition_quadrature_order")
    if not isinstance(order, int) or order < 1:
        raise ValueError("condition quadrature order must be a positive integer")
    scale = configuration.get("theoretical_covariance_scale", 1.0)
    if not isinstance(scale, (int, float)) or not np.isfinite(scale) or scale < 0.0:
        raise ValueError("theoretical covariance scale must be finite and non-negative")
    group_key = configuration.get("theoretical_covariance_scale_metadata_key")
    group_scales = configuration.get("theoretical_covariance_scales", {})
    if group_key is not None and (not isinstance(group_key, str) or not group_key):
        raise ValueError(
            "theoretical covariance scale metadata key must be a non-empty string"
        )
    if not isinstance(group_scales, dict) or any(
        not isinstance(key, str)
        or not key
        or not isinstance(value, (int, float))
        or not np.isfinite(value)
        or value < 0.0
        for key, value in group_scales.items()
    ):
        raise ValueError(
            "theoretical covariance group scales must be finite and non-negative"
        )
    if group_key is not None and not group_scales:
        raise ValueError(
            "theoretical covariance group scales are required with a metadata key"
        )


def write_inversion_configuration(
    problem: InversionProblem, output_dir: str | Path
) -> None:
    """Write inversion provenance excluding the runtime forward-model object.

    The archive captures the parameter order that defines the covariance
    matrix, all likelihood uncertainty records, and optimization-relevant
    problem settings. Reconstructing a runnable problem additionally requires
    selecting the original forward-model implementation and database path.
    """
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    parameters = []
    for parameter in problem.parameters:
        if isinstance(parameter, DirichletConcentrationParameter):
            parameters.append(
                {
                    "parameter_kind": "dirichlet_concentration",
                    "name": parameter.name,
                    "phase": parameter.phase,
                    "unit": parameter.unit,
                    "prior_mean": parameter.prior_mean,
                    "prior_sigma": parameter.prior_sigma,
                    "lower_bound": parameter.lower_bound,
                    "upper_bound": parameter.upper_bound,
                    "key": parameter.key,
                }
            )
        elif isinstance(parameter, ThermodynamicParameter):
            parameters.append(
                {
                    "parameter_kind": "thermodynamic",
                    "name": parameter.name,
                    "phase": parameter.phase,
                    "field": parameter.field,
                    "parameter_type": parameter.parameter_type,
                    "prior_mean": parameter.prior_mean,
                    "prior_sigma": parameter.prior_sigma,
                    "lower_bound": parameter.lower_bound,
                    "upper_bound": parameter.upper_bound,
                    "prior_mean_j_mol": parameter.prior_mean_j_mol,
                    "prior_sigma_j_mol": parameter.prior_sigma_j_mol,
                    "lower_bound_j_mol": parameter.lower_bound_j_mol,
                    "upper_bound_j_mol": parameter.upper_bound_j_mol,
                    "key": parameter.key,
                }
            )
        elif isinstance(parameter, Parameter):
            parameters.append(
                {
                    "parameter_kind": "scalar",
                    "name": parameter.name,
                    "unit": parameter.unit,
                    "prior_mean": parameter.prior_mean,
                    "prior_sigma": parameter.prior_sigma,
                    "lower_bound": parameter.lower_bound,
                    "upper_bound": parameter.upper_bound,
                    "key": parameter.key,
                }
            )
        else:
            raise TypeError("inversion parameters must be Parameter instances")
    configuration = {
        "forward_model": _forward_model_provenance(problem.forward_model),
        "parameter_order": [parameter.key for parameter in problem.parameters],
        "parameters": parameters,
        "prior_covariance": problem.prior_covariance.tolist(),
        "prior_covariance_j_mol2": problem.prior_covariance.tolist(),
        "missing_phase_penalty": problem.missing_phase_penalty,
        "extra_phase_penalty": problem.extra_phase_penalty,
        "marginalize_condition_uncertainty": problem.marginalize_condition_uncertainty,
        "condition_quadrature_order": problem.condition_quadrature_order,
        "objective_criterion": problem.objective_criterion,
        "misfit_form": problem.misfit_form,
        "theoretical_covariance_scale": problem.theoretical_covariance_scale,
        "theoretical_covariance_scale_metadata_key": problem.theoretical_covariance_scale_metadata_key,
        "theoretical_covariance_scales": problem.theoretical_covariance_scales,
    }
    (output / "inversion_configuration.json").write_text(
        json.dumps(configuration, indent=2), encoding="utf-8"
    )
    (output / "experiments.json").write_text(
        json.dumps(
            [experiment_to_dict(item) for item in problem.experiments], indent=2
        ),
        encoding="utf-8",
    )
    (output / "brackets.json").write_text(
        json.dumps([bracket_to_dict(item) for item in problem.brackets], indent=2),
        encoding="utf-8",
    )


def read_inversion_configuration(input_dir: str | Path) -> dict[str, object]:
    """Read and validate a database-parameter inversion archive."""
    input_path = Path(input_dir)
    configuration = json.loads(
        (input_path / "inversion_configuration.json").read_text(encoding="utf-8")
    )
    if not isinstance(configuration, dict):
        raise ValueError("inversion configuration must be a JSON object")
    _validate_inversion_settings(configuration)
    parameter_order = configuration.get("parameter_order")
    parameters = configuration.get("parameters")
    if not isinstance(parameter_order, list) or not isinstance(parameters, list):
        raise ValueError("invalid parameter metadata in inversion archive")
    if any(not isinstance(key, str) or not key for key in parameter_order):
        raise ValueError("inversion parameter keys must be non-empty strings")
    if len(parameter_order) == 0 or len(set(parameter_order)) != len(parameter_order):
        raise ValueError("inversion parameter keys must be unique and non-empty")
    if any(
        not isinstance(parameter, dict) or "key" not in parameter
        for parameter in parameters
    ):
        raise ValueError("inversion parameter metadata must contain keys")
    if any(
        parameter["key"] != parameter_order[index]
        for index, parameter in enumerate(parameters)
    ):
        raise ValueError("inversion parameter metadata order is invalid")
    for parameter in parameters:
        parameter_kind = parameter.get("parameter_kind", "thermodynamic")
        if parameter_kind == "dirichlet_concentration":
            for field in ("name", "phase", "unit", "key"):
                if not isinstance(parameter.get(field), str) or not parameter[field]:
                    raise ValueError(
                        f"Dirichlet concentration parameter field {field} must be non-empty"
                    )
            if parameter["unit"] != "1":
                raise ValueError("Dirichlet concentration parameter unit must be '1'")
            if parameter["key"] != f"dirichlet_concentration:{parameter['phase']}":
                raise ValueError("Dirichlet concentration parameter key is invalid")
            for field in ("prior_mean", "prior_sigma"):
                value = parameter.get(field)
                if not isinstance(value, (int, float)) or not np.isfinite(value):
                    raise ValueError(
                        f"Dirichlet concentration parameter {field} must be finite"
                    )
            if parameter["prior_mean"] <= 0.0 or parameter["prior_sigma"] <= 0.0:
                raise ValueError("Dirichlet concentration priors must be positive")
            lower, upper = parameter.get("lower_bound"), parameter.get("upper_bound")
            if (
                not isinstance(lower, (int, float))
                or not np.isfinite(lower)
                or lower <= 0.0
                or not isinstance(upper, (int, float))
                or not np.isfinite(upper)
                or lower >= upper
            ):
                raise ValueError(
                    "Dirichlet concentration bounds must be positive and increasing"
                )
            if not lower <= parameter["prior_mean"] <= upper:
                raise ValueError(
                    "Dirichlet concentration prior mean must lie within bounds"
                )
            continue
        if parameter_kind == "scalar":
            for field in ("name", "unit"):
                if not isinstance(parameter.get(field), str) or not parameter[field]:
                    raise ValueError(
                        f"scalar parameter field {field} must be non-empty"
                    )
            for field in ("prior_mean", "prior_sigma"):
                value = parameter.get(field)
                if not isinstance(value, (int, float)) or not np.isfinite(value):
                    raise ValueError(f"scalar parameter field {field} must be finite")
            if parameter["prior_sigma"] <= 0.0:
                raise ValueError("scalar parameter prior sigma must be positive")
            if parameter.get("key") != parameter["name"]:
                raise ValueError("scalar parameter key must match its name")
            continue
        for field in ("name", "phase", "field", "parameter_type"):
            if not isinstance(parameter.get(field), str) or not parameter[field]:
                raise ValueError(
                    f"inversion parameter field {field} must be a non-empty string"
                )
        expected_key = (
            parameter["phase"]
            if parameter["field"] == "H0"
            else f"{parameter['phase']}:{parameter['field']}"
        )
        if parameter["key"] != expected_key:
            raise ValueError(
                "inversion parameter key is inconsistent with phase and field"
            )
        if parameter["field"] not in {"H0", "G0", "V0", "Hf", "GH"}:
            raise ValueError("inversion parameter field is unsupported")
        for field in ("prior_mean_j_mol", "prior_sigma_j_mol"):
            value = parameter.get(field)
            if not isinstance(value, (int, float)) or not np.isfinite(value):
                raise ValueError(f"inversion parameter field {field} must be finite")
        if parameter["prior_sigma_j_mol"] <= 0.0:
            raise ValueError("inversion parameter prior sigmas must be positive")
        lower = parameter.get("lower_bound_j_mol")
        upper = parameter.get("upper_bound_j_mol")
        if lower is not None and (
            not isinstance(lower, (int, float)) or not np.isfinite(lower)
        ):
            raise ValueError("inversion parameter lower bounds must be finite or null")
        if upper is not None and (
            not isinstance(upper, (int, float)) or not np.isfinite(upper)
        ):
            raise ValueError("inversion parameter upper bounds must be finite or null")
        if lower is not None and upper is not None and lower >= upper:
            raise ValueError("inversion parameter bounds must increase")
    if [parameter.get("key") for parameter in parameters] != parameter_order:
        raise ValueError("inversion parameter order does not match parameter metadata")
    covariance_payload = configuration.get(
        "prior_covariance", configuration.get("prior_covariance_j_mol2")
    )
    if covariance_payload is None:
        raise ValueError("inversion configuration is missing prior_covariance")
    covariance = np.asarray(covariance_payload, dtype=float)
    if covariance.shape != (len(parameter_order), len(parameter_order)):
        raise ValueError("inversion prior covariance shape does not match parameters")
    if not np.all(np.isfinite(covariance)):
        raise ValueError("inversion prior covariance must contain finite values")
    covariance = 0.5 * (covariance + covariance.T)
    if np.any(np.linalg.eigvalsh(covariance) <= 0.0):
        raise ValueError("inversion prior covariance must be positive definite")
    experiments = json.loads(
        (input_path / "experiments.json").read_text(encoding="utf-8")
    )
    brackets = json.loads((input_path / "brackets.json").read_text(encoding="utf-8"))
    _validate_records(experiments, "inversion experiments")
    _validate_records(brackets, "inversion brackets")
    _validate_experiment_records(experiments, "inversion experiments")
    _validate_bracket_records(brackets, "inversion brackets")
    return {
        "configuration": configuration,
        "experiments": experiments,
        "brackets": brackets,
    }


def write_pt_inversion_configuration(
    problem: PTInversionProblem | PTCalibrationBiasProblem, output_dir: str | Path
) -> None:
    """Write P-T inversion provenance, excluding the runtime forward model."""
    base_problem = (
        problem.base_problem
        if isinstance(problem, PTCalibrationBiasProblem)
        else problem
    )
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    parameters = [
        (
            {
                "parameter_kind": "scalar",
                "name": parameter.name,
                "unit": parameter.unit,
                "prior_mean": parameter.prior_mean,
                "prior_sigma": parameter.prior_sigma,
                "lower_bound": parameter.lower_bound,
                "upper_bound": parameter.upper_bound,
            }
            if isinstance(parameter, Parameter)
            else {
                "parameter_kind": "thermodynamic",
                "name": parameter.name,
                "phase": parameter.phase,
                "field": parameter.field,
                "parameter_type": parameter.parameter_type,
                "prior_mean": parameter.prior_mean,
                "prior_sigma": parameter.prior_sigma,
                "lower_bound": parameter.lower_bound,
                "upper_bound": parameter.upper_bound,
                "prior_mean_j_mol": parameter.prior_mean_j_mol,
                "prior_sigma_j_mol": parameter.prior_sigma_j_mol,
                "lower_bound_j_mol": parameter.lower_bound_j_mol,
                "upper_bound_j_mol": parameter.upper_bound_j_mol,
                "key": parameter.key,
            }
        )
        for parameter in problem.parameters
    ]
    configuration = {
        "forward_model": _forward_model_provenance(base_problem.forward_model),
        "parameter_order": [parameter.name for parameter in problem.parameters],
        "parameters": parameters,
        "prior_covariance": problem.prior_covariance.tolist(),
        "prior_covariance_j_mol2": problem.prior_covariance.tolist(),
        "missing_phase_penalty": base_problem.missing_phase_penalty,
        "extra_phase_penalty": base_problem.extra_phase_penalty,
        "misfit_form": base_problem.misfit_form,
        "theoretical_covariance_scale": base_problem.theoretical_covariance_scale,
        "corrections_j_mol": base_problem.corrections,
        "problem_mode": "bracket" if base_problem.brackets else "composition",
        "bracket_ids": [bracket.id for bracket in base_problem.brackets],
        "calibration_bias": isinstance(problem, PTCalibrationBiasProblem),
    }
    (output / "pt_inversion_configuration.json").write_text(
        json.dumps(configuration, indent=2), encoding="utf-8"
    )
    (output / "experiment.json").write_text(
        json.dumps(experiment_to_dict(base_problem.experiment), indent=2),
        encoding="utf-8",
    )
    (output / "brackets.json").write_text(
        json.dumps([bracket_to_dict(item) for item in base_problem.brackets], indent=2),
        encoding="utf-8",
    )


def write_joint_inversion_configuration(
    problem: JointInversionProblem, output_dir: str | Path
) -> None:
    """Write joint thermodynamic/P-T inversion provenance."""
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    parameters = [
        {
            "parameter_kind": "thermodynamic",
            "name": parameter.name,
            "phase": parameter.phase,
            "field": parameter.field,
            "parameter_type": parameter.parameter_type,
            "prior_mean": parameter.prior_mean,
            "prior_sigma": parameter.prior_sigma,
            "lower_bound": parameter.lower_bound,
            "upper_bound": parameter.upper_bound,
            "prior_mean_j_mol": parameter.prior_mean_j_mol,
            "prior_sigma_j_mol": parameter.prior_sigma_j_mol,
            "lower_bound_j_mol": parameter.lower_bound_j_mol,
            "upper_bound_j_mol": parameter.upper_bound_j_mol,
            "key": parameter.key,
        }
        for parameter in problem.thermodynamic_parameters
    ]
    for pressure, temperature in zip(
        problem.pressure_parameters, problem.temperature_parameters
    ):
        parameters.extend(
            {
                "parameter_kind": "scalar",
                "name": parameter.name,
                "unit": parameter.unit,
                "prior_mean": parameter.prior_mean,
                "prior_sigma": parameter.prior_sigma,
                "lower_bound": parameter.lower_bound,
                "upper_bound": parameter.upper_bound,
                "key": parameter.key,
            }
            for parameter in (pressure, temperature)
        )
    configuration = {
        "forward_model": _forward_model_provenance(problem.forward_model),
        "parameter_order": [parameter.name for parameter in problem.parameters],
        "parameters": parameters,
        "experiment_count": len(problem.experiments),
        "experiment_ids": [experiment.id for experiment in problem.experiments],
        "prior_covariance": problem.prior_covariance.tolist(),
        "missing_phase_penalty": problem.missing_phase_penalty,
        "extra_phase_penalty": problem.extra_phase_penalty,
        "misfit_form": problem.misfit_form,
        "theoretical_covariance_scale": problem.theoretical_covariance_scale,
    }
    (output / "joint_inversion_configuration.json").write_text(
        json.dumps(configuration, indent=2), encoding="utf-8"
    )
    if len(problem.experiments) == 1:
        (output / "experiment.json").write_text(
            json.dumps(experiment_to_dict(problem.experiments[0]), indent=2),
            encoding="utf-8",
        )
    else:
        (output / "experiments.json").write_text(
            json.dumps(
                [experiment_to_dict(experiment) for experiment in problem.experiments],
                indent=2,
            ),
            encoding="utf-8",
        )


def read_joint_inversion_configuration(input_dir: str | Path) -> dict[str, object]:
    """Read and validate a joint thermodynamic/P-T configuration archive."""
    input_path = Path(input_dir)
    configuration = json.loads(
        (input_path / "joint_inversion_configuration.json").read_text(encoding="utf-8")
    )
    if not isinstance(configuration, dict):
        raise TypeError("joint inversion configuration must be a JSON object")
    order = configuration.get("parameter_order")
    parameters = configuration.get("parameters")
    if not isinstance(order, list) or not isinstance(parameters, list):
        raise ValueError("joint configuration must contain parameter arrays")
    if len(order) != len(parameters) or len(set(order)) != len(order):
        raise ValueError("joint parameter order must be unique and aligned")
    if [parameter.get("name") for parameter in parameters] != order:
        raise ValueError("joint parameter metadata order is invalid")
    covariance = np.asarray(configuration.get("prior_covariance"), dtype=float)
    if covariance.shape != (len(order), len(order)) or not np.all(
        np.isfinite(covariance)
    ):
        raise ValueError(
            "joint prior covariance must match parameter count and be finite"
        )
    experiment_count = configuration.get("experiment_count", 1)
    if experiment_count == 1:
        experiment = json.loads(
            (input_path / "experiment.json").read_text(encoding="utf-8")
        )
        return {"configuration": configuration, "experiment": experiment}
    experiments = json.loads(
        (input_path / "experiments.json").read_text(encoding="utf-8")
    )
    if not isinstance(experiments, list) or len(experiments) != experiment_count:
        raise ValueError("joint experiment archive count does not match configuration")
    return {"configuration": configuration, "experiments": experiments}


def read_pt_inversion_configuration(input_dir: str | Path) -> dict[str, object]:
    """Read and validate a P-T configuration archive without a forward model."""
    input_path = Path(input_dir)
    configuration = json.loads(
        (input_path / "pt_inversion_configuration.json").read_text(encoding="utf-8")
    )
    if not isinstance(configuration, dict):
        raise ValueError("P-T configuration must be a JSON object")
    for field in ("missing_phase_penalty", "extra_phase_penalty"):
        value = configuration.get(field)
        if not isinstance(value, (int, float)) or not np.isfinite(value) or value < 0.0:
            raise ValueError(f"P-T setting {field} must be finite and non-negative")
    if configuration.get("misfit_form") not in {
        "weighted_lsq",
        "logit_lsq",
        "alr_lsq",
        "dirichlet_nll",
        "chi_square",
        "raw_lsq",
    }:
        raise ValueError("invalid P-T misfit form")
    scale = configuration.get("theoretical_covariance_scale", 1.0)
    if not isinstance(scale, (int, float)) or not np.isfinite(scale) or scale < 0.0:
        raise ValueError(
            "P-T theoretical covariance scale must be finite and non-negative"
        )
    corrections = configuration.get("corrections_j_mol", {})
    if not isinstance(corrections, dict):
        raise ValueError("P-T corrections must be a JSON object")
    if any(
        not isinstance(key, str)
        or not key
        or not isinstance(value, (int, float))
        or not np.isfinite(value)
        for key, value in corrections.items()
    ):
        raise ValueError(
            "P-T corrections must have finite numeric values and string keys"
        )
    parameter_order = configuration.get("parameter_order")
    has_calibration_bias = configuration.get("calibration_bias", False)
    if not isinstance(has_calibration_bias, bool):
        raise ValueError("P-T calibration_bias flag must be boolean")
    expected_count = 4 if has_calibration_bias else 2
    parameters = configuration.get("parameters")
    if (
        not isinstance(parameter_order, list)
        or len(parameter_order) != expected_count
        or len(set(parameter_order)) != expected_count
    ):
        raise ValueError("invalid P-T parameter order in configuration archive")
    if not isinstance(parameters, list) or len(parameters) != expected_count:
        raise ValueError("P-T parameter metadata count does not match mode")
    if any(not isinstance(parameter, dict) for parameter in parameters):
        raise ValueError("P-T parameter metadata must be JSON objects")
    if [parameter.get("name") for parameter in parameters] != parameter_order:
        raise ValueError("P-T parameter metadata names do not match order")
    if parameter_order[:2] != ["pressure_bar", "temperature_c"]:
        raise ValueError("P-T parameter order must begin with pressure and temperature")
    for index, parameter in enumerate(parameters):
        kind = parameter.get("parameter_kind", "thermodynamic")
        if kind == "scalar":
            if not isinstance(parameter.get("unit"), str) or not parameter["unit"]:
                raise ValueError("P-T scalar parameter units must be non-empty strings")
            mean_field, sigma_field = "prior_mean", "prior_sigma"
            lower_field, upper_field = "lower_bound", "upper_bound"
        elif kind == "thermodynamic":
            mean_field, sigma_field = "prior_mean_j_mol", "prior_sigma_j_mol"
            lower_field, upper_field = "lower_bound_j_mol", "upper_bound_j_mol"
        else:
            raise ValueError("unsupported P-T parameter kind")
        if has_calibration_bias and index >= 2:
            if kind != "scalar":
                raise ValueError("P-T calibration bias parameters must be scalar")
            expected_unit = "bar" if index == 2 else {"degC", "C", "celsius"}
            if (
                parameter.get("unit") != expected_unit
                if isinstance(expected_unit, str)
                else parameter.get("unit") not in expected_unit
            ):
                raise ValueError("P-T calibration bias parameter has an invalid unit")
        for field in (mean_field, sigma_field):
            value = parameter.get(field)
            if not isinstance(value, (int, float)) or not np.isfinite(value):
                raise ValueError(f"P-T parameter field {field} must be finite")
        if parameter[sigma_field] <= 0.0:
            raise ValueError("P-T parameter prior sigmas must be positive")
        lower = parameter.get(lower_field)
        upper = parameter.get(upper_field)
        if lower is not None and not isinstance(lower, (int, float)):
            raise ValueError("P-T parameter lower bound must be numeric or null")
        if upper is not None and not isinstance(upper, (int, float)):
            raise ValueError("P-T parameter upper bound must be numeric or null")
        if lower is not None and upper is not None and lower >= upper:
            raise ValueError("P-T parameter bounds must increase")
    mode = configuration.get("problem_mode")
    if mode not in {"composition", "bracket"}:
        raise ValueError("invalid P-T problem mode in configuration archive")
    bracket_ids = configuration.get("bracket_ids", [])
    if not isinstance(bracket_ids, list) or len(set(bracket_ids)) != len(bracket_ids):
        raise ValueError("P-T bracket IDs must be unique")
    if (mode == "composition" and bracket_ids) or (
        mode == "bracket" and not bracket_ids
    ):
        raise ValueError("P-T problem mode is inconsistent with bracket IDs")
    covariance_payload = configuration.get(
        "prior_covariance", configuration.get("prior_covariance_j_mol2")
    )
    if covariance_payload is None:
        raise ValueError("P-T configuration is missing prior_covariance")
    covariance = np.asarray(covariance_payload, dtype=float)
    if covariance.shape != (expected_count, expected_count) or not np.all(
        np.isfinite(covariance)
    ):
        raise ValueError("invalid P-T prior covariance in configuration archive")
    covariance = 0.5 * (covariance + covariance.T)
    if np.any(np.linalg.eigvalsh(covariance) <= 0.0):
        raise ValueError("P-T prior covariance must be positive definite")
    experiments = json.loads(
        (input_path / "experiment.json").read_text(encoding="utf-8")
    )
    brackets = json.loads((input_path / "brackets.json").read_text(encoding="utf-8"))
    if not isinstance(experiments, dict):
        raise ValueError("P-T experiment must be a JSON object")
    if not isinstance(experiments.get("id"), str) or not experiments["id"]:
        raise ValueError("P-T experiment must have a non-empty string id")
    bracket_records = _validate_records(brackets, "P-T brackets")
    _validate_bracket_records(bracket_records, "P-T brackets")
    if [record["id"] for record in bracket_records] != bracket_ids:
        raise ValueError("P-T bracket IDs do not match archived bracket records")
    return {
        "configuration": configuration,
        "experiment": experiments,
        "brackets": brackets,
    }


def write_grouped_pt_calibration_configuration(
    problem: GroupedPTCalibrationProblem, output_dir: str | Path
) -> None:
    """Archive the grouped calibration hierarchy and its experiment metadata."""
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    parameter_records = [
        {
            "parameter_kind": "thermodynamic",
            "name": parameter.name,
            "phase": parameter.phase,
            "field": parameter.field,
            "parameter_type": parameter.parameter_type,
            "key": parameter.key,
            "unit": parameter.unit,
            "prior_mean": parameter.prior_mean,
            "prior_sigma": parameter.prior_sigma,
            "lower_bound": parameter.lower_bound,
            "upper_bound": parameter.upper_bound,
        }
        for parameter in problem.thermodynamic_parameters
    ]
    concentration_records = [
        {
            "parameter_kind": "dirichlet_concentration",
            "name": parameter.name,
            "phase": parameter.phase,
            "key": parameter.key,
            "unit": parameter.unit,
            "prior_mean": parameter.prior_mean,
            "prior_sigma": parameter.prior_sigma,
            "lower_bound": parameter.lower_bound,
            "upper_bound": parameter.upper_bound,
        }
        for parameter in problem.dirichlet_concentration_parameters
    ]
    parameter_records.extend(concentration_records)
    parameter_records.extend(
        {
            "parameter_kind": "grouped_calibration",
            "name": parameter.name,
            "unit": parameter.unit,
            "prior_mean": parameter.prior_mean,
            "prior_sigma": parameter.prior_sigma,
            "lower_bound": parameter.lower_bound,
            "upper_bound": parameter.upper_bound,
        }
        for parameter in problem.parameters[
            len(problem.thermodynamic_parameters)
            + len(problem.dirichlet_concentration_parameters) :
        ]
    )
    configuration = {
        "forward_model": _forward_model_provenance(problem.forward_model),
        "metadata_key": problem.metadata_key,
        "groups": list(problem.groups),
        "parameter_order": [parameter.name for parameter in problem.parameters],
        "parameter_units": [parameter.unit for parameter in problem.parameters],
        "parameters": parameter_records,
        "thermodynamic_parameter_count": len(problem.thermodynamic_parameters),
        "dirichlet_concentration_parameter_count": len(
            problem.dirichlet_concentration_parameters
        ),
        "thermodynamic_prior_covariance": problem.thermodynamic_prior_covariance.tolist(),
        "observation_bias_keys": list(problem.observation_bias_keys),
        "observation_bias_prior_sigma": problem.observation_bias_prior_sigma,
        "observation_bias_bounds": list(problem.observation_bias_bounds),
        "model_error_scale_groups": list(problem.model_error_scale_groups),
        "model_error_scale_prior_sigma": problem.model_error_scale_prior_sigma,
        "model_error_scale_bounds": list(problem.model_error_scale_bounds),
        "pressure_scale_hyperprior": problem.pressure_scale_hyperprior,
        "temperature_scale_hyperprior": problem.temperature_scale_hyperprior,
        "missing_phase_penalty": problem.missing_phase_penalty,
        "extra_phase_penalty": problem.extra_phase_penalty,
        "misfit_form": problem.misfit_form,
        "theoretical_covariance_scale": problem.theoretical_covariance_scale,
        "corrections": problem.corrections,
    }
    (output / "grouped_pt_calibration_configuration.json").write_text(
        json.dumps(configuration, indent=2), encoding="utf-8"
    )
    (output / "experiments.json").write_text(
        json.dumps(
            [experiment_to_dict(experiment) for experiment in problem.experiments],
            indent=2,
        ),
        encoding="utf-8",
    )


def read_grouped_pt_calibration_configuration(
    input_dir: str | Path,
) -> dict[str, object]:
    """Read and validate grouped P-T calibration provenance."""
    input_path = Path(input_dir)
    configuration = json.loads(
        (input_path / "grouped_pt_calibration_configuration.json").read_text(
            encoding="utf-8"
        )
    )
    if not isinstance(configuration, dict):
        raise ValueError("grouped P-T calibration configuration must be an object")
    metadata_key = configuration.get("metadata_key")
    groups = configuration.get("groups")
    order = configuration.get("parameter_order")
    units = configuration.get("parameter_units")
    parameters = configuration.get("parameters")
    if not isinstance(metadata_key, str) or not metadata_key:
        raise ValueError("grouped P-T calibration metadata_key must be non-empty")
    if (
        not isinstance(groups, list)
        or not groups
        or any(not isinstance(group, str) or not group for group in groups)
        or len(set(groups)) != len(groups)
    ):
        raise ValueError("grouped P-T calibration groups must be unique strings")
    thermo_count = configuration.get("thermodynamic_parameter_count", 0)
    if not isinstance(thermo_count, int) or thermo_count < 0:
        raise ValueError("thermodynamic_parameter_count must be non-negative")
    if parameters is None and thermo_count == 0:
        parameters = [
            {
                "parameter_kind": "grouped_calibration",
                "name": name,
                "unit": unit,
            }
            for name, unit in zip(order or (), units or ())
        ]
    bias_keys = configuration.get("observation_bias_keys", [])
    model_scale_groups = configuration.get("model_error_scale_groups", [])
    concentration_count = configuration.get(
        "dirichlet_concentration_parameter_count", 0
    )
    if not isinstance(concentration_count, int) or concentration_count < 0:
        raise ValueError("dirichlet_concentration_parameter_count must be non-negative")
    if not isinstance(bias_keys, list) or not isinstance(model_scale_groups, list):
        raise ValueError("grouped nuisance group/key metadata must be arrays")
    if any(not isinstance(key, str) or not key for key in bias_keys):
        raise ValueError("grouped observation-bias keys must be non-empty strings")
    if any(
        not isinstance(group, str) or group not in groups
        for group in model_scale_groups
    ):
        raise ValueError(
            "grouped model-error scale groups must match experiment groups"
        )
    expected_parameter_count = (
        thermo_count
        + concentration_count
        + 2 * len(groups)
        + len(bias_keys) * len(groups)
        + len(model_scale_groups)
        + 2
    )
    if (
        not isinstance(parameters, list)
        or len(set(bias_keys)) != len(bias_keys)
        or len(set(model_scale_groups)) != len(model_scale_groups)
        or len(parameters) != expected_parameter_count
        or any(not isinstance(parameter, dict) for parameter in parameters)
    ):
        raise ValueError("grouped P-T calibration parameter records are invalid")
    expected_order = [
        *(parameter["name"] for parameter in parameters[:thermo_count]),
        *(
            parameter["name"]
            for parameter in parameters[
                thermo_count : thermo_count + concentration_count
            ]
        ),
        *(f"pressure_bias[{group}]" for group in groups),
        *(f"temperature_bias[{group}]" for group in groups),
        *(f"bias[{key}|{group}]" for key in bias_keys for group in groups),
        *(f"scale_C_T[{group}]" for group in model_scale_groups),
        "pressure_bias_scale",
        "temperature_bias_scale",
    ]
    expected_units = [
        *(parameter.get("unit") for parameter in parameters[:thermo_count]),
        *(
            parameter.get("unit")
            for parameter in parameters[
                thermo_count : thermo_count + concentration_count
            ]
        ),
        *("bar" for _ in groups),
        *("degC" for _ in groups),
        *(
            parameter.get("unit")
            for parameter in parameters[
                thermo_count + concentration_count + 2 * len(groups) : thermo_count
                + concentration_count
                + 2 * len(groups)
                + len(bias_keys) * len(groups)
            ]
        ),
        *("1" for _ in model_scale_groups),
        "bar",
        "degC",
    ]
    if order != expected_order or units != expected_units:
        raise ValueError("grouped P-T calibration parameter metadata is inconsistent")
    concentration_records = parameters[
        thermo_count : thermo_count + concentration_count
    ]
    for parameter in concentration_records:
        phase = parameter.get("phase")
        lower = parameter.get("lower_bound")
        upper = parameter.get("upper_bound")
        if (
            parameter.get("parameter_kind") != "dirichlet_concentration"
            or not isinstance(phase, str)
            or not phase
            or parameter.get("key") != f"dirichlet_concentration:{phase}"
            or parameter.get("unit") != "1"
        ):
            raise ValueError("invalid grouped Dirichlet concentration parameter")
        for field in ("prior_mean", "prior_sigma", "lower_bound", "upper_bound"):
            value = parameter.get(field)
            if not isinstance(value, (int, float)) or not np.isfinite(value):
                raise ValueError(
                    f"grouped Dirichlet concentration {field} must be finite"
                )
        if (
            parameter["prior_mean"] <= 0.0
            or parameter["prior_sigma"] <= 0.0
            or lower <= 0.0
            or lower >= upper
            or not lower <= parameter["prior_mean"] <= upper
        ):
            raise ValueError("invalid grouped Dirichlet concentration prior bounds")
    for parameter in parameters[thermo_count : thermo_count + concentration_count]:
        phase = parameter.get("phase")
        if (
            parameter.get("parameter_kind") != "dirichlet_concentration"
            or not isinstance(phase, str)
            or not phase
            or parameter.get("key") != f"dirichlet_concentration:{phase}"
            or parameter.get("unit") != "1"
        ):
            raise ValueError(
                "invalid grouped Dirichlet concentration parameter metadata"
            )
        for field in ("prior_mean", "prior_sigma", "lower_bound", "upper_bound"):
            value = parameter.get(field)
            if not isinstance(value, (int, float)) or not np.isfinite(value):
                raise ValueError(
                    f"grouped Dirichlet concentration {field} must be finite"
                )
        if (
            parameter["prior_mean"] <= 0.0
            or parameter["prior_sigma"] <= 0.0
            or parameter["lower_bound"] <= 0.0
            or parameter["lower_bound"] >= parameter["upper_bound"]
            or not parameter["lower_bound"]
            <= parameter["prior_mean"]
            <= parameter["upper_bound"]
        ):
            raise ValueError("invalid grouped Dirichlet concentration prior or bounds")
    thermo_covariance = np.asarray(
        configuration.get("thermodynamic_prior_covariance", []), dtype=float
    ).reshape(thermo_count, thermo_count)
    if thermo_covariance.shape != (thermo_count, thermo_count) or not np.all(
        np.isfinite(thermo_covariance)
    ):
        raise ValueError(
            "thermodynamic prior covariance does not match parameter count"
        )
    if thermo_count and np.any(np.linalg.eigvalsh(thermo_covariance) <= 0.0):
        raise ValueError("thermodynamic prior covariance must be positive definite")
    for key in ("pressure_scale_hyperprior", "temperature_scale_hyperprior"):
        value = configuration.get(key)
        if (
            not isinstance(value, (int, float))
            or not np.isfinite(value)
            or value <= 0.0
        ):
            raise ValueError(f"{key} must be finite and positive")
    bias_sigma = configuration.get("observation_bias_prior_sigma", 0.1)
    bias_bounds = configuration.get("observation_bias_bounds", [-1.0, 1.0])
    if (
        not isinstance(bias_sigma, (int, float))
        or not np.isfinite(bias_sigma)
        or bias_sigma <= 0.0
        or not isinstance(bias_bounds, list)
        or len(bias_bounds) != 2
        or any(
            not isinstance(value, (int, float)) or not np.isfinite(value)
            for value in bias_bounds
        )
        or bias_bounds[0] >= bias_bounds[1]
    ):
        raise ValueError("grouped observation-bias prior metadata is invalid")
    scale_sigma = configuration.get("model_error_scale_prior_sigma")
    scale_bounds = configuration.get("model_error_scale_bounds", [1.0e-3, 5.0])
    if model_scale_groups and (
        not isinstance(scale_sigma, (int, float))
        or not np.isfinite(scale_sigma)
        or scale_sigma <= 0.0
    ):
        raise ValueError("grouped model-error scale prior sigma must be positive")
    if (
        not isinstance(scale_bounds, list)
        or len(scale_bounds) != 2
        or any(
            not isinstance(value, (int, float)) or not np.isfinite(value)
            for value in scale_bounds
        )
        or scale_bounds[0] <= 0.0
        or scale_bounds[0] >= scale_bounds[1]
    ):
        raise ValueError("grouped model-error scale bounds are invalid")
    experiments = json.loads(
        (input_path / "experiments.json").read_text(encoding="utf-8")
    )
    _validate_records(experiments, "grouped P-T calibration experiments")
    _validate_experiment_records(experiments, "grouped P-T calibration experiments")
    observed_groups = {
        item.get("metadata", {}).get(metadata_key) for item in experiments
    }
    if observed_groups != set(groups):
        raise ValueError("archived experiment groups do not match configuration groups")
    return {"configuration": configuration, "experiments": experiments}


def write_input_perturbation_result(
    result: InputPerturbationResult, output_dir: str | Path
) -> None:
    """Write perturbation samples/objectives and covariance provenance."""
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    np.savetxt(
        output / "corrections_samples.csv",
        result.corrections_samples,
        delimiter=",",
        header=",".join(result.parameter_names),
        comments="",
    )
    np.savetxt(output / "objectives.csv", result.objectives, delimiter=",")
    summary = {
        "parameter_names": list(result.parameter_names),
        "source": result.source,
        "correction_mean": result.correction_mean,
        "correction_std": result.correction_std,
        "success_count": result.success_count,
        "thermodynamic_covariance_keys": list(result.thermodynamic_covariance_keys),
        "thermodynamic_covariance_diagnostics": result.thermodynamic_covariance_diagnostics,
    }
    (output / "summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )


def read_input_perturbation_result(input_dir: str | Path) -> InputPerturbationResult:
    """Reload an `InputPerturbationResult` written by its archive helper."""
    input_path = Path(input_dir)
    summary = json.loads((input_path / "summary.json").read_text(encoding="utf-8"))
    parameter_names = tuple(summary["parameter_names"])
    samples = np.loadtxt(
        input_path / "corrections_samples.csv", delimiter=",", skiprows=1
    ).reshape(-1, len(parameter_names))
    objectives = np.atleast_1d(np.loadtxt(input_path / "objectives.csv", delimiter=","))
    return InputPerturbationResult(
        parameter_names=parameter_names,
        source=summary["source"],
        corrections_samples=samples,
        objectives=objectives,
        correction_mean=summary["correction_mean"],
        correction_std=summary["correction_std"],
        success_count=summary["success_count"],
        thermodynamic_covariance_keys=tuple(
            summary.get("thermodynamic_covariance_keys", ())
        ),
        thermodynamic_covariance_diagnostics=summary.get(
            "thermodynamic_covariance_diagnostics"
        ),
    )


def write_mcmc_result(result: MCMCResult, output_dir: str | Path) -> None:
    """Write MCMC samples, log-posterior values, and summary statistics."""
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    np.savetxt(
        output / "samples.csv",
        result.samples,
        delimiter=",",
        header=",".join(result.parameter_names),
        comments="",
    )
    np.savetxt(output / "log_posterior.csv", result.log_posterior, delimiter=",")
    summary = {
        "parameter_names": list(result.parameter_names),
        "acceptance_rate": result.acceptance_rate,
        "posterior_mean": result.posterior_mean,
        "posterior_std": result.posterior_std,
        "credible_interval_95": {
            name: list(interval)
            for name, interval in result.credible_interval_95.items()
        },
        "uncertainty_source": result.uncertainty_source,
        "proposal_covariance_j_mol2": result.proposal_covariance_j_mol2.tolist()
        if result.proposal_covariance_j_mol2 is not None
        else None,
        "proposal_covariance": result.proposal_covariance.tolist()
        if result.proposal_covariance is not None
        else None,
        "seed": result.seed,
        "burn_in": result.burn_in,
        "parameter_units": list(result.parameter_units),
    }
    (output / "summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )


def read_mcmc_result(input_dir: str | Path) -> MCMCResult:
    """Reload an `MCMCResult` previously written by `write_mcmc_result`."""
    input_path = Path(input_dir)
    summary = json.loads((input_path / "summary.json").read_text(encoding="utf-8"))
    parameter_names = tuple(summary["parameter_names"])
    if not parameter_names or any(
        not isinstance(name, str) or not name for name in parameter_names
    ):
        raise ValueError("MCMC parameter names must be non-empty strings")
    units = tuple(summary.get("parameter_units", ()))
    if units and (
        len(units) != len(parameter_names)
        or any(not isinstance(unit, str) or not unit for unit in units)
    ):
        raise ValueError("MCMC parameter units must match parameter names")
    samples = np.loadtxt(input_path / "samples.csv", delimiter=",", skiprows=1).reshape(
        -1, len(parameter_names)
    )
    log_posterior = np.atleast_1d(
        np.loadtxt(input_path / "log_posterior.csv", delimiter=",")
    )
    return MCMCResult(
        parameter_names,
        samples,
        log_posterior,
        summary["acceptance_rate"],
        summary["posterior_mean"],
        summary["posterior_std"],
        {
            name: tuple(interval)
            for name, interval in summary["credible_interval_95"].items()
        },
        summary.get("uncertainty_source", "all"),
        np.asarray(
            summary.get(
                "proposal_covariance", summary.get("proposal_covariance_j_mol2")
            ),
            dtype=float,
        )
        if summary.get("proposal_covariance", summary.get("proposal_covariance_j_mol2"))
        is not None
        else None,
        summary.get("seed"),
        summary.get("burn_in", 0),
        units,
    )


def write_hierarchical_mcmc_result(
    result: HierarchicalMCMCResult, output_dir: str | Path
) -> None:
    """Write joint correction/model-error-scale MCMC samples and metadata."""
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    np.savetxt(
        output / "samples.csv",
        result.samples,
        delimiter=",",
        header=",".join(result.parameter_names),
        comments="",
    )
    np.savetxt(output / "log_posterior.csv", result.log_posterior, delimiter=",")
    summary = {
        "parameter_names": list(result.parameter_names),
        "parameter_units": list(result.parameter_units),
        "acceptance_rate": result.acceptance_rate,
        "posterior_mean": result.posterior_mean,
        "posterior_std": result.posterior_std,
        "credible_interval_95": {
            name: list(interval)
            for name, interval in result.credible_interval_95.items()
        },
        "scale_name": result.scale_name,
        "scale_prior_mean": result.scale_prior_mean,
        "scale_prior_sigma": result.scale_prior_sigma,
        "scale_bounds": list(result.scale_bounds),
        "seed": result.seed,
        "burn_in": result.burn_in,
    }
    (output / "summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )


def read_hierarchical_mcmc_result(
    input_dir: str | Path,
) -> HierarchicalMCMCResult:
    """Reload a hierarchical MCMC result written by its archive helper."""
    input_path = Path(input_dir)
    summary = json.loads((input_path / "summary.json").read_text(encoding="utf-8"))
    names = tuple(summary["parameter_names"])
    units = tuple(summary["parameter_units"])
    if (
        not names
        or len(units) != len(names)
        or any(
            not isinstance(name, str)
            or not name
            or not isinstance(unit, str)
            or not unit
            for name, unit in zip(names, units)
        )
    ):
        raise ValueError(
            "hierarchical MCMC names and units must be non-empty and aligned"
        )
    samples = np.loadtxt(input_path / "samples.csv", delimiter=",", skiprows=1).reshape(
        -1, len(names)
    )
    log_posterior = np.atleast_1d(
        np.loadtxt(input_path / "log_posterior.csv", delimiter=",")
    )
    bounds = tuple(float(value) for value in summary["scale_bounds"])
    if len(bounds) != 2 or bounds[0] < 0.0 or bounds[0] >= bounds[1]:
        raise ValueError("invalid archived scale bounds")
    return HierarchicalMCMCResult(
        names,
        units,
        samples,
        log_posterior,
        float(summary["acceptance_rate"]),
        {name: float(value) for name, value in summary["posterior_mean"].items()},
        {name: float(value) for name, value in summary["posterior_std"].items()},
        {
            name: (float(interval[0]), float(interval[1]))
            for name, interval in summary["credible_interval_95"].items()
        },
        str(summary["scale_name"]),
        float(summary["scale_prior_mean"]),
        float(summary["scale_prior_sigma"]),
        bounds,
        int(summary["seed"]),
        int(summary["burn_in"]),
    )


def write_grouped_hierarchical_mcmc_result(
    result: GroupedHierarchicalMCMCResult, output_dir: str | Path
) -> None:
    """Write joint correction/group-scale MCMC samples and metadata."""
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    np.savetxt(
        output / "samples.csv",
        result.samples,
        delimiter=",",
        header=",".join(result.parameter_names),
        comments="",
    )
    np.savetxt(output / "log_posterior.csv", result.log_posterior, delimiter=",")
    summary = {
        "parameter_names": list(result.parameter_names),
        "parameter_units": list(result.parameter_units),
        "acceptance_rate": result.acceptance_rate,
        "posterior_mean": result.posterior_mean,
        "posterior_std": result.posterior_std,
        "credible_interval_95": {
            name: list(interval)
            for name, interval in result.credible_interval_95.items()
        },
        "metadata_key": result.metadata_key,
        "scale_names": list(result.scale_names),
        "scale_prior_means": result.scale_prior_means,
        "scale_prior_sigmas": result.scale_prior_sigmas,
        "scale_bounds": list(result.scale_bounds),
        "seed": result.seed,
        "burn_in": result.burn_in,
    }
    (output / "summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )


def read_grouped_hierarchical_mcmc_result(
    input_dir: str | Path,
) -> GroupedHierarchicalMCMCResult:
    """Reload a grouped hierarchical MCMC archive."""
    input_path = Path(input_dir)
    summary = json.loads((input_path / "summary.json").read_text(encoding="utf-8"))
    names = tuple(summary["parameter_names"])
    units = tuple(summary["parameter_units"])
    if not names or len(names) != len(units):
        raise ValueError("grouped hierarchical names and units must be aligned")
    samples = np.loadtxt(input_path / "samples.csv", delimiter=",", skiprows=1).reshape(
        -1, len(names)
    )
    log_posterior = np.atleast_1d(
        np.loadtxt(input_path / "log_posterior.csv", delimiter=",")
    )
    bounds = tuple(float(value) for value in summary["scale_bounds"])
    if len(bounds) != 2 or bounds[0] < 0.0 or bounds[0] >= bounds[1]:
        raise ValueError("invalid grouped scale bounds")
    return GroupedHierarchicalMCMCResult(
        names,
        units,
        samples,
        log_posterior,
        float(summary["acceptance_rate"]),
        {name: float(value) for name, value in summary["posterior_mean"].items()},
        {name: float(value) for name, value in summary["posterior_std"].items()},
        {
            name: (float(interval[0]), float(interval[1]))
            for name, interval in summary["credible_interval_95"].items()
        },
        str(summary["metadata_key"]),
        tuple(summary["scale_names"]),
        {str(key): float(value) for key, value in summary["scale_prior_means"].items()},
        {
            str(key): float(value)
            for key, value in summary["scale_prior_sigmas"].items()
        },
        bounds,
        int(summary["seed"]),
        int(summary["burn_in"]),
    )


def write_grouped_observation_bias_mcmc_result(
    result: GroupedObservationBiasMCMCResult, output_dir: str | Path
) -> None:
    """Write joint correction/group-bias MCMC samples and metadata."""
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    np.savetxt(
        output / "samples.csv",
        result.samples,
        delimiter=",",
        header=",".join(result.parameter_names),
        comments="",
    )
    np.savetxt(output / "log_posterior.csv", result.log_posterior, delimiter=",")
    summary = {
        "parameter_names": list(result.parameter_names),
        "parameter_units": list(result.parameter_units),
        "acceptance_rate": result.acceptance_rate,
        "posterior_mean": result.posterior_mean,
        "posterior_std": result.posterior_std,
        "credible_interval_95": {
            name: list(interval)
            for name, interval in result.credible_interval_95.items()
        },
        "metadata_key": result.metadata_key,
        "observation_key": result.observation_key,
        "observation_keys": list(result.observation_keys or (result.observation_key,)),
        "bias_names": list(result.bias_names),
        "bias_prior_mean": result.bias_prior_mean,
        "bias_prior_sigma": result.bias_prior_sigma,
        "bias_bounds": list(result.bias_bounds),
        "seed": result.seed,
        "burn_in": result.burn_in,
    }
    (output / "summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )


def read_grouped_observation_bias_mcmc_result(
    input_dir: str | Path,
) -> GroupedObservationBiasMCMCResult:
    """Reload a grouped observation-bias MCMC archive."""
    input_path = Path(input_dir)
    summary = json.loads((input_path / "summary.json").read_text(encoding="utf-8"))
    names = tuple(summary["parameter_names"])
    units = tuple(summary["parameter_units"])
    if not names or len(names) != len(units):
        raise ValueError("grouped bias names and units must be aligned")
    samples = np.loadtxt(input_path / "samples.csv", delimiter=",", skiprows=1).reshape(
        -1, len(names)
    )
    log_posterior = np.atleast_1d(
        np.loadtxt(input_path / "log_posterior.csv", delimiter=",")
    )
    bounds = tuple(float(value) for value in summary["bias_bounds"])
    if len(bounds) != 2 or bounds[0] >= bounds[1]:
        raise ValueError("invalid archived bias bounds")
    observation_keys = tuple(
        summary.get("observation_keys", (summary["observation_key"],))
    )
    if not observation_keys or any(
        not isinstance(key, str) or not key for key in observation_keys
    ):
        raise ValueError("grouped bias observation_keys must be non-empty strings")
    return GroupedObservationBiasMCMCResult(
        names,
        units,
        samples,
        log_posterior,
        float(summary["acceptance_rate"]),
        {name: float(value) for name, value in summary["posterior_mean"].items()},
        {name: float(value) for name, value in summary["posterior_std"].items()},
        {
            name: (float(interval[0]), float(interval[1]))
            for name, interval in summary["credible_interval_95"].items()
        },
        str(summary["metadata_key"]),
        str(summary["observation_key"]),
        tuple(summary["bias_names"]),
        float(summary["bias_prior_mean"]),
        float(summary["bias_prior_sigma"]),
        bounds,
        int(summary["seed"]),
        int(summary["burn_in"]),
        observation_keys,
    )


def write_combined_grouped_uncertainty_mcmc_result(
    result: CombinedGroupedUncertaintyMCMCResult, output_dir: str | Path
) -> None:
    """Write joint group-scale and group-bias MCMC samples."""
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    np.savetxt(
        output / "samples.csv",
        result.samples,
        delimiter=",",
        header=",".join(result.parameter_names),
        comments="",
    )
    np.savetxt(output / "log_posterior.csv", result.log_posterior, delimiter=",")
    summary = {
        "parameter_names": list(result.parameter_names),
        "parameter_units": list(result.parameter_units),
        "acceptance_rate": result.acceptance_rate,
        "posterior_mean": result.posterior_mean,
        "posterior_std": result.posterior_std,
        "credible_interval_95": {
            name: list(interval)
            for name, interval in result.credible_interval_95.items()
        },
        "metadata_key": result.metadata_key,
        "observation_key": result.observation_key,
        "groups": list(result.groups),
        "observation_keys": list(result.observation_keys or (result.observation_key,)),
        "bias_labels": [list(label) for label in result.bias_labels],
        "scale_bounds": list(result.scale_bounds),
        "scale_prior_mean": result.scale_prior_mean,
        "scale_prior_sigma": result.scale_prior_sigma,
        "bias_bounds": list(result.bias_bounds),
        "bias_prior_mean": result.bias_prior_mean,
        "bias_prior_sigma": result.bias_prior_sigma,
        "seed": result.seed,
        "burn_in": result.burn_in,
    }
    (output / "summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )


def read_combined_grouped_uncertainty_mcmc_result(
    input_dir: str | Path,
) -> CombinedGroupedUncertaintyMCMCResult:
    """Reload a combined grouped uncertainty MCMC archive."""
    input_path = Path(input_dir)
    summary = json.loads((input_path / "summary.json").read_text(encoding="utf-8"))
    names = tuple(summary["parameter_names"])
    units = tuple(summary["parameter_units"])
    if not names or len(names) != len(units):
        raise ValueError("combined hierarchy names and units must be aligned")
    samples = np.loadtxt(input_path / "samples.csv", delimiter=",", skiprows=1).reshape(
        -1, len(names)
    )
    log_posterior = np.atleast_1d(
        np.loadtxt(input_path / "log_posterior.csv", delimiter=",")
    )
    return CombinedGroupedUncertaintyMCMCResult(
        names,
        units,
        samples,
        log_posterior,
        float(summary["acceptance_rate"]),
        {name: float(value) for name, value in summary["posterior_mean"].items()},
        {name: float(value) for name, value in summary["posterior_std"].items()},
        {
            name: (float(interval[0]), float(interval[1]))
            for name, interval in summary["credible_interval_95"].items()
        },
        str(summary["metadata_key"]),
        str(summary["observation_key"]),
        tuple(summary["groups"]),
        (float(summary["scale_bounds"][0]), float(summary["scale_bounds"][1])),
        float(summary["scale_prior_mean"]),
        float(summary["scale_prior_sigma"]),
        (float(summary["bias_bounds"][0]), float(summary["bias_bounds"][1])),
        float(summary["bias_prior_mean"]),
        float(summary["bias_prior_sigma"]),
        int(summary["seed"]),
        int(summary["burn_in"]),
        tuple(summary.get("observation_keys", (summary["observation_key"],))),
        tuple(tuple(label) for label in summary.get("bias_labels", ())),
    )


def write_parallel_tempering_result(
    result: ParallelTemperingResult, output_dir: str | Path
) -> None:
    """Write cold-chain samples and parallel-tempering diagnostics."""
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    np.savetxt(
        output / "samples.csv",
        result.samples,
        delimiter=",",
        header=",".join(result.parameter_names),
        comments="",
    )
    np.savetxt(output / "log_posterior.csv", result.log_posterior, delimiter=",")
    np.savetxt(output / "mode_labels.csv", result.mode_labels, delimiter=",")
    if result.all_temperature_samples is not None:
        np.save(output / "all_temperature_samples.npy", result.all_temperature_samples)
    if result.temperature_mode_labels is not None:
        np.save(output / "temperature_mode_labels.npy", result.temperature_mode_labels)
    if result.all_temperature_objectives is not None:
        np.save(
            output / "all_temperature_objectives.npy", result.all_temperature_objectives
        )
    if result.uniform_reference_objectives is not None:
        np.save(
            output / "uniform_reference_objectives.npy",
            result.uniform_reference_objectives,
        )
    summary = {
        "parameter_names": list(result.parameter_names),
        "parameter_units": list(result.parameter_units),
        "temperatures": list(result.temperatures),
        "mode_occupancy": {
            str(key): value for key, value in result.mode_occupancy.items()
        },
        "mode_means": {str(key): value for key, value in result.mode_means.items()},
        "mode_best_objective": {
            str(key): value for key, value in result.mode_best_objective.items()
        },
        "acceptance_rates": list(result.acceptance_rates),
        "swap_acceptance_rates": list(result.swap_acceptance_rates),
        "crossed_modes": result.crossed_modes,
        "seed": result.seed,
        "burn_in": result.burn_in,
        "proposal_sigma": list(result.proposal_sigma),
        "temperature_mode_occupancy": {
            str(temperature): {str(mode): value for mode, value in occupancy.items()}
            for temperature, occupancy in (
                result.temperature_mode_occupancy or {}
            ).items()
        },
        "mode_transition_count": result.mode_transition_count,
        "mode_transition_rate": result.mode_transition_rate,
        "has_global_temperature_diagnostics": result.all_temperature_samples
        is not None,
        "has_thermodynamic_integration_diagnostics": result.all_temperature_objectives
        is not None
        and result.uniform_reference_objectives is not None,
    }
    (output / "summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )


def read_parallel_tempering_result(
    input_dir: str | Path,
) -> ParallelTemperingResult:
    """Reload a result written by `write_parallel_tempering_result`."""
    input_path = Path(input_dir)
    summary = json.loads((input_path / "summary.json").read_text(encoding="utf-8"))
    names = tuple(summary["parameter_names"])
    units = tuple(summary["parameter_units"])
    if not names or len(units) != len(names):
        raise ValueError("parallel-tempering names and units must be aligned")
    samples = np.loadtxt(input_path / "samples.csv", delimiter=",", skiprows=1).reshape(
        -1, len(names)
    )
    log_posterior = np.atleast_1d(
        np.loadtxt(input_path / "log_posterior.csv", delimiter=",")
    )
    raw_labels = np.atleast_1d(
        np.loadtxt(input_path / "mode_labels.csv", delimiter=",")
    )
    if not np.all(np.isfinite(raw_labels)) or not np.all(
        raw_labels == np.floor(raw_labels)
    ):
        raise ValueError("parallel-tempering mode labels must be integers")
    labels = raw_labels.astype(int)
    if len(samples) != len(log_posterior) or len(samples) != len(labels):
        raise ValueError("parallel-tempering sample arrays must have equal length")
    temperatures = tuple(float(value) for value in summary["temperatures"])
    if (
        not temperatures
        or temperatures[0] != 1.0
        or any(value < 1.0 for value in temperatures)
    ):
        raise ValueError("invalid parallel-tempering temperatures")
    all_temperature_samples = None
    temperature_mode_labels = None
    all_temperature_objectives = None
    uniform_reference_objectives = None
    if summary.get("has_global_temperature_diagnostics", False):
        all_temperature_samples = np.load(input_path / "all_temperature_samples.npy")
        temperature_mode_labels = np.load(input_path / "temperature_mode_labels.npy")
    if summary.get("has_thermodynamic_integration_diagnostics", False):
        all_temperature_objectives = np.load(
            input_path / "all_temperature_objectives.npy"
        )
        uniform_reference_objectives = np.load(
            input_path / "uniform_reference_objectives.npy"
        )
    return ParallelTemperingResult(
        names,
        units,
        temperatures,
        samples,
        log_posterior,
        labels,
        {int(key): float(value) for key, value in summary["mode_occupancy"].items()},
        {
            int(key): {name: float(value) for name, value in values.items()}
            for key, values in summary["mode_means"].items()
        },
        {
            int(key): float(value)
            for key, value in summary["mode_best_objective"].items()
        },
        tuple(float(value) for value in summary["acceptance_rates"]),
        tuple(float(value) for value in summary["swap_acceptance_rates"]),
        bool(summary["crossed_modes"]),
        int(summary["seed"]),
        int(summary["burn_in"]),
        tuple(float(value) for value in summary["proposal_sigma"]),
        {
            float(temperature): {
                int(mode): float(value) for mode, value in occupancy.items()
            }
            for temperature, occupancy in summary.get(
                "temperature_mode_occupancy", {}
            ).items()
        },
        int(summary.get("mode_transition_count", 0)),
        float(summary.get("mode_transition_rate", 0.0)),
        all_temperature_samples,
        temperature_mode_labels,
        all_temperature_objectives,
        uniform_reference_objectives,
    )


def write_multi_chain_mcmc_result(
    result: MultiChainMCMCResult, output_dir: str | Path
) -> None:
    """Archive every MCMC chain and its aggregate convergence diagnostics."""
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    for index, chain in enumerate(result.chains):
        write_mcmc_result(chain, output / f"chain_{index}")
    summary = {
        "parameter_names": list(result.parameter_names),
        "chain_count": len(result.chains),
        "r_hat": result.r_hat,
        "effective_sample_size": result.effective_sample_size,
        "initial_points": [list(point) for point in result.initial_points],
        "seeds": list(result.seeds),
        "parameter_units": list(result.parameter_units),
    }
    (output / "multi_chain_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )


def read_multi_chain_mcmc_result(input_dir: str | Path) -> MultiChainMCMCResult:
    """Reload a `MultiChainMCMCResult` written by its archive helper."""
    input_path = Path(input_dir)
    summary = json.loads(
        (input_path / "multi_chain_summary.json").read_text(encoding="utf-8")
    )
    chains = tuple(
        read_mcmc_result(input_path / f"chain_{index}")
        for index in range(int(summary["chain_count"]))
    )
    units = tuple(
        summary.get("parameter_units", chains[0].parameter_units if chains else ())
    )
    if units and len(units) != len(summary["parameter_names"]):
        raise ValueError("multi-chain parameter units must match parameter names")
    if any(not isinstance(unit, str) or not unit for unit in units):
        raise ValueError("multi-chain parameter units must be non-empty strings")
    if any(
        chain.parameter_units and chain.parameter_units != units for chain in chains
    ):
        raise ValueError("multi-chain archive units disagree with chain units")
    return MultiChainMCMCResult(
        tuple(summary["parameter_names"]),
        chains,
        {name: float(value) for name, value in summary["r_hat"].items()},
        {
            name: float(value)
            for name, value in summary["effective_sample_size"].items()
        },
        tuple(
            tuple(float(value) for value in point)
            for point in summary.get("initial_points", ())
        ),
        tuple(int(seed) for seed in summary.get("seeds", ())),
        units,
    )


def write_posterior_predictive_result(
    result: PosteriorPredictiveResult, output_dir: str | Path
) -> None:
    """Write a continuous posterior-predictive result as JSON."""
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    payload = {
        "observable_names": list(result.observable_names),
        "observed_values": result.observed_values.tolist(),
        "predictive_mean": result.predictive_mean.tolist(),
        "predictive_std": result.predictive_std.tolist(),
        "credible_interval_95": {
            name: list(interval)
            for name, interval in result.credible_interval_95.items()
        },
        "within_95": result.within_95,
        "include_observation_uncertainty": result.include_observation_uncertainty,
        "seed": result.seed,
        "sample_selection": result.sample_selection,
    }
    (output / "posterior_predictive.json").write_text(
        json.dumps(payload, indent=2), encoding="utf-8"
    )


def read_posterior_predictive_result(
    input_dir: str | Path,
) -> PosteriorPredictiveResult:
    """Reload a continuous posterior-predictive result archive."""
    payload = json.loads(
        (Path(input_dir) / "posterior_predictive.json").read_text(encoding="utf-8")
    )
    return PosteriorPredictiveResult(
        tuple(payload["observable_names"]),
        np.asarray(payload["observed_values"], dtype=float),
        np.asarray(payload["predictive_mean"], dtype=float),
        np.asarray(payload["predictive_std"], dtype=float),
        {
            name: tuple(interval)
            for name, interval in payload["credible_interval_95"].items()
        },
        {name: bool(value) for name, value in payload["within_95"].items()},
        bool(payload.get("include_observation_uncertainty", False)),
        payload.get("seed"),
        payload.get("sample_selection", "even"),
    )


def write_bracket_posterior_predictive_result(
    result: BracketPosteriorPredictiveResult, output_dir: str | Path
) -> None:
    """Write a bracket posterior-predictive result as JSON."""
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    payload = {
        "bracket_ids": list(result.bracket_ids),
        "samples_evaluated": result.samples_evaluated,
        "satisfied_fraction": result.satisfied_fraction,
        "mean_misfit": result.mean_misfit,
        "credible_interval_95": list(result.credible_interval_95),
        "misfits": list(result.misfits),
        "seed": result.seed,
        "sample_selection": result.sample_selection,
    }
    (output / "bracket_posterior_predictive.json").write_text(
        json.dumps(payload, indent=2), encoding="utf-8"
    )


def read_bracket_posterior_predictive_result(
    input_dir: str | Path,
) -> BracketPosteriorPredictiveResult:
    """Reload a bracket posterior-predictive result archive."""
    payload = json.loads(
        (Path(input_dir) / "bracket_posterior_predictive.json").read_text(
            encoding="utf-8"
        )
    )
    return BracketPosteriorPredictiveResult(
        tuple(payload["bracket_ids"]),
        int(payload["samples_evaluated"]),
        float(payload["satisfied_fraction"]),
        float(payload["mean_misfit"]),
        tuple(payload["credible_interval_95"]),
        tuple(float(value) for value in payload.get("misfits", ())),
        payload.get("seed"),
        payload.get("sample_selection", "even"),
    )


def write_observation_bias_profile_result(
    result: ObservationBiasProfileResult, output_dir: str | Path
) -> None:
    """Write a grouped observation-bias profile result as JSON."""
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    payload = {
        "biases": result.biases,
        "objective": result.objective,
        "corrections_j_mol": result.corrections_j_mol,
        "bias_prior_penalty": result.bias_prior_penalty,
        "metadata_key": result.metadata_key,
        "observation_key": result.observation_key,
        "records": list(result.records),
        "bias_bounds": list(result.bias_bounds),
        "bias_prior_mean": result.bias_prior_mean,
        "bias_prior_sigma": result.bias_prior_sigma,
    }
    (output / "observation_bias_profile.json").write_text(
        json.dumps(payload, indent=2), encoding="utf-8"
    )


def write_grouped_observation_bias_report(
    result: GroupedObservationBiasReport, output_dir: str | Path
) -> None:
    """Write a grouped multi-key observation-bias report as JSON."""
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    payload = {
        "metadata_key": result.metadata_key,
        "observation_keys": list(result.observation_keys),
        "group_biases": result.group_biases,
        "bias_prior_mean": result.bias_prior_mean,
        "bias_prior_sigma": result.bias_prior_sigma,
        "bias_bounds": list(result.bias_bounds),
        "objective": result.objective,
        "bias_prior_penalty": result.bias_prior_penalty,
    }
    (output / "grouped_observation_bias_report.json").write_text(
        json.dumps(payload, indent=2), encoding="utf-8"
    )


def read_grouped_observation_bias_report(
    input_dir: str | Path,
) -> GroupedObservationBiasReport:
    """Read and validate a grouped multi-key observation-bias report archive."""
    payload = json.loads(
        (Path(input_dir) / "grouped_observation_bias_report.json").read_text(
            encoding="utf-8"
        )
    )
    metadata_key = payload.get("metadata_key")
    observation_keys = payload.get("observation_keys")
    group_biases = payload.get("group_biases")
    if not isinstance(metadata_key, str) or not metadata_key:
        raise ValueError("group bias report metadata_key must be a non-empty string")
    if not isinstance(observation_keys, list) or not observation_keys:
        raise ValueError("group bias report observation_keys must be a non-empty list")
    if not isinstance(group_biases, dict) or not group_biases:
        raise ValueError("group bias report group_biases must be a non-empty dict")
    for key in observation_keys:
        if not isinstance(key, str) or not key:
            raise ValueError(
                "group bias report observation keys must be non-empty strings"
            )
    for group, values in group_biases.items():
        if not isinstance(group, str) or not group or not isinstance(values, dict):
            raise ValueError(
                "group bias report group entries must have string keys and dict values"
            )
        for key, value in values.items():
            if (
                not isinstance(key, str)
                or not key
                or not isinstance(value, (int, float))
                or not np.isfinite(value)
            ):
                raise ValueError(
                    "group bias report entries must contain finite numeric values"
                )
    bias_bounds = payload.get("bias_bounds", (-1.0, 1.0))
    if (
        not isinstance(bias_bounds, list)
        or len(bias_bounds) != 2
        or any(
            not isinstance(value, (int, float)) or not np.isfinite(value)
            for value in bias_bounds
        )
        or bias_bounds[0] >= bias_bounds[1]
    ):
        raise ValueError("group bias report bounds must be finite and increasing")
    prior_mean = payload.get("bias_prior_mean", 0.0)
    prior_sigma = payload.get("bias_prior_sigma", 0.1)
    if not isinstance(prior_mean, (int, float)) or not np.isfinite(prior_mean):
        raise ValueError("group bias report prior mean must be finite")
    if (
        not isinstance(prior_sigma, (int, float))
        or not np.isfinite(prior_sigma)
        or prior_sigma <= 0.0
    ):
        raise ValueError("group bias report prior sigma must be finite and positive")
    objective = payload.get("objective")
    prior_penalty = payload.get("bias_prior_penalty")
    if not isinstance(objective, (int, float)) or not np.isfinite(objective):
        raise ValueError("group bias report objective must be finite")
    if not isinstance(prior_penalty, (int, float)) or not np.isfinite(prior_penalty):
        raise ValueError("group bias report prior penalty must be finite")
    return GroupedObservationBiasReport(
        metadata_key,
        tuple(observation_keys),
        {
            group: {str(key): float(value) for key, value in values.items()}
            for group, values in group_biases.items()
        },
        float(prior_mean),
        float(prior_sigma),
        tuple(float(value) for value in bias_bounds),
        float(objective),
        float(prior_penalty),
    )


def write_grouped_observation_bias_summary(
    result: GroupedObservationBiasSummary, output_dir: str | Path
) -> None:
    """Write a grouped observation-bias summary as JSON."""
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    payload = {
        "metadata_key": result.metadata_key,
        "observation_keys": list(result.observation_keys),
        "group_biases": result.group_biases,
        "bias_prior_mean": result.bias_prior_mean,
        "bias_prior_sigma": result.bias_prior_sigma,
        "bias_bounds": list(result.bias_bounds),
        "objective": result.objective,
        "bias_prior_penalty": result.bias_prior_penalty,
    }
    (output / "grouped_observation_bias_summary.json").write_text(
        json.dumps(payload, indent=2), encoding="utf-8"
    )


def read_grouped_observation_bias_summary(
    input_dir: str | Path,
) -> GroupedObservationBiasSummary:
    """Read and validate a grouped observation-bias summary archive."""
    payload = json.loads(
        (Path(input_dir) / "grouped_observation_bias_summary.json").read_text(
            encoding="utf-8"
        )
    )
    metadata_key = payload.get("metadata_key")
    observation_keys = payload.get("observation_keys")
    group_biases = payload.get("group_biases")
    if not isinstance(metadata_key, str) or not metadata_key:
        raise ValueError("group bias summary metadata_key must be a non-empty string")
    if not isinstance(observation_keys, list) or not observation_keys:
        raise ValueError("group bias summary observation_keys must be a non-empty list")
    if not isinstance(group_biases, dict) or not group_biases:
        raise ValueError("group bias summary group_biases must be a non-empty dict")
    bias_bounds = payload.get("bias_bounds", (-1.0, 1.0))
    if (
        not isinstance(bias_bounds, list)
        or len(bias_bounds) != 2
        or any(
            not isinstance(value, (int, float)) or not np.isfinite(value)
            for value in bias_bounds
        )
        or bias_bounds[0] >= bias_bounds[1]
    ):
        raise ValueError("group bias summary bounds must be finite and increasing")
    prior_mean = payload.get("bias_prior_mean", 0.0)
    prior_sigma = payload.get("bias_prior_sigma", 0.1)
    if not isinstance(prior_mean, (int, float)) or not np.isfinite(prior_mean):
        raise ValueError("group bias summary prior mean must be finite")
    if (
        not isinstance(prior_sigma, (int, float))
        or not np.isfinite(prior_sigma)
        or prior_sigma <= 0.0
    ):
        raise ValueError("group bias summary prior sigma must be finite and positive")
    objective = payload.get("objective")
    prior_penalty = payload.get("bias_prior_penalty")
    if not isinstance(objective, (int, float)) or not np.isfinite(objective):
        raise ValueError("group bias summary objective must be finite")
    if not isinstance(prior_penalty, (int, float)) or not np.isfinite(prior_penalty):
        raise ValueError("group bias summary prior penalty must be finite")
    return GroupedObservationBiasSummary(
        metadata_key,
        tuple(observation_keys),
        {
            group: {str(key): float(value) for key, value in values.items()}
            for group, values in group_biases.items()
        },
        float(prior_mean),
        float(prior_sigma),
        tuple(float(value) for value in bias_bounds),
        float(objective),
        float(prior_penalty),
    )


def read_observation_bias_profile_result(
    input_dir: str | Path,
) -> ObservationBiasProfileResult:
    """Read and validate a grouped observation-bias profile archive."""
    payload = json.loads(
        (Path(input_dir) / "observation_bias_profile.json").read_text(encoding="utf-8")
    )
    biases = payload.get("biases")
    if (
        not isinstance(biases, dict)
        or not biases
        or any(
            not isinstance(key, str)
            or not key
            or not isinstance(value, (int, float))
            or not np.isfinite(value)
            for key, value in biases.items()
        )
    ):
        raise ValueError(
            "bias profile biases must be finite numeric values keyed by strings"
        )
    metadata_key = payload.get("metadata_key")
    observation_key = payload.get("observation_key")
    if not isinstance(metadata_key, str) or not metadata_key:
        raise ValueError("bias profile metadata_key must be a non-empty string")
    if not isinstance(observation_key, str) or not observation_key:
        raise ValueError("bias profile observation_key must be a non-empty string")
    for field in ("objective", "bias_prior_penalty"):
        value = payload.get(field)
        if not isinstance(value, (int, float)) or not np.isfinite(value):
            raise ValueError(f"bias profile {field} must be finite")
    corrections = payload.get("corrections_j_mol")
    if not isinstance(corrections, dict):
        raise ValueError("bias profile corrections must be a JSON object")
    records = payload.get("records", [])
    if not isinstance(records, list) or any(
        not isinstance(record, dict) for record in records
    ):
        raise ValueError("bias profile records must be an array of JSON objects")
    bounds = payload.get("bias_bounds", (-1.0, 1.0))
    if (
        not isinstance(bounds, list)
        or len(bounds) != 2
        or any(
            not isinstance(value, (int, float)) or not np.isfinite(value)
            for value in bounds
        )
        or bounds[0] >= bounds[1]
    ):
        raise ValueError("bias profile bounds must be finite and increasing")
    prior_mean = payload.get("bias_prior_mean", 0.0)
    prior_sigma = payload.get("bias_prior_sigma", 0.1)
    if not isinstance(prior_mean, (int, float)) or not np.isfinite(prior_mean):
        raise ValueError("bias profile prior mean must be finite")
    if (
        not isinstance(prior_sigma, (int, float))
        or not np.isfinite(prior_sigma)
        or prior_sigma <= 0.0
    ):
        raise ValueError("bias profile prior sigma must be finite and positive")
    return ObservationBiasProfileResult(
        {key: float(value) for key, value in biases.items()},
        float(payload["objective"]),
        {key: float(value) for key, value in corrections.items()},
        float(payload["bias_prior_penalty"]),
        metadata_key,
        observation_key,
        tuple(records),
        tuple(float(value) for value in bounds),
        float(prior_mean),
        float(prior_sigma),
    )


def write_loocv_result(result: LOOCVResult, output_dir: str | Path) -> None:
    """Write leave-one-out folds (corrections as embedded JSON) and mean misfit."""
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    with (output / "folds.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "held_out_experiment_id",
                "held_out_experiment_ids",
                "held_out_misfits",
                "held_out_log_predictive_densities",
                "misfit",
                "corrections_json",
            ]
        )
        for fold in result.folds:
            writer.writerow(
                [
                    fold.held_out_experiment_id,
                    json.dumps(fold.held_out_experiment_ids),
                    json.dumps(fold.held_out_misfits),
                    json.dumps(fold.held_out_log_predictive_densities),
                    fold.misfit,
                    json.dumps(fold.corrections_j_mol),
                ]
            )
    (output / "summary.json").write_text(
        json.dumps(
            {
                "mean_misfit": result.mean_misfit,
                "group_metadata_key": result.group_metadata_key,
                "mean_log_predictive_density": result.mean_log_predictive_density,
            },
            indent=2,
        ),
        encoding="utf-8",
    )


def read_loocv_result(input_dir: str | Path) -> LOOCVResult:
    """Reload an `LOOCVResult` previously written by `write_loocv_result`."""
    input_path = Path(input_dir)
    folds: list[LOOCVFold] = []
    with (input_path / "folds.csv").open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            held_out_ids = tuple(
                json.loads(row["held_out_experiment_ids"])
                if row.get("held_out_experiment_ids")
                else ()
            )
            held_out_misfits = tuple(
                float(value)
                for value in (
                    json.loads(row["held_out_misfits"])
                    if row.get("held_out_misfits")
                    else ()
                )
            )
            held_out_log_scores = tuple(
                float(value)
                for value in (
                    json.loads(row["held_out_log_predictive_densities"])
                    if row.get("held_out_log_predictive_densities")
                    else ()
                )
            )
            if held_out_ids or held_out_misfits:
                if len(held_out_ids) != len(held_out_misfits):
                    raise ValueError(
                        "LOOCV held-out IDs and misfits must have equal lengths"
                    )
                if any(
                    not isinstance(identifier, str) or not identifier
                    for identifier in held_out_ids
                ):
                    raise ValueError("LOOCV held-out IDs must be non-empty strings")
                if any(not np.isfinite(value) for value in held_out_misfits):
                    raise ValueError("LOOCV held-out misfits must be finite")
            if held_out_log_scores and len(held_out_ids) != len(held_out_log_scores):
                raise ValueError(
                    "LOOCV held-out IDs and log predictive scores must have equal lengths"
                )
            elif not row.get("held_out_experiment_ids") and not row.get(
                "held_out_misfits"
            ):
                held_out_ids = ()
                held_out_misfits = ()
            folds.append(
                LOOCVFold(
                    row["held_out_experiment_id"],
                    json.loads(row["corrections_json"]),
                    float(row["misfit"]),
                    held_out_ids,
                    held_out_misfits,
                    held_out_log_scores,
                )
            )
    summary = json.loads((input_path / "summary.json").read_text(encoding="utf-8"))
    return LOOCVResult(
        tuple(folds),
        summary["mean_misfit"],
        summary.get("group_metadata_key"),
        summary.get("mean_log_predictive_density"),
    )


def write_grid_posterior(result: GridPosteriorResult, output_dir: str | Path) -> None:
    """Write grid axes, flattened density/log-posterior arrays, and summary stats."""
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    for index, axis in enumerate(result.axes):
        np.savetxt(output / f"axis_{index}.csv", axis, delimiter=",")
    np.savetxt(
        output / "posterior_density_flat.csv",
        result.posterior_density.ravel(),
        delimiter=",",
    )
    np.savetxt(
        output / "log_posterior_flat.csv", result.log_posterior.ravel(), delimiter=","
    )
    summary = {
        "parameter_names": list(result.parameter_names),
        "shape": list(result.posterior_density.shape),
        "marginal_mean": result.marginal_mean,
        "marginal_std": result.marginal_std,
    }
    (output / "summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )


def read_grid_posterior(input_dir: str | Path) -> GridPosteriorResult:
    """Reload a `GridPosteriorResult` previously written by `write_grid_posterior`."""
    input_path = Path(input_dir)
    summary = json.loads((input_path / "summary.json").read_text(encoding="utf-8"))
    shape = tuple(summary["shape"])
    axes = tuple(
        np.atleast_1d(np.loadtxt(input_path / f"axis_{index}.csv", delimiter=","))
        for index in range(len(summary["parameter_names"]))
    )
    posterior_density = np.loadtxt(
        input_path / "posterior_density_flat.csv", delimiter=","
    ).reshape(shape)
    log_posterior = np.loadtxt(
        input_path / "log_posterior_flat.csv", delimiter=","
    ).reshape(shape)
    return GridPosteriorResult(
        tuple(summary["parameter_names"]),
        axes,
        log_posterior,
        posterior_density,
        summary["marginal_mean"],
        summary["marginal_std"],
    )


def write_posterior_resolution_result(
    result: PosteriorResolutionResult, output_dir: str | Path
) -> None:
    """Persist posterior-averaged local-resolution matrices as JSON."""
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    payload = {
        "parameter_names": list(result.parameter_names),
        "sample_count": result.sample_count,
        "mean_matrix": np.asarray(result.mean_matrix, dtype=float).tolist(),
        "lower_matrix": np.asarray(result.lower_matrix, dtype=float).tolist(),
        "upper_matrix": np.asarray(result.upper_matrix, dtype=float).tolist(),
        "mean_diagonal": np.asarray(result.mean_diagonal, dtype=float).tolist(),
        "lower_diagonal": np.asarray(result.lower_diagonal, dtype=float).tolist(),
        "upper_diagonal": np.asarray(result.upper_diagonal, dtype=float).tolist(),
        "selected_chain_indices": list(result.selected_chain_indices),
        "selected_draw_indices": list(result.selected_draw_indices),
        "selection_strategy": result.selection_strategy,
        "observable_names": list(result.observable_names),
        "credible_interval": list(result.credible_interval),
        "step_fraction": result.step_fraction,
        "observation_covariance_source": result.observation_covariance_source,
        "theoretical_covariance_source": result.theoretical_covariance_source,
        "observation_covariance": (
            np.asarray(result.observation_covariance, dtype=float).tolist()
            if result.observation_covariance is not None
            else None
        ),
        "theoretical_covariance": (
            np.asarray(result.theoretical_covariance, dtype=float).tolist()
            if result.theoretical_covariance is not None
            else None
        ),
        "prior_covariance": (
            np.asarray(result.prior_covariance, dtype=float).tolist()
            if result.prior_covariance is not None
            else None
        ),
    }
    (output / "posterior_resolution.json").write_text(
        json.dumps(payload, indent=2, allow_nan=False), encoding="utf-8"
    )


def read_posterior_resolution_result(
    input_dir: str | Path,
) -> PosteriorResolutionResult:
    """Load and validate posterior-averaged local-resolution summaries."""
    path = Path(input_dir) / "posterior_resolution.json"
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise TypeError("posterior resolution archive must be a JSON object")
    names = payload.get("parameter_names")
    if (
        not isinstance(names, list)
        or not names
        or any(not isinstance(name, str) or not name for name in names)
        or len(set(names)) != len(names)
    ):
        raise ValueError("parameter_names must be unique non-empty strings")
    sample_count = payload.get("sample_count")
    if (
        isinstance(sample_count, bool)
        or not isinstance(sample_count, int)
        or sample_count <= 0
    ):
        raise ValueError("sample_count must be a positive integer")
    count = len(names)
    matrix_names = ("mean_matrix", "lower_matrix", "upper_matrix")
    diagonal_names = ("mean_diagonal", "lower_diagonal", "upper_diagonal")
    matrices: dict[str, np.ndarray] = {}
    diagonals: dict[str, np.ndarray] = {}
    for name in matrix_names:
        try:
            value = np.asarray(payload[name], dtype=float)
        except (KeyError, TypeError, ValueError) as error:
            raise ValueError(f"{name} must be a finite square matrix") from error
        if value.shape != (count, count) or not np.all(np.isfinite(value)):
            raise ValueError(f"{name} must be a finite square matrix")
        matrices[name] = value
    for name in diagonal_names:
        try:
            value = np.asarray(payload[name], dtype=float)
        except (KeyError, TypeError, ValueError) as error:
            raise ValueError(f"{name} must be a finite parameter vector") from error
        if value.shape != (count,) or not np.all(np.isfinite(value)):
            raise ValueError(f"{name} must be a finite parameter vector")
        diagonals[name] = value
    chain_indices = payload.get("selected_chain_indices", [])
    draw_indices = payload.get("selected_draw_indices", [])
    strategy = payload.get("selection_strategy", "unspecified")
    if not isinstance(chain_indices, list) or not isinstance(draw_indices, list):
        raise TypeError("selected sample indices must be arrays")
    if bool(chain_indices) != bool(draw_indices):
        raise ValueError("chain and draw index provenance must be supplied together")
    if chain_indices and (
        len(chain_indices) != sample_count
        or len(draw_indices) != sample_count
        or any(
            isinstance(index, bool) or not isinstance(index, int) or index < 0
            for index in (*chain_indices, *draw_indices)
        )
    ):
        raise ValueError(
            "selected sample indices must be non-negative integers matching sample_count"
        )
    if not isinstance(strategy, str) or not strategy:
        raise ValueError("selection_strategy must be a non-empty string")
    observable_names = payload.get("observable_names", [])
    if (
        not isinstance(observable_names, list)
        or any(not isinstance(name, str) or not name for name in observable_names)
        or len(set(observable_names)) != len(observable_names)
    ):
        raise ValueError("observable_names must be unique non-empty strings")
    credible_interval = payload.get("credible_interval", [0.05, 0.95])
    if (
        not isinstance(credible_interval, list)
        or len(credible_interval) != 2
        or any(
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not np.isfinite(value)
            for value in credible_interval
        )
        or not 0.0 <= credible_interval[0] < credible_interval[1] <= 1.0
    ):
        raise ValueError("credible_interval must contain ordered probabilities")
    step_fraction = payload.get("step_fraction")
    if step_fraction is not None and (
        isinstance(step_fraction, bool)
        or not isinstance(step_fraction, (int, float))
        or not np.isfinite(step_fraction)
        or step_fraction <= 0.0
    ):
        raise ValueError("step_fraction must be finite and positive or null")
    observation_source = payload.get("observation_covariance_source", "unspecified")
    theoretical_source = payload.get("theoretical_covariance_source", "unspecified")
    valid_sources = {"unspecified", "explicit", "problem_assembled"}
    if (
        observation_source not in valid_sources
        or theoretical_source not in valid_sources
    ):
        raise ValueError("covariance provenance labels are invalid")
    covariance_matrices: dict[str, np.ndarray | None] = {}
    for name in (
        "observation_covariance",
        "theoretical_covariance",
        "prior_covariance",
    ):
        raw_covariance = payload.get(name)
        if raw_covariance is None:
            covariance_matrices[name] = None
            continue
        try:
            covariance = np.asarray(raw_covariance, dtype=float)
        except (TypeError, ValueError) as error:
            raise ValueError(
                f"{name} must be a finite square matrix or null"
            ) from error
        if (
            covariance.ndim != 2
            or covariance.shape[0] != covariance.shape[1]
            or not np.all(np.isfinite(covariance))
        ):
            raise ValueError(f"{name} must be a finite square matrix or null")
        if name == "prior_covariance":
            if covariance.shape != (len(names), len(names)):
                raise ValueError(f"{name} must match parameter_names length")
        elif observable_names and covariance.shape != (
            len(observable_names),
            len(observable_names),
        ):
            raise ValueError(f"{name} must match observable_names length")
        covariance_matrices[name] = covariance
    return PosteriorResolutionResult(
        tuple(names),
        sample_count,
        matrices["mean_matrix"],
        matrices["lower_matrix"],
        matrices["upper_matrix"],
        diagonals["mean_diagonal"],
        diagonals["lower_diagonal"],
        diagonals["upper_diagonal"],
        tuple(chain_indices),
        tuple(draw_indices),
        strategy,
        tuple(observable_names),
        (float(credible_interval[0]), float(credible_interval[1])),
        float(step_fraction) if step_fraction is not None else None,
        observation_source,
        theoretical_source,
        covariance_matrices["observation_covariance"],
        covariance_matrices["theoretical_covariance"],
        covariance_matrices["prior_covariance"],
    )
