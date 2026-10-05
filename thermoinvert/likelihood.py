"""Likelihood terms for sparse phase-equilibrium observations."""

from math import inf, lgamma

import numpy as np

from .data import EquilibriumBracket, Experiment
from .forward_model import EquilibriumResult

MISFIT_FORMS = (
    "weighted_lsq",
    "logit_lsq",
    "alr_lsq",
    "logistic_normal_nll",
    "dirichlet_nll",
    "chi_square",
    "raw_lsq",
)


def _residual_term(
    calculated: float, observed: float, sigma: float | None, misfit_form: str
) -> float:
    """One squared-residual contribution, per `misfit_form`:

    - ``weighted_lsq`` (default): ((calculated - observed) / sigma) ** 2 --
      mc_fit's LSQ/weighted-chi-square forms (identical formula in mc_fit's
      own documentation).
    - ``chi_square``: (calculated - observed) ** 2 / calculated -- mc_fit's
      classic chi-square statistic, documented there as "dubious" and rarely
      used; included for completeness, undefined contribution when
      ``calculated`` is zero (skipped in that case).
    - ``raw_lsq``: (calculated - observed) ** 2, unweighted by any sigma.
    """
    if misfit_form == "raw_lsq":
        return (calculated - observed) ** 2
    if misfit_form == "chi_square":
        return (calculated - observed) ** 2 / calculated if calculated != 0.0 else 0.0
    if not sigma or sigma <= 0.0:
        return 0.0
    return ((calculated - observed) / sigma) ** 2


def _logit_residual_term(
    calculated: float, observed: float, sigma: float | None
) -> float:
    """Squared residual in logit space with a delta-method sigma."""
    if not 0.0 < calculated < 1.0 or not 0.0 < observed < 1.0:
        raise ValueError("logit_lsq requires calculated and observed values in (0, 1)")
    if not sigma or sigma <= 0.0:
        return 0.0
    logit = lambda value: np.log(value / (1.0 - value))
    transformed_sigma = sigma / (observed * (1.0 - observed))
    return ((logit(calculated) - logit(observed)) / transformed_sigma) ** 2


def _continuous_residuals(
    experiment: Experiment, prediction: EquilibriumResult
) -> tuple[dict[str, float], dict[str, tuple[float, float, float | None]]]:
    residuals: dict[str, float] = {}
    observations: dict[str, tuple[float, float, float | None]] = {}
    if experiment.phase_fractions:
        for phase, observed in experiment.phase_fractions.items():
            key = f"fraction:{phase}"
            calculated = prediction.phase_fractions.get(phase, 0.0)
            residuals[key] = calculated - observed
            observations[key] = (
                calculated,
                observed,
                (
                    experiment.phase_fraction_sigma.get(phase)
                    if experiment.phase_fraction_sigma
                    else None
                ),
            )
    if experiment.phase_compositions:
        for phase, observed_composition in experiment.phase_compositions.items():
            calculated_composition = prediction.phase_compositions.get(phase, {})
            phase_sigmas = (
                experiment.phase_composition_sigma.get(phase, {})
                if experiment.phase_composition_sigma
                else {}
            )
            for component, observed in observed_composition.items():
                key = f"composition:{phase}:{component}"
                calculated = calculated_composition.get(component, 0.0)
                residuals[key] = calculated - observed
                observations[key] = (calculated, observed, phase_sigmas.get(component))
    biases = experiment.observation_bias or {}
    if any(
        not isinstance(key, str)
        or key not in residuals
        or not isinstance(value, (int, float))
        or not np.isfinite(value)
        for key, value in biases.items()
    ):
        raise ValueError("observation_bias must use continuous keys and finite values")
    for key, bias in biases.items():
        residuals[key] += float(bias)
        calculated, observed, sigma = observations[key]
        observations[key] = (calculated, observed - float(bias), sigma)
    return residuals, observations


def _validated_covariance(
    covariance_values: list[list[float]],
    keys: tuple[str, ...],
    residuals: dict[str, float],
    name: str,
    canonical_keys: tuple[str, ...],
) -> np.ndarray:
    if set(keys) != set(residuals) or len(keys) != len(residuals):
        raise ValueError(f"{name}_keys must identify every continuous observation")
    covariance = np.asarray(covariance_values, dtype=float)
    if covariance.shape != (len(keys), len(keys)):
        raise ValueError(f"{name} shape must match its keys")
    covariance = 0.5 * (covariance + covariance.T)
    if np.any(np.linalg.eigvalsh(covariance) <= 0.0):
        raise ValueError(f"{name} must be positive definite")
    indices = [keys.index(key) for key in canonical_keys]
    return covariance[np.ix_(indices, indices)]


def _correlated_residual_term(
    experiment: Experiment,
    residuals: dict[str, float],
    observations: dict[str, tuple[float, float, float | None]],
    theoretical_covariance_scale: float = 1.0,
) -> float:
    keys = (
        experiment.observation_covariance_keys or experiment.theoretical_covariance_keys
    )
    covariance = np.zeros((len(keys), len(keys)))
    if experiment.observation_covariance is not None:
        covariance += _validated_covariance(
            experiment.observation_covariance,
            experiment.observation_covariance_keys,
            residuals,
            "observation_covariance",
            keys,
        )
    elif experiment.theoretical_covariance is not None:
        sigmas = [observations[key][2] for key in keys]
        if any(sigma is None or sigma <= 0.0 for sigma in sigmas):
            raise ValueError(
                "independent analytical sigmas are required without observation_covariance"
            )
        validated_sigmas = [sigma for sigma in sigmas if sigma is not None]
        covariance += np.diag(np.square(validated_sigmas))
    if experiment.theoretical_covariance is not None:
        covariance += theoretical_covariance_scale**2 * _validated_covariance(
            experiment.theoretical_covariance,
            experiment.theoretical_covariance_keys,
            residuals,
            "theoretical_covariance",
            keys,
        )
    if np.any(np.linalg.eigvalsh(covariance) <= 0.0):
        raise ValueError(
            "combined analytical and theoretical covariance must be positive definite"
        )
    vector = np.asarray([residuals[key] for key in keys], dtype=float)
    return float(vector @ np.linalg.solve(covariance, vector))


def _alr_residual_term(
    experiment: Experiment,
    residuals: dict[str, float],
    observations: dict[str, tuple[float, float, float | None]],
    theoretical_covariance_scale: float,
    include_normalization: bool = False,
) -> float:
    """Return a covariance-aware additive-log-ratio composition misfit."""
    if experiment.phase_fractions:
        raise ValueError("alr_lsq does not support phase-fraction observations")
    composition_keys = tuple(key for key in residuals if key.startswith("composition:"))
    if not composition_keys or len(composition_keys) != len(residuals):
        raise ValueError("alr_lsq requires composition observations only")

    transformed_residuals: list[float] = []
    transform_rows: list[np.ndarray] = []
    raw_keys: list[str] = []
    compositions = experiment.phase_compositions or {}
    for phase in sorted(compositions):
        keys = tuple(
            f"composition:{phase}:{component}"
            for component in sorted(compositions[phase])
        )
        if len(keys) < 2:
            raise ValueError("alr_lsq requires at least two components per phase")
        reference = keys[-1]
        calculated = np.asarray([observations[key][0] for key in keys], dtype=float)
        observed = np.asarray([observations[key][1] for key in keys], dtype=float)
        if np.any(calculated <= 0.0) or np.any(observed <= 0.0):
            raise ValueError("alr_lsq requires strictly positive compositions")
        calculated = calculated / np.sum(calculated)
        observed = observed / np.sum(observed)
        for index, key in enumerate(keys[:-1]):
            transformed_residuals.append(
                float(
                    np.log(calculated[index] / calculated[-1])
                    - np.log(observed[index] / observed[-1])
                )
            )
            row = np.zeros(len(composition_keys), dtype=float)
            row[composition_keys.index(key)] = 1.0 / observed[index]
            row[composition_keys.index(reference)] = -1.0 / observed[-1]
            transform_rows.append(row)
        raw_keys.extend(keys)

    transform = np.asarray(transform_rows)
    covariance = np.zeros((len(composition_keys), len(composition_keys)))
    if experiment.observation_covariance is not None:
        covariance += _validated_covariance(
            experiment.observation_covariance,
            experiment.observation_covariance_keys,
            residuals,
            "observation_covariance",
            composition_keys,
        )
    else:
        sigmas = [observations[key][2] for key in composition_keys]
        if any(sigma is None or sigma <= 0.0 for sigma in sigmas):
            raise ValueError("alr_lsq requires positive composition sigmas")
        validated_sigmas = [sigma for sigma in sigmas if sigma is not None]
        covariance += np.diag(np.square(validated_sigmas))
    if experiment.theoretical_covariance is not None:
        covariance += theoretical_covariance_scale**2 * _validated_covariance(
            experiment.theoretical_covariance,
            experiment.theoretical_covariance_keys,
            residuals,
            "theoretical_covariance",
            composition_keys,
        )
    transformed_covariance = transform @ covariance @ transform.T
    transformed_covariance = 0.5 * (transformed_covariance + transformed_covariance.T)
    if np.any(np.linalg.eigvalsh(transformed_covariance) <= 0.0):
        raise ValueError("transformed composition covariance must be positive definite")
    vector = np.asarray(transformed_residuals, dtype=float)
    quadratic = float(vector @ np.linalg.solve(transformed_covariance, vector))
    if not include_normalization:
        return quadratic
    sign, logdet = np.linalg.slogdet(transformed_covariance)
    if sign <= 0.0:
        raise ValueError("transformed composition covariance must be positive definite")
    return 0.5 * (quadratic + logdet + len(vector) * np.log(2.0 * np.pi))


def _dirichlet_nll(
    experiment: Experiment,
    residuals: dict[str, float],
    observations: dict[str, tuple[float, float, float | None]],
    concentration_overrides: dict[str, float] | None = None,
) -> float:
    """Evaluate closed or explicitly censored phase compositions."""
    if experiment.phase_fractions:
        raise ValueError("dirichlet_nll does not support phase-fraction observations")
    compositions = experiment.phase_compositions or {}
    if not compositions or any(
        not isinstance(values, dict) or len(values) < 2
        for values in compositions.values()
    ):
        raise ValueError("dirichlet_nll requires multi-component compositions")
    if (
        experiment.observation_covariance is not None
        or experiment.theoretical_covariance is not None
    ):
        raise ValueError("dirichlet_nll does not support covariance matrices")

    total = 0.0
    for phase in sorted(compositions):
        keys = tuple(
            f"composition:{phase}:{component}"
            for component in sorted(compositions[phase])
        )
        calculated = np.asarray([observations[key][0] for key in keys], dtype=float)
        observed = np.asarray([observations[key][1] for key in keys], dtype=float)
        sigmas = [observations[key][2] for key in keys]
        if np.any(calculated <= 0.0) or np.any(observed < 0.0):
            raise ValueError(
                "Dirichlet predictions must be positive and observations non-negative"
            )
        calculated /= np.sum(calculated)
        censored = observed == 0.0
        detection_limit = (experiment.phase_composition_detection_limit or {}).get(
            phase
        )
        if np.any(censored) and detection_limit is None:
            raise ValueError(
                "zero Dirichlet components require an explicit phase detection limit"
            )
        if detection_limit is not None and not np.any(censored):
            raise ValueError("phase detection limit requires a zero-marked component")
        explicit_concentration = (concentration_overrides or {}).get(phase)
        if explicit_concentration is None:
            explicit_concentration = (
                experiment.phase_composition_concentration or {}
            ).get(phase)
        if explicit_concentration is not None:
            concentration = float(explicit_concentration)
        else:
            if any(sigma is None or sigma <= 0.0 for sigma in sigmas):
                raise ValueError(
                    "dirichlet_nll requires positive composition sigmas or an "
                    "explicit phase concentration"
                )
            validated_sigmas = [sigma for sigma in sigmas if sigma is not None]
            variance_ratios = np.asarray(
                [
                    calculated[index] * (1.0 - calculated[index]) / float(sigma) ** 2
                    for index, sigma in enumerate(validated_sigmas)
                ],
                dtype=float,
            )
            concentration = max(float(np.mean(variance_ratios) - 1.0), 1.0e-6)
        alpha = concentration * calculated
        if not np.any(censored):
            observed /= np.sum(observed)
            total += lgamma(concentration) - sum(
                lgamma(float(value)) for value in alpha
            )
            total -= float(np.sum((alpha - 1.0) * np.log(observed)))
            continue

        censored_alpha = float(np.sum(alpha[censored]))
        active_alpha = alpha[~censored]
        active_observed = observed[~censored]
        active_total = float(np.sum(active_observed))
        if active_total <= 0.0:
            raise ValueError("censored Dirichlet composition needs detected components")
        active_observed = active_observed / active_total
        active_concentration = float(np.sum(active_alpha))
        total += lgamma(active_concentration) - sum(
            lgamma(float(value)) for value in active_alpha
        )
        total -= float(np.sum((active_alpha - 1.0) * np.log(active_observed)))
        probability = _regularized_beta_cdf(
            float(detection_limit), censored_alpha, active_concentration
        )
        if probability <= 0.0:
            return inf
        total -= float(np.log(probability))
    return float(total)


def _beta_continued_fraction(a: float, b: float, x: float) -> float:
    """Evaluate the incomplete-beta continued fraction without SciPy."""
    maximum_iterations = 200
    epsilon = 3.0e-14
    minimum = 1.0e-300
    qab = a + b
    qap = a + 1.0
    qam = a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < minimum:
        d = minimum
    d = 1.0 / d
    fraction = d
    for iteration in range(1, maximum_iterations + 1):
        twice = 2 * iteration
        coefficient = iteration * (b - iteration) * x / ((qam + twice) * (a + twice))
        d = 1.0 + coefficient * d
        if abs(d) < minimum:
            d = minimum
        c = 1.0 + coefficient / c
        if abs(c) < minimum:
            c = minimum
        d = 1.0 / d
        fraction *= d * c
        coefficient = (
            -(a + iteration) * (qab + iteration) * x / ((a + twice) * (qap + twice))
        )
        d = 1.0 + coefficient * d
        if abs(d) < minimum:
            d = minimum
        c = 1.0 + coefficient / c
        if abs(c) < minimum:
            c = minimum
        d = 1.0 / d
        delta = d * c
        fraction *= delta
        if abs(delta - 1.0) < epsilon:
            return fraction
    raise ArithmeticError("incomplete beta continued fraction did not converge")


def _regularized_beta_cdf(x: float, a: float, b: float) -> float:
    """Regularized incomplete beta CDF for positive shape parameters."""
    if not 0.0 <= x <= 1.0 or a <= 0.0 or b <= 0.0:
        raise ValueError("beta CDF requires x in [0,1] and positive shapes")
    if x == 0.0:
        return 0.0
    if x == 1.0:
        return 1.0
    log_front = lgamma(a + b) - lgamma(a) - lgamma(b) + a * np.log(x) + b * np.log1p(-x)
    front = float(np.exp(log_front))
    if x < (a + 1.0) / (a + b + 2.0):
        return float(front * _beta_continued_fraction(a, b, x) / a)
    complement = front * _beta_continued_fraction(b, a, 1.0 - x) / b
    return float(np.clip(1.0 - complement, 0.0, 1.0))


def experiment_misfit(
    experiment: Experiment,
    prediction: EquilibriumResult,
    missing_phase_penalty: float = 100.0,
    extra_phase_penalty: float = 10.0,
    misfit_form: str = "weighted_lsq",
    theoretical_covariance_scale: float = 1.0,
    dirichlet_concentrations: dict[str, float] | None = None,
) -> float:
    """Return a transparent phase/fraction misfit for the Phase 1 model."""
    if misfit_form not in MISFIT_FORMS:
        raise ValueError(f"misfit_form must be one of {MISFIT_FORMS}")
    if (
        not np.isfinite(theoretical_covariance_scale)
        or theoretical_covariance_scale < 0.0
    ):
        raise ValueError("theoretical_covariance_scale must be finite and non-negative")
    if not prediction.valid:
        return inf
    observed = set(experiment.observed_phases)
    predicted = set(prediction.stable_phases)
    missing = observed - predicted
    extra = (predicted - observed) - set(experiment.absent_phases)
    score = missing_phase_penalty * len(missing)
    score += extra_phase_penalty * len(extra)
    score += missing_phase_penalty * len(
        predicted.intersection(experiment.absent_phases)
    )

    residuals, observations = _continuous_residuals(experiment, prediction)
    if misfit_form == "alr_lsq":
        score += _alr_residual_term(
            experiment, residuals, observations, theoretical_covariance_scale
        )
        return score
    if misfit_form == "logistic_normal_nll":
        score += _alr_residual_term(
            experiment,
            residuals,
            observations,
            theoretical_covariance_scale,
            include_normalization=True,
        )
        return score
    if misfit_form == "dirichlet_nll":
        if dirichlet_concentrations is not None and any(
            not isinstance(phase, str)
            or not phase
            or not isinstance(concentration, (int, float))
            or not np.isfinite(concentration)
            or concentration <= 0.0
            for phase, concentration in dirichlet_concentrations.items()
        ):
            raise ValueError("Dirichlet concentration overrides must be positive")
        score += _dirichlet_nll(
            experiment, residuals, observations, dirichlet_concentrations
        )
        return score
    if (
        experiment.observation_covariance is not None
        or experiment.theoretical_covariance is not None
    ):
        if misfit_form != "weighted_lsq":
            raise ValueError(
                "covariance likelihood requires misfit_form='weighted_lsq'"
            )
        score += _correlated_residual_term(
            experiment, residuals, observations, theoretical_covariance_scale
        )
    else:
        for calculated, observed_value, sigma in observations.values():
            if misfit_form == "logit_lsq":
                score += _logit_residual_term(calculated, observed_value, sigma)
            else:
                score += _residual_term(calculated, observed_value, sigma, misfit_form)
    return score


def experiment_log_likelihood(
    experiment: Experiment,
    prediction: EquilibriumResult,
    misfit_form: str = "weighted_lsq",
    theoretical_covariance_scale: float = 1.0,
    dirichlet_concentrations: dict[str, float] | None = None,
    phase_error_probability: float = 0.01,
) -> float:
    """Return a normalized observation log likelihood for posterior scoring.

    Supported continuous models are weighted Gaussian, logit-normal,
    covariance-aware ALR logistic-normal, and Dirichlet. Discrete observed and
    absent phases use independent Bernoulli error probability. Raw LSQ,
    chi-square, and unnormalized ALR objectives are intentionally rejected.
    """
    if not 0.0 < phase_error_probability < 0.5 or not np.isfinite(
        phase_error_probability
    ):
        raise ValueError("phase_error_probability must lie strictly between 0 and 0.5")
    if misfit_form not in {
        "weighted_lsq",
        "logit_lsq",
        "logistic_normal_nll",
        "dirichlet_nll",
    }:
        raise ValueError(
            "normalized predictive log likelihood is unavailable for this misfit_form"
        )
    if (
        not np.isfinite(theoretical_covariance_scale)
        or theoretical_covariance_scale < 0.0
    ):
        raise ValueError("theoretical_covariance_scale must be finite and non-negative")
    if not prediction.valid:
        return -inf

    log_likelihood = 0.0
    predicted_phases = set(prediction.stable_phases)
    for phase, present in (
        *((phase, True) for phase in experiment.observed_phases),
        *((phase, False) for phase in experiment.absent_phases),
    ):
        matches = (phase in predicted_phases) == present
        log_likelihood += np.log(
            1.0 - phase_error_probability if matches else phase_error_probability
        )

    residuals, observations = _continuous_residuals(experiment, prediction)
    if not observations:
        if not experiment.observed_phases and not experiment.absent_phases:
            raise ValueError("experiment contains no predictive observables")
        return float(log_likelihood)

    if misfit_form == "dirichlet_nll":
        return float(
            log_likelihood
            - _dirichlet_nll(
                experiment, residuals, observations, dirichlet_concentrations
            )
        )
    if misfit_form == "logistic_normal_nll":
        negative_log_density = _alr_residual_term(
            experiment,
            residuals,
            observations,
            theoretical_covariance_scale,
            include_normalization=True,
        )
        jacobian_log = 0.0
        compositions = experiment.phase_compositions or {}
        for phase in sorted(compositions):
            values = np.asarray(
                [
                    observations[f"composition:{phase}:{component}"][1]
                    for component in sorted(compositions[phase])
                ],
                dtype=float,
            )
            if np.any(values <= 0.0):
                raise ValueError(
                    "logistic_normal_nll requires positive observed compositions"
                )
            values /= np.sum(values)
            jacobian_log += float(np.sum(np.log(values)))
        return float(log_likelihood - negative_log_density - jacobian_log)

    if misfit_form == "logit_lsq":
        if (
            experiment.observation_covariance is not None
            or experiment.theoretical_covariance is not None
        ):
            raise ValueError(
                "logit-normal log scores do not support covariance matrices"
            )
        for calculated, observed, sigma in observations.values():
            if not 0.0 < calculated < 1.0 or not 0.0 < observed < 1.0:
                raise ValueError("logit_lsq log scores require values in (0, 1)")
            if sigma is None or sigma <= 0.0:
                raise ValueError("logit_lsq log scores require positive sigmas")
            transformed_sigma = sigma / (observed * (1.0 - observed))
            transformed_residual = np.log(calculated / (1.0 - calculated)) - np.log(
                observed / (1.0 - observed)
            )
            log_likelihood -= 0.5 * (transformed_residual / transformed_sigma) ** 2
            log_likelihood -= np.log(transformed_sigma)
            log_likelihood -= 0.5 * np.log(2.0 * np.pi)
            log_likelihood -= np.log(observed * (1.0 - observed))
        return float(log_likelihood)

    keys = tuple(
        experiment.observation_covariance_keys or experiment.theoretical_covariance_keys
    )
    if (
        experiment.observation_covariance is not None
        or experiment.theoretical_covariance is not None
    ):
        if not keys:
            raise ValueError(
                "covariance log scores require covariance observation keys"
            )
        covariance = np.zeros((len(keys), len(keys)), dtype=float)
        if experiment.observation_covariance is not None:
            covariance += _validated_covariance(
                experiment.observation_covariance,
                experiment.observation_covariance_keys,
                residuals,
                "observation_covariance",
                keys,
            )
        else:
            sigmas = [observations[key][2] for key in keys]
            if any(sigma is None or sigma <= 0.0 for sigma in sigmas):
                raise ValueError(
                    "positive analytical sigmas are required with theoretical covariance"
                )
            covariance += np.diag(
                [float(sigma) ** 2 for sigma in sigmas if sigma is not None]
            )
        if experiment.theoretical_covariance is not None:
            covariance += theoretical_covariance_scale**2 * _validated_covariance(
                experiment.theoretical_covariance,
                experiment.theoretical_covariance_keys,
                residuals,
                "theoretical_covariance",
                keys,
            )
    else:
        keys = tuple(observations)
        sigmas = [observations[key][2] for key in keys]
        if any(sigma is None or sigma <= 0.0 for sigma in sigmas):
            raise ValueError("weighted Gaussian log scores require positive sigmas")
        covariance = np.diag(
            [float(sigma) ** 2 for sigma in sigmas if sigma is not None]
        )
    sign, log_determinant = np.linalg.slogdet(covariance)
    if sign <= 0.0:
        raise ValueError("predictive observation covariance must be positive definite")
    vector = np.asarray([residuals[key] for key in keys], dtype=float)
    quadratic = float(vector @ np.linalg.solve(covariance, vector))
    dimension = len(vector)
    log_likelihood -= 0.5 * (
        quadratic + log_determinant + dimension * np.log(2.0 * np.pi)
    )
    return float(log_likelihood)


def bracket_misfit(
    bracket: EquilibriumBracket,
    lower_prediction: EquilibriumResult,
    upper_prediction: EquilibriumResult,
    phase_penalty: float = 100.0,
) -> float:
    """Score violations of phase-presence inequalities at bracket endpoints."""
    if not lower_prediction.valid or not upper_prediction.valid:
        return inf

    score = 0.0
    lower_phases = set(lower_prediction.stable_phases)
    upper_phases = set(upper_prediction.stable_phases)
    score += phase_penalty * len(set(bracket.lower_required_phases) - lower_phases)
    score += phase_penalty * len(
        lower_phases.intersection(bracket.lower_forbidden_phases)
    )
    score += phase_penalty * len(set(bracket.upper_required_phases) - upper_phases)
    score += phase_penalty * len(
        upper_phases.intersection(bracket.upper_forbidden_phases)
    )
    return score
