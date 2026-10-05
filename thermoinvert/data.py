"""Data structures for thermodynamic parameter inversion."""

from dataclasses import dataclass, field
from math import isfinite
from typing import Any


@dataclass(frozen=True)
class Experiment:
    """One phase-equilibrium observation at reported P-T conditions."""

    id: str
    pressure_bar: float
    temperature_c: float
    bulk_composition: dict[str, float]
    observed_phases: tuple[str, ...]
    pressure_sigma_bar: float | None = None
    temperature_sigma_c: float | None = None
    absent_phases: tuple[str, ...] = ()
    phase_compositions: dict[str, dict[str, float]] | None = None
    phase_composition_sigma: dict[str, dict[str, float]] | None = None
    phase_fractions: dict[str, float] | None = None
    phase_fraction_sigma: dict[str, float] | None = None
    observation_covariance: list[list[float]] | None = None
    observation_covariance_keys: tuple[str, ...] = ()
    theoretical_covariance: list[list[float]] | None = None
    theoretical_covariance_keys: tuple[str, ...] = ()
    observation_bias: dict[str, float] | None = None
    metadata: dict[str, Any] = field(default_factory=dict)
    phase_composition_concentration: dict[str, float] | None = None
    phase_composition_detection_limit: dict[str, float] | None = None

    def __post_init__(self) -> None:
        concentrations = self.phase_composition_concentration
        if concentrations is not None and not isinstance(concentrations, dict):
            raise ValueError("phase_composition_concentration must be a mapping")
        compositions = self.phase_compositions or {}
        for phase, concentration in (concentrations or {}).items():
            if phase not in compositions:
                raise ValueError(
                    f"Dirichlet concentration for {phase!r} requires a phase composition"
                )
            if (
                not isinstance(concentration, (int, float))
                or not isfinite(concentration)
                or concentration <= 0.0
            ):
                raise ValueError(
                    "phase composition concentrations must be finite and positive"
                )
        detection_limits = self.phase_composition_detection_limit
        if detection_limits is not None and not isinstance(detection_limits, dict):
            raise ValueError("phase_composition_detection_limit must be a mapping")
        for phase, limit in (detection_limits or {}).items():
            observed = compositions.get(phase)
            if not isinstance(observed, dict) or not any(
                value == 0.0 for value in observed.values()
            ):
                raise ValueError(
                    f"A detection limit for {phase!r} requires a zero-marked component"
                )
            if (
                not isinstance(limit, (int, float))
                or not isfinite(limit)
                or not 0.0 < limit < 1.0
            ):
                raise ValueError(
                    "phase composition detection limits must lie strictly between zero and one"
                )


@dataclass(frozen=True)
class EquilibriumBracket:
    """Two endpoint observations that bound a phase transition."""

    id: str
    lower: Experiment
    upper: Experiment
    lower_required_phases: tuple[str, ...] = ()
    upper_required_phases: tuple[str, ...] = ()
    lower_forbidden_phases: tuple[str, ...] = ()
    upper_forbidden_phases: tuple[str, ...] = ()


@dataclass(frozen=True)
class Parameter:
    """A bounded scalar inversion parameter with explicit physical units."""

    name: str
    prior_mean: float = 0.0
    prior_sigma: float = 1.0
    lower_bound: float | None = None
    upper_bound: float | None = None
    unit: str = "1"

    def __post_init__(self) -> None:
        if not self.name:
            raise ValueError("Parameter name must be non-empty")
        if not self.unit:
            raise ValueError("Parameter unit must be non-empty")
        if (
            not isinstance(self.prior_mean, (int, float))
            or not isinstance(self.prior_sigma, (int, float))
            or not self._finite(self.prior_mean)
            or not self._finite(self.prior_sigma)
        ):
            raise ValueError("Parameter prior values must be finite numbers")
        if self.prior_sigma <= 0.0:
            raise ValueError("Parameter prior_sigma must be positive")
        if self.lower_bound is not None and not self._finite(self.lower_bound):
            raise ValueError("Parameter lower_bound must be finite or null")
        if self.upper_bound is not None and not self._finite(self.upper_bound):
            raise ValueError("Parameter upper_bound must be finite or null")
        if (
            self.lower_bound is not None
            and self.upper_bound is not None
            and self.lower_bound >= self.upper_bound
        ):
            raise ValueError("Parameter bounds must increase")

    @staticmethod
    def _finite(value: float) -> bool:
        return bool(isfinite(value))

    @property
    def prior_mean_j_mol(self) -> float:
        return self.prior_mean

    @property
    def prior_sigma_j_mol(self) -> float:
        return self.prior_sigma

    def bounds(self) -> tuple[float, float]:
        lower = self.lower_bound
        upper = self.upper_bound
        if lower is None:
            lower = self.prior_mean - 5.0 * self.prior_sigma
        if upper is None:
            upper = self.prior_mean + 5.0 * self.prior_sigma
        if lower >= upper:
            raise ValueError(f"Invalid bounds for {self.name}: {lower} >= {upper}")
        return lower, upper

    @property
    def key(self) -> str:
        """Correction-map key for generic scalar inversion parameters."""
        return self.name


@dataclass(frozen=True)
class ScalarParameter(Parameter):
    """Backward-compatible name for the generic native-unit parameter."""


@dataclass(frozen=True, init=False)
class DirichletConcentrationParameter(Parameter):
    """Bounded positive concentration parameter for one phase composition."""

    phase: str

    def __init__(
        self,
        phase: str,
        prior_mean: float = 20.0,
        prior_sigma: float = 10.0,
        lower_bound: float = 1.0e-3,
        upper_bound: float = 1.0e4,
    ) -> None:
        if not phase:
            raise ValueError("Dirichlet concentration phase must be non-empty")
        object.__setattr__(self, "name", f"dirichlet_kappa[{phase}]")
        object.__setattr__(self, "prior_mean", prior_mean)
        object.__setattr__(self, "prior_sigma", prior_sigma)
        object.__setattr__(self, "lower_bound", lower_bound)
        object.__setattr__(self, "upper_bound", upper_bound)
        object.__setattr__(self, "unit", "1")
        object.__setattr__(self, "phase", phase)
        Parameter.__post_init__(self)
        if lower_bound <= 0.0:
            raise ValueError("Dirichlet concentration lower_bound must be positive")
        if prior_mean <= 0.0 or not lower_bound <= prior_mean <= upper_bound:
            raise ValueError(
                "Dirichlet concentration prior_mean must be positive and within bounds"
            )

    @property
    def key(self) -> str:
        return f"dirichlet_concentration:{self.phase}"


@dataclass(frozen=True, init=False)
class ThermodynamicParameter(Parameter):
    """A bounded correction to one database field of one species.

    ``field`` selects which additive database record is perturbed. Only
    ``StandardThermoModel.Constant`` entries (``H0``, ``G0``, ``V0``) are
    actually consumed by Reaktoro's equilibrium solver for this database
    format; ``Hf`` and ``GH`` are informational/reference copies retained
    for schema completeness but have no effect on computed equilibria.
    """

    name: str
    phase: str
    field: str = "H0"
    parameter_type: str = "dHf"
    prior_mean_j_mol: float = 0.0
    prior_sigma_j_mol: float = 1000.0
    lower_bound_j_mol: float | None = None
    upper_bound_j_mol: float | None = None

    def __init__(
        self,
        name: str,
        phase: str,
        field: str = "H0",
        parameter_type: str = "dHf",
        prior_mean_j_mol: float = 0.0,
        prior_sigma_j_mol: float = 1000.0,
        lower_bound_j_mol: float | None = None,
        upper_bound_j_mol: float | None = None,
    ) -> None:
        if not phase or not field or not parameter_type:
            raise ValueError(
                "ThermodynamicParameter phase, field, and parameter_type must be non-empty"
            )
        object.__setattr__(self, "name", name)
        object.__setattr__(self, "prior_mean", prior_mean_j_mol)
        object.__setattr__(self, "prior_sigma", prior_sigma_j_mol)
        object.__setattr__(self, "lower_bound", lower_bound_j_mol)
        object.__setattr__(self, "upper_bound", upper_bound_j_mol)
        object.__setattr__(self, "unit", "J/mol")
        object.__setattr__(self, "phase", phase)
        object.__setattr__(self, "field", field)
        object.__setattr__(self, "parameter_type", parameter_type)
        object.__setattr__(self, "prior_mean_j_mol", prior_mean_j_mol)
        object.__setattr__(self, "prior_sigma_j_mol", prior_sigma_j_mol)
        object.__setattr__(self, "lower_bound_j_mol", lower_bound_j_mol)
        object.__setattr__(self, "upper_bound_j_mol", upper_bound_j_mol)
        Parameter.__post_init__(self)

    @property
    def key(self) -> str:
        """Correction-dict key: bare phase name for the default H0 field,
        otherwise ``phase:field`` so multiple fields on one species coexist."""
        return self.phase if self.field == "H0" else f"{self.phase}:{self.field}"

    def bounds(self) -> tuple[float, float]:
        lower = self.lower_bound_j_mol
        upper = self.upper_bound_j_mol
        if lower is None:
            lower = self.prior_mean_j_mol - 5.0 * self.prior_sigma_j_mol
        if upper is None:
            upper = self.prior_mean_j_mol + 5.0 * self.prior_sigma_j_mol
        if lower >= upper:
            raise ValueError(f"Invalid bounds for {self.name}: {lower} >= {upper}")
        return lower, upper
