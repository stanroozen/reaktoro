"""Reaktoro-backed equilibrium forward model."""

from dataclasses import dataclass
import hashlib
from pathlib import Path
import sys
from typing import Mapping

import numpy as np

from .data import Experiment
from .thermo import DatabaseOverlay

_repo_root = Path(__file__).resolve().parents[1]
_local_release = _repo_root / "build" / "Reaktoro" / "Release"
if _local_release.is_dir() and str(_local_release) not in sys.path:
    sys.path.insert(0, str(_local_release))

try:
    import reaktoro4py as rkt
except ImportError:  # pragma: no cover - depends on the active environment
    import reaktoro as rkt


@dataclass
class EquilibriumResult:
    valid: bool
    stable_phases: tuple[str, ...]
    phase_fractions: dict[str, float]
    phase_compositions: dict[str, dict[str, float]]
    diagnostics: dict[str, object]


def database_file_provenance(database_path: str | Path) -> dict[str, object]:
    """Return reproducible identity metadata for a thermodynamic database."""
    path = Path(database_path)
    if not path.is_file():
        raise FileNotFoundError(f"Database file does not exist: {path}")
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return {
        "database_path": str(path.resolve()),
        "database_sha256": digest.hexdigest(),
        "database_size_bytes": path.stat().st_size,
    }


class ReaktoroForwardModel:
    """Build an isolated Reaktoro system for each enthalpy parameter proposal."""

    def __init__(
        self,
        database_path: str | Path,
        aqueous_species: list[str],
        mineral_names: list[str],
        bulk_species: Mapping[str, str] | None = None,
        stability_threshold_mol: float = 1.0e-10,
    ):
        self.database_path = Path(database_path)
        self.aqueous_species = aqueous_species
        self.mineral_names = mineral_names
        self.bulk_species = dict(
            bulk_species
            or {
                "Al": "Al+3",
                "Ca": "Ca+2",
                "Fe": "Fe+2",
                "K": "K+",
                "Mg": "Mg+2",
                "Mn": "Mn+2",
                "Si": "SiO2,aq",
                "Zn": "Zn2+",
                "C": "HCO3-",
                "S": "HS-",
            }
        )
        self.stability_threshold_mol = stability_threshold_mol

    def provenance(self) -> dict[str, object]:
        """Describe all forward-model settings needed to identify its behavior."""
        version = getattr(rkt, "__version__", None)
        return {
            **database_file_provenance(self.database_path),
            "engine_module": rkt.__name__,
            "engine_version": str(version) if version is not None else None,
            "aqueous_species": list(self.aqueous_species),
            "mineral_names": list(self.mineral_names),
            "bulk_species": dict(self.bulk_species),
            "aqueous_activity_model": "ideal",
            "stability_threshold_mol": self.stability_threshold_mol,
        }

    def _build_system(self, database_path: Path):
        database = rkt.Database.fromFile(str(database_path))
        aqueous = rkt.AqueousPhase(" ".join(self.aqueous_species))
        aqueous.setActivityModel(rkt.ActivityModelIdealAqueous())
        minerals = rkt.MineralPhases(rkt.StringList(self.mineral_names))
        return rkt.ChemicalSystem(database, aqueous, minerals)

    def _initial_state(self, system, experiment: Experiment):
        state = rkt.ChemicalState(system)
        state.setSpeciesAmounts(np.full(system.species().size(), 1.0e-16))

        def set_if_present(name: str, amount: float) -> None:
            try:
                state.set(name, float(amount), "mol")
            except RuntimeError:
                pass

        for element, amount in experiment.bulk_composition.items():
            species_name = self.bulk_species.get(element)
            if species_name:
                set_if_present(species_name, amount)
        set_if_present("H2O", float(experiment.metadata.get("water_moles", 1.0)))
        set_if_present("H+", 1.0e-7)
        set_if_present("OH-", 1.0e-7)
        return state

    @staticmethod
    def _phase_composition(system, state, phase_name: str) -> dict[str, float]:
        """Return normalized species mole fractions for one Reaktoro phase."""
        phase_index = system.phases().index(phase_name)
        phase = system.phase(phase_index)
        amounts: dict[str, float] = {}
        for species in phase.species():
            species_name = species.name()
            try:
                amount = max(float(state.speciesAmount(species_name)), 0.0)
            except Exception:
                amount = 0.0
            if amount > 0.0:
                amounts[species_name] = amount
        total = sum(amounts.values())
        return (
            {name: amount / total for name, amount in amounts.items()}
            if total > 0.0
            else {}
        )

    def predict(
        self,
        experiment: Experiment,
        corrections_j_mol: Mapping[str, float] | None = None,
    ) -> EquilibriumResult:
        corrections = dict(corrections_j_mol or {})
        temporary_database: Path | None = None
        database_path = self.database_path
        try:
            if corrections:
                temporary_database = DatabaseOverlay(self.database_path).create(
                    corrections
                )
                database_path = temporary_database
            system = self._build_system(database_path)
            specs = rkt.EquilibriumSpecs(system)
            specs.temperature()
            specs.pressure()
            solver = rkt.EquilibriumSolver(specs)
            conditions = rkt.EquilibriumConditions(specs)
            conditions.temperature(experiment.temperature_c, "celsius")
            conditions.pressure(experiment.pressure_bar, "bar")
            state = self._initial_state(system, experiment)
            result = solver.solve(state, conditions)
            if not result.succeeded():
                return EquilibriumResult(
                    False, (), {}, {}, {"solver": str(result.optima)}
                )

            amounts = {}
            for phase in self.mineral_names:
                try:
                    amount = max(float(state.speciesAmount(phase)), 0.0)
                except Exception:
                    amount = 0.0
                if amount > self.stability_threshold_mol:
                    amounts[phase] = amount
            total = sum(amounts.values())
            fractions = (
                {phase: amount / total for phase, amount in amounts.items()}
                if total
                else {}
            )
            compositions = {
                phase: self._phase_composition(system, state, phase)
                for phase in amounts
            }
            return EquilibriumResult(
                True,
                tuple(sorted(amounts)),
                fractions,
                compositions,
                {
                    "temperature_c": experiment.temperature_c,
                    "pressure_bar": experiment.pressure_bar,
                },
            )
        except Exception as exc:
            return EquilibriumResult(False, (), {}, {}, {"exception": str(exc)})
        finally:
            if temporary_database is not None:
                temporary_database.unlink(missing_ok=True)
