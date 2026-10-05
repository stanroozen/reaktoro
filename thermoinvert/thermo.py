"""Safe thermodynamic database overlays for Reaktoro."""

import json
import os
import tempfile
from pathlib import Path
from typing import Mapping


class DatabaseOverlay:
    """Create a temporary JSON database with selected dHf corrections."""

    def __init__(self, source_path: str | os.PathLike[str]):
        self.source_path = Path(source_path)
        if not self.source_path.is_file():
            raise FileNotFoundError(self.source_path)

    def create(
        self,
        corrections_j_mol: Mapping[str, float],
    ) -> Path:
        """Create an overlay applying each correction exactly once.

        Keys are either a bare species name (perturbs ``H0``, the default
        field consumed by Reaktoro's equilibrium solver) or a compound
        ``"species:field"`` key targeting one of ``H0``, ``G0``, ``V0``
        (all read from ``StandardThermoModel.Constant``), ``Hf``
        (``ThermoReference``, informational only), or ``GH``
        (``Metadata.PerpleX_Params``, informational only).
        """
        with self.source_path.open("r", encoding="utf-8") as handle:
            database = json.load(handle)
        species = database.get("Species", {})
        grouped: dict[str, dict[str, float]] = {}
        for key, correction in corrections_j_mol.items():
            name, _, field_name = key.partition(":")
            field_name = field_name or "H0"
            fields = grouped.setdefault(name, {})
            fields[field_name] = fields.get(field_name, 0.0) + float(correction)
        for name, field_deltas in grouped.items():
            if name not in species:
                raise KeyError(
                    f"Species {name!r} is absent from {self.source_path.name}"
                )
            entry = species[name]
            thermo_reference = entry.get("ThermoReference")
            constant = entry.get("StandardThermoModel", {}).get("Constant")
            if not isinstance(thermo_reference, dict) or not isinstance(constant, dict):
                raise ValueError(
                    f"Species {name!r} has no editable constant thermodynamic record"
                )
            for field_name, delta in field_deltas.items():
                if field_name == "Hf":
                    if "Hf" not in thermo_reference:
                        raise ValueError(f"Species {name!r} has no ThermoReference Hf")
                    thermo_reference["Hf"] += delta
                elif field_name in {"H0", "G0", "V0"}:
                    if field_name not in constant:
                        raise ValueError(
                            f"Species {name!r} has no constant {field_name}"
                        )
                    constant[field_name] += delta
                elif field_name == "GH":
                    metadata_params = entry.get("Metadata", {}).get(
                        "PerpleX_Params", {}
                    )
                    if "GH" not in metadata_params:
                        raise ValueError(f"Species {name!r} has no Perple_X GH")
                    metadata_params["GH"] += delta
                else:
                    raise ValueError(
                        f"Unsupported thermodynamic overlay field: {field_name}"
                    )
        file_descriptor, output_name = tempfile.mkstemp(
            prefix="thermoinvert_", suffix=".json"
        )
        os.close(file_descriptor)
        output = Path(output_name)
        with output.open("w", encoding="utf-8") as handle:
            json.dump(database, handle, indent=2)
        return output
