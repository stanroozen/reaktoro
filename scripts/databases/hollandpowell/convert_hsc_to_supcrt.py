#!/usr/bin/env python3
"""
Convert Holland-Powell species Gibbs energies from HSC to SUPCRT convention
in a Reaktoro-style JSON database.

Per Perple_X documentation (perplex.ethz.ch/perplex_thermodynamic_data_file.html):
    G_SUPCRT = G_HSC + Tr * S_elements
where S_elements is the stoichiometric sum of standard elemental entropies at Tr, Pr.

HSC convention:   G0_HSC    = Hf - Tr * Sr
SUPCRT convention: G0_SUPCRT = Hf - Tr * (Sr - S_elements) = G0_HSC + Tr * S_elements

IMPORTANT — DO NOT DOUBLE-CONVERT
-----------------------------------
The companion script convert_tc_ds_to_reaktoro_json.py reads raw tc-ds*.txt files,
which store Hf and Sr. It already computes Gf in SUPCRT convention:
    Gf = Hf - Tr * (Sr - S_elements)
The output JSON files (e.g. tc-ds62-reaktoro.json) are therefore already in
SUPCRT convention. Applying THIS script on those files produces incorrect,
doubly-shifted Gf values and must be avoided.

Use this script ONLY when the input JSON contains Gf values that were stored in
the raw HSC convention: G0 = Hf - Tr * Sr  (no elemental entropy correction).

This script assumes species are in the form:
{
  "Species": {
    "name": {
      "Elements": "2:Mg 1:Si 4:O",
      "StandardThermoModel": {
        "HollandPowell": {
          "Gf": ...,  # must be raw HSC-convention input (not already SUPCRT)
          ...
        }
      }
    }
  }
}
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, Tuple

TR_DEFAULT = 298.15

# J/mol/K values consistent with Perple_X HSC_conversion usage and the
# tc-ds converter constants currently used in this repository.
ELEMENT_ENTROPY_J_PER_MOL_K = {
    "Na": 51.30,
    "Mg": 32.68,
    "Al": 28.35,
    "Si": 18.81,
    "K": 64.68,
    "Ca": 41.63,
    "Ti": 30.63,
    "Mn": 32.01,
    "Fe": 27.28,
    "Ni": 4.184 * 7.14,
    "Zr": 4.184 * 9.32,
    "Cl": 4.184 * 53.288 / 2.0,
    "O": 205.20 / 2.0,
    "H": 130.70 / 2.0,
    "C": 5.74,
    # Not present in the legacy constant set above. Kept at 0.0 to avoid
    # introducing unverified values silently; warnings are emitted when used.
    "Cu": 0.0,
    "Cr": 0.0,
    "S": 0.0,
    "N": 0.0,
}


def parse_elements_string(elements: str) -> Dict[str, float]:
    """
    Parse strings like "2:Mg 4:O 1:Si" into {"Mg": 2.0, "O": 4.0, "Si": 1.0}.
    """
    result: Dict[str, float] = {}
    if not elements.strip():
        return result

    for token in elements.split():
        if ":" not in token:
            raise ValueError(f"Invalid Elements token '{token}' (missing ':')")
        coeff_str, symbol = token.split(":", 1)
        if not symbol:
            raise ValueError(f"Invalid Elements token '{token}' (missing symbol)")

        coeff = float(coeff_str)
        result[symbol] = result.get(symbol, 0.0) + coeff

    return result


def compute_selements(
    elements_map: Dict[str, float],
) -> Tuple[float, Dict[str, float], Dict[str, float]]:
    """
    Return (Selements, known_contrib, unknown_contrib).
    """
    known_contrib: Dict[str, float] = {}
    unknown_contrib: Dict[str, float] = {}
    total = 0.0

    for symbol, coeff in elements_map.items():
        if symbol in ELEMENT_ENTROPY_J_PER_MOL_K:
            contrib = coeff * ELEMENT_ENTROPY_J_PER_MOL_K[symbol]
            known_contrib[symbol] = contrib
            total += contrib
        else:
            unknown_contrib[symbol] = coeff

    return total, known_contrib, unknown_contrib


def convert_database(data: Dict, tr: float) -> Tuple[int, int]:
    species = data.get("Species")
    if not isinstance(species, dict):
        raise ValueError("Input JSON must contain a 'Species' object.")

    converted = 0
    warned = 0

    for name, spec in species.items():
        if not isinstance(spec, dict):
            continue

        stm = spec.get("StandardThermoModel", {})
        hp = stm.get("HollandPowell", {}) if isinstance(stm, dict) else {}

        if not isinstance(hp, dict) or "Gf" not in hp:
            continue

        elements_str = spec.get("Elements", "")
        if not isinstance(elements_str, str) or not elements_str.strip():
            continue

        elems = parse_elements_string(elements_str)
        selements, _, unknown = compute_selements(elems)

        gf_hsc = float(hp["Gf"])
        shift = tr * selements
        gf_supcrt = gf_hsc + shift

        hp["Gf"] = gf_supcrt

        metadata = spec.setdefault("Metadata", {})
        if isinstance(metadata, dict):
            metadata["GfConventionInput"] = "HSC"
            metadata["GfConventionOutput"] = "SUPCRT"
            metadata["GfHSCOriginal"] = gf_hsc
            metadata["HSCtoSUPCRTShift_J_per_mol"] = shift
            metadata["ReferenceTemperature_K"] = tr

            if unknown:
                warned += 1
                metadata["HSCtoSUPCRTUnknownElements"] = unknown

        converted += 1

    return converted, warned


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Convert Reaktoro HP Gibbs energies from HSC to SUPCRT convention"
    )
    parser.add_argument("-i", "--input", required=True, help="Input JSON path")
    parser.add_argument("-o", "--output", default=None, help="Output JSON path")
    parser.add_argument(
        "--tr", type=float, default=TR_DEFAULT, help="Reference temperature in K"
    )
    parser.add_argument("--indent", type=int, default=2, help="JSON indentation")
    args = parser.parse_args()

    inpath = Path(args.input)
    if not inpath.exists():
        raise FileNotFoundError(f"Input file not found: {inpath}")

    if args.output:
        outpath = Path(args.output)
    else:
        outpath = inpath.with_name(f"{inpath.stem}-supcrt{inpath.suffix}")

    data = json.loads(inpath.read_text(encoding="utf-8"))
    converted, warned = convert_database(data, tr=args.tr)

    outpath.write_text(json.dumps(data, indent=args.indent), encoding="utf-8")

    print(f"Converted species: {converted}")
    print(f"Output file: {outpath}")
    if warned:
        print(f"Warnings: {warned} species had unknown elemental entropy symbols.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
