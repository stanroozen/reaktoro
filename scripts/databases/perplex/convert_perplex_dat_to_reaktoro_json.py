#!/usr/bin/env python3
"""
Convert Perple_X .dat databases in embedded/databases/perplex to Reaktoro-style JSON.

Scope and behavior:
- Parses endmember entries only (name/EoS/composition/properties ... end).
- Skips make definitions in header (begin_makes ... end_makes).
- Skips pseudo entries with null composition.
- Handles apparent Gibbs convention by reading header tags:
  - If HSC_conversion is enabled and entry uses GH, convert to SUPCRT:
      Gf = GH + Tr * S_elements
    where S_elements is built from component elemental entropy values in header.
  - G0 is taken as SUPCRT convention directly.

Notes:
- This script keeps all parsed Perple_X parameters in Metadata.
- For EoS 8/9 entries with the required keys, a HollandPowell model block is emitted.
"""

from __future__ import annotations

import argparse
import json
import re
from dataclasses import dataclass, field
from pathlib import Path
import sys
from typing import Dict, List, Optional, Tuple

TR_DEFAULT = 298.15

PAIR_RE = re.compile(r"([A-Za-z][A-Za-z0-9_\-]*)\s*=\s*([^\s|]+)")
COMP_RE = re.compile(r"([A-Za-z0-9_\-]+)\(([^\)]+)\)")
ELEM_RE = re.compile(r"([A-Z][a-z]?)(\d*(?:\.\d+)?)")

UPPERCASE_DIGRAPH_MAP = {
    "NA": "Na",
    "MG": "Mg",
    "AL": "Al",
    "SI": "Si",
    "CA": "Ca",
    "TI": "Ti",
    "MN": "Mn",
    "FE": "Fe",
    "NI": "Ni",
    "ZR": "Zr",
    "CL": "Cl",
    "CU": "Cu",
    "CR": "Cr",
}

# Fallback component entropies (J/mol/K) for HSC->SUPCRT conversion when
# header component entropies are absent (e.g., hp11ver.dat).
COMPONENT_ENTROPY_FALLBACK = {
    "NA2O": 205.175,
    "MGO": 135.255,
    "AL2O3": 364.425,
    "SIO2": 223.96,
    "K2O": 231.935,
    "CAO": 144.205,
    "TIO2": 235.87,
    "MNO": 134.795,
    "FEO": 129.855,
    "NIO": 132.375,
    "ZRO2": 244.33,
    "CL2": 223.08,
    "O2": 205.15,
    "H2O": 233.255,
    "CO2": 210.89,
    "CUO": 135.725,
    "CR2O3": 358.811,
    "S2": 64.1,
    "F2": 202.79,
    "N2": 191.610,
}


@dataclass
class HeaderInfo:
    title: str = ""
    tr: float = TR_DEFAULT
    hsc_conversion: bool = False
    components_weight: Dict[str, float] = field(default_factory=dict)
    components_entropy: Dict[str, float] = field(default_factory=dict)


@dataclass
class Entry:
    name: str
    eos: int
    composition: Dict[str, float]
    params: Dict[str, float]
    source_lines: List[str]


def strip_comment(line: str) -> str:
    return line.split("|", 1)[0].strip()


def to_float(token: str) -> float:
    return float(token.replace("D", "E").replace("d", "e"))


def parse_formula_to_elements(formula: str) -> Dict[str, float]:
    # Perple_X component names may be all-uppercase (e.g., AL2O3, SIO2, NA2O),
    # so a standard CamelCase element parser is insufficient.
    elems: Dict[str, float] = {}

    if any(ch.islower() for ch in formula):
        for sym, num in ELEM_RE.findall(formula):
            coeff = float(num) if num else 1.0
            elems[sym] = elems.get(sym, 0.0) + coeff
        return elems

    i = 0
    text = formula.upper()
    while i < len(text):
        ch = text[i]
        if not ch.isalpha():
            i += 1
            continue

        sym = None
        if i + 1 < len(text):
            digraph = text[i : i + 2]
            sym = UPPERCASE_DIGRAPH_MAP.get(digraph)
            if sym is not None:
                i += 2
        if sym is None:
            sym = ch
            i += 1

        j = i
        while j < len(text) and (text[j].isdigit() or text[j] == "."):
            j += 1
        coeff = float(text[i:j]) if j > i else 1.0
        i = j

        elems[sym] = elems.get(sym, 0.0) + coeff

    return elems


def parse_composition(comp_line: str) -> Dict[str, float]:
    s = strip_comment(comp_line)
    if s.lower() == "null":
        return {}

    comp: Dict[str, float] = {}
    for cname, coeff in COMP_RE.findall(s):
        comp[cname] = comp.get(cname, 0.0) + to_float(coeff)
    return comp


def parse_header(lines: List[str]) -> Tuple[HeaderInfo, int]:
    info = HeaderInfo()
    idx = 0

    # Title = first non-empty, non-comment line
    while idx < len(lines):
        s = strip_comment(lines[idx])
        if s:
            info.title = s
            idx += 1
            break
        idx += 1

    in_components = False
    in_makes = False

    while idx < len(lines):
        raw = lines[idx]
        s = strip_comment(raw)
        low = s.lower()

        if not s:
            idx += 1
            continue

        if low.startswith("t("):
            # fallback for malformed files; handled below if standard var parsing exists
            pass

        if low == "hsc_conversion":
            info.hsc_conversion = True
            idx += 1
            continue

        if low == "begin_components":
            in_components = True
            idx += 1
            continue

        if low == "end_components":
            in_components = False
            idx += 1
            continue

        if low == "begin_makes":
            in_makes = True
            idx += 1
            continue

        if low == "end_makes":
            in_makes = False
            idx += 1
            continue

        if in_makes:
            idx += 1
            continue

        if in_components:
            toks = s.split()
            if len(toks) >= 2:
                name = toks[0]
                try:
                    info.components_weight[name] = to_float(toks[1])
                except Exception:
                    pass
                if len(toks) >= 3:
                    try:
                        info.components_entropy[name] = to_float(toks[2])
                    except Exception:
                        pass
            idx += 1
            continue

        # standard variables: capture T reference
        if low.startswith("t("):
            toks = s.split()
            if len(toks) >= 2:
                try:
                    info.tr = to_float(toks[1])
                except Exception:
                    pass
            idx += 1
            continue

        if low == "end":
            return info, idx + 1

        idx += 1

    return info, idx


def parse_entry(lines: List[str], start: int) -> Tuple[Optional[Entry], int]:
    i = start

    while i < len(lines):
        s = strip_comment(lines[i])
        if s:
            break
        i += 1

    if i >= len(lines):
        return None, i

    first_raw = lines[i]
    first = strip_comment(first_raw)
    m = re.match(r"^([A-Za-z0-9_\-\(\)\.]+)\s+EoS\s*=\s*([0-9]+)", first)
    if not m:
        return None, i + 1

    name = m.group(1)
    eos = int(m.group(2))
    source_lines = [lines[i]]
    i += 1

    # composition line: first non-empty/comment line after entry header
    while i < len(lines) and not strip_comment(lines[i]):
        source_lines.append(lines[i])
        i += 1
    if i >= len(lines):
        return None, i

    comp_line = lines[i]
    source_lines.append(comp_line)
    composition = parse_composition(comp_line)
    i += 1

    params: Dict[str, float] = {}

    # H values are often provided in the comment of the EoS line: "| H= ..."
    h_match = re.search(r"\bH\s*=\s*([-+]?\d+(?:\.\d+)?(?:[eEdD][-+]?\d+)?)", first_raw)
    if h_match:
        try:
            params["H"] = to_float(h_match.group(1))
        except Exception:
            pass

    while i < len(lines):
        raw = lines[i]
        s = strip_comment(raw)
        source_lines.append(raw)
        if not s:
            i += 1
            continue
        if s.lower() == "end":
            i += 1
            break

        for k, v in PAIR_RE.findall(s):
            try:
                params[k] = to_float(v)
            except Exception:
                continue

        i += 1

    return Entry(
        name=name,
        eos=eos,
        composition=composition,
        params=params,
        source_lines=source_lines,
    ), i


def composition_to_elements(composition: Dict[str, float]) -> Dict[str, float]:
    elements: Dict[str, float] = {}
    for comp_name, coeff in composition.items():
        part = parse_formula_to_elements(comp_name)
        for sym, n in part.items():
            elements[sym] = elements.get(sym, 0.0) + coeff * n
    return elements


def fmt_coeff(x: float) -> str:
    if abs(x - round(x)) < 1e-12:
        return str(int(round(x)))
    return f"{x:.8f}".rstrip("0").rstrip(".")


def elements_to_formula(elements: Dict[str, float]) -> str:
    order = [
        "H",
        "C",
        "N",
        "Na",
        "Mg",
        "Al",
        "Si",
        "K",
        "Ca",
        "Ti",
        "Mn",
        "Fe",
        "Ni",
        "Zr",
        "Cu",
        "Cr",
        "S",
        "Cl",
        "F",
        "O",
    ]
    parts: List[str] = []
    used = set()
    for sym in order:
        v = elements.get(sym, 0.0)
        if abs(v) > 1e-12:
            c = fmt_coeff(v)
            parts.append(sym if c == "1" else f"{sym}{c}")
            used.add(sym)
    for sym in sorted(elements.keys()):
        if sym in used:
            continue
        v = elements[sym]
        if abs(v) > 1e-12:
            c = fmt_coeff(v)
            parts.append(sym if c == "1" else f"{sym}{c}")
    return "".join(parts)


def elements_to_string(elements: Dict[str, float]) -> str:
    return " ".join(
        f"{fmt_coeff(v)}:{k}" for k, v in elements.items() if abs(v) > 1e-12
    )


def _component_entropy(comp: str, header: HeaderInfo) -> Optional[float]:
    if comp in header.components_entropy:
        return header.components_entropy[comp]
    return COMPONENT_ENTROPY_FALLBACK.get(comp.upper())


def _selements_for_entry(entry: Entry, header: HeaderInfo) -> float:
    s_elems = 0.0
    for comp, coeff in entry.composition.items():
        s_comp = _component_entropy(comp, header)
        if s_comp is not None:
            s_elems += coeff * s_comp
    return s_elems


def compute_gf_supcrt(entry: Entry, header: HeaderInfo) -> Tuple[Optional[float], str]:
    params = entry.params
    if "GH" in params:
        gh = params["GH"]
        if header.hsc_conversion:
            s_elems = _selements_for_entry(entry, header)
            return gh + header.tr * s_elems, "GH_HSC_to_SUPCRT"
        return gh, "GH_passthrough"

    if "G0" in params:
        g0 = params["G0"]
        # hp11-style cards can store HSC values in G0 without HSC_conversion tag.
        # Detect via the defining HSC relation: G0 ~= H - Tr*S0.
        if "H" in params and "S0" in params:
            hsc_like = abs((params["H"] - header.tr * params["S0"]) - g0) < 500.0
            if hsc_like:
                s_elems = _selements_for_entry(entry, header)
                return g0 + header.tr * s_elems, "G0_detected_HSC_to_SUPCRT"
        return g0, "G0_passthrough"

    return None, "missing_G"


def try_build_holland_powell(entry: Entry, gf: Optional[float]) -> Optional[Dict]:
    p = entry.params
    # Works for H&P style EoS 8/9 cards in perplex files.
    if entry.eos not in (8, 9):
        return None
    required = ["S0", "V0", "c1", "c2", "c3", "c5", "b1", "b6", "b8"]
    if not all(k in p for k in required):
        return None

    # Perplex stores K0 as bar-scale values in many files; tc-ds parsers use kbar.
    k0_raw = p["b6"]
    kappa0 = k0_raw * 1.0e8 if abs(k0_raw) < 1.0e5 else k0_raw * 1.0e5

    k0pp_raw = p.get("b7", 0.0)
    kappa0pp = k0pp_raw / 1.0e5 if abs(k0pp_raw) < 1.0e-4 else k0pp_raw / 1.0e8

    model = {
        "Gf": gf if gf is not None else p.get("G0", 0.0),
        "Hf": p.get("H", p.get("Hf", 0.0)),
        "Sr": p["S0"],
        "Vr": p["V0"] / 1.0e5,
        "a": p["c1"],
        "b": p["c2"],
        "c": p["c3"],
        "d": p["c5"],
        "alpha0": p["b1"],
        "kappa0": kappa0,
        "kappa0p": p["b8"],
        "kappa0pp": kappa0pp,
        "numatoms": 0.0,
        "Tmax": 9999.0,
    }
    return model


def build_constant_fallback(
    entry: Entry, gf: Optional[float], header: HeaderInfo
) -> Optional[Dict[str, Dict]]:
    if gf is None:
        return None

    p = entry.params
    s0 = p.get("S0")
    h0 = p.get("H", p.get("Hf"))
    if h0 is None and s0 is not None:
        h0 = gf + header.tr * s0

    v0_j_per_bar = p.get("V0")
    v0 = v0_j_per_bar / 1.0e5 if v0_j_per_bar is not None else None

    thermo_reference: Dict[str, object] = {
        "Gf": gf,
        "S0": s0,
        "V0": v0_j_per_bar,
    }
    if h0 is not None:
        thermo_reference["Hf"] = h0

    constant_model: Dict[str, object] = {
        "G0": gf,
    }
    if h0 is not None:
        constant_model["H0"] = h0
    if v0 is not None:
        constant_model["V0"] = v0

    return {
        "ThermoReference": thermo_reference,
        "StandardThermoModel": {"Constant": constant_model},
    }


def convert_file(dat_path: Path, out_path: Optional[Path] = None) -> Path:
    lines = dat_path.read_text(encoding="utf-8", errors="replace").splitlines()
    header, idx = parse_header(lines)

    species: Dict[str, Dict] = {}

    i = idx
    while i < len(lines):
        entry, j = parse_entry(lines, i)
        if entry is None:
            i = j
            continue
        i = j

        # Endmembers only: skip null composition pseudo entries.
        if not entry.composition:
            continue

        elems = composition_to_elements(entry.composition)
        formula = elements_to_formula(elems)
        elements_str = elements_to_string(elems)
        gf, gf_mode = compute_gf_supcrt(entry, header)

        item: Dict[str, object] = {
            "Name": entry.name,
            "Formula": formula,
            "Elements": elements_str,
            "AggregateState": "Solid",
            "Source": f"Perple_X {dat_path.name}",
            "Metadata": {
                "PerpleX_EoS": entry.eos,
                "PerpleX_Composition": entry.composition,
                "PerpleX_Params": entry.params,
                "HSC_conversion_enabled": header.hsc_conversion,
                "GfConversionMode": gf_mode,
                "ReferenceTemperature_K": header.tr,
            },
        }

        hp_model = try_build_holland_powell(entry, gf)
        if hp_model is not None:
            hp_model["numatoms"] = sum(elems.values())
            item["StandardThermoModel"] = {"HollandPowell": hp_model}
        else:
            fallback = build_constant_fallback(entry, gf, header)
            if fallback is not None:
                if gf_mode == "G0_passthrough" and entry.params.get("S0") is None:
                    print(
                        f"WARNING: {entry.name}: G0_passthrough mode with no S0 — "
                        f"H0 cannot be derived (H0 = G0 + Tr*S0 requires S0). "
                        f"StandardThermoModel.Constant will contain G0 only; "
                        f"H/S/V queries off-reference will return degenerate values. "
                        f"Typical cause: old-style PerpleX database (e.g. b89/b92/ba96) "
                        f"that stores SUPCRT G0 without elemental entropies or per-species S0.",
                        file=sys.stderr,
                    )
                item.update(fallback)

        species[entry.name] = item

    db = {
        "Database": dat_path.name,
        "Title": header.title,
        "Species": species,
    }

    if out_path is None:
        out_path = dat_path.with_name(f"{dat_path.stem}-reaktoro.json")

    out_path.write_text(json.dumps(db, indent=2), encoding="utf-8")
    return out_path


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Convert Perple_X .dat endmembers to Reaktoro-style JSON"
    )
    parser.add_argument("-i", "--input", help="Input .dat file path")
    parser.add_argument(
        "--all", action="store_true", help="Convert all *.dat files in this folder"
    )
    parser.add_argument("-o", "--output", help="Output JSON path (single-file mode)")
    args = parser.parse_args()

    here = Path(__file__).resolve().parent
    db_dir = (here / "../../../embedded/databases/perplex").resolve()

    if not args.input and not args.all:
        raise SystemExit("Use --input <file> or --all")

    if args.all:
        paths = sorted(p for p in db_dir.glob("*.dat") if p.is_file())
        for p in paths:
            out = convert_file(p)
            print(f"Converted {p.name} -> {out.name}")
        return 0

    inpath = Path(args.input)
    if not inpath.is_absolute():
        inpath = (Path.cwd() / inpath).resolve()

    outpath = Path(args.output).resolve() if args.output else None
    out = convert_file(inpath, outpath)
    print(f"Converted {inpath.name} -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
