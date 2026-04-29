#!/usr/bin/env python3
"""
Convert THERMOCALC Holland & Powell DS6 text databases (e.g., tc-ds62.txt)
into Reaktoro JSON database format, and export covariance data to JSON.

This script mirrors the key parsing/conversion logic used in Perple_X:
- Perple_X/src/hp2010tover.f
- Perple_X/src/read_hp_covariance.f

Outputs:
1) Reaktoro species database JSON with HollandPowell standard thermo parameters.
2) Covariance JSON with packed upper-triangle covariance data and metadata.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple


# Mapping used by Perple_X read_hp_covariance.f / hp2010tover.f
# hp-index (1-based) -> ver-index (1-based)
HP_TO_VER = [4, 7, 3, 9, 2, 8, 6, 1, 5, 13, 14, 15, 12, 19, 10, 11, 18, 16, 17]

# 1-based component atom counts used in Perple_X sanity check
COMP_ATOMS = [
    3.0,
    2.0,
    5.0,
    3.0,
    3.0,
    2.0,
    3.0,
    2.0,
    2.0,
    2.0,
    3.0,
    2.0,
    2.0,
    3.0,
    3.0,
    2.0,
    5.0,
    2.0,
    0.0,
]

# Element entropy constants used in Perple_X conversion
SNA = 51.30
SMG = 32.68
SAL = 28.35
SSI = 18.81
SK = 64.68
SCA = 41.63
STI = 30.63
SMN = 32.01
SFE = 27.28
SNI = 4.184 * 7.14
SZR = 4.184 * 9.32
SCL = 4.184 * 53.288 / 2.0
SO = 205.20 / 2.0
SH = 130.70 / 2.0
SC = 5.74

TR = 298.15

# After Perple_X "convert to oxide stoichiometry", components are interpreted as:
# 1 NaO0.5, 2 MgO, 3 AlO1.5, 4 SiO2, 5 KO0.5, 6 CaO, 7 TiO2, 8 MnO, 9 FeO,
# 10 NiO, 11 ZrO2, 12 Cl, 13 O, 14 H, 15 CO2, 16 CuO, 17 CrO1.5, 18 S, 19 e-
VER_COMPONENT_ELEMENTAL = {
    1: {"Na": 1.0, "O": 0.5},
    2: {"Mg": 1.0, "O": 1.0},
    3: {"Al": 1.0, "O": 1.5},
    4: {"Si": 1.0, "O": 2.0},
    5: {"K": 1.0, "O": 0.5},
    6: {"Ca": 1.0, "O": 1.0},
    7: {"Ti": 1.0, "O": 2.0},
    8: {"Mn": 1.0, "O": 1.0},
    9: {"Fe": 1.0, "O": 1.0},
    10: {"Ni": 1.0, "O": 1.0},
    11: {"Zr": 1.0, "O": 2.0},
    12: {"Cl": 1.0},
    13: {"O": 1.0},
    14: {"H": 1.0},
    15: {"C": 1.0, "O": 2.0},
    16: {"Cu": 1.0, "O": 1.0},
    17: {"Cr": 1.0, "O": 1.5},
    18: {"S": 1.0},
}

# Formula ordering tuned for geochemical readability
FORMULA_ORDER = [
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
    "O",
]


@dataclass
class Entry:
    code: str
    comp_ver: List[float]  # length 19, 1-based semantics in comments
    atoms_raw: float
    hf_j: float
    sr_j_per_mol_k: float
    vr_j_per_bar: float
    a_j_per_mol_k: float
    b_j_per_mol_k2: float
    c_jk_per_mol: float
    d_j_per_mol_khalf: float
    alpha0_per_k: float
    kappa0_pa: float
    kappa0p: float
    kappa0pp_per_pa: float
    lam: float
    extras: List[float]
    aggregate_state: str
    charge: float

    @property
    def numatoms(self) -> float:
        return self.atoms_raw

    @property
    def gf_j(self) -> float:
        dsf = (
            self.sr_j_per_mol_k
            - self.comp_ver[0] * SNA
            - self.comp_ver[1] * SMG
            - self.comp_ver[2] * SAL
            - self.comp_ver[3] * SSI
            - self.comp_ver[4] * SK
            - self.comp_ver[5] * SCA
            - self.comp_ver[6] * STI
            - self.comp_ver[7] * SMN
            - self.comp_ver[8] * SFE
            - self.comp_ver[9] * SNI
            - self.comp_ver[10] * SZR
            - self.comp_ver[11] * SCL
            - self.comp_ver[12] * SO
            - self.comp_ver[13] * SH
            - self.comp_ver[14] * SC
        )
        return self.hf_j - TR * dsf


def parse_floats(line: str) -> List[float]:
    return [float(tok) for tok in line.split()]


def parse_component_card(line: str) -> Tuple[str, List[float], float]:
    # First try token-based parsing (works for DS636 cards that include a long
    # species name between code and composition pairs).
    toks = line.split()
    if toks and toks[0] and toks[0][0].isalpha():
        code = toks[0]
        comp = [0.0] * 19
        atoms_raw = 0.0
        started = False
        j = 1
        while j + 1 < len(toks):
            try:
                idx = int(toks[j])
            except ValueError:
                j += 1
                continue

            if idx == 0:
                if started:
                    break
                j += 1
                continue

            try:
                val = float(toks[j + 1])
            except ValueError:
                j += 1
                continue

            if not (1 <= idx <= 19):
                if started:
                    break
                j += 1
                continue

            started = True
            ver_idx = HP_TO_VER[idx - 1]
            comp[ver_idx - 1] = val
            atoms_raw += val
            j += 2

        if started:
            return code, comp, atoms_raw

    # Fallback for older DS6 fixed-column cards.
    text = line.rstrip("\n").ljust(132)
    code = text[:8].strip()
    if not code:
        raise ValueError("Empty phase name in component card")

    comp = [0.0] * 19
    atoms_raw = 0.0
    started = False

    for i in range(19):
        nst = 14 + i * 12  # Fortran nst=15+(i-1)*12 (1-based)
        if nst + 9 >= len(text):
            break

        if text[nst] == "0" and text[nst - 1] == " ":
            break

        idx_str = text[nst - 1 : nst + 1].strip()
        if not idx_str:
            continue

        idx = int(idx_str)
        if not (1 <= idx <= 19):
            raise ValueError(f"Unexpected HP component index {idx} in entry {code}")

        val_str = text[nst + 1 : nst + 10].strip()
        if not val_str:
            continue
        val = float(val_str)

        started = True
        ver_idx = HP_TO_VER[idx - 1]
        comp[ver_idx - 1] = val
        atoms_raw += val

    if not started:
        raise ValueError("No composition pairs found in component card")

    return code, comp, atoms_raw


def almost_int(x: float, tol: float = 1e-10) -> bool:
    return abs(x - round(x)) <= tol


def fmt_coeff(x: float) -> str:
    if almost_int(x):
        i = int(round(x))
        return "" if i == 1 else str(i)
    s = f"{x:.8f}".rstrip("0").rstrip(".")
    return s


def fmt_elements_coeff(x: float) -> str:
    if almost_int(x):
        return str(int(round(x)))
    return f"{x:.8f}".rstrip("0").rstrip(".")


def formula_from_elements(elements: Dict[str, float]) -> str:
    parts: List[str] = []
    used = set()
    for sym in FORMULA_ORDER:
        val = elements.get(sym, 0.0)
        if abs(val) > 1e-12:
            parts.append(f"{sym}{fmt_coeff(val)}")
            used.add(sym)
    for sym in sorted(elements.keys()):
        if sym in used:
            continue
        val = elements[sym]
        if abs(val) > 1e-12:
            parts.append(f"{sym}{fmt_coeff(val)}")
    return "".join(parts)


def convert_to_oxide_stoich(comp: List[float]) -> List[float]:
    # comp is 0-based list length 19; formulas below follow Fortran 1-based indexing
    c = comp[:]
    c[0] /= 2.0  # Na2O -> NaO0.5 basis
    c[2] /= 2.0  # Al2O3 -> AlO1.5 basis
    c[4] /= 2.0  # K2O -> KO0.5 basis
    c[13] /= 2.0  # H2O -> H basis for this representation

    c[12] = (
        c[12]
        - c[0]
        - c[1]
        - 3.0 * c[2]
        - 2.0 * c[3]
        - c[4]
        - c[5]
        - 2.0 * c[6]
        - c[7]
        - c[8]
        - c[9]
        - 2.0 * c[10]
        - c[13]
        - 2.0 * c[14]
        - c[15]
        - 1.5 * c[16]
    )

    c[12] /= 2.0  # O2 -> O basis
    c[11] /= 2.0  # Cl2 -> Cl basis
    c[16] /= 2.0  # Cr2O3 -> CrO1.5 basis
    c[17] /= 2.0  # S2 -> S basis
    return c


def elements_from_comp_oxide_basis(comp_oxide: List[float]) -> Dict[str, float]:
    elements: Dict[str, float] = {}
    for i in range(1, 19):
        coeff = comp_oxide[i - 1]
        if abs(coeff) < 1e-14:
            continue
        for sym, sto in VER_COMPONENT_ELEMENTAL.get(i, {}).items():
            elements[sym] = elements.get(sym, 0.0) + coeff * sto

    # Tiny numerical cleanup
    for key in list(elements.keys()):
        if abs(elements[key]) < 1e-12:
            elements.pop(key)
    return elements


def parse_tc_ds(filepath: Path) -> Tuple[List[Entry], List[float]]:
    lines = filepath.read_text(encoding="utf-8", errors="replace").splitlines()

    # Find first line that parses as a component card (supports DS62/DS636 headers).
    start = 0
    for i, line in enumerate(lines):
        try:
            parse_component_card(line)
            start = i
            break
        except Exception:
            continue

    entries: List[Entry] = []
    i = start
    cov_start = None

    while i < len(lines):
        line = lines[i]
        s = line.strip()

        # blank before covariance block
        if s == "":
            cov_start = i + 1
            break

        # Some DS6 files start covariance immediately with numeric rows and no
        # blank separator. Detect that transition explicitly.
        if s[0] in "+-.0123456789":
            cov_start = i
            break

        code, comp, atoms_raw = parse_component_card(line)

        # Next three lines are numeric property records
        if i + 3 >= len(lines):
            raise ValueError(f"Unexpected EOF while parsing entry {code}")

        l2 = parse_floats(lines[i + 1])
        l3 = parse_floats(lines[i + 2])
        l4 = parse_floats(lines[i + 3])

        if len(l2) < 3 or len(l3) < 4 or len(l4) < 3:
            raise ValueError(f"Malformed thermo records for entry {code}")

        # DS636 stores lam/transition code as first value of l2; DS62 stores
        # only (Hf, Sr, Vr) in l2 and puts lam on l4.
        if len(l2) >= 4 and abs(l2[0] - round(l2[0])) < 1.0e-9:
            lam = l2[0]
            hf_j = l2[1] * 1000.0
            sr_j_per_mol_k = l2[2] * 1000.0
            vr_j_per_bar = l2[3]
            lam_from_l2 = True
        else:
            hf_j = l2[0] * 1000.0
            sr_j_per_mol_k = l2[1] * 1000.0
            vr_j_per_bar = l2[2]
            lam_from_l2 = False

        a_j_per_mol_k = l3[0] * 1000.0
        b_j_per_mol_k2 = l3[1] * 1000.0
        c_jk_per_mol = l3[2] * 1000.0
        d_j_per_mol_khalf = l3[3] * 1000.0

        alpha0_per_k = l4[0]
        kappa0_pa = l4[1] * 1.0e8
        kappa0p = l4[2]

        # DS62 includes kappa0pp and lam on l4; DS636 frequently omits kappa0pp.
        if len(l4) >= 4 and abs(l4[3]) < 1.0:
            kappa0pp_per_pa = l4[3] / 1.0e8
            extras_offset = 4
        else:
            kappa0pp_per_pa = 0.0
            extras_offset = 3

        if lam_from_l2:
            extras = l4[extras_offset:]
        else:
            lam = l4[4] if len(l4) >= 5 else 0.0
            extras = l4[5:] if len(l4) >= 6 else []

        if lam == -1.0:
            aggregate_state = "Aqueous"
        elif abs(vr_j_per_bar) < 1e-14:
            aggregate_state = "Gas"
        else:
            aggregate_state = "Solid"

        # charge convention follows Perple_X handling of comp(19) electrons
        charge = -comp[18]

        entries.append(
            Entry(
                code=code,
                comp_ver=comp,
                atoms_raw=atoms_raw,
                hf_j=hf_j,
                sr_j_per_mol_k=sr_j_per_mol_k,
                vr_j_per_bar=vr_j_per_bar,
                a_j_per_mol_k=a_j_per_mol_k,
                b_j_per_mol_k2=b_j_per_mol_k2,
                c_jk_per_mol=c_jk_per_mol,
                d_j_per_mol_khalf=d_j_per_mol_khalf,
                alpha0_per_k=alpha0_per_k,
                kappa0_pa=kappa0_pa,
                kappa0p=kappa0p,
                kappa0pp_per_pa=kappa0pp_per_pa,
                lam=lam,
                extras=extras,
                aggregate_state=aggregate_state,
                charge=charge,
            )
        )

        i += 4
        # DS636+ style empty separator line
        if i < len(lines) and lines[i].strip() == "":
            i += 1

    cov_values: List[float] = []
    if cov_start is not None and cov_start < len(lines):
        tail_tokens: List[float] = []
        for line in lines[cov_start:]:
            for tok in line.split():
                try:
                    tail_tokens.append(float(tok))
                except ValueError:
                    pass

        n = len(entries)
        expected = n * (n + 1) // 2 + 1  # +1 because Perple_X skips first value
        if len(tail_tokens) >= expected:
            cov_values = tail_tokens[1:expected]
        elif len(tail_tokens) > 1:
            # Keep best-effort packed vector if truncated/unexpected
            cov_values = tail_tokens[1:]

    return entries, cov_values


def build_reaktoro_database(entries: List[Entry]) -> Dict:
    species: Dict[str, Dict] = {}

    for e in entries:
        comp_oxide = convert_to_oxide_stoich(e.comp_ver)
        elements = elements_from_comp_oxide_basis(comp_oxide)
        formula = formula_from_elements(elements)
        elements_str = " ".join(f"{fmt_elements_coeff(v)}:{k}" for k, v in elements.items())

        attrs: Dict[str, object] = {
            "Name": e.code,
            "Formula": formula,
            "Elements": elements_str,
            "AggregateState": e.aggregate_state,
            "StandardThermoModel": {
                "HollandPowell": {
                    "Gf": e.gf_j,
                    "Hf": e.hf_j,
                    "Sr": e.sr_j_per_mol_k,
                    "Vr": e.vr_j_per_bar / 1.0e5,  # J/bar -> J/Pa
                    "a": e.a_j_per_mol_k,
                    "b": e.b_j_per_mol_k2,
                    "c": e.c_jk_per_mol,
                    "d": e.d_j_per_mol_khalf,
                    "alpha0": e.alpha0_per_k,
                    "kappa0": e.kappa0_pa,
                    "kappa0p": e.kappa0p,
                    "kappa0pp": e.kappa0pp_per_pa,
                    "numatoms": e.numatoms,
                    "Tmax": 9999.0,
                }
            },
            "Source": "THERMOCALC DS6 (tc-ds)",
            "Metadata": {
                "OriginalCode": e.code,
                "TransitionCode_lam": e.lam,
                "TransitionExtras": e.extras,
            },
        }

        if abs(e.charge) > 1e-12:
            attrs["Charge"] = e.charge

        species[e.code] = attrs

    return {"Species": species}


def covariance_payload(
    entries: List[Entry], packed_upper: List[float], include_matrix: bool
) -> Dict:
    names = [e.code for e in entries]
    n = len(names)

    payload: Dict[str, object] = {
        "Schema": "Reaktoro.HollandPowell.Covariance.v1",
        "Entities": names,
        "NumEntities": n,
        "PackedUpperTriangle": packed_upper,
        "PackedLength": len(packed_upper),
    }

    diag = []
    idx = 0
    for i in range(n):
        if idx < len(packed_upper):
            diag.append(
                packed_upper[idx] * 1.0e3
            )  # matches read_hp_covariance .dia scaling
        idx += n - i
    payload["DiagonalTimes1e3"] = diag

    if include_matrix and len(packed_upper) >= n * (n + 1) // 2:
        mat = [[0.0 for _ in range(n)] for _ in range(n)]
        k = 0
        for i in range(n):
            for j in range(i, n):
                v = packed_upper[k]
                mat[i][j] = v
                mat[j][i] = v
                k += 1
        payload["MatrixSymmetric"] = mat

    return payload


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Convert tc-ds62 style HP database to Reaktoro JSON + covariance JSON"
    )
    parser.add_argument("-i", "--input", required=True, help="Path to tc-ds*.txt file")
    parser.add_argument(
        "--out-db", default=None, help="Output Reaktoro JSON database path"
    )
    parser.add_argument("--out-cov", default=None, help="Output covariance JSON path")
    parser.add_argument(
        "--include-cov-matrix",
        action="store_true",
        help="Include full symmetric covariance matrix",
    )
    parser.add_argument("--indent", type=int, default=2, help="JSON indentation")
    args = parser.parse_args()

    infile = Path(args.input)
    if not infile.exists():
        raise FileNotFoundError(f"Input file not found: {infile}")

    stem = infile.stem
    out_db = (
        Path(args.out_db) if args.out_db else infile.with_name(f"{stem}-reaktoro.json")
    )
    out_cov = (
        Path(args.out_cov)
        if args.out_cov
        else infile.with_name(f"{stem}-covariance.json")
    )

    entries, packed_cov = parse_tc_ds(infile)
    if not entries:
        raise RuntimeError("No entries parsed from input file.")

    db = build_reaktoro_database(entries)
    cov = covariance_payload(entries, packed_cov, args.include_cov_matrix)

    out_db.write_text(json.dumps(db, indent=args.indent), encoding="utf-8")
    out_cov.write_text(json.dumps(cov, indent=args.indent), encoding="utf-8")

    print(f"Parsed entries: {len(entries)}")
    print(f"Wrote Reaktoro JSON database: {out_db}")
    print(f"Wrote covariance JSON: {out_cov}")
    if packed_cov:
        print(f"Covariance packed values: {len(packed_cov)}")
    else:
        print("Covariance data not detected or could not be parsed.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
