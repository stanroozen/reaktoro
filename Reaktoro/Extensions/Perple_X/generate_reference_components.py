#!/usr/bin/env python3
"""
Generate Perple_X-based reference CSVs for electrostatic and HKF components.

Inputs:
  test/epsh2o_matrix.csv
  test/gfunc_matrix.csv
  test/dh_matrix.csv
  test/born_matrix.csv
  test/hkf_matrix.csv

Outputs:
  test/epsh2o_reference.csv
  test/gfunc_reference.csv
  test/dh_reference.csv
  test/born_reference.csv
  test/hkf_reference.csv
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

ETA = 694656.968
THETA = 228.0
PSI = 2600.0
EPSILON0 = 78.47
NEUTRAL_RADIUS = 3.082
CDH = -42182668.74


def epsh2o(v_jbar: float, temperatureK: float) -> float:
    # Sverjenski 2014 / Fernandez et al. 1997 dielectric constant for pure water.
    # v_jbar: molar volume of water in J/bar (= cm³/mol / 10).
    # Fortran rlib.f: epsh2o = exp(...) * (0.1801526833D1 / v) ** (...)
    # where v is in J/bar. Equivalently: (18.01526833 / vcm3) when using vcm3 = v_jbar*10.
    # NOTE: the constant is 18.01526833 (Mw of water in g/mol), NOT 1.801526833.
    vcm3 = v_jbar * 10.0
    sqrtt = math.sqrt(temperatureK - 273.15) if temperatureK >= 273.15 else 0.0
    eps = math.exp(-8.016651e-5 * temperatureK + 4.769870482 - 0.06871618 * sqrtt)
    eps *= (18.01526833 / vcm3) ** (
        -1.576377e-3 * temperatureK + 1.185462878 + 0.06810288 * sqrtt
    )
    return eps


def gfunc(rho: float, pressureBar: float, temperatureK: float) -> float:
    if rho > 1.0:
        return 0.0
    g = (-6.557892e-6 * temperatureK + 9.3295764e-3) * temperatureK - 4.096745422
    g *= (1.0 - rho) ** (
        (1.268348e-5 * temperatureK - 1.767275512e-2) * temperatureK + 9.98834792
    )
    if temperatureK > 428.15 and pressureBar < 1000.0:
        tf = temperatureK / 300.0 - 1.427166667
        g -= (tf**4.8 + 0.366666e-15 * tf**16) * (
            (
                ((5.01799e-14 * pressureBar - 5.0224e-11) * pressureBar - 1.504074e-7)
                * pressureBar
                + 2.507672e-4
            )
            * pressureBar
            - 0.1003157
        )
    return g


def debye_huckel(
    msol: float, vsolv_cm3: float, epsilon: float, temperatureK: float
) -> float:
    denom = (epsilon * temperatureK) ** 3
    return CDH * math.sqrt(10.0 * msol / vsolv_cm3 / denom)


def born_omega(z: float, omega0: float, born_radius: float, gf: float) -> float:
    if abs(z) < 1e-10:
        return omega0
    absz = abs(z)
    return ETA * z * (z / (born_radius + absz * gf) - 1.0 / (NEUTRAL_RADIUS + gf))


def preprocess_hkf(params: dict, Tr: float = 298.15, Pr: float = 1.0) -> dict:
    # yr is the Born Y reference constant (dZ/dT at Tr, Pr, where Z=-1/epsilon).
    # From tlib.f: data psi, theta, yr, eta/2600d0, 228d0, -5.79865d-5, 694656.968d0/
    # It is NOT 1/(theta-Tr) = -0.01426 — that would cause ~Megajoule errors at high T.
    yr = -5.79865e-5  # Born Y reference constant [K^-1], from tlib.f
    fp_ref = math.log(PSI + Pr)

    b8 = (
        -params["S0"]
        + params["c1"] * math.log(Tr)
        + params["c1"]
        + params["omega0"] * yr
    )
    b8 += math.log(Tr / (Tr - THETA)) * params["c2"] / (THETA * THETA)

    b9 = (-params["omega0"] * yr - params["c1"] + params["S0"]) * Tr
    b9 += params["omega0"] - params["a1"] * Pr - params["a2"] * fp_ref + params["G0"]
    b9 += params["c2"] / THETA

    b10 = -params["a3"] * Pr - params["a4"] * fp_ref
    b11 = -params["c2"] / ((Tr - THETA) * THETA)
    b12 = params["c2"] / (THETA * THETA)
    b13 = -params["c1"] - params["c2"] / (THETA * THETA)

    params.update({"b8": b8, "b9": b9, "b10": b10, "b11": b11, "b12": b12, "b13": b13})

    z = params["charge"]
    if abs(z) > 1e-10:
        q2 = z * z
        # tlib.f line ~11286: b9 = 5d9*eta*z^2 / (1.622323167d9*eta*z + 5d9*omega0)
        # Equivalent to: re = eta*z^2 / (eta*z/3.082 + omega0)
        numerator = 5e9 * ETA * q2
        denominator = 1.622323167e9 * ETA * z + 5e9 * params["omega0"]
        params["bornRadius"] = numerator / denominator
    else:
        params["bornRadius"] = NEUTRAL_RADIUS

    return params


def hkf_gibbs(params: dict, P: float, T: float, epsilon: float, gf: float) -> float:
    ft = T - THETA
    fp = math.log(PSI + P)
    omega = born_omega(params["charge"], params["omega0"], params["bornRadius"], gf)
    G = (
        params["b9"]
        + (params["b8"] + params["b12"] * math.log(ft) + params["b13"] * math.log(T))
        * T
    )
    G += params["b11"] * ft
    G += params["a1"] * P + params["a2"] * fp
    G += (params["a3"] * P + params["a4"] * fp + params["b10"]) / ft
    G += omega * (1.0 / epsilon - 1.0) - params["omega0"] / EPSILON0
    return G


def read_csv(path: Path):
    with path.open() as f:
        return list(csv.DictReader(f))


def write_csv(path: Path, header: list[str], rows: list[list]):
    with path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)


def main() -> None:
    base = Path(__file__).parent / "test"

    # epsh2o
    eps_rows = []
    for r in read_csv(base / "epsh2o_matrix.csv"):
        if r["enabled"].strip() != "1":
            continue
        T = float(r["T_K"])
        vol_cm3 = float(r["vol_cm3"])
        eps = epsh2o(vol_cm3 / 10.0, T)
        eps_rows.append([T, vol_cm3, eps])
    write_csv(
        base / "epsh2o_reference.csv", ["T_K", "vol_cm3", "epsilon_ref"], eps_rows
    )

    # gfunc
    g_rows = []
    for r in read_csv(base / "gfunc_matrix.csv"):
        if r["enabled"].strip() != "1":
            continue
        P = float(r["P_bar"])
        T = float(r["T_K"])
        rho = float(r["rho_gcc"])
        g_rows.append([P, T, rho, gfunc(rho, P, T)])
    write_csv(
        base / "gfunc_reference.csv", ["P_bar", "T_K", "rho_gcc", "g_ref"], g_rows
    )

    # Debye-Hückel
    dh_rows = []
    for r in read_csv(base / "dh_matrix.csv"):
        if r["enabled"].strip() != "1":
            continue
        msol = float(r["msol"])
        vsolv = float(r["vsolv_cm3"])
        eps = float(r["epsilon"])
        T = float(r["T_K"])
        dh_rows.append([msol, vsolv, eps, T, debye_huckel(msol, vsolv, eps, T)])
    write_csv(
        base / "dh_reference.csv",
        ["msol", "vsolv_cm3", "epsilon", "T_K", "adh_ref"],
        dh_rows,
    )

    # Born omega
    born_rows = []
    for r in read_csv(base / "born_matrix.csv"):
        if r["enabled"].strip() != "1":
            continue
        z = float(r["z"])
        omega0 = float(r["omega0"])
        born_radius = float(r["born_radius"])
        P = float(r["P_bar"])
        T = float(r["T_K"])
        rho = float(r["rho_gcc"])
        g = gfunc(rho, P, T)
        born_rows.append(
            [
                z,
                omega0,
                born_radius,
                P,
                T,
                rho,
                g,
                born_omega(z, omega0, born_radius, g),
            ]
        )
    write_csv(
        base / "born_reference.csv",
        ["z", "omega0", "born_radius", "P_bar", "T_K", "rho_gcc", "g", "omega_ref"],
        born_rows,
    )

    # HKF Gibbs
    # The hkf_matrix.csv stores parameters in DEW database SI units (J, Pa).
    # The PerplexHKF engine uses J/bar (same as Perple_X datafiles): only a1 and a3
    # need conversion from J/(mol·Pa) and J·K/(mol·Pa) to J/(mol·bar) and J·K/(mol·bar).
    BAR = 1e5  # Pa/bar

    hkf_rows = []
    for r in read_csv(base / "hkf_matrix.csv"):
        if r["enabled"].strip() != "1":
            continue
        # Read DEW database parameters; apply only the Pa→bar conversion for a1 and a3
        params = {
            "G0": float(r["G0"]),
            "S0": float(r["S0"]),
            "omega0": float(r["omega0"]),
            "charge": float(r["charge"]),
            "a1": float(r["a1"]) * BAR,  # J/(mol·Pa) → J/(mol·bar)
            "a2": float(r["a2"]),
            "a3": float(r["a3"]) * BAR,  # J·K/(mol·Pa) → J·K/(mol·bar)
            "a4": float(r["a4"]),
            "c1": float(r["c1"]),
            "c2": float(r["c2"]),
        }
        P = float(r["P_bar"])
        T = float(r["T_K"])
        epsilon = float(r["epsilon"])
        gf = float(r["gf"])
        params = preprocess_hkf(params)
        G = hkf_gibbs(params, P, T, epsilon, gf)
        hkf_rows.append([r["species"], P, T, epsilon, gf, G])
    write_csv(
        base / "hkf_reference.csv",
        ["species", "P_bar", "T_K", "epsilon", "gf", "G_ref"],
        hkf_rows,
    )

    print(
        "Wrote epsh2o_reference.csv, gfunc_reference.csv, dh_reference.csv, born_reference.csv, hkf_reference.csv"
    )


if __name__ == "__main__":
    main()
