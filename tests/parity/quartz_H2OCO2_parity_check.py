"""
quartz_H2OCO2_parity_check.py

Parity benchmark between:
- hydratedDatabase-style reconstruction
- full minimization sweep

This script compares both curves at 800 C for 9 and 10 kbar and exits non-zero
if configured release gates are not satisfied.
"""

import argparse
import importlib
import os
import sys

if os.name == "nt":
    ep = sys.prefix
    env_paths = [
        ep,
        os.path.join(ep, "Library", "mingw-w64", "bin"),
        os.path.join(ep, "Library", "usr", "bin"),
        os.path.join(ep, "Library", "bin"),
        os.path.join(ep, "Scripts"),
        os.path.join(ep, "bin"),
    ]
    sr = os.environ.get("SystemRoot", r"C:\Windows")
    os.environ["PATH"] = ";".join(
        [
            p
            for p in env_paths
            + [os.path.join(sr, "System32"), sr, os.path.join(sr, "System32", "Wbem")]
            if os.path.isdir(p)
        ]
    )

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BENCHMARK_DIR = os.path.dirname(SCRIPT_DIR)
ROOT_DIR = os.path.dirname(BENCHMARK_DIR)

for _build_pkg in [
    os.path.join(ROOT_DIR, "build", "python", "package"),
    os.path.join(ROOT_DIR, "build", "python", "package"),
    os.path.join(ROOT_DIR, "build", "python", "package"),
]:
    _rkt_inner = os.path.join(_build_pkg, "reaktoro")
    if os.path.isdir(_rkt_inner):
        if _build_pkg not in sys.path:
            sys.path.insert(0, _build_pkg)
        os.environ["PATH"] = _rkt_inner + os.pathsep + os.environ.get("PATH", "")
        break

try:
    from reaktoro import *
except (ModuleNotFoundError, ImportError):
    _pyd_dir = None
    for _d in [
        os.path.join(ROOT_DIR, "build", "Reaktoro", "Release"),
        os.path.join(ROOT_DIR, "build", "Reaktoro", "Release"),
        os.path.join(ROOT_DIR, "build", "Reaktoro", "Debug"),
        os.path.join(ROOT_DIR, "build", "Reaktoro", "Release"),
    ]:
        if not os.path.isdir(_d):
            continue
        sys.path.insert(0, _d)
        sys.modules.pop("reaktoro4py", None)
        try:
            _m = importlib.import_module("reaktoro4py")
            globals().update(
                {k: getattr(_m, k) for k in dir(_m) if not k.startswith("_")}
            )
            _pyd_dir = _d
            break
        except (ModuleNotFoundError, ImportError):
            continue
    if _pyd_dir is None:
        raise ModuleNotFoundError(
            "Reaktoro import failed for all local build candidates."
        )

try:
    Warnings.disable(906)
except Exception:
    pass

if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)

OUT_PLOT = os.path.join(SCRIPT_DIR, "quartz_H2OCO2_parity_delta.png")

T_FIXED_C = 800.0
P_KBAR_LIST = [9.0, 10.0]
XCO2_GRID = np.concatenate([[0.0], np.linspace(0.005, 0.85, 60)])
N_GAS_MOLES = 1000.0

SI_SPECIES = ["SiO2_aq", "HSiO3-", "Si2O4_aq", "Si3O6_aq"]
SI_COUNT = {"SiO2_aq": 1, "HSiO3-": 1, "Si2O4_aq": 2, "Si3O6_aq": 3}
HYDRATION = {"SiO2_aq": 2, "HSiO3-": 1, "Si2O4_aq": 4, "Si3O6_aq": 6}


def build_coupled_system(dew_db, supcrt_db):
    quartz = supcrt_db.species("Quartz")
    h2og = supcrt_db.species("H2O(g)")
    co2g = supcrt_db.species("CO2(g)")

    db = Database(dew_db.species())
    db.addSpecies(quartz)
    db.addSpecies(h2og)
    db.addSpecies(co2g)

    aq = AqueousPhase("WATER,AQ H+ OH- SiO2_aq HSiO3- Si2O4_aq Si3O6_aq")
    aq.setActivityModel(ActivityModelPerplexDEW(ActivityDHModel.ExtendedDH))

    gas = GaseousPhase("H2O(g) CO2(g)")
    params = ActivityModelParamsPerplexGFSM()
    opts = PerpleXHybridEosOptions()
    opts.water = PerpleXWaterEos.ZhangDuan09
    opts.co2 = PerpleXCO2Eos.ZhangDuan09
    params.hybridEosOptions = opts
    gas.setActivityModel(ActivityModelPerplexGFSM(params))

    mineral = MineralPhase("Quartz")
    return ChemicalSystem(db, aq, gas, mineral)


def total_si_molality(aqp):
    total = 0.0
    for sp in SI_SPECIES:
        total += SI_COUNT[sp] * float(aqp.speciesMolality(sp))
    return total


def solve_full_sweep(system, xco2_array, t_c, p_bar):
    solver = EquilibriumSolver(system)
    cond = EquilibriumConditions(system)
    cond.temperature(t_c, "celsius")
    cond.pressure(p_bar, "bar")

    xs = []
    ms = []
    for xco2 in xco2_array:
        state = ChemicalState(system)
        state.set("WATER,AQ", 1.0, "kg")
        state.set("Quartz", 10.0, "mol")
        state.set("SiO2_aq", 1e-6, "mol")

        n_co2 = max(xco2, 1e-12) * N_GAS_MOLES
        n_h2o = max(1.0 - xco2, 1e-12) * N_GAS_MOLES
        state.set("CO2(g)", n_co2, "mol")
        state.set("H2O(g)", n_h2o, "mol")

        res = solver.solve(state, cond)
        if not res.succeeded():
            continue

        aqp = AqueousProps(state)
        m_si = total_si_molality(aqp)
        if np.isfinite(m_si) and m_si > 0.0:
            xs.append(xco2)
            ms.append(m_si)

    return np.array(xs), np.array(ms)


def build_dew_system(dew_db, supcrt_db):
    mineral_sp = supcrt_db.species("Quartz")
    db2 = Database(dew_db.species())
    db2.addSpecies(mineral_sp)
    aq = AqueousPhase("WATER,AQ H+ OH- SiO2_aq HSiO3- Si2O4_aq Si3O6_aq")
    aq.setActivityModel(ActivityModelPerplexDEW(ActivityDHModel.ExtendedDH))
    mineral = MineralPhase("Quartz")
    return ChemicalSystem(db2, aq, mineral)


def build_gfsm_system(supcrt_db):
    h2og = supcrt_db.species("H2O(g)")
    co2g = supcrt_db.species("CO2(g)")
    gas_db = Database([h2og, co2g])
    gas = GaseousPhase("H2O(g) CO2(g)")

    params = ActivityModelParamsPerplexGFSM()
    opts = PerpleXHybridEosOptions()
    opts.water = PerpleXWaterEos.ZhangDuan09
    opts.co2 = PerpleXCO2Eos.ZhangDuan09
    params.hybridEosOptions = opts
    gas.setActivityModel(ActivityModelPerplexGFSM(params))
    return ChemicalSystem(gas_db, gas)


def solve_dew_baseline(dew_system, t_c, p_bar):
    solver = EquilibriumSolver(dew_system)
    cond = EquilibriumConditions(dew_system)
    cond.temperature(t_c, "celsius")
    cond.pressure(p_bar, "bar")

    state = ChemicalState(dew_system)
    state.set("WATER,AQ", 1.0, "kg")
    state.set("SiO2_aq", 1e-6, "mol")
    state.set("Quartz", 10.0, "mol")

    res = solver.solve(state, cond)
    if not res.succeeded():
        raise RuntimeError(f"DEW baseline failed at T={t_c} C, P={p_bar} bar")

    aqp = AqueousProps(state)
    return {sp: float(aqp.speciesMolality(sp)) for sp in SI_SPECIES}


def gfsm_a_h2o_sweep(gfsm_system, xco2_array, t_c, p_bar):
    ln_ref = None
    out = {}
    for xco2 in xco2_array:
        y_h2o = 1.0 - xco2
        state = ChemicalState(gfsm_system)
        state.setTemperature(t_c, "celsius")
        state.setPressure(p_bar, "bar")
        state.set("H2O(g)", y_h2o * N_GAS_MOLES, "mol")
        if xco2 > 0.0:
            state.set("CO2(g)", xco2 * N_GAS_MOLES, "mol")

        cp = ChemicalProps(state)
        ln_a = float(cp.speciesActivityLn("H2O(g)"))
        if ln_ref is None:
            ln_ref = ln_a
        out[xco2] = float(np.exp(ln_a - ln_ref))
    return out


def hydrated_curve(m0, a_h2o, xco2_array):
    xs = []
    ms = []
    for xco2 in xco2_array:
        a = a_h2o.get(xco2, np.nan)
        if not np.isfinite(a):
            continue
        m_si = 0.0
        for sp in SI_SPECIES:
            m_si += SI_COUNT[sp] * m0[sp] * (a ** HYDRATION[sp])
        if np.isfinite(m_si) and m_si > 0.0:
            xs.append(xco2)
            ms.append(m_si)
    return np.array(xs), np.array(ms)


def log10_interp(x_src, y_src, x_target):
    mask = np.isfinite(x_src) & np.isfinite(y_src) & (y_src > 0.0)
    x = np.asarray(x_src[mask])
    y = np.asarray(y_src[mask])
    if x.size < 2:
        return np.full_like(x_target, np.nan, dtype=float)
    order = np.argsort(x)
    x = x[order]
    y = y[order]
    lx = np.log10(y)
    vals = np.interp(x_target, x, lx, left=np.nan, right=np.nan)
    vals[x_target < x[0]] = np.nan
    vals[x_target > x[-1]] = np.nan
    return vals


def compute_metrics(x_ref, y_ref, x_cmp, y_cmp):
    mask_ref = np.isfinite(x_ref) & np.isfinite(y_ref) & (y_ref > 0.0)
    xr = np.asarray(x_ref[mask_ref])
    yr = np.asarray(y_ref[mask_ref])
    if xr.size < 2:
        return None

    lref = np.log10(yr)
    lcmp = log10_interp(np.asarray(x_cmp), np.asarray(y_cmp), xr)
    valid = np.isfinite(lcmp)
    if np.count_nonzero(valid) < 2:
        return None

    delta = lcmp[valid] - lref[valid]
    abs_delta = np.abs(delta)
    return {
        "x": xr[valid],
        "delta_log10": delta,
        "median_abs_log10": float(np.median(abs_delta)),
        "p95_abs_log10": float(np.percentile(abs_delta, 95.0)),
        "max_abs_log10": float(np.max(abs_delta)),
        "n": int(np.count_nonzero(valid)),
    }


def run_parity(median_gate, p95_gate):
    dew_db = DEWDatabase("dew2024-aqueous")
    supcrt_db = SupcrtDatabase("supcrtbl")

    dew_sys = build_dew_system(dew_db, supcrt_db)
    gfsm_sys = build_gfsm_system(supcrt_db)
    coupled_sys = build_coupled_system(dew_db, supcrt_db)

    summary = {}
    all_ok = True

    for p_kbar in P_KBAR_LIST:
        p_bar = p_kbar * 1000.0

        m0 = solve_dew_baseline(dew_sys, T_FIXED_C, p_bar)
        a_h2o = gfsm_a_h2o_sweep(gfsm_sys, XCO2_GRID, T_FIXED_C, p_bar)
        x_ref, y_ref = hydrated_curve(m0, a_h2o, XCO2_GRID)

        x_cmp, y_cmp = solve_full_sweep(coupled_sys, XCO2_GRID, T_FIXED_C, p_bar)

        metrics = compute_metrics(x_ref, y_ref, x_cmp, y_cmp)
        if metrics is None:
            summary[p_kbar] = {"ok": False, "error": "insufficient overlap"}
            all_ok = False
            continue

        ok = (
            metrics["median_abs_log10"] <= median_gate
            and metrics["p95_abs_log10"] <= p95_gate
        )
        all_ok = all_ok and ok
        summary[p_kbar] = {
            "ok": ok,
            "median_abs_log10": metrics["median_abs_log10"],
            "p95_abs_log10": metrics["p95_abs_log10"],
            "max_abs_log10": metrics["max_abs_log10"],
            "n": metrics["n"],
            "x": metrics["x"],
            "delta": metrics["delta_log10"],
        }

    return all_ok, summary


def plot_delta(summary):
    fig, ax = plt.subplots(figsize=(10, 6))
    for p_kbar, data in sorted(summary.items()):
        if "x" not in data:
            continue
        color = "#1f77b4" if np.isclose(p_kbar, 9.0) else "#d62728"
        ls = "-" if data.get("ok", False) else "--"
        ax.plot(
            data["x"],
            data["delta"],
            color=color,
            linestyle=ls,
            linewidth=2.0,
            label=f"{p_kbar:.0f} kbar",
        )

    ax.axhline(0.0, color="black", linewidth=1.0, alpha=0.7)
    ax.set_xlabel("X_CO2")
    ax.set_ylabel("Delta log10(m): full minimization - hydrated curve")
    ax.set_title("Quartz H2O-CO2 parity residuals")
    ax.grid(True, alpha=0.3)
    ax.legend()
    plt.tight_layout()
    plt.savefig(OUT_PLOT, dpi=150, bbox_inches="tight")
    plt.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--median-gate", type=float, default=0.20)
    parser.add_argument("--p95-gate", type=float, default=0.45)
    args = parser.parse_args()

    ok, summary = run_parity(args.median_gate, args.p95_gate)
    plot_delta(summary)

    print("Quartz H2O-CO2 parity check")
    print(
        f"Gates: median_abs_log10 <= {args.median_gate:.3f}, p95_abs_log10 <= {args.p95_gate:.3f}"
    )
    for p_kbar in sorted(summary.keys()):
        data = summary[p_kbar]
        if "error" in data:
            print(f"  P={p_kbar:.1f} kbar: FAILED ({data['error']})")
            continue
        print(
            f"  P={p_kbar:.1f} kbar: median={data['median_abs_log10']:.4f}, "
            f"p95={data['p95_abs_log10']:.4f}, max={data['max_abs_log10']:.4f}, "
            f"n={data['n']}, ok={data['ok']}"
        )

    print(f"Saved: {OUT_PLOT}")
    raise SystemExit(0 if ok else 1)


def test_quartz_h2o_co2_parity():
    """Pytest entry point: quartz H2O-CO2 full-minimization vs hydrated-curve parity."""
    ok, summary = run_parity(median_gate=0.20, p95_gate=0.45)
    failing = [
        f"P={p:.1f} kbar: median={d.get('median_abs_log10', float('nan')):.4f}, "
        f"p95={d.get('p95_abs_log10', float('nan')):.4f}"
        for p, d in sorted(summary.items())
        if not d.get("ok", False)
    ]
    assert ok, "Quartz H2O-CO2 parity failed:\n" + "\n".join(failing)


if __name__ == "__main__":
    main()
