"""
quartz_H2OCO2_parity_check_nondeps.py

No-extra-dependency parity check for quartz H2O-CO2 benchmark.
Uses only Python stdlib + reaktoro module already built in workspace.
"""

import argparse
import importlib
import math
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
    from reaktoro import *  # noqa: F401,F403
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

T_FIXED_C = 800.0
P_KBAR_LIST = [9.0, 10.0]
XCO2_GRID = [0.0] + [0.005 + i * (0.845 / 59.0) for i in range(60)]
N_GAS_MOLES = 1000.0

SI_SPECIES = ["SiO2_aq", "HSiO3-", "Si2O4_aq", "Si3O6_aq"]
SI_COUNT = {"SiO2_aq": 1, "HSiO3-": 1, "Si2O4_aq": 2, "Si3O6_aq": 3}
HYDRATION = {"SiO2_aq": 2, "HSiO3-": 1, "Si2O4_aq": 4, "Si3O6_aq": 6}


def percentile(sorted_vals, q):
    if not sorted_vals:
        return float("nan")
    if len(sorted_vals) == 1:
        return sorted_vals[0]
    pos = q * (len(sorted_vals) - 1)
    lo = int(math.floor(pos))
    hi = int(math.ceil(pos))
    if lo == hi:
        return sorted_vals[lo]
    frac = pos - lo
    return sorted_vals[lo] * (1.0 - frac) + sorted_vals[hi] * frac


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
        out[xco2] = math.exp(ln_a - ln_ref)
    return out


def hydrated_curve(m0, a_h2o, xco2_array):
    xs = []
    ms = []
    for xco2 in xco2_array:
        a = a_h2o.get(xco2)
        if a is None or not math.isfinite(a):
            continue
        m_si = 0.0
        for sp in SI_SPECIES:
            m_si += SI_COUNT[sp] * m0[sp] * (a ** HYDRATION[sp])
        if math.isfinite(m_si) and m_si > 0.0:
            xs.append(xco2)
            ms.append(m_si)
    return xs, ms


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
        if math.isfinite(m_si) and m_si > 0.0:
            xs.append(xco2)
            ms.append(m_si)

    return xs, ms


def interp_log10(xs, ys, x):
    if len(xs) < 2:
        return None
    if x < xs[0] or x > xs[-1]:
        return None
    for i in range(len(xs) - 1):
        x0, x1 = xs[i], xs[i + 1]
        if x0 <= x <= x1:
            if ys[i] <= 0.0 or ys[i + 1] <= 0.0:
                return None
            t = 0.0 if x1 == x0 else (x - x0) / (x1 - x0)
            y0 = math.log10(ys[i])
            y1 = math.log10(ys[i + 1])
            return y0 * (1.0 - t) + y1 * t
    return None


def compute_metrics(x_ref, y_ref, x_cmp, y_cmp):
    pairs = sorted(zip(x_cmp, y_cmp), key=lambda p: p[0])
    xc = [p[0] for p in pairs]
    yc = [p[1] for p in pairs]

    deltas = []
    for xr, yr in zip(x_ref, y_ref):
        if yr <= 0.0:
            continue
        lc = interp_log10(xc, yc, xr)
        if lc is None:
            continue
        deltas.append(lc - math.log10(yr))

    if len(deltas) < 2:
        return None

    absd = sorted(abs(v) for v in deltas)
    return {
        "median_abs_log10": percentile(absd, 0.5),
        "p95_abs_log10": percentile(absd, 0.95),
        "max_abs_log10": absd[-1],
        "n": len(deltas),
    }


def run_parity(median_gate: float = 0.20, p95_gate: float = 0.45):
    """Run the quartz H2O-CO2 parity check.

    Returns ``(ok, messages)`` where *ok* is True if all pressure points pass and
    *messages* is a list of per-pressure summary strings suitable for reporting.
    """
    dew_db = DEWDatabase("dew2024-aqueous")
    supcrt_db = SupcrtDatabase("supcrtbl")

    dew_sys = build_dew_system(dew_db, supcrt_db)
    gfsm_sys = build_gfsm_system(supcrt_db)
    coupled_sys = build_coupled_system(dew_db, supcrt_db)

    ok_all = True
    messages = []

    for p_kbar in P_KBAR_LIST:
        p_bar = p_kbar * 1000.0

        m0 = solve_dew_baseline(dew_sys, T_FIXED_C, p_bar)
        a_h2o = gfsm_a_h2o_sweep(gfsm_sys, XCO2_GRID, T_FIXED_C, p_bar)
        x_ref, y_ref = hydrated_curve(m0, a_h2o, XCO2_GRID)

        x_cmp, y_cmp = solve_full_sweep(coupled_sys, XCO2_GRID, T_FIXED_C, p_bar)

        metrics = compute_metrics(x_ref, y_ref, x_cmp, y_cmp)
        if metrics is None:
            ok_all = False
            messages.append(f"P={p_kbar:.1f} kbar: FAILED (insufficient overlap)")
            continue

        ok = (
            metrics["median_abs_log10"] <= median_gate
            and metrics["p95_abs_log10"] <= p95_gate
        )
        ok_all = ok_all and ok
        messages.append(
            f"P={p_kbar:.1f} kbar: "
            f"median={metrics['median_abs_log10']:.4f}, "
            f"p95={metrics['p95_abs_log10']:.4f}, "
            f"max={metrics['max_abs_log10']:.4f}, "
            f"n={metrics['n']}, ok={ok}"
        )

    return ok_all, messages


def test_quartz_h2o_co2_parity_nondeps():
    """Pytest entry point: quartz H2O-CO2 full-minimization vs hydrated-curve parity (no extra deps)."""
    ok, messages = run_parity(median_gate=0.20, p95_gate=0.45)
    assert ok, "Quartz H2O-CO2 parity failed:\n" + "\n".join(messages)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--median-gate", type=float, default=0.20)
    parser.add_argument("--p95-gate", type=float, default=0.45)
    args = parser.parse_args()

    print("Quartz H2O-CO2 parity check (no extra Python deps)")
    print(f"Gates: median <= {args.median_gate:.3f}, p95 <= {args.p95_gate:.3f}")

    ok_all, messages = run_parity(args.median_gate, args.p95_gate)
    for msg in messages:
        print(f"  {msg}")
    raise SystemExit(0 if ok_all else 1)


if __name__ == "__main__":
    main()

