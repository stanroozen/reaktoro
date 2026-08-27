"""
quartz_H2OCO2_solubility_fullMinimization.py

Quartz solubility in H2O-CO2 at T=800 C, P=9 and 10 kbar.

This benchmark performs a full Gibbs free energy minimization at every XCO2
point in a coupled system containing:
- Aqueous phase (PerplexDEW)
- Gas phase (PerplexGFSM Zhang-Duan 2009)
- Quartz mineral phase

Unlike hydratedDatabase/postprocessing workflows, this script does not apply an
analytic a_H2O scaling after the solve. The plotted curve is taken directly
from minimizer outputs at each composition.

Output:
- quartz_H2OCO2_fullMinimization.png
"""

import os
import sys
import importlib

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

# CRITICAL: import autodiff BEFORE reaktoro4py to register pybind11 type casters
# for autodiff::real. Without this, plain Python floats fail in cond.temperature() etc.
try:
    import autodiff  # noqa: F401
except ModuleNotFoundError:
    pass

import numpy as np
import pandas as pd
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
    print(f"Using local reaktoro4py from {_pyd_dir}.")


def _ensure_required_symbols():
    required = [
        "DEWDatabase",
        "SupcrtDatabase",
        "ActivityModelPerplexDEW",
        "ActivityModelPerplexGFSM",
        "ActivityDHModel",
        "PerpleXWaterEos",
        "PerpleXCO2Eos",
    ]
    return all(name in globals() for name in required)


def _load_local_reaktoro4py():
    _pyd_dir = None
    for _d in [
        os.path.join(ROOT_DIR, "build", "Reaktoro", "Release"),
        os.path.join(ROOT_DIR, "build", "Reaktoro", "Release"),
        os.path.join(ROOT_DIR, "build", "Reaktoro", "Release"),
    ]:
        if not os.path.isdir(_d):
            continue
        if _d in sys.path:
            sys.path.remove(_d)
        sys.path.insert(0, _d)
        sys.modules.pop("reaktoro4py", None)
        try:
            _m = importlib.import_module("reaktoro4py")
            globals().update({k: getattr(_m, k) for k in dir(_m) if not k.startswith("_")})
            _pyd_dir = _d
            if _ensure_required_symbols():
                break
        except (ModuleNotFoundError, ImportError):
            continue
    if _pyd_dir is None or not _ensure_required_symbols():
        raise ModuleNotFoundError(
            "Required Reaktoro symbols unavailable after local reaktoro4py fallback."
        )
    print(f"Using local reaktoro4py from {_pyd_dir}.")


if not _ensure_required_symbols():
    _load_local_reaktoro4py()

try:
    Warnings.disable(906)
except Exception:
    pass


def to_real(value):
    try:
        return autodiff.real(value)
    except Exception:
        return value

T_FIXED_C = 800.0
P_KBAR_LIST = [9.0, 10.0]
XCO2_GRID = np.concatenate([[0.0], np.linspace(0.005, 0.85, 60)])
N_GAS_MOLES = 1000.0

SI_SPECIES = ["SiO2(aq)", "HSiO3-(aq)", "Si2O4(aq)", "Si3O6(aq)"]
SI_COUNT = {"SiO2(aq)": 1, "HSiO3-(aq)": 1, "Si2O4(aq)": 2, "Si3O6(aq)": 3}

CSV_FILE = os.path.join(SCRIPT_DIR, "quartz_H2OCO2_DEW_testset.csv")
OUT_PATH = os.path.join(SCRIPT_DIR, "quartz_H2OCO2_fullMinimization.png")

PRESSURE_STYLE = {
    9.0: {"color": "#1f77b4", "label": "9 kbar"},
    10.0: {"color": "#d62728", "label": "10 kbar"},
}


def build_coupled_system(dew_db, supcrt_db):
    quartz = supcrt_db.species("Quartz")
    h2og = supcrt_db.species("H2O(g)")
    co2g = supcrt_db.species("CO2(g)")

    db = Database(dew_db.species())
    db.addSpecies(quartz)
    db.addSpecies(h2og)
    db.addSpecies(co2g)

    aq = AqueousPhase("H2O(aq) H+(aq) OH-(aq) SiO2(aq) HSiO3-(aq) Si2O4(aq) Si3O6(aq)")
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


def solve_full_sweep(system, xco2_array, T_C, P_bar):
    solver = EquilibriumSolver(system)
    cond = EquilibriumConditions(system)
    cond.temperature(to_real(T_C), "celsius")
    cond.pressure(to_real(P_bar), "bar")

    xs = []
    ms = []

    for xco2 in xco2_array:
        state = ChemicalState(system)
        state.set("H2O(aq)", 1.0, "kg")
        state.set("Quartz", 10.0, "mol")
        state.set("SiO2(aq)", 1e-6, "mol")

        n_co2 = max(xco2, 1e-12) * N_GAS_MOLES
        n_h2o = max(1.0 - xco2, 1e-12) * N_GAS_MOLES
        state.set("CO2(g)", n_co2, "mol")
        state.set("H2O(g)", n_h2o, "mol")

        res = solver.solve(state, cond)
        if not res.succeeded():
            print(f"    Skipped XCO2={xco2:.4f} (no convergence)")
            continue

        aqp = AqueousProps(state)
        m_si = total_si_molality(aqp)
        if np.isfinite(m_si) and m_si > 0.0:
            xs.append(xco2)
            ms.append(m_si)

    return np.array(xs), np.array(ms)


def load_exp():
    df = pd.read_csv(CSV_FILE)
    for c in ("xco2", "molality_m", "P_kbar"):
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df.dropna(subset=["xco2", "molality_m"]).reset_index(drop=True)


def plot_results(model, exp_df):
    fig, ax = plt.subplots(figsize=(10, 6))

    for P_kbar in sorted(model.keys()):
        style = PRESSURE_STYLE.get(P_kbar, {"color": "gray", "label": f"{P_kbar} kbar"})
        x, y = model[P_kbar]
        ax.plot(
            x,
            y,
            color=style["color"],
            linewidth=2.2,
            linestyle="-",
            label=f"Full minimization ({style['label']})",
        )

        sub = exp_df[np.isclose(exp_df["P_kbar"], P_kbar)]
        if not sub.empty:
            ax.scatter(
                sub["xco2"],
                sub["molality_m"],
                color=style["color"],
                s=70,
                marker="o",
                edgecolors="black",
                linewidths=0.7,
                label=f"Experiment ({style['label']})",
            )

    ax.set_yscale("log")
    ax.set_xlabel("X_CO2 (mole fraction)")
    ax.set_ylabel("Quartz solubility (mol/kg-H2O)")
    ax.set_title("Quartz Solubility in H2O-CO2 at 800 C\nFull Gibbs minimization sweep")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=9)
    ax.set_xlim(0.0, 0.85)
    plt.tight_layout()
    plt.savefig(OUT_PATH, dpi=150, bbox_inches="tight")
    plt.close()


def main():
    print("=" * 72)
    print("Quartz H2O-CO2 full minimization benchmark")
    print("=" * 72)

    print("[1] Loading databases...")
    dew_db = DEWDatabase("dew2024-aqueous")
    supcrt_db = SupcrtDatabase("supcrtbl")

    print("[2] Building coupled system...")
    system = build_coupled_system(dew_db, supcrt_db)

    print("[3] Loading experimental data...")
    exp_df = load_exp()
    print(f"    {len(exp_df)} experimental points loaded.")

    results = {}
    for P_kbar in P_KBAR_LIST:
        P_bar = P_kbar * 1000.0
        print(f"[4] Full sweep at P={P_kbar:.1f} kbar...")
        results[P_kbar] = solve_full_sweep(system, XCO2_GRID, T_FIXED_C, P_bar)

    print("[5] Plotting...")
    plot_results(results, exp_df)
    print(f"Saved: {OUT_PATH}")


if __name__ == "__main__":
    main()

