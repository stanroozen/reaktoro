"""
quartz_H2OCO2_solubility_hydratedDatabase.py

Quartz solubility in H2O-CO2 at T=800 C, P=9 and 10 kbar.

This script is the hydrated-database style benchmark counterpart:
- DEW baseline Si species molalities are obtained once at pure H2O.
- GFSM provides a_H2O(XCO2).
- The hydrated-equivalent Si molality is computed as:
    m_Si = sum_i nu_i * m0_i * a_H2O^(n_i)
  where n_i is the hydration number of each anhydrous DEW Si species.

Output:
- quartz_H2OCO2_hydratedDatabase.png
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

import numpy as np
import autodiff
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BENCHMARK_DIR = os.path.dirname(SCRIPT_DIR)
ROOT_DIR = os.path.dirname(BENCHMARK_DIR)

# Prefer build package first.
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

try:
    Warnings.disable(906)
except Exception:
    pass

T_FIXED_C = 800.0
P_KBAR_LIST = [9.0, 10.0]
XCO2_GRID = np.concatenate([[0.0], np.linspace(0.005, 0.85, 60)])
N_GAS_MOLES = 1000.0

SI_SPECIES = ["SiO2_aq", "HSiO3-", "Si2O4_aq", "Si3O6_aq"]
SI_COUNT = {"SiO2_aq": 1, "HSiO3-": 1, "Si2O4_aq": 2, "Si3O6_aq": 3}
HYDRATION = {"SiO2_aq": 2, "HSiO3-": 1, "Si2O4_aq": 4, "Si3O6_aq": 6}

CSV_FILE = os.path.join(SCRIPT_DIR, "quartz_H2OCO2_DEW_testset.csv")
OUT_PATH = os.path.join(SCRIPT_DIR, "quartz_H2OCO2_hydratedDatabase.png")

PRESSURE_STYLE = {
    9.0: {"color": "#1f77b4", "label": "9 kbar"},
    10.0: {"color": "#d62728", "label": "10 kbar"},
}


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


def solve_dew_baseline(dew_system, T_C, P_bar):
    solver = EquilibriumSolver(dew_system)
    cond = EquilibriumConditions(dew_system)
    cond.temperature(T_C, "celsius")
    cond.pressure(P_bar, "bar")

    state = ChemicalState(dew_system)
    state.set("WATER,AQ", 1.0, "kg")
    state.set("SiO2_aq", 1e-6, "mol")
    state.set("Quartz", 10.0, "mol")

    res = solver.solve(state, cond)
    if not res.succeeded():
        raise RuntimeError(f"DEW baseline failed at T={T_C} C, P={P_bar} bar")

    aqp = AqueousProps(state)
    return {sp: float(aqp.speciesMolality(sp)) for sp in SI_SPECIES}


def gfsm_aH2O_sweep(gfsm_system, xco2_array, T_C, P_bar):
    ln_ref = None
    out = {}
    for xco2 in xco2_array:
        y_h2o = 1.0 - xco2
        state = ChemicalState(gfsm_system)
        state.setTemperature(T_C, "celsius")
        state.setPressure(P_bar, "bar")
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
        m = 0.0
        for sp in SI_SPECIES:
            m += SI_COUNT[sp] * m0[sp] * (a ** HYDRATION[sp])
        if np.isfinite(m) and m > 0.0:
            xs.append(xco2)
            ms.append(m)
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
            linewidth=2.0,
            linestyle="-",
            label=f"Hydrated database equivalent ({style['label']})",
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
    ax.set_title("Quartz Solubility in H2O-CO2 at 800 C\nHydrated-database benchmark")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=9)
    ax.set_xlim(0.0, 0.85)
    plt.tight_layout()
    plt.savefig(OUT_PATH, dpi=150, bbox_inches="tight")
    plt.close()


def main():
    print("=" * 72)
    print("Quartz H2O-CO2 hydratedDatabase benchmark")
    print("=" * 72)

    print("[1] Loading databases...")
    dew_db = DEWDatabase("dew2024-aqueous")
    supcrt_db = SupcrtDatabase("supcrtbl")

    print("[2] Building systems...")
    dew_sys = build_dew_system(dew_db, supcrt_db)
    gfsm_sys = build_gfsm_system(supcrt_db)

    print("[3] Loading experimental data...")
    exp_df = load_exp()
    print(f"    {len(exp_df)} experimental points loaded.")

    results = {}

    for P_kbar in P_KBAR_LIST:
        P_bar = P_kbar * 1000.0
        print(f"[4] P={P_kbar:.1f} kbar")

        m0 = solve_dew_baseline(dew_sys, T_FIXED_C, P_bar)
        a_h2o = gfsm_aH2O_sweep(gfsm_sys, XCO2_GRID, T_FIXED_C, P_bar)
        results[P_kbar] = hydrated_curve(m0, a_h2o, XCO2_GRID)

        x4 = 0.4
        if x4 in a_h2o:
            print(f"    a_H2O(XCO2=0.4) = {a_h2o[x4]:.4f}")

    print("[5] Plotting...")
    plot_results(results, exp_df)
    print(f"Saved: {OUT_PATH}")


if __name__ == "__main__":
    main()

