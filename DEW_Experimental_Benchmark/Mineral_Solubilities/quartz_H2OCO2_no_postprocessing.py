"""
quartz_H2OCO2_no_postprocessing.py

Same setup as quartz_H2OCO2_solubility_coupled.py, but WITHOUT the a_H2O^n
post-processing step.

Three curves are shown:
  "Raw DEW (flat)"   -- baseline molality from ONE DEW solve at pure H2O,
                        replicated as a constant across all XCO2.
                        This is what comes out of the database with no
                        hydration correction whatsoever.
  "DEW+GFSM corrected"  -- the coupled result from the parent script
                        (a_H2O^n applied, shown for reference).

The figure makes clear that without the correction, the DEW database predicts
a completely flat quartz solubility vs XCO2, because the anhydrous SiO2_aq
formula contains zero H2O units.
"""

import os
import sys
import importlib

# PATH must be set BEFORE any numpy/autodiff imports on Windows.
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
import autodiff  # must be imported before reaktoro4py to register autodiff::Real converters
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BENCHMARK_DIR = os.path.dirname(SCRIPT_DIR)
ROOT_DIR = os.path.dirname(BENCHMARK_DIR)

# Prepend DEW build package dir so "from reaktoro import *" finds the DEW build first.
# Also add the inner reaktoro dir to PATH so Reaktoro.dll is resolved.
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
            "reaktoro4py not found. Check build/build/build/Reaktoro/Release."
        )
    print(f"Using local reaktoro4py from {_pyd_dir}.")

try:
    Warnings.disable(906)
except Exception:
    pass

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
T_FIXED_C = 800.0
P_KBAR_LIST = [9.0, 10.0]
XCO2_GRID = np.concatenate([[0.0], np.linspace(0.005, 0.85, 60)])
N_GAS_MOLES = 1000.0

SI_SPECIES = ["SiO2_aq", "HSiO3-", "Si2O4_aq", "Si3O6_aq"]
SI_COUNT = {"SiO2_aq": 1, "HSiO3-": 1, "Si2O4_aq": 2, "Si3O6_aq": 3}
HYDRATION = {"SiO2_aq": 2, "HSiO3-": 1, "Si2O4_aq": 4, "Si3O6_aq": 6}

CSV_FILE = os.path.join(SCRIPT_DIR, "quartz_H2OCO2_DEW_testset.csv")
OUT_PATH = os.path.join(SCRIPT_DIR, "quartz_H2OCO2_no_postprocessing.png")

PRESSURE_PALETTE = {
    10.0: {"color": "#d62728", "label": "10 kbar"},
    9.0: {"color": "#1f77b4", "label": "9 kbar"},
}


# ---------------------------------------------------------------------------
# Systems
# ---------------------------------------------------------------------------
def build_dew_system(dew_db, supcrt_db, dh_model="ExtendedDH"):
    mineral_sp = supcrt_db.species("Quartz")
    db2 = Database(dew_db.species())
    db2.addSpecies(mineral_sp)
    aq = AqueousPhase("WATER,AQ H+ OH- SiO2_aq HSiO3- Si2O4_aq Si3O6_aq")
    _dh = (
        ActivityDHModel.ExtendedDH
        if dh_model == "ExtendedDH"
        else ActivityDHModel.Davies
    )
    aq.setActivityModel(ActivityModelPerplexDEW(_dh))
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


# ---------------------------------------------------------------------------
# Solves
# ---------------------------------------------------------------------------
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
        raise RuntimeError(f"DEW baseline failed at T={T_C}C, P={P_bar}bar.")
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


# ---------------------------------------------------------------------------
# Corrections
# ---------------------------------------------------------------------------
def no_correction(m0, xco2_array):
    """Return the flat DEW baseline across all XCO2 (no a_H2O^n applied)."""
    m_flat = sum(SI_COUNT[sp] * m0[sp] for sp in SI_SPECIES)
    return np.array(xco2_array), np.full(len(xco2_array), m_flat)


def with_correction(m0, a_h2o_map, xco2_array):
    """Apply the a_H2O^n hydration correction."""
    xs, ms = [], []
    for xco2 in xco2_array:
        a = a_h2o_map.get(xco2, np.nan)
        if not np.isfinite(a):
            continue
        m = sum(
            SI_COUNT[sp] * m0[sp] * a ** HYDRATION[sp]
            for sp in SI_SPECIES
            if m0[sp] > 0.0
        )
        if np.isfinite(m) and m > 0.0:
            xs.append(xco2)
            ms.append(m)
    return np.array(xs), np.array(ms)


# ---------------------------------------------------------------------------
# Experimental data
# ---------------------------------------------------------------------------
def load_exp():
    df = pd.read_csv(CSV_FILE)
    for col in ("molality_m", "xco2", "P_kbar", "T_C"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df.dropna(subset=["xco2", "molality_m"]).reset_index(drop=True)


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
REF_MARKERS = {
    "Newton_Manning_2000": ("o", "Newton & Manning (2000)  10 kbar"),
    "Shmulovich_Graham_Yardley_2001": ("s", "Shmulovich et al. (2001)  9 kbar"),
}


def make_plot(flat_results, corr_results, exp_df):
    fig, ax = plt.subplots(figsize=(10, 6))

    for P_kbar in sorted(flat_results.keys()):
        sty = PRESSURE_PALETTE.get(P_kbar, {"color": "gray", "label": f"{P_kbar} kbar"})

        # --- flat (no correction) ---
        xf, mf = flat_results[P_kbar]
        ax.plot(
            xf,
            mf,
            color=sty["color"],
            linewidth=2,
            linestyle="--",
            label=f"Raw DEW â€” no correction ({sty['label']})",
            zorder=3,
        )

        # --- with a_H2O^n correction ---
        xc, mc = corr_results[P_kbar]
        ax.plot(
            xc,
            mc,
            color=sty["color"],
            linewidth=2,
            linestyle="-",
            label=f"DEW + a_Hâ‚‚O^n correction ({sty['label']})",
            zorder=4,
        )

    for ref, (marker, label) in REF_MARKERS.items():
        sub = exp_df[exp_df["reference"] == ref]
        if sub.empty:
            continue
        p_vals = sub["P_kbar"].unique()
        color = (
            PRESSURE_PALETTE.get(float(p_vals[0]), {"color": "gray"})["color"]
            if len(p_vals) == 1
            else "gray"
        )
        ax.scatter(
            sub["xco2"],
            sub["molality_m"],
            color=color,
            marker=marker,
            s=80,
            zorder=5,
            label=label,
            edgecolors="black",
            linewidths=0.7,
        )

    ax.set_yscale("log")
    ax.set_xlabel("X$_{CO_2}$ (mole fraction)", fontsize=13)
    ax.set_ylabel("Quartz solubility (mol/kg-Hâ‚‚O)", fontsize=13)
    ax.set_title(
        "Quartz Solubility in Hâ‚‚O-COâ‚‚  â€”  T = 800 Â°C\n"
        "Dashed: raw DEW (flat, no correction)  |  Solid: DEW + a$_{H_2O}^{n}$ correction",
        fontsize=11,
        fontweight="bold",
    )
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(True, which="both", alpha=0.3)
    ax.set_xlim(left=0.0)
    plt.tight_layout()
    plt.savefig(OUT_PATH, dpi=150, bbox_inches="tight")
    print(f"Saved: {OUT_PATH}")
    plt.close()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print("=" * 70)
    print("Quartz solubility -- raw DEW (flat) vs DEW + a_H2O^n correction")
    print("=" * 70)

    print("\n[1] Loading databases...")
    dew_db = DEWDatabase("dew2024-aqueous")
    supcrt_db = SupcrtDatabase("supcrtbl")

    print("[2] Building systems (DEW ExtendedDH + GFSM ZhangDuan09)...")
    dew_sys = build_dew_system(dew_db, supcrt_db)
    gfsm_sys = build_gfsm_system(supcrt_db)

    print("[3] Loading experimental data...")
    exp_df = load_exp()
    print(f"    {len(exp_df)} data points.")

    flat_results = {}
    corr_results = {}

    for P_kbar in P_KBAR_LIST:
        P_bar = P_kbar * 1000.0
        print(f"\n[4] P = {P_kbar} kbar")

        print("    DEW baseline solve (pure H2O)...")
        m0 = solve_dew_baseline(dew_sys, T_FIXED_C, P_bar)
        m_total = sum(SI_COUNT[sp] * m0[sp] for sp in SI_SPECIES)
        print(f"    m0 total Si = {m_total:.5f} mol/kg")
        for sp in SI_SPECIES:
            print(f"      {sp}: {m0[sp]:.6f} mol/kg")

        print("    GFSM a_H2O sweep...")
        a_map = gfsm_aH2O_sweep(gfsm_sys, XCO2_GRID, T_FIXED_C, P_bar)
        print(
            f"    a_H2O at XCO2=0.4: {
                a_map[min(XCO2_GRID, key=lambda x: abs(x - 0.4))]:.4f}"
        )

        # No correction: flat line
        flat_results[P_kbar] = no_correction(m0, XCO2_GRID)
        print(
            f"    Flat (no correction) = {flat_results[P_kbar][1][0]:.5f} mol/kg everywhere"
        )

        # With correction
        corr_results[P_kbar] = with_correction(m0, a_map, XCO2_GRID)
        mc = corr_results[P_kbar][1]
        if len(mc):
            print(f"    Corrected range: {mc.min():.5f} â€“ {mc.max():.5f} mol/kg")

    print("\n[5] Plotting...")
    make_plot(flat_results, corr_results, exp_df)
    print("\nDone.")


if __name__ == "__main__":
    main()

