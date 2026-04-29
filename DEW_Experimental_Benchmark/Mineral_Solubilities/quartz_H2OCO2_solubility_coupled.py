"""
quartz_H2OCO2_solubility_coupled.py

Quartz (SiO2) solubility in H2O-CO2 mixtures as a function of XCO2.

Fixed conditions: T = 800 C, P = 9 and 10 kbar
Aqueous model  : PerplexDEW (ExtendedDH or Davies)
Gas EOS        : H2O Zhang-Duan 2005 + CO2 Zhang-Duan 2009 (GFSM)
Experimental   : Newton & Manning (2000), Shmulovich et al. (2001)

Redesign vs quartz_H2OCO2_solubility.py
-----------------------------------------
The original script builds THREE separate systems and does 51 redundant DEW
solves (Si solubility is flat vs XCO2 due to anhydrous DEW formulas).

This script uses TWO systems with a minimal, unified loop:

  System 1 -- DEW aqueous (WATER,AQ + Si species + Quartz)
              Solved ONCE at pure H2O to get the per-species baseline
              molalities m0.  No loop over XCO2 needed.

  System 2 -- GFSM gas phase (H2O(g) + CO2(g))
              Evaluated at FIXED composition for each XCO2 via
              ChemicalProps (no equilibrium solve needed -- gas
              composition is the control variable, not computed).
              Provides ln_a_H2O(XCO2) and the pure-H2O reference.

  Inline correction (one line per point):
              m_Si = m0(SiO2_aq)*a^2 + m0(HSiO3-)*a^1
                   + 2*m0(Si2O4)*a^4 + 3*m0(Si3O6)*a^6

Net result: the 51 duplicate DEW solves are eliminated; a_H2O comes
entirely from C++ GFSM; the explicit "post-processing" Python functions
(compute_aH2O_gfsm, compute_activity_corrected_curve, solve_xco2_sweep,
get_baseline_species_molalities) from the original are gone.

Why the inline a^n step is still needed
-----------------------------------------
The DEW database stores Si species with anhydrous formulas (SiO2_aq =
SiO2, not H4SiO4), so H2O does not appear in the dissolution
stoichiometry.  The Gibbs minimiser therefore cannot propagate water
activity changes to Si solubility automatically.  The a_H2O^n term is
the mandatory thermodynamic correction encoding:
  SiO2(s) + 2 H2O <==> H4SiO4_aq  (= SiO2_aq in anhydrous convention)
"""

import os
import sys
import argparse
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
import matplotlib.pyplot as plt

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BENCHMARK_DIR = os.path.dirname(SCRIPT_DIR)
ROOT_DIR = os.path.dirname(BENCHMARK_DIR)

try:
    from reaktoro import *  # noqa: F401,F403
except ModuleNotFoundError:
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
        except ModuleNotFoundError:
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

# =============================================================================
# Configuration
# =============================================================================

T_FIXED_C = 800.0
P_KBAR_LIST = [9.0, 10.0]
XCO2_MODEL = np.concatenate([[0.0], np.linspace(0.005, 0.85, 60)])
N_GAS_MOLES = 1000.0  # moles of total gas for GFSM composition evaluation

# Si species: hydration numbers n and Si counts per mol species.
#   Physical basis: SiO2_aq = SiO2+2H2O (H4SiO4) -> n=2, Si=1
#                   HSiO3-  = SiO2+H2O            -> n=1, Si=1
#                   Si2O4   = 2SiO2+4H2O           -> n=4, Si=2
#                   Si3O6   = 3SiO2+6H2O           -> n=6, Si=3
SI_SPECIES = ["SiO2_aq", "HSiO3-", "Si2O4_aq", "Si3O6_aq"]
SI_COUNT = {"SiO2_aq": 1, "HSiO3-": 1, "Si2O4_aq": 2, "Si3O6_aq": 3}
HYDRATION = {"SiO2_aq": 2, "HSiO3-": 1, "Si2O4_aq": 4, "Si3O6_aq": 6}

PERPLEXDEW_SYMBOLS = (
    "ActivityModelPerplexDEW",
    "ActivityDHModel",
    "ActivityModelPerplexGFSM",
    "ActivityModelParamsPerplexGFSM",
    "PerpleXHybridEosOptions",
    "PerpleXCO2Eos",
    "PerpleXWaterEos",
)

CSV_FILE = os.path.join(SCRIPT_DIR, "quartz_H2OCO2_DEW_testset.csv")
PLOT_TITLE = "Quartz Solubility in H\u2082O-CO\u2082 (T\u2009=\u2009800\u00b0C)"
Y_LABEL = "Quartz Solubility (mol/kg-H\u2082O)"
OUTPUT_PREFIX = "quartz_H2OCO2"


# =============================================================================
# CLI
# =============================================================================


def parse_args():
    p = argparse.ArgumentParser(
        description="Quartz solubility in H2O-CO2: DEW baseline + GFSM a_H2O (no 51 flat solves)."
    )
    p.add_argument(
        "--dh-model",
        default="ExtendedDH",
        choices=["Davies", "ExtendedDH"],
        help="Debye-Huckel variant (default: ExtendedDH).",
    )
    return p.parse_args()


def output_paths(dh_model):
    tag = f"coupled_{dh_model}"
    return {
        "solubility": os.path.join(
            SCRIPT_DIR, f"{OUTPUT_PREFIX}_solubility_vs_xco2_{tag}.png"
        ),
        "residuals": os.path.join(SCRIPT_DIR, f"{OUTPUT_PREFIX}_residuals_{tag}.png"),
    }


# =============================================================================
# Helpers
# =============================================================================


def ensure_perplexdew():
    missing = [n for n in PERPLEXDEW_SYMBOLS if n not in globals()]
    if not missing:
        return
    for _d in [
        os.path.join(ROOT_DIR, "build", "Reaktoro", "Release"),
        os.path.join(ROOT_DIR, "build", "Reaktoro", "Release"),
    ]:
        if not os.path.isdir(_d):
            continue
        if _d in sys.path:
            sys.path.remove(_d)
        sys.path.insert(0, _d)
        sys.modules.pop("reaktoro4py", None)
        importlib.invalidate_caches()
        try:
            _m = importlib.import_module("reaktoro4py")
        except ModuleNotFoundError:
            continue
        for name in PERPLEXDEW_SYMBOLS:
            if hasattr(_m, name):
                globals()[name] = getattr(_m, name)
        if not [n for n in PERPLEXDEW_SYMBOLS if n not in globals()]:
            return
    raise RuntimeError("PerplexDEW symbols missing after search.")


# =============================================================================
# System 1: DEW aqueous (Quartz solubility baseline)
# =============================================================================


def build_dew_system(dew_db, supcrt_db, dh_model="ExtendedDH"):
    ensure_perplexdew()
    mineral_sp = supcrt_db.species("Quartz")
    db2 = Database(dew_db.species())
    db2.addSpecies(mineral_sp)

    aq_str = "WATER,AQ H+ OH- SiO2_aq HSiO3- Si2O4_aq Si3O6_aq"
    aq = AqueousPhase(aq_str)
    _dh = (
        ActivityDHModel.ExtendedDH
        if dh_model == "ExtendedDH"
        else ActivityDHModel.Davies
    )
    aq.setActivityModel(ActivityModelPerplexDEW(_dh))
    mineral = MineralPhase("Quartz")
    system = ChemicalSystem(db2, aq, mineral)
    print(f"  DEW system  : [{aq_str}] ({dh_model}) + Quartz")
    return system


def solve_dew_baseline(dew_system, T_C, P_bar):
    """Solve pure-H2O DEW system once; return per-species Si molalities."""
    solver = EquilibriumSolver(dew_system)
    conditions = EquilibriumConditions(dew_system)
    conditions.temperature(T_C, "celsius")
    conditions.pressure(P_bar, "bar")

    state = ChemicalState(dew_system)
    state.set("WATER,AQ", 1.0, "kg")
    state.set("SiO2_aq", 1e-6, "mol")
    state.set("Quartz", 10.0, "mol")

    result = solver.solve(state, conditions)
    if not result.succeeded():
        raise RuntimeError(f"DEW baseline solve failed at T={T_C}C, P={P_bar}bar.")

    aqprops = AqueousProps(state)
    m0 = {sp: float(aqprops.speciesMolality(sp)) for sp in SI_SPECIES}
    m_total = sum(SI_COUNT[sp] * m0[sp] for sp in SI_SPECIES)
    print(
        f"      m0(SiO2_aq)={m0['SiO2_aq']:.4f}  m0(HSiO3-)={m0['HSiO3-']:.5f}"
        f"  m0(Si2O4)={m0['Si2O4_aq']:.5f}  total_Si={m_total:.4f} mol/kg"
    )
    return m0


# =============================================================================
# System 2: GFSM gas phase (a_H2O at fixed XCO2 composition)
# =============================================================================


def build_gfsm_system(supcrt_db):
    ensure_perplexdew()
    h2og_sp = supcrt_db.species("H2O(g)")
    co2g_sp = supcrt_db.species("CO2(g)")
    gas_db = Database([h2og_sp, co2g_sp])

    gas_phase = GaseousPhase("H2O(g) CO2(g)")
    gfsm_params = ActivityModelParamsPerplexGFSM()
    hybrid_opts = PerpleXHybridEosOptions()
    hybrid_opts.water = PerpleXWaterEos.ZhangDuan05
    hybrid_opts.co2 = PerpleXCO2Eos.ZhangDuan09
    gfsm_params.hybridEosOptions = hybrid_opts
    gas_phase.setActivityModel(ActivityModelPerplexGFSM(gfsm_params))
    system = ChemicalSystem(gas_db, gas_phase)
    print("  GFSM system : [H2O(g) CO2(g)] H2O=ZhangDuan05, CO2=ZhangDuan09")
    return system


def compute_aH2O_gfsm_sweep(gfsm_system, xco2_array, T_C, P_bar, n_total=N_GAS_MOLES):
    """Evaluate GFSM gas-phase activity of H2O at fixed XCO2 compositions.

    Uses ChemicalProps (no equilibrium solve) to evaluate the activity model
    at the given T, P, and gas-phase mole fractions.  Returns a_H2O(XCO2)
    normalised to the pure-H2O reference (XCO2=0).
    """
    # Re-use the same state object to evaluate ChemicalProps at different compositions
    # without running the equilibrium solver -- pure activity model evaluation.
    ln_a_pure = None
    a_h2o_out = {}

    for xco2 in xco2_array:
        y_h2o = 1.0 - xco2
        n_h2o = y_h2o * n_total
        n_co2 = xco2 * n_total

        state = ChemicalState(gfsm_system)
        state.setTemperature(T_C, "celsius")
        state.setPressure(P_bar, "bar")
        state.set("H2O(g)", n_h2o, "mol")
        if n_co2 > 0.0:
            state.set("CO2(g)", n_co2, "mol")

        cprops = ChemicalProps(state)
        ln_a = float(cprops.speciesActivityLn("H2O(g)"))

        if ln_a_pure is None:
            ln_a_pure = ln_a  # reference at XCO2=0

        a_h2o_out[xco2] = float(np.exp(ln_a - ln_a_pure))

    return a_h2o_out


# =============================================================================
# Main solve loop
# =============================================================================


def compute_corrected_solubility(m0, a_h2o_map, xco2_array):
    """Apply inline a_H2O^n correction to baseline molalities.

    m0       : dict {species_name: baseline_molality}  (from DEW at pure H2O)
    a_h2o_map: dict {xco2: a_H2O}                      (from GFSM)
    """
    xco2_out = []
    m_out = []
    for xco2 in xco2_array:
        a = a_h2o_map.get(xco2, np.nan)
        if not np.isfinite(a):
            continue
        m_corr = sum(
            SI_COUNT[sp] * m0[sp] * a ** HYDRATION[sp]
            for sp in SI_SPECIES
            if m0[sp] > 0.0
        )
        if np.isfinite(m_corr) and m_corr > 0.0:
            xco2_out.append(xco2)
            m_out.append(m_corr)
    return np.array(xco2_out), np.array(m_out)


# =============================================================================
# Experimental Data
# =============================================================================


def load_experimental_data():
    df = pd.read_csv(CSV_FILE)
    for col in ("molality_m", "xco2", "P_kbar", "T_C"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df.dropna(subset=["xco2", "molality_m"]).reset_index(drop=True)


# =============================================================================
# Plotting
# =============================================================================

PRESSURE_PALETTE = {
    10.0: {"color": "#d62728", "label": "10 kbar"},
    9.0: {"color": "#1f77b4", "label": "9 kbar"},
}
REF_MARKERS = {
    "Newton_Manning_2000": ("o", "Newton & Manning (2000)  10 kbar"),
    "Shmulovich_Graham_Yardley_2001": ("s", "Shmulovich et al. (2001)  9 kbar"),
}


def plot_solubility(model_results, exp_df, output_path, dh_model):
    fig, ax = plt.subplots(figsize=(9, 6))
    for P_kbar in sorted(model_results.keys()):
        xco2_arr, m_arr = model_results[P_kbar]
        if len(xco2_arr) == 0:
            continue
        sty = PRESSURE_PALETTE.get(P_kbar, {"color": "gray", "label": f"{P_kbar} kbar"})
        ax.plot(
            xco2_arr,
            m_arr,
            color=sty["color"],
            linewidth=2,
            zorder=4,
            label=f"DEW+GFSM ({dh_model}) -- {sty['label']}",
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
        if "Name" in sub.columns:
            for _, row in sub.iterrows():
                ax.annotate(
                    str(row["Name"]),
                    (row["xco2"], row["molality_m"]),
                    textcoords="offset points",
                    xytext=(5, 3),
                    fontsize=7.5,
                    color=color,
                )
    ax.set_yscale("log")
    ax.set_xlabel("X_CO2 (mole fraction)", fontsize=13)
    ax.set_ylabel(Y_LABEL, fontsize=13)
    ax.set_title(
        PLOT_TITLE + f"\nDEW+GFSM ({dh_model}): 1 DEW baseline + GFSM a_H2O sweep",
        fontsize=11,
        fontweight="bold",
    )
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(True, which="both", alpha=0.3)
    ax.set_xlim(left=0.0)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"    Saved: {output_path}")
    plt.close()


def plot_residuals(model_results, exp_df, output_path, dh_model):
    all_rows = []
    print(f"\n    Residuals ({dh_model}):")
    for _, row in exp_df.iterrows():
        P = float(row["P_kbar"])
        if P not in model_results:
            continue
        xco2_arr, m_arr = model_results[P]
        if len(xco2_arr) == 0:
            continue
        m_mod = float(
            np.interp(row["xco2"], xco2_arr, m_arr, left=np.nan, right=np.nan)
        )
        if not np.isfinite(m_mod) or row["molality_m"] <= 0:
            continue
        res = (m_mod - row["molality_m"]) / row["molality_m"] * 100.0
        name = row.get("Name", "") if "Name" in exp_df.columns else ""
        print(
            f"      {name} ({P:.0f} kbar): XCO2={row['xco2']:.3f}  "
            f"exp={row['molality_m']:.4f}  model={m_mod:.4f}  res={res:+.1f}%"
        )
        all_rows.append(
            {
                "xco2": row["xco2"],
                "res": res,
                "name": name,
                "color": PRESSURE_PALETTE.get(P, {"color": "gray"})["color"],
            }
        )
    if not all_rows:
        print("    No valid residuals.")
        return
    rdf = pd.DataFrame(all_rows).sort_values("xco2")
    colors = ["steelblue" if r >= 0 else "tomato" for r in rdf["res"]]
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.bar(
        rdf["xco2"],
        rdf["res"],
        width=0.018,
        color=colors,
        edgecolor="black",
        linewidth=0.5,
    )
    ax.axhline(0, color="black", linewidth=1)
    for level in (50, -50, 100, -100):
        ax.axhline(level, color="gray", linestyle="--", linewidth=0.7, alpha=0.6)
    for _, r in rdf.iterrows():
        va = "bottom" if r["res"] >= 0 else "top"
        ax.text(r["xco2"], r["res"], f"  {r['name']}", ha="center", va=va, fontsize=7.5)
    ax.set_xlabel("X_CO2", fontsize=12)
    ax.set_ylabel("Residual (%)", fontsize=12)
    ax.set_title(
        f"Quartz Solubility -- Model vs Experiment\nDEW+GFSM ({dh_model}), T=800 C",
        fontsize=11,
    )
    ax.set_xlim(left=0.0)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"    Saved: {output_path}")
    plt.close()


# =============================================================================
# Main
# =============================================================================


def main():
    args = parse_args()
    dh_model = args.dh_model
    paths = output_paths(dh_model)

    print("=" * 80)
    print("Quartz Solubility in H2O-CO2  --  DEW baseline + GFSM a_H2O")
    print(f"  DEW model : PerplexDEW ({dh_model})")
    print("  Gas EOS   : GFSM (H2O=ZhangDuan05, CO2=ZhangDuan09)")
    print(f"  a_H2O     : GFSM evaluated at FIXED XCO2 (no aqueous coupling)")
    print(f"  P list    : {P_KBAR_LIST} kbar,  T = {T_FIXED_C} C")
    print(f"  XCO2 grid : {len(XCO2_MODEL)} pts, 0 to {XCO2_MODEL[-1]:.3f}")
    print(
        f"  Improvement over original: 1 DEW solve (not 51) + {len(XCO2_MODEL)} GFSM evals"
    )
    print("=" * 80)

    print("\n[1] Loading databases...")
    dew_db = DEWDatabase("dew2024-aqueous")
    supcrt_db = SupcrtDatabase("supcrtbl")
    print("    DEW2024-aqueous + SUPCRTBL loaded.")

    print("\n[2] Building systems...")
    ensure_perplexdew()
    dew_system = build_dew_system(dew_db, supcrt_db, dh_model=dh_model)
    gfsm_system = build_gfsm_system(supcrt_db)

    print("\n[3] Loading experimental data...")
    exp_df = load_experimental_data()
    print(f"    {len(exp_df)} data points.")
    for ref, grp in exp_df.groupby("reference"):
        print(
            f"      {ref}: {len(grp)} pts, "
            f"P={sorted(grp['P_kbar'].unique().tolist())} kbar, "
            f"XCO2=[{grp['xco2'].min():.3f}-{grp['xco2'].max():.3f}]"
        )

    model_results = {}

    for P_kbar in P_KBAR_LIST:
        P_bar = P_kbar * 1000.0
        print(f"\n[4] P = {P_kbar} kbar ({P_bar:.0f} bar)")

        # Step 1: ONE DEW solve at pure H2O for baseline molalities
        print(f"  DEW baseline (pure H2O):")
        m0 = solve_dew_baseline(dew_system, T_FIXED_C, P_bar)

        # Step 2: GFSM a_H2O sweep -- no equilibrium solve, just activity evaluation
        print(f"  GFSM a_H2O sweep ({len(XCO2_MODEL)} pts, no solver call)...")
        a_h2o_map = compute_aH2O_gfsm_sweep(gfsm_system, XCO2_MODEL, T_FIXED_C, P_bar)
        a_vals = [a_h2o_map[x] for x in XCO2_MODEL if x in a_h2o_map]
        if a_vals:
            print(f"    a_H2O range: {min(a_vals):.4f} to {max(a_vals):.4f}")
            # Print a few diagnostic points
            for xco2 in [0.0, 0.2, 0.4, 0.6, 0.8]:
                closest = min(XCO2_MODEL, key=lambda x: abs(x - xco2))
                print(f"    XCO2={closest:.3f}: a_H2O={a_h2o_map[closest]:.4f}")

        # Step 3: Inline correction
        xco2_out, m_out = compute_corrected_solubility(m0, a_h2o_map, XCO2_MODEL)
        model_results[P_kbar] = (xco2_out, m_out)
        if len(m_out) > 0:
            print(f"  Corrected m range: {m_out.min():.4f} - {m_out.max():.4f} mol/kg")

    print("\n[5] Plotting...")
    plot_solubility(model_results, exp_df, paths["solubility"], dh_model)
    plot_residuals(model_results, exp_df, paths["residuals"], dh_model)
    print("\n  Done.")


if __name__ == "__main__":
    main()

