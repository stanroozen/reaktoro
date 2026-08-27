import os
import argparse
import subprocess
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DEW_DIR = os.path.normpath(os.path.join(SCRIPT_DIR, "..", "quartz_DEW"))
PERPLEX_DIR = os.path.normpath(os.path.join(SCRIPT_DIR, "..", "quartz_PerplexDEW"))
QUARTZ_SCRIPT = os.path.join(DEW_DIR, "quartz_solubility_analysis_v2_dew24.py")
OUTPUT_PREFIX = "quartz"


# Human-readable display names for the two backends
DEW_LABEL = "DEW (ZD2005 EOS + Sverjensky2014 ε + DEW HKF integral)"
PERPLEX_LABEL_FMT = "{tag} (ZD2005 EOS + Sverjensky2014 ε + Perple_X HKF path)"
PERPLEX_DEW_ACTIVITY_LABEL_FMT = (
    "{tag} (ZD2005 EOS + Sverjensky2014 ε + DEW activity path on PerplexDEW run)"
)


def perplex_display_label(perplex_tag):
    if str(perplex_tag).endswith("_DEWActivity"):
        return PERPLEX_DEW_ACTIVITY_LABEL_FMT.format(tag=perplex_tag)
    return PERPLEX_LABEL_FMT.format(tag=perplex_tag)


def assumptions_text(perplex_tag, plot_min_pred, plot_max_molality, plot_max_abs_delta):
    return (
        "Molality shown: total dissolved Si (element basis, mol/kg H2O).\n"
        "System basis: pure-water aqueous phase (no H2O-CO2 mixture), Quartz mineral.\n"
        "Species basis: WATER,AQ + H+ OH- + SiO2_aq + selected Si complexes.\n"
        f"Backend A: {DEW_LABEL}\n"
        f"Backend B: {perplex_display_label(perplex_tag)}\n"
        "Shared models: ZD2005 water EOS, Sverjensky 2014 dielectric constant.\n"
        "Difference: HKF Born-solvation integral path (DEW C++ vs Perple_X Fortran).\n"
        "Hydrated-mixture workflow: OFF in this benchmark script.\n"
        f"Plot filters: min molality={plot_min_pred:g}, max molality={plot_max_molality:g}, max |delta|={plot_max_abs_delta:g}."
    )


def backend_file(kind, backend, perplex_tag="PerplexDEW_Davies_DEWActivity"):
    backend_name = backend
    if backend == "PerplexDEW":
        backend_name = perplex_tag

    filename = f"{OUTPUT_PREFIX}_{kind}_dew24_{backend_name}.csv"
    primary_dir = DEW_DIR if backend == "DEW" else PERPLEX_DIR
    primary = os.path.join(primary_dir, filename)
    fallback = os.path.join(DEW_DIR, filename)
    return primary if os.path.exists(primary) else fallback


def _perplex_tag_settings(perplex_tag):
    tag = str(perplex_tag or "PerplexDEW_Davies").strip()
    dh = "ExtendedDH" if "ExtendedDH" in tag else "Davies"
    activity_model = "DEW" if tag.endswith("_DEWActivity") else "PerplexDEW"
    return dh, activity_model


def run_backend(backend, perplex_tag="PerplexDEW_Davies_DEWActivity"):
    cmd = [sys.executable, QUARTZ_SCRIPT, "--backend", backend]
    if backend == "PerplexDEW":
        dh, activity_model = _perplex_tag_settings(perplex_tag)
        cmd.extend(["--dh-model", dh, "--perplex-activity-model", activity_model])
    print(f"Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True, cwd=DEW_DIR)


def load_outputs(backend, perplex_tag="PerplexDEW_Davies_DEWActivity"):
    residuals = pd.read_csv(backend_file("residuals", backend, perplex_tag))
    curves = pd.read_csv(backend_file("curves", backend, perplex_tag))
    return residuals, curves


def summarize_residuals(df, label):
    valid = df.dropna(subset=["abs_diff", "rel_diff_pct"])
    mae = valid["abs_diff"].abs().mean()
    rmse = np.sqrt((valid["abs_diff"] ** 2).mean())
    mape = valid["rel_diff_pct"].abs().mean()
    print(f"{label}: N={len(valid)}, MAE={mae:.4e}, RMSE={rmse:.4e}, MAPE={mape:.2f}%")
    return {"N": len(valid), "MAE": mae, "RMSE": rmse, "MAPE": mape}


def compare_residuals(
    dew,
    perplex,
    plot_min_pred=1e-4,
    perplex_tag="PerplexDEW_Davies_DEWActivity",
    plot_max_molality=100.0,
    plot_max_abs_delta=10.0,
):
    key_cols = ["T_C", "P_kbar", "reference", "experiment_type", "molality_m"]
    cols = key_cols + ["predicted_m", "abs_diff", "rel_diff_pct"]
    merged = dew[cols].merge(
        perplex[cols],
        on=key_cols,
        how="inner",
        suffixes=("_dew", "_perplex"),
    )

    merged["predicted_delta"] = (
        merged["predicted_m_perplex"] - merged["predicted_m_dew"]
    )
    merged["abs_diff_delta"] = merged["abs_diff_perplex"] - merged["abs_diff_dew"]

    out_csv = os.path.join(
        SCRIPT_DIR, f"{OUTPUT_PREFIX}_backend_residuals_comparison.csv"
    )
    merged.to_csv(out_csv, index=False)
    print(f"Saved residual comparison: {out_csv}")

    # Remove near-zero predicted points from plotting because they dominate
    # the axis and hide meaningful trends (critical-region spike around ~350 C).
    plot_df = merged[
        (merged["predicted_m_dew"] >= plot_min_pred)
        & (merged["predicted_m_perplex"] >= plot_min_pred)
    ].copy()
    removed = len(merged) - len(plot_df)
    print(
        f"Plot filter: excluded {removed} rows with predicted_m < {plot_min_pred:g} in either backend"
    )

    fig, axes = plt.subplots(2, 1, figsize=(12, 9), sharex=True)
    axes[0].scatter(plot_df["T_C"], plot_df["predicted_delta"], s=18, alpha=0.7)
    axes[0].axhline(0.0, color="gray", linestyle="--", linewidth=1)
    axes[0].set_ylabel(
        f"{perplex_display_label(perplex_tag)}\n− {DEW_LABEL}\npredicted molality (mol/kg)"
    )
    axes[0].grid(True, alpha=0.3)

    axes[1].scatter(plot_df["T_C"], plot_df["abs_diff_delta"], s=18, alpha=0.7)
    axes[1].axhline(0.0, color="gray", linestyle="--", linewidth=1)
    axes[1].set_ylabel(
        f"{perplex_display_label(perplex_tag)}\n− {DEW_LABEL}\n|residual| difference"
    )
    axes[1].set_xlabel("Temperature (°C)")
    axes[1].grid(True, alpha=0.3)

    axes[0].text(
        0.01,
        0.99,
        assumptions_text(
            perplex_tag,
            plot_min_pred,
            plot_max_molality,
            plot_max_abs_delta,
        ),
        transform=axes[0].transAxes,
        fontsize=8,
        verticalalignment="top",
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
    )

    out_png = os.path.join(
        SCRIPT_DIR, f"{OUTPUT_PREFIX}_backend_residuals_comparison.png"
    )
    plt.tight_layout()
    plt.savefig(out_png, dpi=250)
    print(f"Saved residual comparison plot: {out_png}")


def compare_curves(
    dew,
    perplex,
    plot_min_pred=1e-4,
    plot_max_molality=100.0,
    plot_max_abs_delta=10.0,
    perplex_tag="PerplexDEW_Davies_DEWActivity",
):
    key_cols = ["curve_type", "P_kbar", "T_C"]
    merged = dew.merge(
        perplex,
        on=key_cols,
        how="inner",
        suffixes=("_dew", "_perplex"),
    )
    merged["molality_delta"] = merged["molality_perplex"] - merged["molality_dew"]

    outliers_raw = merged[
        (merged["molality_dew"] > plot_max_molality)
        | (merged["molality_perplex"] > plot_max_molality)
        | (np.abs(merged["molality_delta"]) > plot_max_abs_delta)
    ].copy()
    if len(outliers_raw) > 0:
        print(
            f"Curve blowout diagnostic: {len(outliers_raw)} raw outlier rows detected before filtering"
        )
        by_p = (
            outliers_raw.groupby("P_kbar", dropna=False)
            .size()
            .sort_values(ascending=False)
            .head(5)
        )
        for p, n in by_p.items():
            print(f"  P={p}: N={int(n)} outlier rows")
        worst = outliers_raw.loc[np.abs(outliers_raw["molality_delta"]).idxmax()]
        print(
            "  Worst row: "
            f"T={worst['T_C']}, P={worst['P_kbar']}, "
            f"DEW={worst['molality_dew']}, Perplex={worst['molality_perplex']}, "
            f"delta={worst['molality_delta']}"
        )

    # Filter near-zero molality rows from curve plotting/comparison output to
    # avoid critical-region spikes dominating the plot scale.
    merged_plot = merged[
        (merged["molality_dew"] >= plot_min_pred)
        & (merged["molality_perplex"] >= plot_min_pred)
        & np.isfinite(merged["molality_dew"])
        & np.isfinite(merged["molality_perplex"])
        & np.isfinite(merged["molality_delta"])
        & (merged["molality_dew"] <= plot_max_molality)
        & (merged["molality_perplex"] <= plot_max_molality)
        & (np.abs(merged["molality_delta"]) <= plot_max_abs_delta)
    ].copy()

    out_csv = os.path.join(SCRIPT_DIR, f"{OUTPUT_PREFIX}_backend_curves_comparison.csv")
    merged_plot.to_csv(out_csv, index=False)
    print(f"Saved curve comparison: {out_csv}")
    print(
        "Curve plot filter: excluded "
        f"{len(merged) - len(merged_plot)} rows "
        f"(min molality={plot_min_pred:g}, max molality={plot_max_molality:g}, max |delta|={plot_max_abs_delta:g})"
    )

    isobar = merged_plot[merged_plot["curve_type"] == "isobar"].copy()
    pressures = sorted(isobar["P_kbar"].dropna().unique())

    fig, ax = plt.subplots(figsize=(13.5, 8))
    for p in pressures:
        subset = isobar[isobar["P_kbar"] == p]
        if len(subset) == 0:
            continue
        ax.plot(
            subset["T_C"],
            subset["molality_delta"],
            linewidth=1.3,
            label=f"{p:.3g} kbar",
        )

    ax.axhline(0.0, color="black", linestyle="--", linewidth=1)
    ax.set_xlabel("Temperature (°C)")
    ax.set_ylabel(
        f"{perplex_display_label(perplex_tag)}\n− {DEW_LABEL}\nmolality (mol/kg H₂O)"
    )
    ax.set_title("Isobar curve molality difference: PerplexDEW − DEW")
    ax.grid(True, alpha=0.3)
    if len(pressures) <= 12:
        ax.legend(ncol=2, fontsize=8)

    ax.text(
        0.01,
        0.99,
        assumptions_text(
            perplex_tag,
            plot_min_pred,
            plot_max_molality,
            plot_max_abs_delta,
        ),
        transform=ax.transAxes,
        fontsize=8,
        verticalalignment="top",
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
    )

    out_png = os.path.join(SCRIPT_DIR, f"{OUTPUT_PREFIX}_backend_curves_comparison.png")
    plt.tight_layout()
    plt.savefig(out_png, dpi=250)
    print(f"Saved curve comparison plot: {out_png}")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run quartz benchmark with DEW and PerplexDEW and compare outputs."
    )
    parser.add_argument(
        "--skip-run",
        action="store_true",
        help="Skip running simulations and only compare existing CSV outputs.",
    )
    parser.add_argument(
        "--perplex-tag",
        default="PerplexDEW_Davies_DEWActivity",
        choices=[
            "PerplexDEW_Davies_DEWActivity",
            "PerplexDEW_ExtendedDH_DEWActivity",
            "PerplexDEW_Davies",
            "PerplexDEW_ExtendedDH",
        ],
        help="Exact PerplexDEW output tag to compare against DEW.",
    )
    parser.add_argument(
        "--plot-min-pred",
        type=float,
        default=1e-4,
        help="Exclude rows from comparison plots when predicted_m is below this threshold in either backend.",
    )
    parser.add_argument(
        "--plot-max-molality",
        type=float,
        default=100.0,
        help="Exclude curve rows when either backend molality exceeds this threshold.",
    )
    parser.add_argument(
        "--plot-max-abs-delta",
        type=float,
        default=10.0,
        help="Exclude curve rows when |PerplexDEW-DEW| exceeds this threshold.",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    if not args.skip_run:
        run_backend("DEW")
        run_backend("PerplexDEW", args.perplex_tag)

    dew_resid, dew_curves = load_outputs("DEW", args.perplex_tag)
    perplex_resid, perplex_curves = load_outputs("PerplexDEW", args.perplex_tag)

    print("Residual metrics:")
    summarize_residuals(dew_resid, DEW_LABEL)
    summarize_residuals(perplex_resid, perplex_display_label(args.perplex_tag))

    compare_residuals(
        dew_resid,
        perplex_resid,
        plot_min_pred=args.plot_min_pred,
        perplex_tag=args.perplex_tag,
        plot_max_molality=args.plot_max_molality,
        plot_max_abs_delta=args.plot_max_abs_delta,
    )
    compare_curves(
        dew_curves,
        perplex_curves,
        plot_min_pred=args.plot_min_pred,
        plot_max_molality=args.plot_max_molality,
        plot_max_abs_delta=args.plot_max_abs_delta,
        perplex_tag=args.perplex_tag,
    )


if __name__ == "__main__":
    main()
