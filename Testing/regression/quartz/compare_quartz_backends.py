import os
import argparse
import subprocess
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
QUARTZ_SCRIPT = os.path.join(SCRIPT_DIR, "quartz_solubility_analysis_v2_dew24.py")
OUTPUT_PREFIX = "quartz"


def backend_file(kind, backend):
    return os.path.join(SCRIPT_DIR, f"{OUTPUT_PREFIX}_{kind}_dew24_{backend}.csv")


def run_backend(backend):
    cmd = [sys.executable, QUARTZ_SCRIPT, "--backend", backend]
    print(f"Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True, cwd=SCRIPT_DIR)


def load_outputs(backend):
    residuals = pd.read_csv(backend_file("residuals", backend))
    curves = pd.read_csv(backend_file("curves", backend))
    return residuals, curves


def summarize_residuals(df, label):
    valid = df.dropna(subset=["abs_diff", "rel_diff_pct"])
    mae = valid["abs_diff"].abs().mean()
    rmse = np.sqrt((valid["abs_diff"] ** 2).mean())
    mape = valid["rel_diff_pct"].abs().mean()
    print(f"{label}: N={len(valid)}, MAE={mae:.4e}, RMSE={rmse:.4e}, MAPE={mape:.2f}%")
    return {"N": len(valid), "MAE": mae, "RMSE": rmse, "MAPE": mape}


def compare_residuals(dew, perplex):
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

    fig, axes = plt.subplots(2, 1, figsize=(12, 9), sharex=True)
    axes[0].scatter(merged["T_C"], merged["predicted_delta"], s=18, alpha=0.7)
    axes[0].axhline(0.0, color="gray", linestyle="--", linewidth=1)
    axes[0].set_ylabel("PerplexDEW - DEW predicted m")
    axes[0].grid(True, alpha=0.3)

    axes[1].scatter(merged["T_C"], merged["abs_diff_delta"], s=18, alpha=0.7)
    axes[1].axhline(0.0, color="gray", linestyle="--", linewidth=1)
    axes[1].set_ylabel("PerplexDEW - DEW residual")
    axes[1].set_xlabel("Temperature (C)")
    axes[1].grid(True, alpha=0.3)

    out_png = os.path.join(
        SCRIPT_DIR, f"{OUTPUT_PREFIX}_backend_residuals_comparison.png"
    )
    plt.tight_layout()
    plt.savefig(out_png, dpi=250)
    print(f"Saved residual comparison plot: {out_png}")


def compare_curves(dew, perplex):
    key_cols = ["curve_type", "P_kbar", "T_C"]
    merged = dew.merge(
        perplex,
        on=key_cols,
        how="inner",
        suffixes=("_dew", "_perplex"),
    )
    merged["molality_delta"] = merged["molality_perplex"] - merged["molality_dew"]

    out_csv = os.path.join(SCRIPT_DIR, f"{OUTPUT_PREFIX}_backend_curves_comparison.csv")
    merged.to_csv(out_csv, index=False)
    print(f"Saved curve comparison: {out_csv}")

    isobar = merged[merged["curve_type"] == "isobar"].copy()
    pressures = sorted(isobar["P_kbar"].dropna().unique())

    fig, ax = plt.subplots(figsize=(12, 7))
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
    ax.set_xlabel("Temperature (C)")
    ax.set_ylabel("PerplexDEW - DEW molality")
    ax.set_title("Curve Delta by Pressure")
    ax.grid(True, alpha=0.3)
    if len(pressures) <= 12:
        ax.legend(ncol=2, fontsize=8)

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
    return parser.parse_args()


def main():
    args = parse_args()

    if not args.skip_run:
        run_backend("DEW")
        run_backend("PerplexDEW")

    dew_resid, dew_curves = load_outputs("DEW")
    perplex_resid, perplex_curves = load_outputs("PerplexDEW")

    print("Residual metrics:")
    summarize_residuals(dew_resid, "DEW")
    summarize_residuals(perplex_resid, "PerplexDEW")

    compare_residuals(dew_resid, perplex_resid)
    compare_curves(dew_curves, perplex_curves)


if __name__ == "__main__":
    main()
