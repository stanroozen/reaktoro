import os
import pandas as pd

BASE = os.path.dirname(os.path.abspath(__file__))

PAIRS = [
    (
        os.path.join(BASE, "diagnostics_chnosz_pourbaix_points.csv"),
        os.path.join(BASE, "diagnostics_reaktoro_pourbaix_points.csv"),
        "Fe_Pourbaix",
    ),
    (
        os.path.join(BASE, "diagnostics_chnosz_mosaic_points.csv"),
        os.path.join(BASE, "diagnostics_reaktoro_mosaic_points.csv"),
        "Fe_Mosaic",
    ),
]

all_rows = []
summary_rows = []

for ch_path, rk_path, case_name in PAIRS:
    ch = pd.read_csv(ch_path)
    rk = pd.read_csv(rk_path)

    merged = ch.merge(rk, on=["case", "pH", "Eh_V"], how="inner")

    merged["pred_match"] = (
        merged["pred_species_chnosz"].astype(str).str.lower()
        == merged["pred_species_reaktoro"].astype(str).str.lower()
    )

    merged["top_match"] = (
        merged["top1_species_chnosz"].astype(str).str.lower()
        == merged["top1_species_reaktoro"].astype(str).str.lower()
    )

    all_rows.append(merged)

    n = len(merged)
    pred_m = int(merged["pred_match"].sum())
    top_m = int(merged["top_match"].sum())
    summary_rows.append(
        {
            "case": case_name,
            "points": n,
            "pred_match_count": pred_m,
            "pred_match_fraction": pred_m / n if n else 0.0,
            "top_match_count": top_m,
            "top_match_fraction": top_m / n if n else 0.0,
        }
    )

all_df = pd.concat(all_rows, ignore_index=True)
summary_df = pd.DataFrame(summary_rows)

all_csv = os.path.join(BASE, "diagnostics_points_merged.csv")
sum_csv = os.path.join(BASE, "diagnostics_points_summary.csv")

all_df.to_csv(all_csv, index=False)
summary_df.to_csv(sum_csv, index=False)

print("Wrote:", all_csv)
print("Wrote:", sum_csv)
print(summary_df.to_string(index=False))
