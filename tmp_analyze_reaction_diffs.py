from pathlib import Path
import sys
import pandas as pd


def main():
    p = Path("reaction_column_diffs.csv")
    if not p.exists():
        print("missing file", file=sys.stderr)
        sys.exit(1)
    df = pd.read_csv(p)
    diff_cols = [c for c in df.columns if c.startswith("d")]
    rows = sorted(
        [(c, df[c].abs().max(), df[c].abs().mean()) for c in diff_cols],
        key=lambda x: x[1],
        reverse=True,
    )
    print("Top diffs by max:")
    for c, m, mean in rows[:10]:
        idx = df[c].abs().idxmax()
        print(
            f"{c}: max={m:.6g}, mean={mean:.6g}, at T={df.loc[idx, 'T_C']} C, P={df.loc[idx, 'P_kb']} kb"
        )
    for key in [
        "dDeltaGro_cal",
        "dG0_H2O_cal",
        "drho_gcm3",
        "deps",
        "dlogK",
        "dDeltaVr_cm3",
    ]:
        if key in df.columns:
            s = df[key].abs()
            idxs = s.nlargest(5).index
            print(f"\nTop 5 {key}:")
            print(df.loc[idxs, ["T_C", "P_kb", key]].to_string(index=False))


if __name__ == "__main__":
    main()
