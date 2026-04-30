import pandas as pd
import numpy as np

res_path = r"DEW_Experimental_Benchmark/quartz_residuals_dew24_PerplexDEW_Davies_uncertainty.csv"
cur_path = r"DEW_Experimental_Benchmark/quartz_curves_dew24_PerplexDEW_Davies_uncertainty.csv"

res = pd.read_csv(res_path)
cur = pd.read_csv(cur_path)

res_cols = {c.lower(): c for c in res.columns}
cur_cols = {c.lower(): c for c in cur.columns}

T_res_col = res_cols.get('t_c', 'T_C')
P_res_col = res_cols.get('p_kbar', 'P_kbar')
exp_col = res_cols.get('exp_molality_m', 'exp_molality_m')

T_cur_col = cur_cols.get('t_c', 'T_C')
P_cur_col = cur_cols.get('p_kbar', 'P_kbar')
ctype_col = cur_cols.get('curve_type', 'curve_type')
med_col = cur_cols.get('pred_molality_med', 'pred_molality_med')
lo_col = cur_cols.get('pred_molality_lo', 'pred_molality_lo')
hi_col = cur_cols.get('pred_molality_hi', 'pred_molality_hi')

def pick_row_by_T(dfsub, T_target, tcol):
    if dfsub.empty:
        return None
    dT = (dfsub[tcol] - T_target).abs()
    i = dT.idxmin()
    if float(dT.loc[i]) <= 1e-8:
        return dfsub.loc[i]
    return None

rows = []
for idx, r in res.iterrows():
    T = float(r[T_res_col])
    P = r[P_res_col]

    if pd.isna(P):
        subset = cur[cur[ctype_col] == 'Psat']
        prow = pick_row_by_T(subset, T, T_cur_col)
        label = 'Psat'
    else:
        P = float(P)
        pvals = pd.to_numeric(cur[P_cur_col], errors='coerce')
        dP = (pvals - P).abs()
        iP = dP.idxmin()
        if pd.isna(iP) or float(dP.loc[iP]) > 1e-6:
            prow = None
            label = f'P_{P:g}kbar'
        else:
            p_near = float(pvals.loc[iP])
            ctype = f'P_{p_near:g}kbar'
            subset = cur[cur[ctype_col] == ctype]
            if subset.empty:
                subset = cur[(pd.to_numeric(cur[P_cur_col], errors='coerce') - p_near).abs() <= 1e-6]
            prow = pick_row_by_T(subset, T, T_cur_col)
            label = ctype

    if prow is None:
        med = lo = hi = np.nan
    else:
        med = float(prow[med_col])
        lo = float(prow[lo_col])
        hi = float(prow[hi_col])

    rows.append({
        'idx': int(idx),
        'T_C': T,
        'P_kbar(or Psat)': label if pd.isna(P) else (f"{float(P):g}" if prow is not None else label),
        'exp_molality_m': float(r[exp_col]) if not pd.isna(r[exp_col]) else np.nan,
        'pred_molality_med': med,
        'pred_molality_lo': lo,
        'pred_molality_hi': hi,
        'pred_range_width': (hi - lo) if (not np.isnan(hi) and not np.isnan(lo)) else np.nan,
    })

out = pd.DataFrame(rows)
fmt = out.copy()
for c in ['T_C','exp_molality_m','pred_molality_med','pred_molality_lo','pred_molality_hi','pred_range_width']:
    fmt[c] = fmt[c].map(lambda x: f"{x:.8g}" if pd.notna(x) else "NaN")

print(fmt.to_string(index=False))
exp_vals = pd.to_numeric(out['exp_molality_m'], errors='coerce')
print(f"\nexperimental range over the 4 points: min={exp_vals.min():.8g}, max={exp_vals.max():.8g}")
