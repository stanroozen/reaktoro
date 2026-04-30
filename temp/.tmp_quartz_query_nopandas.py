import csv, math
from pathlib import Path

res_path = Path('DEW_Experimental_Benchmark/quartz_residuals_dew24_PerplexDEW_Davies_uncertainty.csv')
cur_path = Path('DEW_Experimental_Benchmark/quartz_curves_dew24_PerplexDEW_Davies_uncertainty.csv')

def to_float(x):
    if x is None:
        return float('nan')
    s = str(x).strip()
    if s == '' or s.lower() == 'nan':
        return float('nan')
    return float(s)

with res_path.open(newline='', encoding='utf-8-sig') as f:
    res_rows = list(csv.DictReader(f))
with cur_path.open(newline='', encoding='utf-8-sig') as f:
    cur_rows = list(csv.DictReader(f))

# normalize keys access by expected exact names
for r in res_rows:
    r['T_C_f'] = to_float(r.get('T_C'))
    r['P_kbar_f'] = to_float(r.get('P_kbar'))
    r['exp_molality_m_f'] = to_float(r.get('exp_molality_m'))
for r in cur_rows:
    r['T_C_f'] = to_float(r.get('T_C'))
    r['P_kbar_f'] = to_float(r.get('P_kbar'))
    r['pred_molality_med_f'] = to_float(r.get('pred_molality_med'))
    r['pred_molality_lo_f'] = to_float(r.get('pred_molality_lo'))
    r['pred_molality_hi_f'] = to_float(r.get('pred_molality_hi'))


def nearest_T_match(rows, T):
    if not rows:
        return None
    best = min(rows, key=lambda rr: abs(rr['T_C_f'] - T))
    if abs(best['T_C_f'] - T) <= 1e-8:
        return best
    return None

out = []
for idx, rr in enumerate(res_rows):
    T = rr['T_C_f']
    P = rr['P_kbar_f']
    expv = rr['exp_molality_m_f']
    prow = None
    label = 'Psat'

    if math.isnan(P):
        subset = [c for c in cur_rows if str(c.get('curve_type','')).strip() == 'Psat']
        prow = nearest_T_match(subset, T)
        pdisp = 'Psat'
    else:
        valid_p = [c for c in cur_rows if not math.isnan(c['P_kbar_f'])]
        if valid_p:
            p_near = min(valid_p, key=lambda c: abs(c['P_kbar_f'] - P))['P_kbar_f']
            if abs(p_near - P) <= 1e-6:
                ctype = f"P_{p_near:g}kbar"
                subset = [c for c in cur_rows if str(c.get('curve_type','')).strip() == ctype]
                if not subset:
                    subset = [c for c in cur_rows if (not math.isnan(c['P_kbar_f']) and abs(c['P_kbar_f']-p_near)<=1e-6)]
                prow = nearest_T_match(subset, T)
                label = ctype
            else:
                label = f"P_{P:g}kbar"
        pdisp = f"{P:g}" if prow is not None else label

    if prow is None:
        med = lo = hi = float('nan')
    else:
        med = prow['pred_molality_med_f']
        lo = prow['pred_molality_lo_f']
        hi = prow['pred_molality_hi_f']

    width = hi - lo if (not math.isnan(hi) and not math.isnan(lo)) else float('nan')

    out.append((idx, T, pdisp, expv, med, lo, hi, width))

def fmt(x):
    return 'NaN' if (isinstance(x,float) and math.isnan(x)) else f"{x:.8g}" if isinstance(x,float) else str(x)

header = ['idx','T_C','P_kbar(or Psat)','exp_molality_m','pred_molality_med','pred_molality_lo','pred_molality_hi','pred_range_width']
print(' | '.join(header))
for row in out:
    print(' | '.join([str(row[0]), fmt(row[1]), row[2], fmt(row[3]), fmt(row[4]), fmt(row[5]), fmt(row[6]), fmt(row[7])]))

expvals = [r[3] for r in out if not (isinstance(r[3],float) and math.isnan(r[3]))]
print(f"experimental range over the 4 points: min={min(expvals):.8g}, max={max(expvals):.8g}")
