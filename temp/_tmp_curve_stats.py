import pandas as pd
p = r"c:\Users\stanroozen\Documents\Projects\reaktoro-dev\reaktoro\DEW_Experimental_Benchmark\quartz_curves_dew24_PerplexDEW.csv"
df = pd.read_csv(p)
# numeric
for c in ['P_kbar','T_C','molality']:
    df[c]=pd.to_numeric(df[c], errors='coerce')
print('rows', len(df))
print('non_nan_molality', df['molality'].notna().sum())
print('nan_molality', df['molality'].isna().sum())
print('isobar rows', (df['curve_type']=='isobar').sum())
print('psat rows', (df['curve_type']=='psat').sum())
print('isobar non_nan', df.loc[df.curve_type=='isobar','molality'].notna().sum())
print('psat non_nan', df.loc[df.curve_type=='psat','molality'].notna().sum())
# by pressure for isobar
iso = df[df.curve_type=='isobar'].copy()
agg = iso.groupby('P_kbar')['molality'].agg(total='size', ok=lambda s: s.notna().sum())
agg['fail']=agg['total']-agg['ok']
agg['ok_pct']=100*agg['ok']/agg['total']
print('\nby pressure (worst 20 by ok_pct):')
print(agg.sort_values(['ok_pct','P_kbar']).head(20).to_string())
print('\nfully failed pressures:', (agg['ok']==0).sum())
print('fully successful pressures:', (agg['fail']==0).sum())
print('constant-value check (nunique on non-nan, first 10):')
nu = iso.groupby('P_kbar')['molality'].nunique(dropna=True).sort_values().head(10)
print(nu.to_string())
