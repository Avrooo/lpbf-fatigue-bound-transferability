"""Internal-consistency screens of FatigueData-AM2022 (raw files only). Usage: python3 -I data_screens.py XLSX JSON"""
import sys, numpy as np, pandas as pd
sys.argv = ['x', sys.argv[1], sys.argv[2], 'all']; exec(open('reproduce.py').read().split("# ---------------- independent model code")[0])
P = par.copy(); P['id'] = P['id'].astype(int); col = lambda s: [c for c in P.columns if str(c).replace('\n', ' ').strip().startswith(s)][0]
print('S-N rows', len(sn), 'curves', sn.id.nunique(), '| parameter rows', len(P)); print(P[col('method of extraction')].value_counts().to_dict())
print('non-positive S or N:', int(((sn.S <= 0) | (sn.N <= 0)).sum()), '| missing:', int(sn[['N', 'S']].isna().sum().sum()), '| duplicate points within a curve:', int(sn.duplicated(['id', 'N', 'S']).sum()))
Rn = pd.to_numeric(P[col('load ratio')], errors='coerce'); print('R unparsable (e.g. "0.2, 0.6"):', int(Rn.isna().sum()), 'of', len(P), '| R>=1 or R<-5:', int(((Rn >= 1) | (Rn < -5)).sum()))
uts = pd.to_numeric(P.set_index('id')[col('ultimate tensile')], errors='coerce'); ys = pd.to_numeric(P.set_index('id')[col('yield strength')], errors='coerce')
m = sn.assign(uts=sn.id.map(uts), ys=sn.id.map(ys)); print('UTS recorded for', f'{100 * m.uts.notna().mean():.0f}% of rows | S>UTS: rows', int((m.S > m.uts).sum()), 'curves', m[m.S > m.uts].id.nunique(), '| S>YS: curves', m[m.S > m.ys].id.nunique())
sl = []
for i, g in sn[sn.ro == 0].groupby('id'):
    if len(g) >= 6 and g.S.nunique() >= 3: sl.append(np.polyfit(np.log10(g.S), np.log10(g.N), 1)[0])
sl = np.array(sl); print('curves with >=6 failures:', len(sl), '| slope > -1 (life not falling with stress):', int((sl > -1).sum()), '| median slope', round(float(np.median(sl)), 2))
