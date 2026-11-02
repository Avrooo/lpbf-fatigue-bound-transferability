"""Robustness: extreme-study influence (leave-one-paper-out jackknife of the headline estimates), source-study size, and stress ratio.
Uses the record file written by coverage_conservatism.py (strict external predictions); same weighting as the headline analysis.
Usage: python3 influence_checks.py RECORDS.csv XLSX JSON"""
import sys, numpy as np, pandas as pd
rec = pd.read_csv(sys.argv[1]); X_, J_ = sys.argv[2], sys.argv[3]; sys.argv = ['x', X_, J_, 'primary']; exec(open('reproduce.py').read().split("# ---------------- pair records")[0])
PG = [0.00001, 0.0001, 0.0003, 0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.10, 0.20, 0.30, 0.50]; K = PG.index(0.05)
nf = D[D.ro == 0].groupby('id').size(); rec['src_nfail'] = rec.sid.map(nf)
rec['r_single'] = rec[f'br_single_{K}'] / rec.n; rec['r_pooled'] = rec[f'br_pooled_{K}'] / rec.n; rec['r_hier'] = rec[f'br_hier_{K}'] / rec.n
def head(R_, drop_src=None, drop_tgt=None):
    if drop_src is not None: R_ = R_[R_.A != drop_src]
    if drop_tgt is not None: R_ = R_[R_.B != drop_tgt]
    pv = R_.groupby(['mat', 'R', 'A', 'B'])[['r_single', 'r_pooled', 'r_hier']].mean().groupby(['A', 'B']).mean()
    perB = pv.reset_index().groupby('B')[['r_single', 'r_pooled', 'r_hier']]
    # source-paper mean within each target, then target mean (same two-stage structure as the headline)
    t = pv.groupby(level='B').mean(); return 100 * t.mean()
base = head(rec); print('headline (single, pooled, hier):', base.round(2).tolist())
rows = []
for a in sorted(rec.A.unique()): h = head(rec, drop_src=a); rows.append(('source', a, *h.values))
for b in sorted(rec.B.unique()): h = head(rec, drop_tgt=b); rows.append(('target', b, *h.values))
J = pd.DataFrame(rows, columns=['role', 'paper', 'single', 'pooled', 'hier']); J['gap'] = J.single - J.pooled; J.to_csv('influence_jackknife.csv', index=False, float_format='%.3f')
for role in ('source', 'target'):
    s = J[J.role == role]; print(f'{role} jackknife ({len(s)} papers): single {s.single.min():.2f}-{s.single.max():.2f}, pooled {s.pooled.min():.2f}-{s.pooled.max():.2f}, hier {s.hier.min():.2f}-{s.hier.max():.2f}, gap {s.gap.min():.2f}-{s.gap.max():.2f}')
    print(s.assign(dev=(s.gap - (base['r_single'] - base['r_pooled'])).abs()).sort_values('dev', ascending=False).head(3)[['paper', 'single', 'pooled', 'gap']].round(2).to_string(index=False))
bins = [(6, 7), (8, 11), (12, 19), (20, 999)]; sz = []
for lo, hi in bins:
    s = rec[(rec.src_nfail >= lo) & (rec.src_nfail <= hi)]; h = head(s); sz.append(dict(source_failures=f'{lo}-{hi}' if hi < 999 else f'>={lo}', source_curves=s.sid.nunique(), source_papers=s.A.nunique(), single=h['r_single'], pooled=h['r_pooled'], hier=h['r_hier'], gap=h['r_single'] - h['r_pooled']))
for Rv in (-1.0, 0.1):
    s = rec[rec.R == Rv]; h = head(s); sz.append(dict(source_failures=f'R={Rv}', source_curves=s.sid.nunique(), source_papers=s.A.nunique(), single=h['r_single'], pooled=h['r_pooled'], hier=h['r_hier'], gap=h['r_single'] - h['r_pooled']))
S = pd.DataFrame(sz); S.to_csv('influence_by_size_and_R.csv', index=False, float_format='%.2f'); print(S.round(2).to_string(index=False))
