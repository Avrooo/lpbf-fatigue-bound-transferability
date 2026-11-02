"""Learning curves: external breach rate of the pooled censored 5% prediction bound as the number of training papers K grows.
For every target paper (within its alloy and R group) the pooled model is fitted to K randomly chosen OTHER papers of the same group
(10 random draws per K) and also to ALL other papers; the breach rate is measured on all failures of the target paper.
Fixed-target design: every K is evaluated on exactly the same targets (those with at least Kmax other papers in their group).
Weighting: units (alloy, R, target paper) are averaged within each target paper, then equally over target papers; bootstrap over target papers.
Usage: python3 -I learning_curves.py XLSX JSON [SEED]"""
import sys, numpy as np, pandas as pd
X_, J_ = sys.argv[1], sys.argv[2]; LC_SEED = int(sys.argv[3]) if len(sys.argv) > 3 else 51
sys.argv = ['x', X_, J_, 'primary']; exec(open('reproduce.py').read().split("# ---------------- pair records")[0])
rng = np.random.default_rng(LC_SEED); KS = [1, 3, 5, 10, 20]; ND = 10; rec = []
for mat in ALLOYS:
    for R in (-1, 0.1):
        G = D[(D.mat == mat) & (D.R == R)].reset_index(drop=True)
        if G.DOI.nunique() < 6: continue
        y = np.log10(G.N.values); x = np.log10(G.S.values / np.median(G.S.values)); Cm = G[['asb', 'hip', 'hthi']].values.astype(float); cens = G.ro.values.astype(bool); doi = G.DOI.values; P = np.unique(doi)
        def bound(train, te):
            tr = np.isin(doi, train)
            if (~cens[tr]).sum() < 5 or x[tr].std() < 1e-6: return None
            cols = Cm[tr].std(0) > 0; Xtr = np.column_stack([np.ones(tr.sum()), x[tr], Cm[tr][:, cols]]); b, s = censored_lognormal_em(Xtr, y[tr], cens[tr])
            return np.column_stack([np.ones(te.sum()), x[te], Cm[te][:, cols]]) @ b + Z05 * s
        for B in P:
            te = doi == B; others = P[P != B]; nfail = (~cens[te]).sum()
            if nfail == 0: continue
            ref = bound(others, te)
            if ref is not None: rec.append((mat, R, B, 'all', len(others), np.sum(~cens[te] & (y[te] < ref)) / te.sum()))
            for K in KS:
                if K > len(others): continue
                rates = [np.sum(~cens[te] & (y[te] < bd)) / te.sum() for bd in (bound(rng.choice(others, K, replace=False), te) for _ in range(ND)) if bd is not None]
                if rates: rec.append((mat, R, B, K, len(others), float(np.mean(rates))))
        print(mat, R, 'done', flush=True)
r = pd.DataFrame(rec, columns=['mat', 'R', 'B', 'K', 'avail', 'rate']); r['K'] = r.K.astype(str); r.to_csv('learning_curve_units.csv', index=False)
def table(minK, Ks, label, sub=None):
    units = r[r.K == str(minK)][['mat', 'R', 'B']].drop_duplicates()
    if sub: units = units[units.mat == sub]
    cols = {}
    for K in Ks + ['all']:
        q = r[r.K == str(K)].merge(units, on=['mat', 'R', 'B']); cols[K] = q.groupby('B').rate.mean()
    papers = cols['all'].index; M = pd.DataFrame({K: cols[K].reindex(papers) for K in cols}).dropna(); n = len(M); idx = rng.integers(0, n, (4000, n))
    print(f'\n{label}: {len(units)} (alloy,R,paper) units, {n} target papers, identical for every K')
    for K in Ks + ['all']:
        v = M[K].values; bs = v[idx].mean(1); print(f'  K={str(K):>3s}: {100 * v.mean():5.1f}% [{100 * np.percentile(bs, 2.5):4.1f}, {100 * np.percentile(bs, 97.5):5.1f}]')
    for K in Ks:
        dv = (M[K] - M['all']).values; bs = dv[idx].mean(1); print(f'    K={K:>2d} minus all-others: {100 * dv.mean():+5.1f} pts [{100 * np.percentile(bs, 2.5):+5.1f}, {100 * np.percentile(bs, 97.5):+5.1f}]')
table(20, KS, 'FIXED TARGETS A (at least 20 other papers)'); table(10, KS[:4], 'FIXED TARGETS B (at least 10 other papers)')
for mat in ALLOYS: table(10, [1, 5, 10], f'FIXED TARGETS B, {mat}', sub=mat)
