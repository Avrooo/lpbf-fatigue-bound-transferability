"""Breach rates by alloy (main-text Table 3): single-study 5% prediction, single-study 95/95, pooled censored and pooled failures-only,
from the saved pair-curve records results/pairs_primary.csv, with the staged weights of Section 3.5 and a two-way paper bootstrap
(source and target papers resampled independently within the alloy, 2,000 resamples).
Usage: python3 alloy_results.py [seed=11]   (writes results/alloy_results.csv)"""
import sys, numpy as np, pandas as pd
SEED = int(sys.argv[1]) if len(sys.argv) > 1 else 11; NB = 2000
d = pd.read_csv('results/pairs_primary.csv'); M = ['single_pred5', 'single_tol9595', 'pooled_censored_pred5', 'pooled_failures_only_pred5']; rows = []
for a in ['Ti-6Al-4V', 'AlSi10Mg', 'IN718', '316L']:
    g = d[d.mat == a].copy()
    for m in M: g[m] = g[m] / g.n
    pv = g.groupby(['mat', 'R', 'A', 'B'])[M].mean().groupby(['A', 'B']).mean(); Ai = pd.factorize(pv.index.get_level_values('A'))[0]; Bi = pd.factorize(pv.index.get_level_values('B'))[0]; nA, nB = Ai.max() + 1, Bi.max() + 1; P = pv.values
    def est(wA, wB):
        wa = wA[Ai]; den = np.bincount(Bi, wa, nB); ok = (den > 0) & (wB > 0)
        return np.array([(wB[ok] * np.bincount(Bi, wa * P[:, j], nB)[ok] / den[ok]).sum() / wB[ok].sum() for j in range(P.shape[1])]) * 100
    e = est(np.ones(nA), np.ones(nB)); rs = np.random.default_rng(SEED)
    bs = np.array([est(np.bincount(rs.integers(0, nA, nA), minlength=nA).astype(float), np.bincount(rs.integers(0, nB, nB), minlength=nB).astype(float)) for _ in range(NB)])
    for j, m in enumerate(M): rows.append(dict(alloy=a, method=m, papers=g.B.nunique(), estimate=e[j], ci_low=np.nanpercentile(bs[:, j], 2.5), ci_high=np.nanpercentile(bs[:, j], 97.5)))
T = pd.DataFrame(rows); T.to_csv('results/alloy_results.csv', index=False, float_format='%.3f'); print(T.round(1).to_string())
