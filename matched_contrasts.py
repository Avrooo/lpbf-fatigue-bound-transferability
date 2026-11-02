"""Paired coverage-conservatism contrasts: at matched mean bound life H, breach(single) - breach(pooled) (and hierarchical - pooled), with paired publication-level bootstrap intervals.
Both curves are read off the same bootstrap resample (two-way resampling of source and target papers, 2,000 resamples), so the interval is for the difference, not for the two rates separately.
Heights are interpolated along each method's frontier (p-grid); a height outside a method's sampled range gives NaN.
Usage: python3 matched_contrasts.py RECORDS.csv [seed=11]   (writes matched_contrasts.csv)"""
import sys, numpy as np, pandas as pd
rec = pd.read_csv(sys.argv[1]); SEED = int(sys.argv[2]) if len(sys.argv) > 2 else 11
PG = [0.00001, 0.0001, 0.0003, 0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.10, 0.20, 0.30, 0.50]; K5 = PG.index(0.05); METH = ['single', 'pooled', 'hier']
cols = [f'{kind}_{m}_{k}' for m in METH for k in range(len(PG)) for kind in ('r', 'h')]
for m in METH:
    for k in range(len(PG)): rec[f'r_{m}_{k}'] = rec[f'br_{m}_{k}'] / rec.n
PV = rec.groupby(['mat', 'R', 'A', 'B'])[cols].mean().groupby(['A', 'B']).mean(); Aix = {a: i for i, a in enumerate(sorted({a for a, b in PV.index}))}; Bix = {b: i for i, b in enumerate(sorted({b for a, b in PV.index}))}
Ai = np.array([Aix[a] for a, b in PV.index]); Bi = np.array([Bix[b] for a, b in PV.index]); P = PV.values; nA, nB = len(Aix), len(Bix); ix = {c: i for i, c in enumerate(cols)}
def wavg(wA, wB):
    wa = wA[Ai]; den = np.bincount(Bi, wa, nB); ok = (den > 0) & (wB > 0); return np.array([(wB[ok] * np.bincount(Bi, wa * P[:, j], nB)[ok] / den[ok]).sum() / wB[ok].sum() for j in range(P.shape[1])])
est = wavg(np.ones(nA), np.ones(nB)); rs = np.random.default_rng(SEED); NB = 2000
boots = np.array([wavg(np.bincount(rs.integers(0, nA, nA), minlength=nA).astype(float), np.bincount(rs.integers(0, nB, nB), minlength=nB).astype(float)) for _ in range(NB)])
def at_h(v, m, H):
    r = 100 * np.array([v[ix[f'r_{m}_{k}']] for k in range(len(PG))]); h = np.array([v[ix[f'h_{m}_{k}']] for k in range(len(PG))]); o = np.argsort(h)
    return float(np.interp(H, h[o], r[o])) if h.min() <= H <= h.max() else np.nan
# Reference heights are the mean heights of the pooled bound at p = 1%, 2.5%, 5%, 10%, 20%; the height is recomputed inside every bootstrap resample, so the intervals are paired.
rows = []
for p in (0.01, 0.025, 0.05, 0.10, 0.20):
    kk = PG.index(p); Hf = lambda v: v[ix[f'h_pooled_{kk}']]
    for a, b in (('single', 'pooled'), ('hier', 'pooled')):
        f = lambda v: at_h(v, a, Hf(v)) - at_h(v, b, Hf(v)); e = f(est); bs = np.array([f(x) for x in boots])
        ra = np.array([at_h(x, a, Hf(x)) for x in boots]); rb = np.array([at_h(x, b, Hf(x)) for x in boots])
        rows.append(dict(pooled_p=p, ref_height_log10N=round(Hf(est), 3), contrast=f'{a} minus {b} (breach %, pp)', estimate=e, ci_low=np.nanpercentile(bs, 2.5), ci_high=np.nanpercentile(bs, 97.5),
                         method_pct=at_h(est, a, Hf(est)), method_lo=np.nanpercentile(ra, 2.5), method_hi=np.nanpercentile(ra, 97.5), pooled_pct=at_h(est, b, Hf(est)), pooled_lo=np.nanpercentile(rb, 2.5), pooled_hi=np.nanpercentile(rb, 97.5), boot_defined=int(np.isfinite(bs).sum())))
T = pd.DataFrame(rows); T.to_csv('matched_contrasts.csv', index=False, float_format='%.3f'); pd.set_option('display.width', 300); print(T.round(2).to_string())
