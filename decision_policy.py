"""Reliability-informed decision test (simple, pre-specified).
For each source curve A and target paper B the engineer holds: the source curve, a pooled model fitted WITHOUT B, and B's planned stresses/process settings (no fatigue lives of B).
Rule(tau): accept the single-study 5% bound if it lies no more than tau log10 cycles above the pooled 5% bound (mean over B's stress levels); otherwise fall back to the pooled bound.
tau = 0 ('never exceed the pooled bound') is fixed in advance; other tau values are shown only as a descriptive frontier.
Baselines: pooled only, single-study only, single-study bound lowered by a fixed margin (the p-grid frontier of coverage_conservatism.py), hierarchical.
Uses the record file written by coverage_conservatism.py (same pairs, same weighting, same two-way bootstrap).
Usage: python3 decision_policy.py RECORDS.csv FRONTIER.csv [seed=11]"""
import sys, numpy as np, pandas as pd
rec = pd.read_csv(sys.argv[1]); FR = pd.read_csv(sys.argv[2]); SEED = int(sys.argv[3]) if len(sys.argv) > 3 else 11
PG = [0.00001, 0.0001, 0.0003, 0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.10, 0.20, 0.30, 0.50]; K = PG.index(0.05)
TAUS = [-0.5, -0.25, -0.1, 0.0, 0.1, 0.25, 0.5, 1.0, 5.0]
d = rec.A.to_frame().copy(); delta = rec[f'h_single_{K}'] - rec[f'h_pooled_{K}']; cols = []
for t in TAUS:
    use = delta <= t
    d[f'r_{t}'] = np.where(use, rec[f'br_single_{K}'], rec[f'br_pooled_{K}']) / rec.n; d[f'h_{t}'] = np.where(use, rec[f'h_single_{K}'], rec[f'h_pooled_{K}'])
    d[f'a_{t}'] = use.astype(float); d[f'ar_{t}'] = use * rec[f'br_single_{K}'] / rec.n; cols += [f'r_{t}', f'h_{t}', f'a_{t}', f'ar_{t}']
d['r_pooled'] = rec[f'br_pooled_{K}'] / rec.n; d['h_pooled'] = rec[f'h_pooled_{K}']; d['r_single'] = rec[f'br_single_{K}'] / rec.n; d['h_single'] = rec[f'h_single_{K}']; cols += ['r_pooled', 'h_pooled', 'r_single', 'h_single']
d['mat'], d['R'], d['B'], d['sid'] = rec.mat, rec.R, rec.B, rec.sid
PV = d.groupby(['mat', 'R', 'A', 'B'])[cols].mean().groupby(['A', 'B']).mean(); Aix = {a: i for i, a in enumerate(sorted({a for a, b in PV.index}))}; Bix = {b: i for i, b in enumerate(sorted({b for a, b in PV.index}))}
Ai = np.array([Aix[a] for a, b in PV.index]); Bi = np.array([Bix[b] for a, b in PV.index]); P = PV.values; nA, nB = len(Aix), len(Bix); ix = {c: i for i, c in enumerate(cols)}
def wavg(wA, wB):
    wa = wA[Ai]; den = np.bincount(Bi, wa, nB); ok = (den > 0) & (wB > 0); return np.array([(wB[ok] * np.bincount(Bi, wa * P[:, j], nB)[ok] / den[ok]).sum() / wB[ok].sum() for j in range(P.shape[1])])
est = wavg(np.ones(nA), np.ones(nB)); rs = np.random.default_rng(SEED); NB = 2000
boots = np.array([wavg(np.bincount(rs.integers(0, nA, nA), minlength=nA).astype(float), np.bincount(rs.integers(0, nB, nB), minlength=nB).astype(float)) for _ in range(NB)])
ci = np.percentile(boots, [2.5, 97.5], axis=0)
rows = []
for t in TAUS:
    r, h, a, ar = (ix[f'{k}_{t}'] for k in ('r', 'h', 'a', 'ar'))
    acc_rate = lambda v: 100 * v[ar] / v[a] if v[a] > 0 else np.nan; bs = np.array([acc_rate(b) for b in boots])
    rows.append(dict(rule=f'accept single if <= pooled + {t} log10', breach_pct=100 * est[r], breach_lo=100 * ci[0, r], breach_hi=100 * ci[1, r], mean_bound_log10N=est[h], bound_lo=ci[0, h], bound_hi=ci[1, h],
                     share_single_accepted_pct=100 * est[a], breach_among_accepted_pct=acc_rate(est), accepted_lo=np.nanpercentile(bs, 2.5), accepted_hi=np.nanpercentile(bs, 97.5)))
for lab, rr, hh in (('pooled only', 'r_pooled', 'h_pooled'), ('single-study only', 'r_single', 'h_single')):
    rows.append(dict(rule=lab, breach_pct=100 * est[ix[rr]], breach_lo=100 * ci[0, ix[rr]], breach_hi=100 * ci[1, ix[rr]], mean_bound_log10N=est[ix[hh]], bound_lo=ci[0, ix[hh]], bound_hi=ci[1, ix[hh]]))
T = pd.DataFrame(rows); T.to_csv('decision_policy_summary.csv', index=False, float_format='%.3f')
# equal-height comparison with the pooled 5% bound and with the single-study fixed-margin (p-grid) frontier
def r_at(rr, hh, target):
    o = np.argsort(hh); return float(np.interp(target, hh[o], rr[o])) if hh.min() <= target <= hh.max() else np.nan
def policy_curve(v): return np.array([100 * v[ix[f'r_{t}']] for t in TAUS]), np.array([v[ix[f'h_{t}']] for t in TAUS])
Href = est[ix['h_pooled']]; sf = FR[FR.method == 'single'].sort_values('p'); out = []
for lab, fn in (('policy tau-frontier', lambda v: r_at(*policy_curve(v), Href)),):
    e = fn(est); bs = np.array([fn(b) for b in boots]); out.append(dict(quantity=f'breach % of {lab} at the pooled-5% bound height ({Href:.3f})', estimate=e, ci_low=np.nanpercentile(bs, 2.5), ci_high=np.nanpercentile(bs, 97.5)))
out.append(dict(quantity='breach % of pooled-only at its own height', estimate=100 * est[ix['r_pooled']], ci_low=100 * ci[0, ix['r_pooled']], ci_high=100 * ci[1, ix['r_pooled']]))
out.append(dict(quantity='breach % of single-study fixed-margin (p-grid) at the pooled-5% bound height', estimate=r_at(sf.breach_pct.values, sf.mean_bound_log10N.values, Href), ci_low=np.nan, ci_high=np.nan))
O = pd.DataFrame(out); O.to_csv('decision_policy_equal_height.csv', index=False, float_format='%.3f')
pd.set_option('display.width', 250); pd.set_option('display.max_columns', 20); print(T.round(3).to_string()); print(O.round(3).to_string())
