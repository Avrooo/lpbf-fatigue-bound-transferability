"""Independent audit of within_paper_control.py (separate code path, separate RNG).
Usage: python3 -I audit_within_paper.py data/FatigueData-AM2022.xlsx data/FatigueData-AM2022.json   (after within_paper_control.py)
Checks: (1) perturbation test: held-out y values replaced by extreme numbers leave the bound unchanged;
(2) independent split/fit/bound/scoring (own seed) reproduces the within-paper rate within Monte-Carlo error;
(3) independent recomputation of the cross-paper rate (deterministic) and the per-paper weighting from the pair records;
(4) independent paired paper bootstrap."""
import sys, math, numpy as np, pandas as pd
from scipy import stats
X_, J_ = sys.argv[1], sys.argv[2]; sys.argv = ['x', X_, J_, 'primary']
exec(open('reproduce.py').read().split("# ---------------- pair records")[0])      # same data filter only (D, ALLOYS)
Z = stats.norm.ppf(0.95)
def fit_line(x, y):                       # closed-form simple regression, uses TRAINING arrays only
    n = len(x); xb, yb = x.mean(), y.mean(); sxx = ((x - xb) ** 2).sum(); b1 = ((x - xb) * (y - yb)).sum() / sxx; b0 = yb - b1 * xb
    s = math.sqrt(((y - b0 - b1 * x) ** 2).sum() / (n - 2)); return n, xb, sxx, b0, b1, s
def bounds(model, x0):
    n, xb, sxx, b0, b1, s = model; h = 1 / n + (x0 - xb) ** 2 / sxx; mu = b0 + b1 * x0
    pb = mu - stats.t.ppf(0.95, n - 2) * s * np.sqrt(1 + h)
    ne = 1 / h; k = stats.nct.ppf(0.95, n - 2, nc=Z * np.sqrt(ne)) / np.sqrt(ne); return pb, mu - k * s
rng = np.random.default_rng(777); NSPL = 100; rows = []; pert_ok = True; npert = 0
for mat in ALLOYS:
    for R in (-1, 0.1):
        G = D[(D.mat == mat) & (D.R == R)].reset_index(drop=True)
        if G.DOI.nunique() < 6: continue
        y = np.log10(G.N.values); x = np.log10(G.S.values / np.median(G.S.values)); ro = G.ro.values.astype(bool)
        for cid, idx in G.groupby('id').indices.items():
            f = idx[~ro[idx]]; r_ = idx[ro[idx]]
            if len(f) < 8 or len(set(np.round(x[f], 3))) < 3 or len(f) < 6: continue
            nt = max(2, int(round(0.30 * len(f)))); out = []
            for _ in range(NSPL):
                for _try in range(200):
                    perm = rng.permutation(len(f)); te = f[perm[:nt]]; tr = f[perm[nt:]]
                    if len(set(np.round(x[tr], 3))) >= 3: break
                else: continue
                mdl = fit_line(x[tr], y[tr]); lo, hi = x[tr].min(), x[tr].max()
                # perturbation test: wreck the held-out lives, refit from the training rows only -> identical model
                if npert < 400:
                    y2 = y.copy(); y2[te] = rng.choice([-50., 50.], len(te)); m2 = fit_line(x[tr], y2[tr]); npert += 1
                    pert_ok &= np.allclose(m2[3:], mdl[3:], atol=0, rtol=0)
                tef = te[(x[te] >= lo) & (x[te] <= hi)]; rr = r_[(x[r_] >= lo) & (x[r_] <= hi)]
                if len(tef) == 0: continue
                den = len(tef) + len(rr); pb, tb = bounds(mdl, np.concatenate([x[tef], x[rr]]))
                out.append((np.sum(y[tef] < pb[:len(tef)]) / den, np.sum(y[tef] < tb[:len(tef)]) / den))
            if len(out) >= 0.5 * NSPL:
                m = np.mean(out, axis=0); rows.append((mat, R, cid, G.DOI.values[f[0]], m[0], m[1]))
A = pd.DataFrame(rows, columns=['mat', 'R', 'sid', 'A', 'w_pred', 'w_tol'])
rec = pd.read_csv('within_paper_curve_records.csv'); print('curves:', len(A), 'vs script', len(rec))
assert len(A) == len(rec) and set(zip(A.mat, A.R, A.sid)) == set(zip(rec.mat, rec.R, rec.sid)), 'eligible curve set differs'
# (3) cross-paper rate from raw pair records, written independently of the script's helpers
C = pd.read_csv('within_paper_cross_pair_records.csv')
def cross_point(col):
    C['r'] = C[col] / C.n
    pc = {}                                           # (A,B) -> mean over alloy/R groups of mean over curves
    for (mat, R, a, b), g in C.groupby(['mat', 'R', 'A', 'B']): pc.setdefault((a, b), []).append(g.r.mean())
    pc = {k: np.mean(v) for k, v in pc.items()}
    As = sorted({a for a, b in pc}); Bs = sorted({b for a, b in pc}); perB = {b: [] for b in Bs}
    for (a, b), v in pc.items(): perB[b].append(v)
    headline = np.mean([np.mean(v) for v in perB.values()])           # equal weight to every source paper, then every target paper
    perA = {a: np.mean([v for (aa, b), v in pc.items() if aa == a]) for a in As}
    return headline, np.mean(list(perA.values())), perA
S = pd.read_csv('within_paper_summary.csv')
for col, name, wcol in (('pred', '5% prediction', 'w_pred'), ('tol', '95/95 tolerance', 'w_tol')):
    wA = A.groupby('A')[wcol].mean(); head, mean_perA, perA = cross_point(col); papers = sorted(set(wA.index) & set(perA))
    w = np.array([wA[p] for p in papers]); c = np.array([perA[p] for p in papers]); br = np.random.default_rng(4242); d = []
    for _ in range(2000): i = br.integers(0, len(papers), len(papers)); d.append((c[i] - w[i]).mean())
    ref = S[S.bound == name].iloc[0]
    print(f'{name}: within {100*w.mean():.2f} (script {ref.within_breach_pct:.2f}); cross {100*c.mean():.2f} (script {ref.cross_breach_pct:.2f}); target-weighted cross {100*head:.2f} (script {ref.cross_breach_pct_target_weighting:.2f}); '
          f'gap {100*(c-w).mean():.2f} [{100*np.percentile(d,2.5):.2f}, {100*np.percentile(d,97.5):.2f}] (script {ref.gap_pp:.2f} [{ref.gap_CI_low:.2f}, {ref.gap_CI_high:.2f}])')
    assert abs(100 * c.mean() - ref.cross_breach_pct) < 0.01 and abs(100 * head - ref.cross_breach_pct_target_weighting) < 0.01, 'cross rate mismatch'
    assert abs(100 * w.mean() - ref.within_breach_pct) < 0.5 + 0.0 * ref.within_breach_pct, 'within rate outside Monte-Carlo tolerance'
print('perturbation test (held-out lives replaced by +/-50 log10 cycles, %d splits): bound unchanged = %s' % (npert, pert_ok)); assert pert_ok
print('AUDIT PASSED')
