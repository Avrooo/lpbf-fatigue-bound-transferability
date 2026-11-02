"""Explicit test of the two-way cluster-bootstrap interval, by simulation with a known data-generating process (DGP).
Design: the DGP is calibrated to one real alloy-R group (curve slopes, curve-level life/stress-location relation, between-paper and
between-curve scatter, within-curve scatter) and uses the REAL specimen stress designs of that group. For each replication K papers are
drawn (with replacement) from the real designs, new lives are simulated, the SAME estimator and weighting as in reproduce.py are applied
(single-study 5% prediction bound; pooled leave-one-paper-out failures-only 5% prediction bound), and the 95% percentile bootstrap interval
(resampling source papers and target papers independently) is compared with the truth.
Truth = expected value of the point estimate at the same K, from many independent replications (finite-K estimand, so it contains the
pooled-fit training-size effect); a second reference, the K=300 limit, shows the finite-K bias of the point estimate separately.
Simplifications: one alloy-R group, no covariates, no runouts, Gaussian errors. Usage:
  python3 -I coverage_test.py XLSX JSON GROUP K REPS NB SEED   e.g.  ... "Ti-6Al-4V|-1" 60 400 1000 1"""
import sys, numpy as np, pandas as pd
from scipy.stats import t as tdist, norm
X_, J_, GROUP, K, REPS, NB, SEED = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6]), int(sys.argv[7])
sys.argv = ['x', X_, J_, 'primary']; exec(open('reproduce.py').read().split("# ---------------- independent model code")[0])
mat, R = GROUP.split('|'); R = float(R)
G = D[(D.mat == mat) & (D.R == R)].copy(); G['x'] = np.log10(G.S / G.S.median()); G['y'] = np.log10(G.N); F = G[G.ro == 0]
Z05 = norm.ppf(0.05)
# ---------------- calibrate DGP from the real failures ----------------
cur = []
for sid, g in F.groupby('id'):
    if len(g) >= 6 and g.x.round(3).nunique() >= 3:
        A = np.c_[np.ones(len(g)), g.x - g.x.mean()]; c = np.linalg.lstsq(A, g.y, rcond=None)[0]; s2 = ((g.y - A @ c) ** 2).sum() / (len(g) - 2)
        cur.append((g.DOI.iloc[0], g.x.mean(), c[0], c[1], s2, len(g), (s2 * np.linalg.inv(A.T @ A)[1, 1]) ** .5))
C = pd.DataFrame(cur, columns=['doi', 'm', 'a', 'b', 's2', 'n', 'seb'])
sigma = float(np.sqrt(np.average(C.s2, weights=C.n - 2))); bbar = float(C.b.median()); lowse = C[C.seb <= C.seb.median()]; tslope = float((lowse.b.quantile(.75) - lowse.b.quantile(.25)) / 1.349)
gam, mu0 = np.polyfit(C.m, C.a, 1); res = C.a - (mu0 + gam * C.m); v_tot = max(res.var() - (sigma ** 2 / C.n).mean(), 0)
w = [g.a.sub(mu0 + gam * g.m).var() - sigma ** 2 / g.n.mean() for _, g in C.groupby('doi') if len(g) >= 2]; tc = float(np.sqrt(min(max(np.mean(w), 0), v_tot))); tp = float(np.sqrt(v_tot - tc ** 2))
PAR = dict(mu0=mu0, gam=gam, bbar=bbar, tslope=tslope, tp=tp, tc=tc, sigma=sigma)
designs = {d: [(sid, g.x.values) for sid, g in G[G.DOI == d].groupby('id')] for d in G.DOI.unique()}; dois = list(designs)
# ---------------- one synthetic data set and the estimator ----------------
def make(rng, K):
    pap, cu, xs, ys = [], [], [], []; cid = 0
    for k, d in enumerate(rng.choice(len(dois), K)):
        u = rng.normal(0, PAR['tp'])
        for sid, x in designs[dois[d]]:
            m = x.mean(); v = rng.normal(0, PAR['tc']); b = PAR['bbar'] + rng.normal(0, PAR['tslope'])
            ys.append(PAR['mu0'] + PAR['gam'] * m + u + v + b * (x - m) + rng.normal(0, PAR['sigma'], len(x))); xs.append(x); pap += [k] * len(x); cu += [cid] * len(x); cid += 1
    return np.array(pap), np.array(cu), np.concatenate(xs), np.concatenate(ys)
def rate_matrices(pap, cu, x, y, K):
    n = len(y); X = np.c_[np.ones(n), x]
    XtX = np.zeros((K, 2, 2)); Xty = np.zeros((K, 2)); yty = np.zeros(K); cnt = np.bincount(pap, minlength=K)
    for k in range(K):
        m = pap == k; XtX[k] = X[m].T @ X[m]; Xty[k] = X[m].T @ y[m]; yty[k] = y[m] @ y[m]
    T = XtX.sum(0); Tv = Xty.sum(0); Ty = yty.sum(); bound_p = np.empty(n)
    for k in range(K):
        A = T - XtX[k]; bv = np.linalg.solve(A, Tv - Xty[k]); sse = (Ty - yty[k]) - 2 * bv @ (Tv - Xty[k]) + bv @ A @ bv; s = np.sqrt(sse / (n - cnt[k] - 2)); m = pap == k
        bound_p[m] = X[m] @ bv + Z05 * s
    br_p = y < bound_p; curves = np.unique(cu); cp = np.array([pap[cu == c][0] for c in curves])
    Rs = []; Rp = []; src = []
    for c, a in zip(curves, cp):
        m = cu == c
        if m.sum() < 6 or len(np.unique(np.round(x[m], 3))) < 3: continue
        xa, ya = x[m], y[m]; nn = m.sum(); A = np.c_[np.ones(nn), xa]; b = np.linalg.lstsq(A, ya, rcond=None)[0]; sd = np.sqrt(((ya - A @ b) ** 2).sum() / (nn - 2)); xb = xa.mean(); Sxx = ((xa - xb) ** 2).sum()
        h = 1 / nn + (x - xb) ** 2 / Sxx; bd = b[0] + b[1] * x - tdist.ppf(0.95, nn - 2) * sd * np.sqrt(1 + h)
        sel = (pap != a) & (x >= xa.min()) & (x <= xa.max()); nB = np.bincount(pap[sel], minlength=K).astype(float)
        with np.errstate(invalid='ignore', divide='ignore'):
            Rs.append(np.bincount(pap[sel & (y < bd)], minlength=K) / nB); Rp.append(np.bincount(pap[sel & br_p], minlength=K) / nB)
        src.append(a)
    src = np.array(src); I = np.zeros((len(src), K)); I[np.arange(len(src)), src] = 1
    out = []
    for Rc in (np.array(Rs), np.array(Rp)):          # curves x targets -> targets x source papers (mean over curves of a source paper)
        ok = ~np.isnan(Rc); Rz = np.where(ok, Rc, 0.0); num = Rz.T @ I; den = ok.T.astype(float) @ I; out.append((np.where(den > 0, num / np.maximum(den, 1), 0.0), (den > 0).astype(float)))
    return out
def estimate(Rm, M, wa, wb):
    num = wa @ Rm.T; den = wa @ M.T; ok = den > 0; rate = np.where(ok, num / np.maximum(den, 1e-12), 0.0); return (wb * rate * ok).sum(1) / (wb * ok).sum(1)
def one(rng, K, nb):
    pap, cu, x, y = make(rng, K); res = []
    for Rm, M in rate_matrices(pap, cu, x, y, K):
        pt = estimate(Rm, M, np.ones((1, K)), np.ones((1, K)))[0]
        if nb:
            wa = rng.multinomial(K, np.ones(K) / K, nb).astype(float); wb = rng.multinomial(K, np.ones(K) / K, nb).astype(float); bs = estimate(Rm, M, wa, wb)
            bj = estimate(Rm, M, wa, wa)       # variant: the same resampled paper set plays both source and target roles
            res.append((pt, *np.percentile(bs, [2.5, 97.5]), *np.percentile(bj, [2.5, 97.5])))
        else: res.append((pt,) + (np.nan,) * 4)
    return res
if __name__ == '__main__':
    print('group', GROUP, 'papers', len(dois), 'DGP', {k: round(float(v), 3) for k, v in PAR.items()}, 'K', K)
    rng = np.random.default_rng(SEED)
    truth_pts = np.array([[r[0] for r in one(rng, K, 0)] for _ in range(max(1500, 3 * REPS))])     # truth at same K (point estimates only)
    truth = truth_pts.mean(0); big = np.array([[r[0] for r in one(rng, 300, 0)] for _ in range(40)]).mean(0)
    runs = [one(rng, K, NB) for _ in range(REPS)]
    for j, name in enumerate(['single_pred5', 'pooled_failures_only_pred5']):
        a = np.array([r[j] for r in runs]); cover = ((a[:, 1] <= truth[j]) & (truth[j] <= a[:, 2]))
        cj = ((a[:, 3] <= truth[j]) & (truth[j] <= a[:, 4]))
        print(f'  [{name}] sd of the point estimate over {len(truth_pts)} independent replications = {100 * truth_pts[:, j].std():.2f} pts; implied 95% interval width = {100 * 3.92 * truth_pts[:, j].std():.2f} pts')
        print(f'  [{name}] VARIANT same-paper-set resampling: coverage={100 * cj.mean():5.1f}%  mean width={100 * (a[:, 4] - a[:, 3]).mean():5.2f} pts')
        print(f'{name:28s} truth(K={K})={100 * truth[j]:5.2f}%  K=300 limit={100 * big[j]:5.2f}%  mean estimate={100 * a[:, 0].mean():5.2f}%  coverage={100 * cover.mean():5.1f}% (MC s.e. {100 * np.sqrt(cover.mean() * (1 - cover.mean()) / REPS):.1f})  '
              f'miss-low={100 * np.mean(a[:, 2] < truth[j]):4.1f}% miss-high={100 * np.mean(a[:, 1] > truth[j]):4.1f}%  mean width={100 * (a[:, 2] - a[:, 1]).mean():5.2f} pts  sd(estimate)={100 * a[:, 0].std():5.2f} pts')
