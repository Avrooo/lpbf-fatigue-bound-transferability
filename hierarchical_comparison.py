"""Experiment 2: hierarchical (random-intercept) censored model versus the pooled censored model, leave-one-paper-out.
Within each alloy and stress-ratio group and for every target paper B, both models are fitted to the other papers only, with the same covariate columns.
  pooled censored:        bound = X b + z05 * sigma_pooled                       (as in reproduce.py; EM)
  hierarchical censored:  y_ij = X b + u_i + e_ij, u_i ~ N(0,tau^2), e_ij ~ N(0,s^2), right-censored runouts, random effect integrated out by Gauss-Hermite quadrature.
                          For a NEW paper u is not estimated from B: the predictive distribution is N(X b, s^2 + tau^2), bound = X b + z05 * sqrt(s^2 + tau^2).
Both bounds are scored on identical target specimens (inside each eligible source curve's stress range) with the same hierarchical weighting and the
same two-way paper bootstrap as the headline analysis (reproduce.py).
Usage: python3 hierarchical_comparison.py XLSX JSON [seed=11]   (writes hierarchical_comparison.csv, hierarchical_by_alloy.csv, hierarchical_fit_diagnostics.csv)"""
import sys, math, time, numpy as np, pandas as pd
from scipy.optimize import minimize
X_, J_ = sys.argv[1], sys.argv[2]; SEED = int(sys.argv[3]) if len(sys.argv) > 3 else 11
sys.argv = ['x', X_, J_, 'primary']; exec(open('reproduce.py').read().split("# ---------------- pair records")[0])
GHx, GHw = np.polynomial.hermite_e.hermegauss(21); GHw = GHw / GHw.sum()
def lpf(y, mu, s, c):
    z = (y - mu) / s; return np.where(c, norm.logsf(z), norm.logpdf(z) - np.log(s))
def nll(p, y, X, c, u, K):
    k = X.shape[1]; b = p[:k]; s = math.exp(p[k]); t = math.exp(p[k + 1]); l = lpf(y[:, None], (X @ b)[:, None] + t * GHx[None, :], s, c[:, None]); L = np.zeros((K, len(GHx))); np.add.at(L, u, l); mx = L.max(1, keepdims=True)
    return -(mx[:, 0] + np.log((np.exp(L - mx) * GHw).sum(1))).sum()
def fit_re(y, X, c, g):
    u = np.unique(g, return_inverse=True)[1]; K = u.max() + 1; b0, s0 = censored_lognormal_em(X, y, c); k = X.shape[1]
    best = None
    for t0 in (0.5, 0.2):                       # two starting values for tau; keep the better optimum
        r = minimize(nll, np.r_[b0, math.log(s0 * 0.9), math.log(t0)], args=(y, X, c, u, K), method='L-BFGS-B')
        if best is None or r.fun < best.fun: best = r
    r = best; ok = bool(r.success) and np.isfinite(r.fun) and abs(r.x[k]) < 5 and abs(r.x[k + 1]) < 5
    return r.x[:k], math.exp(r.x[k]), math.exp(r.x[k + 1]), ok, int(r.nit), float(r.fun), str(r.message)
METHODS = ['pooled_censored', 'hierarchical_censored']
records = []; diag = []; t0 = time.time()
for mat in ALLOYS:
    for R in (-1, 0.1):
        G = D[(D.mat == mat) & (D.R == R)].reset_index(drop=True)
        if G.DOI.nunique() < 6: continue
        y = np.log10(G.N.values); x = np.log10(G.S.values / np.median(G.S.values)); C = G[['asb', 'hip', 'hthi']].values.astype(float); C = C[:, C.std(0) > 0]
        X = np.column_stack([np.ones(len(G)), x, C]); cens = G.ro.values.astype(bool); doi = G.DOI.values; cid = G.id.values; curves = {}
        for sid in np.unique(cid):
            m = (cid == sid) & ~cens
            if m.sum() < 6 or len(np.unique(np.round(x[m], 3))) < 3: continue
            curves[sid] = (doi[cid == sid][0], x[m].min(), x[m].max())
        bF, sF, tF, okF, nitF, fF, mF = fit_re(y, X, cens, doi); diag.append(dict(mat=mat, R=R, target='ALL PAPERS (no hold-out)', n_train=len(G), n_train_papers=len(np.unique(doi)), sigma_within=sF, tau_between=tF, converged=okF, n_iter=nitF, message=mF))
        for B in np.unique(doi):
            tr = doi != B; te = ~tr
            bc, sc = censored_lognormal_em(X[tr], y[tr], cens[tr]); bh, sh, th, ok, nit, fval, msg = fit_re(y[tr], X[tr], cens[tr], doi[tr])
            diag.append(dict(mat=mat, R=R, target=B, n_train=int(tr.sum()), n_train_papers=len(np.unique(doi[tr])), sigma_within=sh, tau_between=th, sigma_pooled_censored=sc, sd_total_hier=math.sqrt(sh ** 2 + th ** 2), converged=ok, n_iter=nit, message=msg))
            for sid, (A, lo, hi) in curves.items():
                if A == B: continue
                sel = te & (x >= lo) & (x <= hi)
                if not sel.any(): continue
                bd = {'pooled_censored': X[sel] @ bc + Z05 * sc, 'hierarchical_censored': X[sel] @ bh + Z05 * math.sqrt(sh ** 2 + th ** 2)}; f = ~cens[sel]
                records.append((mat, R, A, sid, B, int(sel.sum()), {m_: int(np.sum(y[sel][f] < bd[m_][f])) for m_ in METHODS}, {m_: float(bd[m_].mean()) for m_ in METHODS}))
    print(mat, 'done', round(time.time() - t0), 's', flush=True)
Dg = pd.DataFrame(diag); Dg.to_csv('hierarchical_fit_diagnostics.csv', index=False, float_format='%.4f')
Dh = Dg[Dg.target != 'ALL PAPERS (no hold-out)']
print(f'hierarchical fits (leave-one-paper-out): {len(Dh)}; converged: {Dh.converged.sum()} ({100 * Dh.converged.mean():.1f}%); tau near zero (<0.02): {(Dh.tau_between < 0.02).sum()}')
# ---------------- weighting and bootstrap (same code as reproduce.py; mean bounds use the same weights) ----------------
def pair_vals(recs, kind, method, mat=None):
    g = {}
    for mat_, R, A, sid, B, n, br, mb in recs:
        if mat and mat_ != mat: continue
        g.setdefault((mat_, R, A, B), []).append(br[method] / n if kind == 'rate' else mb[method])
    p = {}
    for (mat_, R, A, B), v in g.items(): p.setdefault((A, B), []).append(sum(v) / len(v))
    return {k: sum(v) / len(v) for k, v in p.items()}
def weighted(pr, wA, wB):
    perB = {}
    for (A, B), r in pr.items():
        w = wA.get(A, 0)
        if w: s = perB.setdefault(B, [0.0, 0.0]); s[0] += w * r; s[1] += w
    num = den = 0.0
    for B, (s, w) in perB.items():
        wb = wB.get(B, 0)
        if wb and w: num += wb * s / w; den += wb
    return num / den if den else float('nan')
def boot(pr, rs, nb=2000):
    As = sorted({a for a, b in pr}); Bs = sorted({b for a, b in pr}); est = weighted(pr, {a: 1 for a in As}, {b: 1 for b in Bs}); bs = []
    for _ in range(nb):
        ca = np.bincount(rs.integers(0, len(As), len(As)), minlength=len(As)); cb = np.bincount(rs.integers(0, len(Bs), len(Bs)), minlength=len(Bs)); bs.append(weighted(pr, dict(zip(As, ca)), dict(zip(Bs, cb))))
    return est, np.percentile(bs, [2.5, 97.5]), len(pr), len(As), len(Bs)
def compare(mat, seed):
    out = []
    pr = {m_: pair_vals(records, 'rate', m_, mat) for m_ in METHODS}; pb = {m_: pair_vals(records, 'bound', m_, mat) for m_ in METHODS}
    for m_ in METHODS:
        e, ci, npair, na, nb_ = boot(pr[m_], np.random.default_rng(seed)); eb, cib, *_ = boot(pb[m_], np.random.default_rng(seed))
        out.append(dict(model=m_, breach_pct=100 * e, CI_low=100 * ci[0], CI_high=100 * ci[1], mean_lower_bound_log10N=eb, bound_CI_low=cib[0], bound_CI_high=cib[1], pairs=npair, source_papers=na, target_papers=nb_))
    dr = {k: pr['hierarchical_censored'][k] - pr['pooled_censored'][k] for k in pr['pooled_censored']}; db = {k: pb['hierarchical_censored'][k] - pb['pooled_censored'][k] for k in pb['pooled_censored']}
    e, ci, *_ = boot(dr, np.random.default_rng(seed + 1)); eb, cib, *_ = boot(db, np.random.default_rng(seed + 1))
    out.append(dict(model='paired difference (hierarchical minus pooled)', breach_pct=100 * e, CI_low=100 * ci[0], CI_high=100 * ci[1], mean_lower_bound_log10N=eb, bound_CI_low=cib[0], bound_CI_high=cib[1]))
    return out
res = compare(None, SEED)
sig = Dh.groupby('mat')[['sigma_within', 'tau_between']].mean()
Sall = pd.DataFrame(res); Sall['sigma_within_mean_over_fits'] = [Dh.sigma_within.mean()] * 2 + [np.nan] * 0 + [np.nan] if False else np.nan
Sall.loc[Sall.model == 'hierarchical_censored', 'sigma_within_mean_over_fits'] = Dh.sigma_within.mean(); Sall.loc[Sall.model == 'hierarchical_censored', 'tau_between_mean_over_fits'] = Dh.tau_between.mean()
Sall.loc[Sall.model == 'pooled_censored', 'sigma_pooled_mean_over_fits'] = Dh.sigma_pooled_censored.mean()
Sall['converged_pct'] = 100 * Dh.converged.mean(); Sall.to_csv('hierarchical_comparison.csv', index=False, float_format='%.4f')
ra = []
for mat in ALLOYS:
    for r_ in compare(mat, SEED + 7): r_['alloy'] = mat; ra.append(r_)
A_ = pd.DataFrame(ra); sg = Dh.groupby('mat').agg(sigma_within_mean=('sigma_within', 'mean'), tau_between_mean=('tau_between', 'mean'), sigma_pooled_mean=('sigma_pooled_censored', 'mean'), converged_pct=('converged', lambda s: 100 * s.mean()), n_fits=('converged', 'size'))
A_ = A_.merge(sg, left_on='alloy', right_index=True, how='left'); A_.to_csv('hierarchical_by_alloy.csv', index=False, float_format='%.4f')
pd.set_option('display.width', 250); pd.set_option('display.max_columns', 30)
print(Sall.round(3).to_string()); print(A_.round(3).to_string())
