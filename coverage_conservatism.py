"""Coverage versus conservatism, strictly external validation, and weighting/bootstrap audits.
For every target paper B (leave-one-paper-out) and every eligible source curve A (from another paper), single-study, pooled censored and hierarchical censored lower bounds are
computed on IDENTICAL target specimens (inside the source curve's stress range) at a grid of lower-tail probabilities p. For each method and p we record the confirmed-breach
rate (failures below the bound; runouts never breaches, kept in the denominator) and the mean bound height (log10 cycles; lower = more conservative).
STRICT=1: the covariate columns of the pooled and hierarchical models are chosen from the TRAINING papers only (target-paper data never enter model specification);
STRICT=0 reproduces the earlier code (columns chosen on the whole alloy/R group) and must give the published headline numbers.
Usage: python3 coverage_conservatism.py XLSX JSON STRICT [seed=11] [which=primary]   (outputs: cc_<tag>_*.csv)"""
import sys, math, time, numpy as np, pandas as pd
from scipy.optimize import minimize
X_, J_ = sys.argv[1], sys.argv[2]; STRICT = int(sys.argv[3]); SEED = int(sys.argv[4]) if len(sys.argv) > 4 else 11; WH = sys.argv[5] if len(sys.argv) > 5 else 'primary'
sys.argv = ['x', X_, J_, WH]; exec(open('reproduce.py').read().split("# ---------------- pair records")[0])
TAG = f'{WH}_strict{STRICT}'
GHx, GHw = np.polynomial.hermite_e.hermegauss(21); GHw = GHw / GHw.sum()
def lpf(y, mu, s, c):
    z = (y - mu) / s; return np.where(c, norm.logsf(z), norm.logpdf(z) - np.log(s))
def nll(p, y, X, c, u, K):
    k = X.shape[1]; b = p[:k]; s = math.exp(p[k]); t = math.exp(p[k + 1]); l = lpf(y[:, None], (X @ b)[:, None] + t * GHx[None, :], s, c[:, None]); L = np.zeros((K, len(GHx))); np.add.at(L, u, l); mx = L.max(1, keepdims=True)
    return -(mx[:, 0] + np.log((np.exp(L - mx) * GHw).sum(1))).sum()
def fit_re(y, X, c, g):
    u = np.unique(g, return_inverse=True)[1]; K = u.max() + 1; b0, s0 = censored_lognormal_em(X, y, c); k = X.shape[1]; best = None
    for t0 in (0.5, 0.2):
        r = minimize(nll, np.r_[b0, math.log(s0 * 0.9), math.log(t0)], args=(y, X, c, u, K), method='L-BFGS-B')
        if best is None or r.fun < best.fun: best = r
    r = best; return r.x[:k], math.exp(r.x[k]), math.exp(r.x[k + 1])
PG = [0.00001, 0.0001, 0.0003, 0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.10, 0.20, 0.30, 0.50]; METH = ['single', 'pooled', 'hier']; P05 = PG.index(0.05)
records = []; t0 = time.time(); ncol_dropped = 0
for mat in ALLOYS:
    for R in (-1, 0.1):
        G = D[(D.mat == mat) & (D.R == R)].reset_index(drop=True)
        if G.DOI.nunique() < 6: continue
        y = np.log10(G.N.values); x = np.log10(G.S.values / np.median(G.S.values)); Call = G[['asb', 'hip', 'hthi']].values.astype(float)
        cens = G.ro.values.astype(bool); doi = G.DOI.values; cid = G.id.values; curves = {}
        for sid in np.unique(cid):
            m = (cid == sid) & ~cens
            if m.sum() < 6 or len(np.unique(np.round(x[m], 3))) < 3: continue
            Xs = np.column_stack([np.ones(m.sum()), x[m]]); bs, sd = ols_fail(Xs, y[m]); curves[sid] = (doi[cid == sid][0], len(y[m]), bs, sd, np.linalg.inv(Xs.T @ Xs), x[m].min(), x[m].max())
        keep_all = Call.std(0) > 0
        for B in np.unique(doi):
            tr = doi != B; te = ~tr
            keep = (Call[tr].std(0) > 0) if STRICT else keep_all
            ncol_dropped += int((keep != keep_all).sum()); C = Call[:, keep]; X = np.column_stack([np.ones(len(G)), x, C])
            bc, sc = censored_lognormal_em(X[tr], y[tr], cens[tr]); bh, sh, th = fit_re(y[tr], X[tr], cens[tr], doi[tr]); sdh = math.sqrt(sh ** 2 + th ** 2)
            for sid, (A, n, bs, sd, XtXi, lo, hi) in curves.items():
                if A == B: continue
                sel = te & (x >= lo) & (x <= hi)
                if not sel.any(): continue
                x0 = np.column_stack([np.ones(sel.sum()), x[sel]]); mu = x0 @ bs; q = np.einsum('ij,jk,ik->i', x0, XtXi, x0); f = ~cens[sel]; ys = y[sel]; mp = X[sel] @ bc; mh = X[sel] @ bh
                rec = dict(mat=mat, R=R, A=A, sid=sid, B=B, n=int(sel.sum()), slope=bc[1])
                for k, p in enumerate(PG):
                    bnd = {'single': mu - tdist.ppf(1 - p, n - 2) * sd * np.sqrt(1 + q), 'pooled': mp + norm.ppf(p) * sc, 'hier': mh + norm.ppf(p) * sdh}
                    for m_ in METH:
                        rec[f'br_{m_}_{k}'] = int(np.sum(ys[f] < bnd[m_][f])); rec[f'h_{m_}_{k}'] = float(bnd[m_].mean())
                        if k == P05: rec[f'amb_{m_}'] = int(np.sum(ys[~f] < bnd[m_][~f]))      # runouts that stopped below the bound (ambiguous)
                tol = mu - np.array([tol_factor(1 / v, n - 2) for v in q]) * sd; rec['br_tol'] = int(np.sum(ys[f] < tol[f])); rec['h_tol'] = float(tol.mean()); rec['nf'] = int(f.sum())
                records.append(rec)
    print(mat, 'done', round(time.time() - t0), 's', flush=True)
Rc = pd.DataFrame(records); print('pair-curve records:', len(Rc), '| covariate columns dropped by the training-only rule:', ncol_dropped)
cols = [f'{kind}_{m_}_{k}' for m_ in METH for k in range(len(PG)) for kind in ('r', 'h')] + ['r_tol', 'h_tol']
for m_ in METH:
    for k in range(len(PG)): Rc[f'r_{m_}_{k}'] = Rc[f'br_{m_}_{k}'] / Rc.n
Rc['r_tol'] = Rc.br_tol / Rc.n
Rc.to_csv(f'cc_{TAG}_records.csv', index=False, float_format='%.6f')
# ---------------- pair values: mean over curves, then over groups, per (A,B) ----------------
def pair_matrix(R_):
    g = R_.groupby(['mat', 'R', 'A', 'B'])[cols].mean().groupby(['A', 'B']).mean(); return g
PV = pair_matrix(Rc); Aix = {a: i for i, a in enumerate(sorted({a for a, b in PV.index}))}; Bix = {b: i for i, b in enumerate(sorted({b for a, b in PV.index}))}
Ai = np.array([Aix[a] for a, b in PV.index]); Bi = np.array([Bix[b] for a, b in PV.index]); P = PV.values; nA, nB = len(Aix), len(Bix)
def wavg(wA, wB, P=P, Ai=Ai, Bi=Bi, nB=nB):
    wa = wA[Ai]; den = np.bincount(Bi, wa, nB); ok = (den > 0) & (wB > 0); out = np.zeros(P.shape[1])
    for j in range(P.shape[1]): out[j] = (wB[ok] * np.bincount(Bi, wa * P[:, j], nB)[ok] / den[ok]).sum() / wB[ok].sum()
    return out
est = wavg(np.ones(nA), np.ones(nB)); rs = np.random.default_rng(SEED); NB = 2000; boots = np.empty((NB, len(cols)))
for b in range(NB): boots[b] = wavg(np.bincount(rs.integers(0, nA, nA), minlength=nA).astype(float), np.bincount(rs.integers(0, nB, nB), minlength=nB).astype(float))
ci = np.percentile(boots, [2.5, 97.5], axis=0); ix = {c: i for i, c in enumerate(cols)}
fr = []
for m_ in METH:
    for k, p in enumerate(PG):
        a, b = ix[f'r_{m_}_{k}'], ix[f'h_{m_}_{k}']
        fr.append(dict(method=m_, p=p, breach_pct=100 * est[a], breach_lo=100 * ci[0, a], breach_hi=100 * ci[1, a], mean_bound_log10N=est[b], bound_lo=ci[0, b], bound_hi=ci[1, b]))
a, b = ix['r_tol'], ix['h_tol']; fr.append(dict(method='single_95/95', p=np.nan, breach_pct=100 * est[a], breach_lo=100 * ci[0, a], breach_hi=100 * ci[1, a], mean_bound_log10N=est[b], bound_lo=ci[0, b], bound_hi=ci[1, b]))
FR = pd.DataFrame(fr); FR.to_csv(f'cc_{TAG}_frontier.csv', index=False, float_format='%.4f')
# ---------------- coverage-conservatism contrasts inside every bootstrap resample ----------------
def curve(v, m_):    # arrays over the p grid: breach rate and bound height
    return np.array([v[ix[f'r_{m_}_{k}']] for k in range(len(PG))]), np.array([v[ix[f'h_{m_}_{k}']] for k in range(len(PG))])
def h_at_rate(r, h, target):      # bound height at which the external breach rate equals target (monotone in p); NaN outside the sampled range
    o = np.argsort(r); r, h = r[o], h[o]
    return float(np.interp(target, r, h)) if r[0] <= target <= r[-1] else np.nan
def r_at_height(r, h, target):
    o = np.argsort(h); r, h = r[o], h[o]
    return float(np.interp(target, h, r)) if h[0] <= target <= h[-1] else np.nan
def contrasts(v):
    out = {}; cv = {m_: curve(v, m_) for m_ in METH}
    for m_ in METH: out[f'height_at_5pct_breach_{m_}'] = h_at_rate(*cv[m_], 0.05)
    for m_ in ('single', 'hier'): out[f'height_gain_at_5pct_pooled_minus_{m_}'] = out['height_at_5pct_breach_pooled'] - out[f'height_at_5pct_breach_{m_}']
    Href = cv['pooled'][1][P05]                                                 # height of the pooled 5% bound
    for m_ in ('single', 'hier'): out[f'breach_at_pooled5_height_{m_}'] = 100 * r_at_height(*cv[m_], Href); out[f'breach_at_pooled5_height_pooled'] = 100 * cv['pooled'][0][P05]
    return out
c0 = contrasts(est); cb = pd.DataFrame([contrasts(boots[b]) for b in range(NB)])
CT = pd.DataFrame([dict(quantity=k, estimate=c0[k], ci_low=np.nanpercentile(cb[k], 2.5), ci_high=np.nanpercentile(cb[k], 97.5), n_boot_defined=int(cb[k].notna().sum())) for k in c0])
slope = Rc.slope.mean(); CT.loc[len(CT)] = dict(quantity='mean pooled Basquin slope d(log10 N)/d(log10 S)', estimate=slope, ci_low=np.nan, ci_high=np.nan, n_boot_defined=len(Rc))
CT.to_csv(f'cc_{TAG}_contrasts.csv', index=False, float_format='%.4f')
# ---------------- audit: alternative weightings, bootstrap variants, ambiguous runouts (p = 5% bounds) ----------------
aud = []; heads = [('single_5%', 'r_single_%d' % P05), ('single_95/95', 'r_tol'), ('pooled_5%', 'r_pooled_%d' % P05), ('hier_5%', 'r_hier_%d' % P05)]
Rc['pair'] = list(zip(Rc.A, Rc.B))
for lab, c in heads:
    base = 100 * est[ix[c]]
    spec = 100 * Rc['b' + c].sum() / Rc.n.sum()                                       # specimen-weighted (all scored specimens in the denominator)
    pairs_eq = 100 * PV[c].mean()                                                                      # every (A,B) pair equal
    srcw = 100 * PV[c].groupby('A').mean().mean(); tgtw = 100 * PV[c].groupby('B').mean().mean()    # source-paper-first, target-paper-first (headline)
    aud.append(dict(metric=lab, headline_target_then_source=base, source_paper_weighted=srcw, target_paper_weighted=tgtw, pair_weighted=pairs_eq, specimen_weighted=spec))
AUD = pd.DataFrame(aud)
# ambiguous runouts counted as breaches
amb = []
for m_, lab in (('single', 'single_5%'), ('pooled', 'pooled_5%'), ('hier', 'hier_5%')):
    tot = 100 * (Rc[f'br_{m_}_{P05}'].sum() + Rc[f'amb_{m_}'].sum()) / Rc.n.sum(); base = 100 * Rc[f'br_{m_}_{P05}'].sum() / Rc.n.sum()
    amb.append(dict(metric=lab, specimen_weighted_confirmed=base, specimen_weighted_if_ambiguous_runouts_counted=tot, ambiguous_runouts=int(Rc[f'amb_{m_}'].sum()), scored_specimens=int(Rc.n.sum())))
AMB = pd.DataFrame(amb)
# bootstrap variants for the headline metrics: independent roles (primary), same paper set for both roles, source papers only, target papers only
dois = sorted(set(Aix) | set(Bix)); di = {d: i for i, d in enumerate(dois)}; Ad = np.array([di[a] for a in Aix]); Bd = np.array([di[b] for b in Bix]); rs2 = np.random.default_rng(SEED + 5)
hc = [ix[c] for _, c in heads]; Ph = P[:, hc]; var = {k: [] for k in ('independent', 'same_paper_set', 'source_only', 'target_only')}
for b in range(NB):
    wi_ = np.bincount(rs2.integers(0, len(dois), len(dois)), minlength=len(dois)).astype(float)
    var['independent'].append(wavg(np.bincount(rs2.integers(0, nA, nA), minlength=nA).astype(float), np.bincount(rs2.integers(0, nB, nB), minlength=nB).astype(float), Ph))
    var['same_paper_set'].append(wavg(wi_[Ad], wi_[Bd], Ph)); var['source_only'].append(wavg(np.bincount(rs2.integers(0, nA, nA), minlength=nA).astype(float), np.ones(nB), Ph)); var['target_only'].append(wavg(np.ones(nA), np.bincount(rs2.integers(0, nB, nB), minlength=nB).astype(float), Ph))
bv = []
for k, v in var.items():
    v = np.array(v); lo, hi = np.percentile(v, [2.5, 97.5], axis=0)
    for j, (lab, _) in enumerate(heads): bv.append(dict(bootstrap=k, metric=lab, ci_low=100 * lo[j], ci_high=100 * hi[j], width=100 * (hi[j] - lo[j])))
BV = pd.DataFrame(bv)
AUD.to_csv(f'cc_{TAG}_audit_weighting.csv', index=False, float_format='%.3f'); AMB.to_csv(f'cc_{TAG}_audit_ambiguous_runouts.csv', index=False, float_format='%.3f'); BV.to_csv(f'cc_{TAG}_audit_bootstrap.csv', index=False, float_format='%.3f')
pd.set_option('display.width', 250); pd.set_option('display.max_columns', 30)
print(f'pairs={len(PV)} sources={nA} targets={nB}'); print(FR.round(3).to_string()); print(CT.round(3).to_string()); print(AUD.round(2).to_string()); print(AMB.round(3).to_string()); print(BV.round(2).to_string())
