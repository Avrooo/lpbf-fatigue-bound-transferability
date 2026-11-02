"""Standalone reproduction of the cross-paper coverage analysis (primary and all-data sets).
Independent of the exploratory scripts: loads the raw FatigueData-AM2022 files, uses (i) matrix-form regression bounds,
(ii) an EM algorithm for right-censored lognormal regression (not the L-BFGS used before), (iii) dictionary/loop weighting and bootstrap.
Usage: python3 reproduce.py XLSX JSON [primary|all|primary_noflag|primary_nopaperflag|primary_highrated|primary_lowrated] [seed] [slope_min=0]"""
import sys, json, math, numpy as np, pandas as pd
from scipy.stats import norm, t as tdist, nct
XLSX, JSONF = sys.argv[1], sys.argv[2]; WHICH = sys.argv[3] if len(sys.argv) > 3 else 'primary'; SEED = int(sys.argv[4]) if len(sys.argv) > 4 else 7
SLOPE_MIN = float(sys.argv[5]) if len(sys.argv) > 5 else 0.0   # flag a curve if the fitted slope of log10 N on log10 S is >= SLOPE_MIN (0 = non-negative slope, as in the manuscript)
ALLOYS = ['Ti-6Al-4V', 'AlSi10Mg', 'IN718', '316L']; Z05 = norm.ppf(0.05); Z95 = norm.ppf(0.95)
# ---------------- load ----------------
sn = pd.read_excel(XLSX, 'S-N'); sn.columns = ['id', 'N', 'S', 'ro']
par = pd.read_excel(XLSX, 'parameter', header=None); hdr = [str(h) for h in par.iloc[1]]; par = par.iloc[2:].copy(); par.columns = hdr; par['id'] = par['dataset id'].astype(int)
meta = par[['id', 'DOI', 'name of the material', 'types of AM', 'types of fatigue tests', 'load ratio', 'rating score']].copy(); meta.columns = ['id', 'DOI', 'mat', 'am', 'test', 'R', 'rating']; meta['rating'] = pd.to_numeric(meta.rating, errors='coerce'); meta['R'] = pd.to_numeric(meta.R, errors='coerce')
arts = json.load(open(JSONF))['fatiguedata_am2022']['articles']; proc = {}; i = 0
for a in arts:
    for ds in a['scidata']['datasets']:
        if ds['fatigue']['fdata_type'] != 'sn': continue
        i += 1; proc[i] = ds.get('processing', {}).get('proc_para', []) or []
assert i == sn.id.nunique() == 1610, (i, sn.id.nunique())
def features(pp):
    surf = [(p.get('method', '') or '').lower() for p in pp if p.get('type') == 'SURF']; last = surf[-1] if surf else ''
    asb = int(last.startswith('as-built') or 'channel as-built' in last); hip = int(any(p.get('type') == 'HIP' for p in pp))
    hthi = int(any(max(p.get('temperature', [0]) or [0]) >= 700 for p in pp if p.get('type') == 'HT')); return asb, hip, hthi
meta = meta[meta.id.isin(sn.id)].copy()
meta[['asb', 'hip', 'hthi']] = [features(proc[k]) for k in meta.id]
D = sn.merge(meta, on='id'); D = D[(D.am == 'L-PBF') & (D.test == 'uniaxial') & D.mat.isin(ALLOYS) & D.R.isin([-1, 0.1])].copy()
rng_ = D.groupby('DOI').N.agg(['min', 'max']); out_of_window = set(rng_[(rng_['min'] < 1e3) | (rng_['max'] > 2e7)].index)
MEDR = D.rating.median()          # median documentation-quality rating over all analysed rows (before the life-range restriction)
if WHICH.startswith('primary'): D = D[~D.DOI.isin(out_of_window)]
if WHICH == 'primary_highrated': D = D[D.rating >= MEDR]
if WHICH == 'primary_lowrated': D = D[D.rating < MEDR]
if WHICH in ('primary_noflag', 'primary_nopaperflag'):      # data-quality sensitivity: drop curves with S>UTS or non-negative-slope fits, and duplicate rows (or whole papers containing such curves)
    ucol = [c for c in par.columns if c.startswith('ultimate tensile')][0]; uts = pd.to_numeric(par.set_index('id')[ucol], errors='coerce'); uts = uts[~uts.index.duplicated()]
    flag = set(D[D.S > D.id.map(uts)].id)
    for i_, g_ in D[D.ro == 0].groupby('id'):
        if len(g_) >= 6 and g_.S.nunique() >= 3 and np.polyfit(np.log10(g_.S), np.log10(g_.N), 1)[0] >= SLOPE_MIN: flag.add(i_)
    dup_ = D.duplicated(['id', 'N', 'S'])
    D = D[~D.id.isin(flag) & ~dup_] if WHICH == 'primary_noflag' else D[~D.DOI.isin(set(D[D.id.isin(flag)].DOI))]
TAG = ''
import os
FIXED = bool(os.environ.get('DIGIT_NOISE')) and os.environ.get('DIGIT_FIXED_ELIGIBILITY') == '1'   # fixed-eligibility mode: same source curves, same target specimens as the unperturbed analysis; only the measured values change
if os.environ.get('DIGIT_NOISE'):      # digitization-error sensitivity: multiplicative reading error on every stress and life value (worst case: all points treated as digitized)
    _sS, _sN, _sd = [float(v) for v in os.environ['DIGIT_NOISE'].split(',')]; _r = np.random.default_rng(int(_sd))
    D['S0'], D['N0'] = D.S.values.copy(), D.N.values.copy()      # unperturbed values (used for eligibility and stress-range masks in fixed-eligibility mode)
    D['S'] = D.S.values * np.exp(_r.normal(0, _sS, len(D))); D['N'] = D.N.values * 10 ** _r.normal(0, _sN, len(D)); TAG = f'_digit_S{_sS:g}_N{_sN:g}_seed{int(_sd)}' + ('_fixed' if os.environ.get('DIGIT_FIXED_ELIGIBILITY') == '1' else '')
if SLOPE_MIN != 0.0: TAG += f'_slope{SLOPE_MIN:g}'      # keeps the pair file of a slope-threshold sensitivity run from overwriting the default one
print(f'set={WHICH}: papers={D.DOI.nunique()} curves={D.id.nunique()} specimens={len(D)} failures={(D.ro == 0).sum()} runouts={D.ro.sum()}')
# ---------------- independent model code ----------------
def censored_lognormal_em(X, y, cens, iters=2000, tol=1e-10):
    """EM for right-censored normal regression; returns beta, sigma."""
    n = len(y); unc = ~cens; beta = np.linalg.lstsq(X[unc], y[unc], rcond=None)[0]; sig = max(np.std(y[unc] - X[unc] @ beta), 0.05)
    for _ in range(iters):
        mu = X @ beta; a = (y - mu) / sig; lam = np.where(cens, norm.pdf(a) / np.maximum(norm.sf(a), 1e-300), 0.0)
        m1 = np.where(cens, mu + sig * lam, y); m2 = np.where(cens, mu ** 2 + sig ** 2 + sig * (y + mu) * lam, y ** 2)
        beta_new = np.linalg.lstsq(X, m1, rcond=None)[0]; mun = X @ beta_new
        sig_new = math.sqrt(max(np.sum(m2 - 2 * mun * m1 + mun ** 2) / n, 1e-12))
        if abs(sig_new - sig) < tol and np.max(np.abs(beta_new - beta)) < tol: beta, sig = beta_new, sig_new; break
        beta, sig = beta_new, sig_new
    return beta, sig
def ols_fail(X, y):
    b = np.linalg.lstsq(X, y, rcond=None)[0]; r = y - X @ b; return b, math.sqrt(r @ r / (len(y) - X.shape[1]))
def tol_factor(n_eff, df): return nct.ppf(0.95, df, nc=Z95 * math.sqrt(n_eff)) / math.sqrt(n_eff)
# ---------------- pair records ----------------
METHODS = ['single_pred5', 'single_tol9595', 'pooled_censored_pred5', 'pooled_failures_only_pred5']
records = []   # (alloy,R,A,sid,B,n_inrange,{method:breaches})
for mat in ALLOYS:
    for R in (-1, 0.1):
        G = D[(D.mat == mat) & (D.R == R)].reset_index(drop=True)
        if G.DOI.nunique() < 6: continue
        y = np.log10(G.N.values); x = np.log10(G.S.values / np.median(G.S0.values if FIXED else G.S.values)); xb = np.log10(G.S0.values / np.median(G.S0.values)) if FIXED else x; C = G[['asb', 'hip', 'hthi']].values.astype(float); C = C[:, C.std(0) > 0]
        X = np.column_stack([np.ones(len(G)), x, C]); cens = G.ro.values.astype(bool); doi = G.DOI.values; cid = G.id.values
        curves = {}
        for sid in np.unique(cid):
            m = (cid == sid) & ~cens
            if m.sum() < 6 or len(np.unique(np.round(xb[m], 3))) < 3: continue
            Xs = np.column_stack([np.ones(m.sum()), x[m]]); curves[sid] = (doi[cid == sid][0], Xs, y[m], xb[m].min(), xb[m].max())
        for B in np.unique(doi):
            tr = doi != B; te = ~tr; Ktr = len(np.unique(doi[tr]))
            bc, sc = censored_lognormal_em(X[tr], y[tr], cens[tr]); bf, sf = ols_fail(X[tr][~cens[tr]], y[tr][~cens[tr]])
            for sid, (A, Xs, ys, lo, hi) in curves.items():
                if A == B: continue
                sel = te & (xb >= lo) & (xb <= hi)
                if not sel.any(): continue
                n = len(ys); bs, sd = ols_fail(Xs, ys); XtXi = np.linalg.inv(Xs.T @ Xs); x0 = np.column_stack([np.ones(sel.sum()), x[sel]])
                mu = x0 @ bs; q = np.einsum('ij,jk,ik->i', x0, XtXi, x0)             # x0' (X'X)^-1 x0 for each target specimen
                pred = mu - tdist.ppf(0.95, n - 2) * sd * np.sqrt(1 + q); tol = mu - np.array([tol_factor(1 / v, n - 2) for v in q]) * sd
                bound = {'single_pred5': pred, 'single_tol9595': tol, 'pooled_censored_pred5': X[sel] @ bc + Z05 * sc, 'pooled_failures_only_pred5': X[sel] @ bf + Z05 * sf}
                br = {m_: int(np.sum(~cens[sel] & (y[sel] < bound[m_]))) for m_ in METHODS}
                records.append((mat, R, A, sid, B, int(sel.sum()), br))
print('pair-curve records:', len(records))
pd.DataFrame([(r[0], r[1], r[2], r[3], r[4], r[5]) + tuple(r[6][m] for m in METHODS) for r in records], columns=['mat', 'R', 'A', 'sid', 'B', 'n'] + METHODS).to_csv(f'pairs_{WHICH}{TAG}.csv', index=False)
# ---------------- weighting (explicit loops) ----------------
def pair_rates(recs, method):
    """equal weight per source curve within (group,A,B); then per (A,B) averaged over groups"""
    by_group = {}
    for mat, R, A, sid, B, n, br in recs: by_group.setdefault((mat, R, A, B), []).append(br[method] / n)
    by_pair = {}
    for (mat, R, A, B), v in by_group.items(): by_pair.setdefault((A, B), []).append(sum(v) / len(v))
    return {k: sum(v) / len(v) for k, v in by_pair.items()}
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
rs = np.random.default_rng(SEED)
for method in METHODS:
    pr = pair_rates(records, method); As = sorted({a for a, b in pr}); Bs = sorted({b for a, b in pr}); est = weighted(pr, {a: 1 for a in As}, {b: 1 for b in Bs}); bs = []; bj = []
    allp = sorted(set(As) | set(Bs))
    for _ in range(2000):
        ca = np.bincount(rs.integers(0, len(As), len(As)), minlength=len(As)); cb = np.bincount(rs.integers(0, len(Bs), len(Bs)), minlength=len(Bs))
        bs.append(weighted(pr, dict(zip(As, ca)), dict(zip(Bs, cb))))
        cp = np.bincount(rs.integers(0, len(allp), len(allp)), minlength=len(allp)); wp = dict(zip(allp, cp)); bj.append(weighted(pr, wp, wp))   # same resampled paper set in both roles
    lo, hi = np.percentile(bs, [2.5, 97.5]); lj, hj = np.percentile(bj, [2.5, 97.5])
    print(f'{method:28s} {100 * est:5.2f}%  independent-roles [{100 * lo:5.2f}, {100 * hi:5.2f}]  same-paper-set [{100 * lj:5.2f}, {100 * hj:5.2f}]   pairs={len(pr)} sources={len(As)} targets={len(Bs)}')

# ---------------- paired differences against the pooled censored bound (same pairs; independent-role bootstrap) ----------------
rs = np.random.default_rng(SEED + 1000); prC = pair_rates(records, 'pooled_censored_pred5'); As = sorted({a for a, b in prC}); Bs = sorted({b for a, b in prC})
for method in ['single_pred5', 'single_tol9595', 'pooled_failures_only_pred5']:
    pr = pair_rates(records, method); dd = {k: pr[k] - prC[k] for k in prC}; est = weighted(dd, {a: 1 for a in As}, {b: 1 for b in Bs}); bs = []
    for _ in range(2000):
        ca = np.bincount(rs.integers(0, len(As), len(As)), minlength=len(As)); cb = np.bincount(rs.integers(0, len(Bs), len(Bs)), minlength=len(Bs)); bs.append(weighted(dd, dict(zip(As, ca)), dict(zip(Bs, cb))))
    lo, hi = np.percentile(bs, [2.5, 97.5]); print(f'paired difference {method} minus pooled_censored: {100 * est:+5.2f} pts [{100 * lo:+5.2f}, {100 * hi:+5.2f}]')
