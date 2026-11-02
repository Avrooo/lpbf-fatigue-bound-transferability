"""Experiment 1: within-paper held-out validation of single-study bounds.
For every eligible source curve (>=6 failures, >=3 stress levels) the failures are split repeatedly into a training part and a held-out part.
The bound is fitted on the training failures only (training keeps >=6 failures and >=3 stress levels) and scored on the held-out failures
(plus the curve's runouts, which are never fitted and never counted as breaches) inside the training stress range.
The cross-paper rate uses exactly the same source curves (full-curve fit, scored on other papers inside the curve's stress range).
Usage: python3 within_paper_control.py XLSX JSON [n_splits=100] [seed=5]   (run from this folder; outputs CSV files here)"""
import sys, math, numpy as np, pandas as pd
X_, J_ = sys.argv[1], sys.argv[2]; NSPL = int(sys.argv[3]) if len(sys.argv) > 3 else 100; SD = int(sys.argv[4]) if len(sys.argv) > 4 else 5
sys.argv = ['x', X_, J_, 'primary']; exec(open('reproduce.py').read().split("# ---------------- pair records")[0])
HOLD = 0.30; NB = 2000
# ---------------- within-paper repeated splits ----------------
rng = np.random.default_rng(SD); within = {}      # (mat,R,sid) -> dict
cross_rec = []                                    # (mat,R,A,sid,B,n,br_pred,br_tol) for the same source curves
n_all_curves = 0; excl_small = 0; excl_nosplit = 0
def bounds(bs, sd, XtXi, n, x0):
    mu = x0 @ bs; q = np.einsum('ij,jk,ik->i', x0, XtXi, x0)
    return mu - tdist.ppf(0.95, n - 2) * sd * np.sqrt(1 + q), mu - np.array([tol_factor(1 / v, n - 2) for v in q]) * sd
for mat in ALLOYS:
    for R in (-1, 0.1):
        G = D[(D.mat == mat) & (D.R == R)].reset_index(drop=True)
        if G.DOI.nunique() < 6: continue
        y = np.log10(G.N.values); x = np.log10(G.S.values / np.median(G.S.values)); cens = G.ro.values.astype(bool); doi = G.DOI.values; cid = G.id.values
        elig = {}
        for sid in np.unique(cid):
            m = (cid == sid) & ~cens
            if m.sum() < 6 or len(np.unique(np.round(x[m], 3))) < 3: continue
            elig[sid] = m; n_all_curves += 1
        good = {}
        for sid, m in elig.items():
            fi = np.where(m)[0]; n = len(fi)
            if n < 8: excl_small += 1; continue           # cannot keep >=6 training failures and hold out >=2
            nt = max(2, int(round(HOLD * n))); ro = np.where((cid == sid) & cens)[0]; recs = {'pred': [], 'tol': [], 'predF': [], 'tolF': []}; used = 0; seen = set(); nfail_eval = []
            for r in range(NSPL):
                ok = False
                for _ in range(200):
                    te_i = rng.choice(fi, nt, replace=False); tr_i = np.setdiff1d(fi, te_i)
                    if len(np.unique(np.round(x[tr_i], 3))) >= 3: ok = True; break
                if not ok: continue
                assert not (set(te_i.tolist()) & set(tr_i.tolist())) and len(tr_i) >= 6 and not cens[tr_i].any() and not cens[te_i].any()   # held-out failures never enter training; training keeps >=6 failures
                Xs = np.column_stack([np.ones(len(tr_i)), x[tr_i]]); bs, sd = ols_fail(Xs, y[tr_i]); XtXi = np.linalg.inv(Xs.T @ Xs); lo, hi = x[tr_i].min(), x[tr_i].max()
                ev_f = te_i[(x[te_i] >= lo) & (x[te_i] <= hi)]; ev_r = ro[(x[ro] >= lo) & (x[ro] <= hi)]; ev = np.concatenate([ev_f, ev_r])
                if len(ev_f) == 0: continue
                assert set(ev_f.tolist()) <= set(te_i.tolist()) and not (set(ev.tolist()) & set(tr_i.tolist())) and cens[ev_r].all() and not cens[ev_f].any()   # scored specimens are held-out failures + runouts only
                pb, tb = bounds(bs, sd, XtXi, len(tr_i), np.column_stack([np.ones(len(ev)), x[ev]])); nf = len(ev_f)
                bp = int(np.sum(y[ev_f] < pb[:nf])); bt = int(np.sum(y[ev_f] < tb[:nf]))
                recs['pred'].append(bp / len(ev)); recs['tol'].append(bt / len(ev)); recs['predF'].append(bp / nf); recs['tolF'].append(bt / nf)
                used += 1; seen.update(ev_f.tolist()); nfail_eval.append(nf)
            if used < 0.5 * NSPL: excl_nosplit += 1; continue   # not enough valid splits: reported, not relaxed
            doi_A = doi[fi][0]; good[sid] = True
            within[(mat, R, sid)] = dict(mat=mat, R=R, sid=sid, A=doi_A, n_fail=n, n_test=nt, used=used, uniq=len(seen), mean_nfail_eval=np.mean(nfail_eval), **{k: float(np.mean(v)) for k, v in recs.items()})
        # cross-paper records for the same source curves (full fit, scored on other papers inside the curve's range)
        for sid in good:
            m = elig[sid]; Xs = np.column_stack([np.ones(m.sum()), x[m]]); bs, sd = ols_fail(Xs, y[m]); XtXi = np.linalg.inv(Xs.T @ Xs); lo, hi = x[m].min(), x[m].max(); A = doi[m][0]
            for B in np.unique(doi):
                if B == A: continue
                sel = (doi == B) & (x >= lo) & (x <= hi)
                if not sel.any(): continue
                pb, tb = bounds(bs, sd, XtXi, m.sum(), np.column_stack([np.ones(sel.sum()), x[sel]])); f = ~cens[sel]
                cross_rec.append((mat, R, A, sid, B, int(sel.sum()), int(np.sum(y[sel][f] < pb[f])), int(np.sum(y[sel][f] < tb[f]))))
W = pd.DataFrame(within.values()); C = pd.DataFrame(cross_rec, columns=['mat', 'R', 'A', 'sid', 'B', 'n', 'pred', 'tol'])
W.to_csv('within_paper_curve_records.csv', index=False, float_format='%.6f'); C.to_csv('within_paper_cross_pair_records.csv', index=False)   # inputs of the audit (audit_within_paper.py)
print(f'source curves eligible for the cross-paper analysis: {n_all_curves}; excluded for <8 failures: {excl_small}; excluded for too few valid splits: {excl_nosplit}; within-paper eligible curves: {len(W)} from {W.A.nunique()} papers')
# ---------------- aggregation ----------------
def paper_values(Wsub, Csub, col):
    """per source paper A: within = mean over its eligible curves; cross = mean over targets B of the curve-averaged pair rate."""
    wv = Wsub.groupby('A')[col].mean(); cc = Csub.assign(r=Csub[col] / Csub.n); pg = cc.groupby(['mat', 'R', 'A', 'B', 'sid']).r.mean().groupby(['mat', 'R', 'A', 'B']).mean().groupby(['A', 'B']).mean()
    cv = pg.groupby('A').mean(); idx = wv.index.intersection(cv.index); return wv[idx], cv[idx]
def boot_papers(wv, cv, rs):
    n = len(wv); w = wv.values; c = cv.values; ests = (w.mean(), c.mean(), c.mean() - w.mean()); bb = []
    for _ in range(NB):
        i = rs.integers(0, n, n); bb.append((w[i].mean(), c[i].mean(), (c[i] - w[i]).mean()))
    bb = np.array(bb); ci = np.percentile(bb, [2.5, 97.5], axis=0); return ests, ci
def pair_rate_dict(Csub, col):
    g = {}
    for r in Csub.itertuples(): g.setdefault((r.mat, r.R, r.A, r.B), []).append(getattr(r, col) / r.n)
    p = {}
    for (mat, R, A, B), v in g.items(): p.setdefault((A, B), []).append(sum(v) / len(v))
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
def headline_cross(Csub, col, rs):
    pr = pair_rate_dict(Csub, col); As = sorted({a for a, b in pr}); Bs = sorted({b for a, b in pr}); est = weighted(pr, {a: 1 for a in As}, {b: 1 for b in Bs}); bs = []
    for _ in range(NB):
        ca = np.bincount(rs.integers(0, len(As), len(As)), minlength=len(As)); cb = np.bincount(rs.integers(0, len(Bs), len(Bs)), minlength=len(Bs)); bs.append(weighted(pr, dict(zip(As, ca)), dict(zip(Bs, cb))))
    return est, np.percentile(bs, [2.5, 97.5])
def summarize(Wsub, Csub, label, extra=True):
    rows = []
    for name, col, fcol in (('5% prediction', 'pred', 'pred'), ('95/95 tolerance', 'tol', 'tol')):
        rs = np.random.default_rng(SD + 100)
        wv, cv = paper_values(Wsub, Csub, col)
        # within values: held-out rate; use the same column names in W
        ests, ci = boot_papers(wv, cv, rs)
        wF = Wsub.groupby('A')[col + 'F'].mean().loc[wv.index].mean()
        r = dict(subset=label, bound=name, n_papers=len(wv), n_curves=len(Wsub), n_test_failures_per_split=round(Wsub.mean_nfail_eval.sum(), 1), n_test_failures_unique=int(Wsub.uniq.sum()),
                 within_breach_pct=100 * ests[0], within_CI_low=100 * ci[0, 0], within_CI_high=100 * ci[1, 0], cross_breach_pct=100 * ests[1], cross_CI_low=100 * ci[0, 1], cross_CI_high=100 * ci[1, 1],
                 gap_pp=100 * ests[2], gap_CI_low=100 * ci[0, 2], gap_CI_high=100 * ci[1, 2], within_breach_pct_failures_only_denominator=100 * wF)
        if extra:
            e, c2 = headline_cross(Csub, col, np.random.default_rng(SD + 200)); r.update(cross_breach_pct_target_weighting=100 * e, cross_target_weighting_CI_low=100 * c2[0], cross_target_weighting_CI_high=100 * c2[1])
        rows.append(r)
    return rows
S = pd.DataFrame(summarize(W, C, 'all alloys'))
S.insert(1, 'n_splits', NSPL); S.insert(2, 'holdout_fraction', HOLD); S.to_csv('within_paper_summary.csv', index=False, float_format='%.3f')
rows = []
for mat in ALLOYS:
    Ws, Cs = W[W.mat == mat], C[C.mat == mat]
    if len(Ws) >= 5 and Ws.A.nunique() >= 5: rows += summarize(Ws, Cs, mat, extra=False)
    else: print(f'{mat}: only {len(Ws)} curves from {Ws.A.nunique()} papers; not reported')
pd.DataFrame(rows).to_csv('within_paper_by_alloy.csv', index=False, float_format='%.3f')
pd.set_option('display.width', 250); pd.set_option('display.max_columns', 30)
print(S.round(2).T.to_string()); print(pd.DataFrame(rows).round(2).to_string())
