"""Experiment 3: between-paper heterogeneity on ONE fixed set of observations.
Table 3 (tau_explained.py) lets the eligible papers/specimens change as covariates are added (missing values become an 'unknown' level).
Here the six covariate specifications are fitted to exactly the same specimens: only curves with build direction, specimen shape, load control,
test frequency and AM machine all recorded are kept, and alloy-R groups need >=6 such papers (the rule used throughout).
Random-intercept right-censored lognormal model per alloy-R group; weighted summaries use specimen counts as in Table 3.
Usage: python3 heterogeneity_fixed_sample.py XLSX JSON [excluded group, e.g. 'AlSi10Mg|0.1']
   (writes heterogeneity_fixed_sample.csv, or heterogeneity_fixed_sample_excl_<group>.csv when a group is excluded)"""
import sys, re, math, numpy as np, pandas as pd
from scipy.optimize import minimize
X_, J_ = sys.argv[1], sys.argv[2]; EXCL = sys.argv[3] if len(sys.argv) > 3 else None
sys.argv = ['x', X_, J_, 'primary']; exec(open('reproduce.py').read().split("# ---------------- pair records")[0])
def num(v):
    m = re.findall(r'[-+]?\d*\.?\d+', str(v)); return float(m[0]) if m and str(v) != 'nan' else np.nan
P = par.set_index('id')
att = pd.DataFrame({'dir': P['direciton of specimen \n(°)'].map(num), 'freq': P['frequency \n(Hz)'].map(num), 'load': P['load control'].astype(str).str.lower(), 'shape': P['specimens description'].astype(str).str.lower().str.extract(r'^(circular|rectangular|square)')[0], 'amm': P['AM machine'].astype(str).str.lower().str.strip()})
att['load'] = att['load'].replace({'nan': np.nan, '': np.nan}); att['amm'] = att['amm'].replace({'nan': np.nan, '': np.nan}); att['fband'] = pd.cut(att.freq, [0, 30, 100, 1e9], labels=['<=30', '30-100', '>100']).astype(object)
D2 = D.merge(att, left_on='id', right_index=True, how='left')
complete = D2[['dir', 'freq', 'load', 'shape', 'amm']].notna().all(axis=1); F = D2[complete].copy()
keep = F.groupby(['mat', 'R']).DOI.nunique(); dropped = keep[keep < 6]; keep = keep[keep >= 6]
if EXCL: em, er = EXCL.split('|'); keep = keep[[not (m_ == em and float(r_) == float(er)) for m_, r_ in keep.index]]; print('excluded from every specification by request:', EXCL)
F = F[F.set_index(['mat', 'R']).index.isin(keep.index)].reset_index(drop=True)
print(f'primary set: {D2.DOI.nunique()} papers, {len(D2)} specimens. Rows with all five attributes recorded: {int(complete.sum())} specimens, {D2[complete].DOI.nunique()} papers.')
print('groups with <6 such papers (dropped, rule not relaxed):', {f'{a}|{b}': int(v) for (a, b), v in dropped.items()})
print(f'FIXED SAMPLE: {F.DOI.nunique()} publications, {F.id.nunique()} curves, {len(F)} specimens ({int((F.ro == 0).sum())} failures, {int(F.ro.sum())} runouts) in {len(keep)} alloy-R groups')
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
    r = best; bound = abs(r.x[k + 1]) >= 5; ok = bool(r.success) and np.isfinite(r.fun) and abs(r.x[k]) < 5 and not bound
    return math.exp(r.x[k]), math.exp(r.x[k + 1]), ok, bool(r.success), bound
def dummies(col, G, minn=0):
    v = G[col].astype(object).where(G[col].notna(), 'unk')
    if minn: vc = v.value_counts(); v = v.where(v.map(vc) >= minn, 'other')
    d = pd.get_dummies(v, drop_first=True).astype(float); return d.loc[:, d.std() > 0]
SETS = [('Stress only', []), ('+ Surface, HIP, heat treatment', ['proc']), ('+ Build direction', ['proc', 'dir']), ('+ Specimen shape, load control', ['proc', 'dir', 'shape', 'load']), ('+ Frequency band', ['proc', 'dir', 'shape', 'load', 'fband']), ('+ AM machine identity', ['proc', 'dir', 'shape', 'load', 'fband', 'amm'])]
rows = []
for (mat, R) in keep.index:
    G = F[(F.mat == mat) & (F.R == R)].reset_index(drop=True); y = np.log10(G.N.values); x = np.log10(G.S.values / np.median(G.S.values)); cens = G.ro.values.astype(bool); g = G.DOI.values
    for name, cols in SETS:
        parts = [np.ones(len(G)), x]
        for c_ in cols:
            if c_ == 'proc':
                C = G[['asb', 'hip', 'hthi']].values.astype(float); C = C[:, C.std(0) > 0]; parts += [C[:, i] for i in range(C.shape[1])]
            else:
                d_ = dummies(c_, G, minn=(40 if c_ == 'amm' else 0)); parts += [d_.values[:, i] for i in range(d_.shape[1])]
        X = np.column_stack(parts); s, t, ok, osucc, bnd = fit_re(y, X, cens, g)
        rows.append(dict(level='group', specification=name, mat=mat, R=R, n_publications=G.DOI.nunique(), n_specimens=len(G), n_covariate_columns=X.shape[1] - 2, sigma_within=s, tau_between=t, converged=ok, optimizer_success=osucc, tau_at_zero_boundary=bnd))
T = pd.DataFrame(rows); T['ICC'] = T.tau_between ** 2 / (T.tau_between ** 2 + T.sigma_within ** 2)
out = []
t0 = T[T.specification == SETS[0][0]]; tau0 = np.average(t0.tau_between, weights=t0.n_specimens)
for name, _ in SETS:
    q = T[T.specification == name]; tw = np.average(q.tau_between, weights=q.n_specimens)
    out.append(dict(level='all groups (specimen-weighted)', specification=name, mat='all', R='all', n_publications=F.DOI.nunique(), n_specimens=len(F), n_covariate_columns=q.n_covariate_columns.mean(), sigma_within=np.average(q.sigma_within, weights=q.n_specimens), tau_between=tw,
                    ICC=np.average(q.ICC, weights=q.n_specimens), tau_variance_removed_pct_vs_stress_only=100 * (1 - (tw / tau0) ** 2), converged=f'{int(q.converged.sum())}/{len(q)}'))
R_ = pd.concat([pd.DataFrame(out), T], ignore_index=True); R_['converged'] = R_.converged.astype(str)
R_.to_csv('heterogeneity_fixed_sample.csv' if not EXCL else 'heterogeneity_fixed_sample_excl_' + EXCL.replace('|', '_R') + '.csv', index=False, float_format='%.4f')
pd.set_option('display.width', 250); pd.set_option('display.max_columns', 30)
print(pd.DataFrame(out).round(3).to_string()); print(T[T.converged == False].round(4).to_string()); print(T[T.specification == SETS[-1][0]].round(3).to_string())
print('machine dummies vs papers per group (all-covariate model):'); 
for (mat, R) in keep.index:
    G = F[(F.mat == mat) & (F.R == R)]; print(f'  {mat} R={R}: papers={G.DOI.nunique()} machines={G.amm.nunique()}  papers sharing a machine with another paper: {int((G.groupby("amm").DOI.nunique() > 1).sum())} of {G.amm.nunique()} machines')
