"""Coverage-conservatism frontier: 13-point grid, monotonicity, and the matched-height contrast recomputed from the records."""
import numpy as np, pandas as pd
from _common import *
PG = [0.00001, 0.0001, 0.0003, 0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.10, 0.20, 0.30, 0.50]
def test_grid_and_monotonicity():
    f = pd.read_csv(need('results', 'cc_primary_strict1_frontier.csv'))
    for m in ('single', 'pooled', 'hier'):
        g = f[f.method == m]; assert len(g) == 13
        assert (np.diff(g.breach_pct.values) >= -1e-9).all() and (np.diff(g.mean_bound_log10N.values) > 0).all()
def test_nominal_5pct_points_equal_the_headline():
    f = pd.read_csv(need('results', 'cc_primary_strict1_frontier.csv')); e = expected(); get = lambda m: f[f.method == m].iloc[8].breach_pct
    assert abs(get('single') - e['primary']['single_pred5']) <= e['primary']['tol'] and abs(get('pooled') - e['primary']['pooled_censored']) <= e['primary']['tol']
    assert abs(get('hier') - e['hierarchical']['hier']) <= e['hierarchical']['tol']
def test_matched_height_contrast_recomputed_from_records():
    d = pd.read_csv(need('results', 'cc_primary_strict1_records.csv'))
    def est(col):
        r = (d[f'br_{col[0]}_{col[1]}'] / d.n).groupby([d.mat, d.R, d.A, d.B]).mean().groupby(['A', 'B']).mean().groupby('B').mean().mean() * 100; return r
    def h(m, k): return d[f'h_{m}_{k}'].groupby([d.mat, d.R, d.A, d.B]).mean().groupby(['A', 'B']).mean().groupby('B').mean().mean()
    H = h('pooled', 8); hs = np.array([h('single', k) for k in range(13)]); rs = np.array([est(('single', k)) for k in range(13)])
    single_at_H = np.interp(H, hs, rs); e = expected()['matched_5pct']
    assert abs(H - 4.2036) < 5e-4 and abs(single_at_H - e['single']) <= e['tol'] and abs(est(('pooled', 8)) - e['pooled']) <= e['tol']
def test_matched_contrasts_file():
    m = pd.read_csv(need('results', 'matched_contrasts.csv')); r = m[(m.pooled_p == 0.05) & m.contrast.str.startswith('single')].iloc[0]; e = expected()['matched_5pct']
    assert abs(r.estimate - e['diff']) <= e['tol'] and abs(r.ci_low - e['diff_ci'][0]) <= e['ci_tol'] and abs(r.ci_high - e['diff_ci'][1]) <= e['ci_tol']
    assert len(m) == 10 and (m[m.contrast.str.startswith('single')].ci_low > 0).all()
