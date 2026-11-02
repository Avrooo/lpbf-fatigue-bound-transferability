"""The numbers quoted in the manuscript and supplement, compared with the saved outputs (expected values frozen in expected_results.json)."""
import re, pandas as pd
from _common import *
E = expected()
def _line(fname, key): return re.search(r'^' + key + r'\s+([\d.]+)%', open(need('results', fname)).read(), re.M)
def test_primary_headline_rates():
    e = E['primary']
    for key, name in (('single_pred5', 'single_pred5'), ('single_tol9595', 'single_tol9595'), ('pooled_censored_pred5', 'pooled_censored'), ('pooled_failures_only_pred5', 'pooled_failures_only')):
        assert abs(float(_line('primary.out', key).group(1)) - e[name]) <= e['tol'], key
    t = open(need('results', 'primary.out')).read(); assert 'pairs=1562 sources=99 targets=123' in t
    assert abs(float(re.search(r'single_pred5 minus pooled_censored: \+([\d.]+)', t).group(1)) - e['diff_single_minus_pooled']) <= e['tol']
def test_matched_conservatism_values():
    m = pd.read_csv(need('results', 'matched_contrasts.csv')); e = E['matched_5pct']; s = m[(m.pooled_p == 0.05) & m.contrast.str.startswith('single')].iloc[0]; h = m[(m.pooled_p == 0.05) & m.contrast.str.startswith('hier')].iloc[0]
    assert abs(s.method_pct - e['single']) <= e['tol'] and abs(s.pooled_pct - e['pooled']) <= e['tol'] and abs(h.method_pct - e['hier']) <= e['tol'] and abs(h.estimate - e['hier_minus_pooled']) <= e['tol']
def test_within_paper_control():
    w = pd.read_csv(need('results', 'within_paper_summary.csv')); r = w[(w.subset == 'all alloys') & (w.bound == '5% prediction')].iloc[0]; e = E['within_paper']
    assert abs(r.within_breach_pct - e['within']) <= e['tol'] and abs(r.cross_breach_pct - e['cross']) <= e['tol'] and abs(r.gap_pp - e['gap']) <= e['tol']
    assert abs(r.gap_CI_low - e['gap_ci'][0]) <= e['ci_tol'] and abs(r.gap_CI_high - e['gap_ci'][1]) <= e['ci_tol']
def test_hierarchical_comparison():
    h = pd.read_csv(need('results', 'hierarchical_comparison.csv')).set_index('model'); e = E['hierarchical']
    assert abs(h.loc['pooled_censored', 'breach_pct'] - e['pooled']) <= e['tol'] and abs(h.loc['hierarchical_censored', 'breach_pct'] - e['hier']) <= e['tol']
    assert abs(h.iloc[2].mean_lower_bound_log10N - e['bound_diff_log10']) <= e['bound_tol'] and (h.converged_pct == 100).all()
def test_leave_one_paper_out_ranges():
    j = pd.read_csv(need('results', 'influence_jackknife.csv')); e = E['jackknife_gap_range']
    for role in ('source', 'target'):
        g = j[j.role == role].gap; assert abs(g.min() - e[role][0]) <= e['tol'] and abs(g.max() - e[role][1]) <= e['tol'], role
def test_digitization_summary():
    s = pd.read_csv(need('results', 'digitization', 'summary.csv')).set_index(['mode', 'scenario']); e = E['digitization_diff_means']
    for k in ('dynamic_moderate', 'dynamic_severe', 'fixed_moderate', 'fixed_severe'): assert abs(s.loc[tuple(k.split('_')), 'diff_mean'] - e[k]) <= e['tol'], k
def test_alloy_table():
    a = pd.read_csv(need('results', 'alloy_results.csv')); e = E['alloy_point_estimates']
    for alloy, vals in e.items():
        if alloy == 'tol': continue
        for m, v in zip(('single_pred5', 'single_tol9595', 'pooled_censored_pred5', 'pooled_failures_only_pred5'), vals):
            assert abs(a[(a.alloy == alloy) & (a.method == m)].estimate.iloc[0] - v) <= e['tol'], (alloy, m)
