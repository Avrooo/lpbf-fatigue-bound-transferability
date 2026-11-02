"""Re-derive the headline rates from the saved pair-curve records with the staged weights (curve, source paper, target paper)."""
import pandas as pd
from _common import *
def _stage(d, col):
    r = d.assign(rate=d[col] / d.n)
    pv = r.groupby(['mat', 'R', 'A', 'B']).rate.mean().groupby(['A', 'B']).mean()      # curves within source paper; both stress ratios averaged
    return 100 * pv.groupby('B').mean().mean()                                          # source papers of each target, then target papers
def test_headline_rates_from_pair_records():
    d = pd.read_csv(need('results', 'pairs_primary.csv')); e = expected()['primary']
    for col, key in (('single_pred5', 'single_pred5'), ('single_tol9595', 'single_tol9595'), ('pooled_censored_pred5', 'pooled_censored'), ('pooled_failures_only_pred5', 'pooled_failures_only')):
        assert abs(_stage(d, col) - e[key]) <= e['tol'], (col, _stage(d, col))
def test_equal_weights_would_give_a_different_number():
    d = pd.read_csv(need('results', 'pairs_primary.csv')).assign(rate=lambda x: x.single_pred5 / x.n)
    spec = 100 * d.single_pred5.sum() / d.n.sum(); assert abs(spec - _stage(d, 'single_pred5')) > 1      # specimen weighting differs (19.0 vs 27.2)
def test_a_paper_with_more_curves_does_not_get_more_weight():
    toy = pd.DataFrame({'mat': 'X', 'R': -1.0, 'A': ['a1'] * 3 + ['a2'], 'B': 'b', 'sid': [1, 2, 3, 4], 'n': 10, 'c': [10, 10, 10, 0]})
    assert abs(_stage(toy, 'c') - 50.0) < 1e-9      # source paper a1 (3 curves, 100%) and a2 (1 curve, 0%) count equally
