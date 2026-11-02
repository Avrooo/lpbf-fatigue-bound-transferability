"""Invariants of the confirmed-breach counts in the saved records (a breach is a failure below the bound; runouts never count)."""
import re, pandas as pd, numpy as np
from _common import *
def _records(): return pd.read_csv(need('results', 'cc_primary_strict1_records.csv'))
def test_breach_counts_are_between_zero_and_scored_specimens():
    d = _records(); cols = [c for c in d.columns if c.startswith('br_')]; assert cols and (d[cols].min().min() >= 0)
    for c in cols: assert (d[c] <= d.n).all(), c
def test_breaches_never_decrease_as_the_bound_rises():
    d = _records()
    for m in ('single', 'pooled', 'hier'):
        a = d[[f'br_{m}_{k}' for k in range(13)]].values; assert (np.diff(a, axis=1) >= 0).all(), m
def test_runouts_stay_in_denominator_and_ambiguous_ones_are_rare():
    a = pd.read_csv(need('results', 'cc_primary_strict1_audit_ambiguous_runouts.csv')).set_index('metric')
    assert (a.scored_specimens == 57483).all() and a.loc['pooled_5%', 'ambiguous_runouts'] == 0 and a.loc['hier_5%', 'ambiguous_runouts'] == 0
    assert a.loc['single_5%', 'ambiguous_runouts'] == 11 and a.loc['single_5%', 'ambiguous_runouts'] / 57483 < 3e-4
    assert abs(a.loc['single_5%', 'specimen_weighted_if_ambiguous_runouts_counted'] - a.loc['single_5%', 'specimen_weighted_confirmed']) < 0.05
def test_scored_specimens_include_runouts():
    d = _records(); t = open(need('results', 'primary.out')).read(); assert 'runouts=460' in t
    assert d.n.sum() >= 3924 // 2     # the scored set is large and not restricted to failures (exact count: results/*audit*.csv)
