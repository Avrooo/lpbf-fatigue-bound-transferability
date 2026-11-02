"""Source and target papers are separate; pooled and hierarchical fits exclude the target paper."""
import pandas as pd
from _common import *
def test_source_never_equals_target():
    for f in ('pairs_primary.csv', 'cc_primary_strict1_records.csv'):
        d = pd.read_csv(need('results', f)); assert (d.A != d.B).all(), f
def test_pairs_stay_within_alloy_and_stress_ratio_group():
    d = pd.read_csv(need('results', 'pairs_primary.csv')); assert d[['mat', 'R']].drop_duplicates().shape[0] == 8
    e = expected()['counts']; assert d[['A', 'B']].drop_duplicates().shape[0] == e['distinct_pairs']
    assert sum(g[['A', 'B']].drop_duplicates().shape[0] for _, g in d.groupby(['mat', 'R'])) == e['group_pair_occurrences']
    assert d.A.nunique() == e['source_papers'] and d.B.nunique() == expected()['counts']['primary'][0]
def test_hierarchical_training_set_excludes_exactly_the_target_paper():
    d = pd.read_csv(need('results', 'hierarchical_fit_diagnostics.csv')); full = d[d.target == 'ALL PAPERS (no hold-out)'].set_index(['mat', 'R']).n_train_papers
    lo = d[d.target != 'ALL PAPERS (no hold-out)']; assert len(lo) == 128 and len(full) == 8
    for (m, r), g in lo.groupby(['mat', 'R']): assert (g.n_train_papers == full[(m, r)] - 1).all(), (m, r)
def test_every_target_paper_is_scored_by_other_papers_only():
    d = pd.read_csv(need('results', 'cc_primary_strict1_records.csv')); assert d.B.nunique() == expected()['counts']['primary'][0]
