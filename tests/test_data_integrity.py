"""Counts of the selected data and basic validity of the analysis table."""
import pandas as pd
from _common import *
def test_counts_in_build_output():
    t = open(need('results', 'build.out')).read()
    assert 'set=primary: papers=123 curves=424 specimens=3924 failures=3464 runouts=460' in t
    assert 'set=all: papers=147 curves=492 specimens=4609 failures=4083 runouts=526' in t
def test_primary_output_header():
    assert out_counts('primary.out') == tuple(expected()['counts']['primary'])
def test_analysis_table_valid():
    d = pd.read_csv(analysis_table()); e = expected()['counts']['primary']
    assert d.DOI.nunique() == e[0] and d.id.nunique() == e[1] and len(d) == e[2]
    assert (d.ro == 0).sum() == e[3] and (d.ro == 1).sum() == e[4]; assert set(d.ro.unique()) <= {0, 1}
    assert (d.S > 0).all() and (d.N > 0).all() and d.notna()[['S', 'N', 'DOI', 'mat', 'R']].all().all()
    assert set(d.mat) == {'Ti-6Al-4V', 'AlSi10Mg', 'IN718', '316L'} and set(d.R) == {-1.0, 0.1}
def test_life_window_applied_to_whole_papers():
    d = pd.read_csv(analysis_table()); g = d.groupby('DOI').N.agg(['min', 'max']); assert (g['min'] >= 1e3).all() and (g['max'] <= 2e7).all()
def test_curve_id_belongs_to_one_paper():
    d = pd.read_csv(analysis_table()); assert (d.groupby('id').DOI.nunique() == 1).all()
