import os, re, unittest, json
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); RES = os.path.join(ROOT, 'results')
def path(*p): return os.path.join(ROOT, *p)
def need(*p):
    f = path(*p)
    if not os.path.exists(f): raise unittest.SkipTest(f'missing {os.path.join(*p)} (run sh run_all.sh first)')
    return f
def analysis_table():
    """Derived analysis table written by build_intermediates.py (not shipped: it is a subset of the source database)."""
    for f in (os.environ.get('LPBF_ANALYSIS_TABLE', ''), path('analysis_table_primary.csv')):
        if f and os.path.exists(f): return f
    raise unittest.SkipTest('analysis_table_primary.csv not found (run sh run_all.sh, or set LPBF_ANALYSIS_TABLE)')
def out_counts(fname):
    t = open(need('results', fname)).read(); m = re.search(r'papers=(\d+) curves=(\d+) specimens=(\d+) failures=(\d+) runouts=(\d+)', t); return tuple(int(x) for x in m.groups())
def expected(): return json.load(open(path('tests', 'expected_results.json')))
