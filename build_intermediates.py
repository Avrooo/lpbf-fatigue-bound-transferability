"""Note on intermediate files. reproduce.py, comparability.py, tau_explained.py and coverage_test.py read ONLY the two raw files
(the 'processing' field of the JSON is parsed on the fly into the surface / HIP / heat-treatment indicators; there are no hidden .pkl inputs).
This script writes the analysis table that all of them build, so it can be inspected: analysis_table_{primary|all}.csv.
Usage: python3 -I build_intermediates.py XLSX JSON"""
import sys
X_, J_ = sys.argv[1], sys.argv[2]
for w in ('primary', 'all'):
    sys.argv = ['x', X_, J_, w]; exec(open('reproduce.py').read().split("# ---------------- independent model code")[0])
    D.to_csv(f'analysis_table_{w}.csv', index=False); print('wrote', f'analysis_table_{w}.csv', len(D))
