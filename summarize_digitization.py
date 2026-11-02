"""Mean (min-max) of the headline rates over the five seeds for each scenario and mode; writes results/digitization/summary.csv"""
import re, glob, os, numpy as np, pandas as pd
rows = []
for mode in ('dynamic', 'fixed'):
    for sc, lab in (('0.01_0.02', 'moderate'), ('0.03_0.05', 'severe')):
        fs = sorted(glob.glob(f'results/digitization/{mode}/out_{sc}_*.txt'))
        if not fs: continue
        v = {k: [] for k in ('single_pred5', 'single_tol9595', 'pooled_censored', 'diff', 'pairs')}
        for f in fs:
            t = open(f).read()
            for k, pat in (('single_pred5', r'^single_pred5\s+([\d.]+)%'), ('single_tol9595', r'^single_tol9595\s+([\d.]+)%'), ('pooled_censored', r'^pooled_censored_pred5\s+([\d.]+)%|^pooled_censored\s+([\d.]+)%')):
                m = re.search(pat, t, re.M); v[k].append(float([g for g in m.groups() if g][0]))
            v['diff'].append(float(re.search(r'paired difference single_pred5 minus pooled_censored:\s*([+\-\d.]+)', t).group(1))); v['pairs'].append(int(re.search(r'pairs=(\d+)', t).group(1)))
        rows.append(dict(mode=mode, scenario=lab, n_seeds=len(fs), **{f'{k}_mean': round(np.mean(a), 2) for k, a in v.items()}, **{f'{k}_min': round(min(a), 2) for k, a in v.items() if k != 'pairs'}, **{f'{k}_max': round(max(a), 2) for k, a in v.items() if k != 'pairs'}))
S = pd.DataFrame(rows); S.to_csv('results/digitization/summary.csv', index=False); pd.set_option('display.width', 250); print(S.to_string())
