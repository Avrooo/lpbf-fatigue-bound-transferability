"""Supplementary figures S1-S5, drawn only from files in results/.  Usage: python3 make_supplement_figures.py  (writes supplement/FigS*.png and .pdf)"""
import re, glob, os, numpy as np, pandas as pd, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
for f in glob.glob('/usr/share/fonts/truetype/liberation/LiberationSans-*.ttf'): fm.fontManager.addfont(f)
plt.rcParams.update({'font.family': 'Liberation Sans', 'font.size': 7, 'axes.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6, 'axes.spines.top': False, 'axes.spines.right': False, 'pdf.fonttype': 42})
BL, OR, GR, LG = '#2a78d6', '#eb6834', '#5b6770', '#a9b3ba'; R = 'results/'; O = 'results/supplement/'; os.makedirs(O, exist_ok=True)
def save(fig, n): fig.savefig(O + n + '.png', dpi=300, bbox_inches='tight'); fig.savefig(O + n + '.pdf', bbox_inches='tight'); plt.close(fig)
# ---- S1: bootstrap-interval calibration by simulation ----
rows = []
ORDER = ['cov_Ti_m1_K20', 'cov_Ti_m1_K60', 'cov_Ti_m1_K120', 'cov_Ti_01_K30', 'cov_316L_m1_K22']
for f in sorted(glob.glob(R + 'cov_*.out'), key=lambda f: ORDER.index(os.path.basename(f)[:-4])):
    t = open(f).read(); g = re.search(r'group (\S+) papers (\d+).*? K (\d+)', t); lab = f"{g.group(1).replace('|', ', R=')}\nK={g.group(3)}"
    for m in ('single_pred5', 'pooled_failures_only_pred5'):
        v = re.search(r'\[' + m + r'\] VARIANT same-paper-set resampling: coverage=\s*([\d.]+)%', t); i = re.search(r'^' + m + r'\s.*?coverage=\s*([\d.]+)%.*?MC s\.e\. ([\d.]+)', t, re.M)
        rows.append(dict(setting=lab, method=m, indep=float(i.group(1)), indep_se=float(i.group(2)), same=float(v.group(1))))
S = pd.DataFrame(rows); S.to_csv(O + 'FigS1_data.csv', index=False)
fig, ax = plt.subplots(figsize=(5.2, 2.6)); sets = list(dict.fromkeys(S.setting)); x = np.arange(len(sets))
for m, c, off, nm in (('single_pred5', OR, -0.12, 'Single-study 5% prediction bound'), ('pooled_failures_only_pred5', BL, 0.12, 'Pooled failures-only 5% bound')):
    d = S[S.method == m].set_index('setting').loc[sets]
    ax.errorbar(x + off - 0.0, d.indep, yerr=1.96 * d.indep_se, fmt='o', color=c, ms=4, lw=0.8, capsize=1.5, label=nm + ', independent-role bootstrap')
    ax.plot(x + off, d.same, 's', mfc='white', mec=c, ms=4, mew=0.9, label=nm + ', same-paper-set bootstrap')
ax.axhline(95, color=GR, lw=0.7, ls='--'); ax.text(len(sets) - 0.55, 95.2, '95% nominal', color=GR, ha='right', va='bottom')
ax.set_xticks(x); ax.set_xticklabels(sets); ax.set_ylabel('Simulated coverage of 95% interval (%)'); ax.set_ylim(90, 101); ax.legend(frameon=False, fontsize=6, loc='lower left', bbox_to_anchor=(0.0, -0.62), ncol=1)
save(fig, 'FigS1_bootstrap_calibration')
# ---- S2: alternative bootstrap intervals ----
B = pd.read_csv(R + 'cc_primary_strict1_audit_bootstrap.csv'); W = pd.read_csv(R + 'cc_primary_strict1_audit_weighting.csv').set_index('metric')['headline_target_then_source']
mets = ['single_5%', 'single_95/95', 'pooled_5%', 'hier_5%']; names = ['Single-study\n5% prediction', 'Single-study\n95/95 tolerance', 'Pooled\ncensored', 'Hierarchical\ncensored']
kinds = [('independent', 'Independent roles (primary)', BL, 'o'), ('same_paper_set', 'Same paper set for both roles', OR, 's'), ('source_only', 'Source papers only', GR, '^'), ('target_only', 'Target papers only', LG, 'v')]
fig, ax = plt.subplots(figsize=(5.2, 2.9))
for mi, m in enumerate(mets):
    for ki, (k, lab, c, mk) in enumerate(kinds):
        r = B[(B.bootstrap == k) & (B.metric == m)].iloc[0]; y = mi + (ki - 1.5) * 0.17
        ax.plot([r.ci_low, r.ci_high], [y, y], color=c, lw=1.4, solid_capstyle='butt', label=lab if mi == 0 else None); ax.plot(W[m], y, mk, color=c, ms=3.5)
ax.set_yticks(range(4)); ax.set_yticklabels(names); ax.invert_yaxis(); ax.set_xlabel('Confirmed-breach rate on other papers (%), point estimate and 95% interval'); ax.axvline(5, color=GR, lw=0.6, ls='--'); ax.legend(frameon=False, fontsize=6, loc='lower right')
save(fig, 'FigS2_bootstrap_alternatives')
# ---- S3: leave-one-paper-out influence ----
J = pd.read_csv(R + 'influence_jackknife.csv'); full = 20.71; rs = np.random.default_rng(3)
fig, axs = plt.subplots(1, 2, figsize=(5.2, 2.0), sharey=True)
for ax, role in zip(axs, ('source', 'target')):
    d = J[J.role == role]; ax.plot(d.gap, rs.uniform(-0.25, 0.25, len(d)), 'o', color=BL, ms=2.6, alpha=0.7, mew=0); ax.axvline(full, color=OR, lw=0.9); ax.set_ylim(-0.35, 0.4); ax.set_yticks([]); ax.set_xlim(19.5, 21.5); ax.spines['left'].set_visible(False)
    ax.set_title(f'Leave one {role} paper out\nn = {len(d)}, range {d.gap.min():.2f}–{d.gap.max():.2f}', fontsize=7); ax.set_xlabel('Single-study minus pooled (pp)'); ax.text(full + 0.03, 0.38, 'all papers: 20.71', color=OR, fontsize=6, va='top')
save(fig, 'FigS3_leave_one_paper_out')
# ---- S4: alternative weightings ----
Wt = pd.read_csv(R + 'cc_primary_strict1_audit_weighting.csv'); cols = ['headline_target_then_source', 'source_paper_weighted', 'pair_weighted', 'specimen_weighted']; lab = ['Target paper\n(primary)', 'Source paper', 'Source–target\npair', 'Specimen']
fig, ax = plt.subplots(figsize=(4.6, 2.5)); x = np.arange(len(cols))
for (m, nm, c, mk, ls) in (('single_5%', 'Single-study 5% prediction', OR, 'o', '-'), ('single_95/95', 'Single-study 95/95 tolerance', OR, 's', ':'), ('pooled_5%', 'Pooled censored', BL, 'o', '-'), ('hier_5%', 'Hierarchical censored', BL, 's', ':')):
    v = Wt.set_index('metric').loc[m, cols].values.astype(float); ax.plot(x, v, marker=mk, ls=ls, color=c, ms=3.5, lw=1, label=nm)
ax.axhline(5, color=GR, lw=0.6, ls='--'); ax.set_xticks(x); ax.set_xticklabels(lab); ax.set_ylabel('Confirmed-breach rate (%)'); ax.set_ylim(0, 30); ax.set_xlabel('Weighting of the breach rate'); ax.legend(frameon=False, fontsize=6, loc='upper center', bbox_to_anchor=(0.5, -0.32), ncol=2)
save(fig, 'FigS4_alternative_weightings')
# ---- S5: stress ratio and source-curve size ----
Z = pd.read_csv(R + 'influence_by_size_and_R.csv'); fig, axs = plt.subplots(1, 2, figsize=(5.2, 2.3), gridspec_kw={'width_ratios': [1, 2]}, sharey=True)
for ax, sub, ttl in ((axs[0], Z[Z.source_failures.str.startswith('R=')], 'Stress ratio'), (axs[1], Z[~Z.source_failures.str.startswith('R=')], 'Failures in the source curve')):
    x = np.arange(len(sub)); ax.bar(x - 0.18, sub.single, 0.34, color=OR, label='Single-study'); ax.bar(x + 0.18, sub.pooled, 0.34, color=BL, label='Pooled censored')
    ax.set_xticks(x); ax.set_xticklabels([f"{a.replace('R=-1.0','R = −1').replace('R=0.1','R = 0.1').replace('>=','≥')}\n{c} curves" for a, c in zip(sub.source_failures, sub.source_curves)]); ax.set_title(ttl, fontsize=7); ax.axhline(5, color=GR, lw=0.6, ls='--')
    for xi, (s, p) in enumerate(zip(sub.single, sub.pooled)): ax.text(xi - 0.18, s + 0.6, f'{s:.1f}', ha='center', fontsize=6); ax.text(xi + 0.18, p + 0.6, f'{p:.1f}', ha='center', fontsize=6)
axs[0].set_ylabel('Confirmed-breach rate (%)'); axs[1].legend(frameon=False, fontsize=6, loc='upper right'); save(fig, 'FigS5_stress_ratio_and_source_size')
print('figures written'); print(S.round(1).to_string())
