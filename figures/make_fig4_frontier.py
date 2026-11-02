# Main-text Figure 4 (coverage-conservatism frontier), drawn from results/cc_primary_strict1_frontier.csv and results/decision_policy_summary.csv. Run from the repository root: python3 figures/make_fig4_frontier.py
import pandas as pd, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib import font_manager
plt.rcParams.update({'font.family': 'Liberation Sans', 'font.size': 7.5, 'axes.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6})
R = 'results/'; F = pd.read_csv(R + 'cc_primary_strict1_frontier.csv'); T = pd.read_csv(R + 'decision_policy_summary.csv')
BL, OR, GR = '#2a78d6', '#eb6834', '#5b6770'; fig, ax = plt.subplots(figsize=(5.5, 3.6), dpi=300)
for m, c, lab, ls in (('single', OR, 'Single-study', '-'), ('pooled', BL, 'Pooled censored', '-'), ('hier', GR, 'Hierarchical censored', '--')):
    d = F[F.method == m].sort_values('p'); ax.plot(d.mean_bound_log10N, d.breach_pct, ls, color=c, lw=1.3, label=lab, zorder=3)
    d5 = d[np.isclose(d.p, 0.05)].iloc[0]; ax.plot(d5.mean_bound_log10N, d5.breach_pct, 'o', ms=6, mfc=c, mec='white', mew=0.8, zorder=5)
ax.axhline(5, color='#999999', lw=0.7, ls=(0, (3, 3)), zorder=1); ax.text(2.05, 5.9, '5% reference', color='#777777', fontsize=7)
h = F[(F.method == 'pooled') & np.isclose(F.p, 0.05)].mean_bound_log10N.iloc[0]; ax.axvline(h, color='#bbbbbb', lw=0.7, ls=':', zorder=1)
s = F[F.method == 'single'].sort_values('mean_bound_log10N'); y15 = np.interp(h, s.mean_bound_log10N, s.breach_pct)
ax.annotate('', xy=(h, 6.5), xytext=(h, y15), arrowprops=dict(arrowstyle='<->', color='#444444', lw=0.8), zorder=4); ax.annotate(f'{y15:.1f}% (single-study)\nvs 6.5% (pooled)\nat the same mean bound', xy=(h, (y15 + 6.5) / 2), xytext=(4.62, 10.0), fontsize=7, va='center', ha='left', color='#333333', arrowprops=dict(arrowstyle='-', color='#444444', lw=0.6))
t = F[F.method == 'single_95/95'].iloc[0]; ax.plot(t.mean_bound_log10N, t.breach_pct, 'D', ms=5.5, mfc='white', mec=OR, mew=1.2, zorder=5); ax.annotate('Single-study 95/95', (t.mean_bound_log10N, t.breach_pct), xytext=(3.35, 33), fontsize=7, color=OR, arrowprops=dict(arrowstyle='-', color=OR, lw=0.6))
p0 = T[T.rule.str.contains('\\+ 0.0')].iloc[0]; ax.plot(p0.mean_bound_log10N, p0.breach_pct, 's', ms=5, mfc='#222222', mec='white', mew=0.6, zorder=6); ax.annotate('Rule: never exceed\nthe pooled bound', (p0.mean_bound_log10N, p0.breach_pct), xytext=(2.35, 13.0), fontsize=7, color='#222222', arrowprops=dict(arrowstyle='-', color='#222222', lw=0.6))
ax.set_xlabel('Mean lower-bound life on the scored specimens (log$_{10}$ cycles); lower = more conservative'); ax.set_ylabel('Confirmed-breach rate on other papers (%)'); ax.set_xlim(2, 5.5); ax.set_ylim(0, 60)
for sp in ('top', 'right'): ax.spines[sp].set_visible(False)
ax.plot([], [], 'o', color='#888888', ms=5, label='Nominal 5% bound'); ax.legend(frameon=False, loc='upper left', fontsize=7.5, bbox_to_anchor=(0.0, 1.0)); fig.tight_layout(); fig.savefig('figures/Fig4_coverage_conservatism.png', dpi=300); fig.savefig('figures/Fig4_coverage_conservatism.pdf'); print('ok')
