# Main-text Figures 3, 5 and 6, redrawn from saved outputs (no typed values).
#   Fig 3 <- results/primary.out        Fig 5 <- results/learning_curves.out        Fig 6 <- results/comparability.out
# Run from the repository root: python3 figures/make_fig3_fig5_fig6.py
import re, sys
from decimal import Decimal, ROUND_HALF_UP
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.family': 'Liberation Sans', 'font.size': 7.5, 'axes.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6})
BL, OR, GR = '#2a78d6', '#eb6834', '#5b6770'
R = 'results/'; O = 'figures/'
def f1(x): return str(Decimal(str(x)).quantize(Decimal('0.1'), rounding=ROUND_HALF_UP))
def clean(ax):
    for s in ('top', 'right'): ax.spines[s].set_visible(False)
def save(fig, name): fig.savefig(O + name + '.png', dpi=300); fig.savefig(O + name + '.pdf'); plt.close(fig); print('wrote', name)

# ---------- Figure 3 ----------
rows = {}
for ln in open(R + 'primary.out', encoding='utf-8'):
    m = re.match(r'^(single_pred5|single_tol9595|pooled_censored_pred5|pooled_failures_only_pred5)\s+([\d.]+)%\s+independent-roles \[\s*([\d.]+),\s*([\d.]+)\]', ln)
    if m: rows[m.group(1)] = tuple(float(m.group(i)) for i in (2, 3, 4))
assert len(rows) == 4, rows
spec = [('single_pred5', 'Single-study 5% prediction bound', OR, 's', True, 4.0), ('single_tol9595', 'Single-study 95/95 tolerance bound', OR, 's', False, 3.0),
        ('pooled_censored_pred5', 'Pooled 5% bound, runouts censored', BL, 'o', True, 1.0), ('pooled_failures_only_pred5', 'Pooled 5% bound, failures only', BL, 'o', False, 0.0)]
fig, ax = plt.subplots(figsize=(7.0, 3.3), dpi=300)
for key, lab, c, mk, filled, y in spec:
    est, lo, hi = rows[key]; yy = y
    ax.plot([lo, hi], [yy, yy], '-', color=c, lw=1.4, solid_capstyle='butt'); ax.plot([lo, lo], [yy - .08, yy + .08], '-', color=c, lw=1.4); ax.plot([hi, hi], [yy - .08, yy + .08], '-', color=c, lw=1.4)
    ax.plot(est, yy, mk, ms=7, mfc=c if filled else 'white', mec=c, mew=1.4, zorder=5)
    ax.text(-1.0, yy, lab, ha='right', va='center', fontsize=8); ax.text(37, yy, f'{f1(est)}  [{f1(lo)}, {f1(hi)}]', va='center', fontsize=8)
ax.axvline(5, color='#444444', lw=1, ls=(0, (4, 3))); ax.text(5.6, 5.0, 'Nominal 5% reference', color='#555555', fontsize=8, va='center')
ax.text(37, 5.0, 'Breach rate [95% CI]', fontweight='bold', fontsize=8.5, va='center')
ax.text(-1.0 - 22, 5.0, 'Fitted to one paper', color=OR, fontweight='bold', fontsize=8.5, ha='left', va='center'); ax.text(-1.0 - 22, 2.0, 'Pooled model (target paper excluded)', color=BL, fontweight='bold', fontsize=8.5, ha='left', va='center', bbox=dict(fc='white', ec='none', pad=1.5), zorder=6)
ax.axhline(2.5, color='#dddddd', lw=0.8, xmin=0.0, xmax=0.88)
ax.set_xlim(0, 36); ax.set_ylim(-0.6, 5.4); ax.set_yticks([]); ax.set_xticks(range(0, 36, 5)); ax.grid(axis='x', color='#e6e6e6', lw=0.6); ax.set_axisbelow(True)
ax.set_xlabel('Cross-paper confirmed-breach rate (%)', fontsize=8.5)
for s in ('top', 'right', 'left'): ax.spines[s].set_visible(False)
fig.subplots_adjust(left=0.36, right=0.78, top=0.97, bottom=0.2); save(fig, 'Fig3_breach_rates')

# ---------- Figure 5 ----------
txt = open(R + 'learning_curves.out', encoding='utf-8').read().splitlines(); panels = []
for tag, nmin in (('A', 20), ('B', 10)):
    i = next(k for k, l in enumerate(txt) if l.startswith(f'FIXED TARGETS {tag} '))
    nt = int(re.search(r'(\d+) target papers', txt[i]).group(1)); pts = []
    for l in txt[i + 1:]:
        m = re.match(r'\s+K=\s*(\w+):\s+([\d.]+)%\s+\[\s*([\d.]+),\s*([\d.]+)\]', l)
        if m: pts.append((m.group(1), float(m.group(2)), float(m.group(3)), float(m.group(4))))
        elif l.strip() == '': break
    panels.append((nt, nmin, pts))
fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.9), dpi=300, sharey=True)
for ax, (nt, nmin, pts), c, lab in zip(axs, panels, (OR, BL), ('(a)', '(b)')):
    xs = {'1': 0, '3': 1, '5': 2, '10': 3, '20': 4, 'all': 5}
    x = [xs[k] for k, *_ in pts]; y = [p[1] for p in pts]; lo = [p[1] - p[2] for p in pts]; hi = [p[3] - p[1] for p in pts]
    ax.axvline(4.5, color='#e6e6e6', lw=1.2); ax.axhline(5, color='#444444', lw=0.9, ls=(0, (4, 3)), zorder=1); ax.text(-0.35, 5.9, 'Nominal 5% reference', color='#555555', fontsize=7.5)
    ax.errorbar(x, y, yerr=[lo, hi], fmt='-o', color=c, ms=6, lw=1.4, capsize=3, zorder=3)
    for xi, p in zip(x, pts): ax.text(xi, p[3] + 1.2, f1(p[1]), ha='center', fontsize=8)
    ax.set_xticks(range(6)); ax.set_xticklabels(['1', '3', '5', '10', '20', 'All']); ax.set_xlim(-0.5, 5.5); ax.set_ylim(0, 42)
    ax.set_title(f'{nt} target papers (≥{nmin} other papers)', fontsize=8.5, fontweight='bold'); ax.text(-0.12 if lab == '(a)' else -0.04, 1.08, lab, transform=ax.transAxes, fontweight='bold', fontsize=9); ax.grid(axis='y', color='#ececec', lw=0.6); ax.set_axisbelow(True); clean(ax)
axs[0].set_ylabel('Cross-paper confirmed-breach rate (%)', fontsize=8.5); fig.text(0.5, 0.02, 'Number of training papers, K  (All = every other paper in the group)', ha='center', fontsize=8.5); fig.text(0.985, 0.02, 'Bars: 95% bootstrap intervals', ha='right', fontsize=7, color='#555555')
fig.tight_layout(rect=(0, 0.05, 1, 1)); save(fig, 'Fig5_learning_curves')

# ---------- Figure 6 ----------
L = []
for ln in open(R + 'comparability.out', encoding='utf-8'):
    m = re.match(r'^L(\d) .*?\s+([\d.]+)% \[\s*([\d.]+),\s*([\d.]+)\] \|\s+\d+\s+(\d+)\s+\d+\s+\d+\s+\d+ \|\s+([\d.]+)% \[\s*([\d.]+),\s*([\d.]+)\]', ln)
    if m: L.append(tuple(float(v) for v in m.groups()))
assert len(L) == 6, L
labs = ['Same alloy\nand stress ratio', '+ Surface,\nHIP, heat\ntreatment', '+ Build\ndirection', '+ Specimen shape\nand load control', '+ Frequency\nband', '+ Same\nAM machine']
fig, ax = plt.subplots(figsize=(7.0, 4.3), dpi=300); ax.axhline(5, color='#777777', lw=0.9, ls=(0, (4, 3)), zorder=1); ax.text(0.5, 1.6, 'Nominal 5% reference', ha='center', color='#666666', fontsize=7.5)
for off, idx, c, mk, lab in ((-0.0, 1, OR, 's', 'Single-study 5% prediction bound'), (0.0, 5, BL, 'o', 'Pooled 5% bound, runouts censored')):
    xs = list(range(6)); y = [r[idx] for r in L]; lo = [r[idx] - r[idx + 1] for r in L]; hi = [r[idx + 2] - r[idx] for r in L]
    ax.errorbar(xs, y, yerr=[lo, hi], fmt='-' + mk, color=c, ms=6, lw=1.4, capsize=3, label=lab, zorder=3)
    for xi, r in zip(xs, L):
        t = f'{f1(r[idx])}\n[{f1(r[idx + 1])}, {f1(r[idx + 2])}]'
        if xi < 4: ax.text(xi, r[idx + 2] + 1.0, t, ha='center', va='bottom', fontsize=6.5, color=c)
        else: ax.text(xi + 0.1, r[idx] + ((3.0 if xi == 4 else 8.0) if idx == 5 else 0), t, ha='left', va='center', fontsize=6.5, color=c)
ax.set_xticks(range(6)); ax.set_xticklabels([f'{l}\n({int(r[4]):,})' for l, r in zip(labs, L)], fontsize=7); ax.set_xlim(-0.5, 5.8); ax.set_ylim(0, 50)
ax.set_xlabel('Matching (each level adds to the previous)\n(Number of paper pairs)', fontsize=8.5); ax.set_ylabel('Cross-paper confirmed-breach rate (%)', fontsize=8.5); clean(ax); ax.legend(frameon=False, loc='upper right', fontsize=7.5); fig.tight_layout(); save(fig, 'Fig6_matched_conditions')
