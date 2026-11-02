# -*- coding: utf-8 -*-
"""Builds supplement/Supplementary_Information.docx from results/ and the figures of make_supplement_figures.py.
Usage: python3 make_supplement_doc.py [analysis_table_primary.csv]   (default: analysis_table_primary.csv written by build_intermediates.py)"""
import sys, re, glob, platform, numpy as np, pandas as pd, docx, matplotlib
from docx.shared import Pt, Inches, Cm
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml.ns import qn
from docx.oxml import OxmlElement
R = 'results/'; O = 'results/supplement/'; AT = sys.argv[1] if len(sys.argv) > 1 else 'analysis_table_primary.csv'
d = docx.Document(); st = d.styles['Normal']; st.font.name = 'Times New Roman'; st.font.size = Pt(10.5); st.element.rPr.rFonts.set(qn('w:eastAsia'), 'Times New Roman')
for s in d.sections: s.left_margin = s.right_margin = Cm(2.2); s.top_margin = s.bottom_margin = Cm(2.2)
def para(t, bold=False, size=None, align=None, italic=False, after=6):
    p = d.add_paragraph(); r = p.add_run(t); r.bold = bold; r.italic = italic
    if size: r.font.size = Pt(size)
    p.alignment = align if align is not None else WD_ALIGN_PARAGRAPH.JUSTIFY; p.paragraph_format.space_after = Pt(after); return p
def head(t): return para(t, bold=True, size=12, align=WD_ALIGN_PARAGRAPH.LEFT, after=4)
def cap(t): return para(t, size=9.5, after=4)
def borders(tbl):
    t = tbl._tbl; pr = t.tblPr; b = OxmlElement('w:tblBorders')
    for e in ('top', 'bottom'):
        x = OxmlElement('w:' + e); x.set(qn('w:val'), 'single'); x.set(qn('w:sz'), '8'); x.set(qn('w:color'), '000000'); b.append(x)
    pr.append(b)
    for c in tbl.rows[0].cells:
        tcPr = c._tc.get_or_add_tcPr(); bd = OxmlElement('w:tcBorders'); x = OxmlElement('w:bottom'); x.set(qn('w:val'), 'single'); x.set(qn('w:sz'), '4'); x.set(qn('w:color'), '000000'); bd.append(x); tcPr.append(bd)
def table(header, rows, widths=None, size=8.5):
    t = d.add_table(rows=1, cols=len(header)); t.alignment = 1
    for c, h in zip(t.rows[0].cells, header): c.text = ''; r = c.paragraphs[0].add_run(h); r.bold = True; r.font.size = Pt(size)
    for row in rows:
        cs = t.add_row().cells
        for c, v in zip(cs, row): c.text = ''; r = c.paragraphs[0].add_run(str(v)); r.font.size = Pt(size)
    for row in t.rows:
        for c in row.cells:
            for p in c.paragraphs: p.paragraph_format.space_after = Pt(0); p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    if widths:
        t.autofit = False
        for i_, w_ in enumerate(widths): t.columns[i_].width = Cm(w_)
        for row in t.rows:
            for c, w in zip(row.cells, widths): c.width = Cm(w)
    borders(t); d.add_paragraph().paragraph_format.space_after = Pt(2); return t
def fig(name, width, text):
    d.add_picture(O + name + '.png', width=Inches(width)); d.paragraphs[-1].alignment = WD_ALIGN_PARAGRAPH.CENTER; cap(text)
def parse_out(f):
    t = open(f).read(); h = re.search(r'papers=(\d+) curves=(\d+) specimens=(\d+)', t); o = {'papers': int(h.group(1)), 'curves': int(h.group(2)), 'specimens': int(h.group(3))}
    for m in re.finditer(r'^(single_pred5|single_tol9595|pooled_censored_pred5|pooled_failures_only_pred5)\s+([\d.]+)%\s+independent-roles \[\s*([\d.]+),\s*([\d.]+)\].*?pairs=(\d+)', t, re.M): o[m.group(1)] = (float(m.group(2)), float(m.group(3)), float(m.group(4))); o['pairs'] = int(m.group(5))
    return o
f1 = lambda a: f'{a:.1f}'; ci = lambda e, l, h, n=1: f'{e:.{n}f} [{l:.{n}f}, {h:.{n}f}]'
para('Supplementary Information', bold=True, size=14, align=WD_ALIGN_PARAGRAPH.CENTER, after=2)
para('Cross-study transferability of lower-tail S–N bounds for laser powder bed fusion alloys: coverage and conservatism', bold=True, size=11, align=WD_ALIGN_PARAGRAPH.CENTER, after=8)
para('This file supports the main manuscript. All results were written by the scripts of the reproducibility package (Section S7) from the files in results/. The few counts quoted from the database description or the main text are named in the notes. Section, table and figure numbers without the prefix S refer to the main text. The main text contains the principal results, including the coverage–conservatism frontier (Figure 4) and the matched-bound comparison (Table 2). This file adds the full counts, the calibration checks and the robustness and sensitivity results behind them.')
# ================= S1
head('S1. Dataset eligibility')
b = open(R + 'build.out').read(); pr = parse_out(R + 'primary.out'); al = re.search(r'set=all: papers=(\d+) curves=(\d+) specimens=(\d+) failures=(\d+) runouts=(\d+)', b); ds = open(R + 'data_screens.out').read(); obs = re.search(r'S-N rows (\d+) curves (\d+)', ds)
para('The source database is FatigueData-AM2022 (Zhang and Xu, reference 35 of the main text). We kept L-PBF specimens tested uniaxially at R = −1 or 0.1 in Ti-6Al-4V, AlSi10Mg, IN718 or 316L. For the primary analysis we then excluded whole publications that contained any selected observation outside 10³–2 × 10⁷ cycles. Source curves needed at least six failures at three or more stress levels. Table S1a gives the counts at each step and Table S1b the counts by alloy and stress ratio.')
cap('Table S1a. Selection steps and counts.')
prim = pd.read_csv(AT); pairs = pd.read_csv(R + 'pairs_primary.csv'); nsrc = pairs[['A', 'sid']].drop_duplicates().shape[0]; npp = pairs[['A', 'B']].drop_duplicates().shape[0]
pf = int((prim.ro == 0).sum()); pro = int((prim.ro == 1).sum())
table(['Step', 'Papers', 'S–N curves', 'Specimens', 'Failures', 'Runouts'], [
    ['Database, S–N part (all materials and processes)', '459', f'{int(obs.group(2)):,}', f'{int(obs.group(1)):,}', '—', '—'],
    ['L-PBF, uniaxial, four alloys, R = −1 or 0.1 (all-paper set)', al.group(1), al.group(2), f'{int(al.group(3)):,}', f'{int(al.group(4)):,}', al.group(5)],
    ['Primary set (life window 10³–2 × 10⁷ applied to whole papers)', str(pr['papers']), str(pr['curves']), f"{pr['specimens']:,}", f'{pf:,}', str(pro)],
    ['Eligible source curves (≥ 6 failures, ≥ 3 stress levels)', '99', '296', '—', '—', '—']], widths=[7.3, 1.5, 1.9, 1.9, 1.7, 2.0])
cap(f'Notes: the database count of 459 articles and the 99 papers and 296 eligible source curves are quoted from the database description and Table 1 of the main text. The database row counts S–N observations and datasets for all materials, so it is not a subset of the later rows. The primary analysis compares source and target papers in {npp:,} distinct pairs. Papers, curves and counts in the last three rows are from results/build.out, results/primary.out and the analysis table.')
cap('Table S1b. Primary set by alloy and stress ratio. Pairs are source–target paper pairs compared within the group.')
rows = []
for (m, r), g in prim.groupby(['mat', 'R']):
    pp = pairs[(pairs.mat == m) & (pairs.R == r)][['A', 'B']].drop_duplicates().shape[0]
    rows.append([m, f'{r:g}', g.DOI.nunique(), g.id.nunique(), len(g), int((g.ro == 0).sum()), int((g.ro == 1).sum()), pp])
rows.sort(key=lambda x: (['Ti-6Al-4V', '316L', 'AlSi10Mg', 'IN718'].index(x[0]), x[1]))
table(['Alloy', 'R', 'Papers', 'Curves', 'Specimens', 'Failures', 'Runouts', 'Pairs'], rows, widths=[2.6, 1.2, 1.6, 1.6, 2.0, 1.8, 1.8, 1.6])
cap(f'Notes: a paper can contribute to more than one group, so group counts do not add to the totals; a pair compared under both stress ratios appears in each group (group sum {sum(r[-1] for r in rows):,} against {npp:,} distinct pairs). {nsrc} of the 296 eligible source curves appear in at least one pair.')
# ================= S2
head('S2. Statistical validation')
para('Intervals for the breach rates come from the two-way paper bootstrap (Section 3.5). To test them, we simulated data with a known data-generating process calibrated to one real alloy and stress-ratio group (curve slopes, between-paper, between-curve and within-curve scatter), using the real stress designs of that group. In each replication K papers were drawn, new lives were simulated, the same estimator and weights as in the main analysis were applied, and the 95% bootstrap interval was compared with the truth at the same K. Figure S1 gives the coverage in five settings. Each used 300 simulated data sets (200 for K = 120 and 500 for 316L) and 500 bootstrap resamples per interval, and the truth at each K was estimated from at least 1,500 independent replications.')
S = pd.read_csv(O + 'FigS1_data.csv'); sg = S[S.method == 'single_pred5']; pl = S[S.method == 'pooled_failures_only_pred5']
fig('FigS1_bootstrap_calibration', 5.6, 'Figure S1. Simulated coverage of the 95% bootstrap interval. Filled circles: independent-role bootstrap (primary); error bars show ±1.96 Monte Carlo standard errors. Open squares: bootstrap that resamples one paper set for both roles. Dashed line: nominal 95%. K is the number of simulated papers.')
para(f"The independent-role interval contained the true rate in {sg.indep.min():.1f}–{sg.indep.max():.1f}% of replications for the single-study bound and {pl.indep.min():.1f}–{pl.indep.max():.1f}% for the pooled failures-only bound. The same-paper-set interval gave {sg.same.min():.1f}–{sg.same.max():.1f}% and {pl.same.min():.1f}–{pl.same.max():.1f}%. We therefore use the independent-role interval, which is conservative for the single-study rate. The simulation covers the 5% prediction bound and the pooled failures-only bound, with no covariates, no runouts and Gaussian errors. It does not cover the 95/95 tolerance bound, the censored pooled model or the hierarchical model.")
para('Calibration of the approximate 95/95 tolerance bound was assessed only empirically. On held-out specimens from their own papers it was breached by ' + ci(*pd.read_csv(R + 'within_paper_summary.csv').query("bound=='95/95 tolerance'")[['within_breach_pct', 'within_CI_low', 'within_CI_high']].iloc[0].values, 2) + '% of failures (Table 6 of the main text), as expected for a conservative bound. This is not an estimated failure probability.')
fig('FigS2_bootstrap_alternatives', 5.6, 'Figure S2. Point estimates and 95% intervals for four bounds under four bootstrap schemes. Independent roles resamples source and target papers separately (primary). Source only and target only resample one role and hold the other fixed.')
Bt = pd.read_csv(R + 'cc_primary_strict1_audit_bootstrap.csv'); wd = Bt.pivot(index='metric', columns='bootstrap', values='width')
para(f"Intervals for the single-study rates widen from {wd.loc['single_5%','same_paper_set']:.1f} points (same paper set) to {wd.loc['single_5%','independent']:.1f} points (independent roles), because the single-study bound depends on both source and target papers. For the pooled and hierarchical models, resampling source papers alone gives intervals only {wd.loc['pooled_5%','source_only']:.1f} and {wd.loc['hier_5%','source_only']:.1f} points wide, because those models are fitted without reference to a particular source paper. Their uncertainty comes from the target papers.")
# ================= S3
head('S3. Coverage–conservatism analysis')
para('The lower-tail probability p of each bound was varied over the 13-value grid in Table S2. Single-study bounds use the prediction limit at level p, and the pooled and hierarchical bounds use μ̂ + z_p σ̂ and μ̂ + z_p (σ̂² + τ̂²)^½. All methods were scored on identical specimens with the staged weights of Section 3.5, and the covariates of the pooled and hierarchical models were chosen from the training papers only. Mean bound life is the mean of log₁₀ N of the bound over the scored specimens. Intervals are 95% percentile intervals from the independent-role bootstrap with 2,000 resamples.')
PGRID = [0.00001, 0.0001, 0.0003, 0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.10, 0.20, 0.30, 0.50]; Fr = pd.read_csv(R + 'cc_primary_strict1_frontier.csv'); nm = {'single': 'Single-study prediction', 'pooled': 'Pooled censored', 'hier': 'Hierarchical censored', 'single_95/95': 'Single-study 95/95 tolerance'}
rows = []
for m in ('single', 'pooled', 'hier'):
    for k_, (_, r) in enumerate(Fr[Fr.method == m].iterrows()): rows.append([nm[m], f'{PGRID[k_]:g}', ci(r.breach_pct, r.breach_lo, r.breach_hi, 2), ci(r.mean_bound_log10N, r.bound_lo, r.bound_hi, 3)])
r = Fr[Fr.method == 'single_95/95'].iloc[0]; rows.append([nm['single_95/95'], '—', ci(r.breach_pct, r.breach_lo, r.breach_hi, 2), ci(r.mean_bound_log10N, r.bound_lo, r.bound_hi, 3)])
cap('Table S2. Full coverage–conservatism frontier. Breach rate in percent and mean lower-bound life in log₁₀ cycles, with 95% bootstrap intervals.')
table(['Method', 'p', 'Breach rate (%) [95% CI]', 'Mean bound life (log₁₀ N) [95% CI]'], rows, widths=[4.6, 1.5, 4.5, 5.0])
Ct = pd.read_csv(R + 'cc_primary_strict1_contrasts.csv')
def cv(q): return Ct[Ct.quantity == q].iloc[0]
a, bq, c = cv('height_at_5pct_breach_single'), cv('height_at_5pct_breach_pooled'), cv('height_at_5pct_breach_hier')
para(f"Interpolation was linear along each frontier and used only inside the sampled range. The single-study frontier has only two grid points below 3.14 log₁₀ cycles (Table S2), so its 5% crossing at {a.estimate:.2f} log₁₀ cycles (95% CI {a.ci_low:.2f}–{a.ci_high:.2f}) is indicative only and is not used in the main text. The pooled and hierarchical frontiers reach 5% at {bq.estimate:.3f} and {c.estimate:.3f} log₁₀ cycles. Table S3 gives the paired contrasts at the reference heights of Table 2, including the hierarchical-minus-pooled differences that are summarized in the text of Section 4.2.")
M = pd.read_csv(R + 'matched_contrasts.csv'); sp = M[M.contrast.str.startswith('single')]; hh = M[M.contrast.str.startswith('hier')].set_index('pooled_p'); rows = []
for _, r in sp.iterrows():
    h_ = hh.loc[r.pooled_p]; rows.append([f'{r.pooled_p*100:g}%', f'{r.ref_height_log10N:.2f}', ci(r.method_pct, r.method_lo, r.method_hi), ci(r.pooled_pct, r.pooled_lo, r.pooled_hi), ci(r.estimate, r.ci_low, r.ci_high), ci(h_.estimate, h_.ci_low, h_.ci_high)])
cap('Table S3. Matched-bound-life contrasts. The reference height is the mean life of the pooled bound at the stated p, recomputed within each bootstrap resample. Breach rates in percent; differences in percentage points; all intervals are 95% bootstrap intervals (2,000 resamples).')
table(['Pooled level', 'Ref. height (log₁₀ N)', 'Single-study (%)', 'Pooled (%)', 'Single − pooled (pp)', 'Hierarchical − pooled (pp)'], rows, widths=[1.8, 2.2, 3.1, 3.1, 3.1, 3.3])
cap('Note: the marginal intervals for the single-study and pooled rates are read from the same resamples as the differences and can differ slightly from the intervals of Table 2 and Section 4.1, which use the resampling sequence of the main analysis.')
# ================= S4
head('S4. Robustness')
J = pd.read_csv(R + 'influence_jackknife.csv'); js = J[J.role == 'source']; jt = J[J.role == 'target']
para(f"Figure S3 repeats the single-study minus pooled difference after removing one paper at a time. The difference stayed between {js.gap.min():.2f} and {js.gap.max():.2f} points when each source paper was removed and between {jt.gap.min():.2f} and {jt.gap.max():.2f} points when each target paper was removed, against 20.71 points for all papers. No single paper accounts for the contrast.")
fig('FigS3_leave_one_paper_out', 5.4, 'Figure S3. Leave-one-paper-out influence on the difference between the single-study and pooled breach rates. Each point is one omitted paper. The line marks the full-sample value.')
Wt = pd.read_csv(R + 'cc_primary_strict1_audit_weighting.csv').set_index('metric')
para(f"Figure S4 gives the breach rates under four ways of weighting. Target-paper weighting is the primary choice (Section 3.5). The ordering of methods did not change. Specimen weighting gives lower rates for every method ({Wt.loc['single_5%','specimen_weighted']:.1f}% single-study, {Wt.loc['pooled_5%','specimen_weighted']:.1f}% pooled), because it gives large papers more influence.")
fig('FigS4_alternative_weightings', 4.8, 'Figure S4. Confirmed-breach rates under alternative weightings of the source–target comparisons. Lines join points for the same bound.')
Ra = pd.read_csv(R + 'cc_primary_strict1_audit_ambiguous_runouts.csv').set_index('metric')
para(f"Figure S5 splits the single-study and pooled rates by stress ratio and by the number of failures in the source curve. Values are point estimates without intervals. The gap between the methods was larger at R = 0.1 than at R = −1, and it was smaller for the ten source curves with 20 or more failures, which come from only six papers. Counting the {int(Ra.loc['single_5%','ambiguous_runouts'])} ambiguous runouts (runouts stopped below the bound) as breaches changed the specimen-weighted single-study rate from {Ra.loc['single_5%','specimen_weighted_confirmed']:.3f}% to {Ra.loc['single_5%','specimen_weighted_if_ambiguous_runouts_counted']:.3f}%. Choosing the pooled-model covariates from the training papers only (the strict setting) gave the same rates as the earlier selection on whole groups (results/cc_strict0.out and cc_strict1.out).")
fig('FigS5_stress_ratio_and_source_size', 5.4, 'Figure S5. Single-study and pooled breach rates by stress ratio and by number of failures in the source curve (point estimates).')
tx = open(R + 'tau_explained.out').read(); rows = []
for m in re.finditer(r'^(none|\+ [^\d]+?)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s*$', tx, re.M): rows.append([m.group(1), m.group(3), m.group(4), m.group(5)])
cap('Table S4. Between-paper SD (τ̂) and within-paper SD (σ̂) of log₁₀ life as covariates are added, with the sample allowed to change (specimen-weighted means over eight alloy and stress-ratio groups).')
table(['Covariates added', 'τ̂', 'σ̂', 'ICC'], rows, widths=[6.0, 2.2, 2.2, 2.2])
para('The sample changes between rows because not every paper reports every covariate, so these values mix a change in model with a change in papers. This is why the main text (Section 4.6, Table 4) uses a fixed sample, and the variance change in the last row here should not be interpreted.')
# ================= S5
head('S5. Data quality and digitization')
ds_n = int(re.search(r'duplicate points within a curve: (\d+)', ds).group(1)); su = re.search(r'S>UTS: rows (\d+) curves (\d+)', ds); sl = re.search(r'slope > -1 .*?: (\d+)', ds)
para(f"Across the whole database, the screens in results/data_screens.out found no non-positive or missing values, {ds_n} duplicate points within curves, {su.group(1)} rows in {su.group(2)} curves with stress above the recorded ultimate strength, and {sl.group(1)} curves with a slope above −1 among curves with at least six failures. The main text flags curves by the slope rule (slope ≥ 0 by default) and removes duplicates. Table S5a gives the effect of each exclusion and of the documentation-quality subsets, and Table S5b the effect of simulated reading error.")
sets = [('Primary set', 'primary.out'), ('All papers, no life window', 'all.out'), ('Flagged curves and duplicates removed (slope ≥ 0)', 'primary_noflag.out'), ('Papers with flagged curves removed (slope ≥ 0)', 'primary_nopaperflag.out'),
        ('Flagged curves removed (slope > −1)', 'primary_noflag_slope_gt_minus1.out'), ('Papers with flagged curves removed (slope > −1)', 'primary_nopaperflag_slope_gt_minus1.out'), ('Higher-rated subset', 'primary_highrated.out'), ('Lower-rated subset', 'primary_lowrated.out')]
rows = []
for lab, f in sets:
    o = parse_out(R + f); rows.append([lab, o['papers'], o['curves'], f"{o['specimens']:,}", o['pairs'], ci(*o['single_pred5']), ci(*o['single_tol9595']), ci(*o['pooled_censored_pred5']), ci(*o['pooled_failures_only_pred5'])])
cap('Table S5a. Breach rates (%) [95% CI, independent-role bootstrap] under data exclusions and documentation-quality subsets.')
table(['Dataset', 'Papers', 'Curves', 'Specimens', 'Pairs', 'Single 5%', 'Single 95/95', 'Pooled censored', 'Pooled failures-only'], rows, widths=[3.8, 1.1, 1.1, 1.4, 1.1, 2.0, 2.0, 2.0, 2.0], size=7.5)
Dg = pd.read_csv(R + 'digitization/summary.csv'); rows = []
for _, r in Dg.iterrows():
    rows.append([{'dynamic': 'Eligibility re-evaluated', 'fixed': 'Eligibility fixed'}[r['mode']], {'moderate': 'Moderate (1% stress, 5% life)', 'severe': 'Severe (3% stress, 12% life)'}[r.scenario], int(round(r.pairs_mean)),
                 f"{r.single_pred5_mean:.1f} ({r.single_pred5_min:.1f}–{r.single_pred5_max:.1f})", f"{r.pooled_censored_mean:.1f} ({r.pooled_censored_min:.1f}–{r.pooled_censored_max:.1f})", f"{r.diff_mean:.1f} ({r.diff_min:.1f}–{r.diff_max:.1f})"])
cap('Table S5b. Simulated reading error added to every stress and life value. Breach rates (%) are means over five random replicates, with the range in parentheses. In the first mode eligibility and pairing were re-evaluated on the perturbed data, so the number of pairs changed; in the second mode they were fixed at the unperturbed values (1,562 pairs).')
table(['Mode', 'Scenario', 'Pairs', 'Single 5%', 'Pooled censored', 'Difference (pp)'], rows, widths=[3.6, 4.4, 1.3, 2.9, 2.9, 2.9])
para('The perturbation is random and independent between observations. It cannot represent systematic extraction errors such as a misread axis or observations assigned to the wrong curve, which the database does not allow us to identify.')
# ================= S6
head('S6. Hierarchical model and decision rules')
Hc = pd.read_csv(R + 'hierarchical_comparison.csv'); Hd = pd.read_csv(R + 'hierarchical_fit_diagnostics.csv'); lo = Hd[Hd.target != 'ALL PAPERS (no hold-out)']
para(f'The hierarchical model is a right-censored lognormal regression with a paper-level random intercept (Section 3.8), fitted without the target paper. All {len(lo)} leave-one-paper-out fits converged ({int(lo.converged.sum())} of {len(lo)}). Table S6a gives the comparison with the pooled model by alloy. Its lower breach rate came with lower bounds and, at equal mean bound life, it lies on the pooled frontier (Table S3).')
Ha = pd.read_csv(R + 'hierarchical_by_alloy.csv'); rows = []
for a_ in ('Ti-6Al-4V', 'AlSi10Mg', 'IN718', '316L'):
    p = Ha[(Ha.alloy == a_) & (Ha.model == 'pooled_censored')].iloc[0]; h = Ha[(Ha.alloy == a_) & (Ha.model == 'hierarchical_censored')].iloc[0]
    rows.append([a_, int(p.target_papers), ci(p.breach_pct, p.CI_low, p.CI_high), ci(h.breach_pct, h.CI_low, h.CI_high), f'{p.mean_lower_bound_log10N:.3f}', f'{h.mean_lower_bound_log10N:.3f}', f'{h.tau_between_mean:.2f}'])
for lab, mname in (('All alloys', 'hierarchical_censored'),):
    p = Hc[Hc.model == 'pooled_censored'].iloc[0]; h = Hc[Hc.model == mname].iloc[0]; rows.append([lab, int(p.target_papers), ci(p.breach_pct, p.CI_low, p.CI_high), ci(h.breach_pct, h.CI_low, h.CI_high), f'{p.mean_lower_bound_log10N:.3f}', f'{h.mean_lower_bound_log10N:.3f}', f'{h.tau_between_mean_over_fits:.2f}'])
cap('Table S6a. Pooled and hierarchical censored bounds by alloy. Breach rates in percent with 95% bootstrap intervals; mean lower-bound life in log₁₀ cycles; τ̂ is the mean between-paper SD over the fits. The intervals come from the bootstrap of the hierarchical-comparison script and differ slightly from those of Table 3 (Monte Carlo variation between resampling sequences).')
table(['Alloy', 'Target papers', 'Pooled breach (%)', 'Hierarchical breach (%)', 'Pooled bound', 'Hierarchical bound', 'Mean τ̂'], rows, widths=[2.3, 1.8, 3.2, 3.2, 1.9, 2.2, 1.5])
Dp = pd.read_csv(R + 'decision_policy_summary.csv'); rows = []
for _, r in Dp.iterrows():
    rows.append([r.rule.replace('accept single if <= pooled + ', 'Accept single-study if ≤ pooled + ').replace(' log10', ' log₁₀').replace('+ -', '− ').replace('single-study only', 'Single-study only').replace('pooled only', 'Pooled only'), ci(r.breach_pct, r.breach_lo, r.breach_hi), f'{r.mean_bound_log10N:.3f}', '—' if pd.isna(r.share_single_accepted_pct) else f'{r.share_single_accepted_pct:.1f}', '—' if pd.isna(r.breach_among_accepted_pct) else ci(r.breach_among_accepted_pct, r.accepted_lo, r.accepted_hi)])
cap('Table S6b. Margin rule: use the single-study bound when it lies no more than the stated margin above the pooled bound at the same stress, otherwise use the pooled bound. Breach rates in percent with 95% bootstrap intervals; mean bound life in log₁₀ cycles.')
table(['Rule', 'Breach (%) [95% CI]', 'Mean bound (log₁₀ N)', 'Single-study accepted (%)', 'Breach among accepted (%)'], rows, widths=[5.4, 3.0, 2.6, 2.6, 3.2], size=8)
De = pd.read_csv(R + 'decision_policy_equal_height.csv')
para(f"At the mean bound life of the pooled 5% bound (4.204 log₁₀ cycles), the margin rule, traced over its range of margins and read at that height, gave {De.iloc[0].estimate:.1f}% ({De.iloc[0].ci_low:.1f}–{De.iloc[0].ci_high:.1f}), against {De.iloc[1].estimate:.1f}% for the pooled bound at its own height ({De.iloc[1].ci_low:.1f}–{De.iloc[1].ci_high:.1f}). The rule therefore did not improve on pooling alone at equal conservatism.")
# ================= S7
head('S7. Reproducibility')
para('The package (analysis code, the commands below and all outputs in results/) is archived at [ZENODO DOI] and available at [GITHUB URL]. The only external input is the public FatigueData-AM2022 database (figshare DOI 10.6084/m9.figshare.22337629); scripts/fetch_data.py downloads it and checks its SHA-256 checksum. Every random step has a fixed seed. Point estimates are deterministic, and bootstrap intervals move by at most about 0.4 points between seeds. sh run_all.sh reproduces everything except the two optional slow steps (--with-digitization for Table S5b and --with-coverage-test for Figure S1); their archived outputs are in results/.')
X, Jx = 'data/FatigueData-AM2022.xlsx', 'data/FatigueData-AM2022.json'
mp = [('Table 1, Section 4.1, Figure 3', f'reproduce.py {X} {Jx} primary 11', 'primary.out'), ('Section 4.7, Table 5 (all papers; flagged curves; rated subsets)', 'reproduce.py … all | primary_noflag | primary_nopaperflag | primary_highrated | primary_lowrated 11 (slope > −1: final argument −1)', 'all.out, primary_*.out; Tables S5a'),
      ('Section 4.7, reading error; Table S5b', 'sh run_all.sh --with-digitization (digitization_sensitivity.sh both)', 'digitization/summary.csv'), ('Figure 4, Table 2, Sections 3.10 and 4.2; Tables S2, S3', 'coverage_conservatism.py X J 1 11; matched_contrasts.py cc_primary_strict1_records.csv 11', 'cc_primary_strict1_*.csv, matched_contrasts.csv'),
      ('Table 3, Section 4.3', 'alloy_results.py 11', 'alloy_results.csv'), ('Section 4.4 learning curves', 'learning_curves.py X J 51 (and 52)', 'learning_curves*.out'), ('Section 4.5 matched conditions', 'comparability.py X J', 'comparability.out'), ('Section 4.6, Table 4, Table S4', 'heterogeneity_fixed_sample.py X J; tau_explained.py X J', 'heterogeneity_fixed_sample*.csv, tau_explained.out'),
      ('Section 4.8, Tables 6–7', 'within_paper_control.py X J 100 5; audit_within_paper.py X J', 'within_paper_*.csv, audit_within_paper.out'), ('Section 4.9, Tables 8–9, Table S6a', 'hierarchical_comparison.py X J 11', 'hierarchical_*.csv'),
      ('Table S6b, Section 4.9', 'decision_policy.py cc_primary_strict1_records.csv cc_primary_strict1_frontier.csv 11', 'decision_policy_*.csv'), ('Figures S3–S5, Section 4.7', 'influence_checks.py cc_primary_strict1_records.csv X J', 'influence_*.csv'),
      ('Figure S1', 'run_all.sh --with-coverage-test (coverage_test.py)', 'cov_*.out'), ('Section 2, Section S5', 'data_screens.py X J', 'data_screens.out'), ('Figures S1–S5', 'make_supplement_figures.py', 'supplement/FigS*.png/.pdf'), ('This file', 'make_supplement_doc.py', 'supplement/Supplementary_Information.docx')]
cap('Table S7. Script-to-result map. X and J are the database files; all commands use python3 -I.')
table(['Manuscript item', 'Command', 'Output in results/'], [list(m) for m in mp], widths=[5.0, 7.8, 4.0], size=7.5)
para(f"Environment: Python {platform.python_version()}, NumPy {np.__version__}, SciPy {__import__('scipy').__version__}, pandas {pd.__version__}, matplotlib {matplotlib.__version__}, python-docx {docx.__version__ if hasattr(docx, '__version__') else ''}. The versions used for the reported results are pinned in requirements.txt.")
d.core_properties.author = 'Sakib Shahrear Avro'; d.core_properties.title = 'Supplementary Information'; d.core_properties.comments = ''
d.save(O + 'Supplementary_Information.docx'); print('saved', len(d.paragraphs))
