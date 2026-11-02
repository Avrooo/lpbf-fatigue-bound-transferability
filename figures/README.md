# Main-text figures

| Figure | File | How it is produced |
|---|---|---|
| Figure 3 (single-study vs pooled breach rates) | `Fig3_breach_rates.{png,pdf}` | **Regenerated** by `python3 figures/make_fig3_fig5_fig6.py` from `results/primary.out` |
| Figure 4 (coverage-conservatism frontier) | `Fig4_coverage_conservatism.{png,pdf}` | **Regenerated** by `python3 figures/make_fig4_frontier.py` from `results/cc_primary_strict1_frontier.csv` and `results/decision_policy_summary.csv` |
| Figure 5 (training-paper availability) | `Fig5_learning_curves.{png,pdf}` | **Regenerated** by `make_fig3_fig5_fig6.py` from `results/learning_curves.out` |
| Figure 6 (matched test conditions) | `Fig6_matched_conditions.{png,pdf}` | **Regenerated** by `make_fig3_fig5_fig6.py` from `results/comparability.out` |
| Figures 1 and 2 (study design, bound schematics) | not in this repository | Schematics with no data values; supplied with the manuscript |

`run_all.sh` runs both scripts. The numbers printed on Figures 3, 5 and 6 are parsed from the saved outputs, not typed. The versions in the submitted manuscript were drawn earlier with a different font and layout (Figure 3 is kept for reference in `static/`); the plotted values are the same. Supplementary Figures S1-S5 are regenerated into `results/supplement/` by `make_supplement_figures.py`.
