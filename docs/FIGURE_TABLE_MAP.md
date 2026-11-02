# Figure and table map (manuscript numbering as submitted to ACME)

`X` = `data/FatigueData-AM2022.xlsx`, `J` = `data/FatigueData-AM2022.json`; commands run from the repository root with `python3 -I`. `run_all.sh` runs everything except the two optional steps.

| Item | Content | Command | Output |
|---|---|---|---|
| Table 1, Section 2 | Primary data set | `build_intermediates.py X J`; `reproduce.py X J primary 11` | `build.out`, `primary.out` |
| Figure 3, Section 4.1 | Single-study versus pooled breach rates | `reproduce.py X J primary 11` | `primary.out`, `figures/Fig3_breach_rates.*` |
| Figure 4, Table 2, Sections 3.10, 4.2 | Frontier; matched-bound contrasts | `coverage_conservatism.py X J 1 11`; `matched_contrasts.py cc_primary_strict1_records.csv 11`; `figures/make_fig4_frontier.py` | `cc_primary_strict1_*.csv`, `matched_contrasts.csv`, `figures/Fig4_*` |
| Table 3, Section 4.3 | Breach rates by alloy | `alloy_results.py 11` | `alloy_results.csv` |
| Figure 5, Section 4.4 | Training-paper availability | `learning_curves.py X J 51` (and 52) | `learning_curves*.out`, `learning_curve_units_seed51.csv` |
| Figure 6, Section 4.5 | Matched test conditions | `comparability.py X J` | `comparability.out` |
| Table 4, Section 4.6 | Fixed-sample between-paper variation | `heterogeneity_fixed_sample.py X J` (and `'AlSi10Mg|0.1'`) | `heterogeneity_fixed_sample*.csv` |
| Table 5, Section 4.7 | Data exclusions, rated subsets, reading error | `reproduce.py X J all\|primary_noflag\|primary_nopaperflag\|primary_highrated\|primary_lowrated 11` (`-1` for slope > −1); `digitization_sensitivity.sh both` | `all.out`, `primary_*.out`, `digitization/summary.csv` |
| Tables 6–7, Section 4.8 | Within-paper control | `within_paper_control.py X J 100 5`; `audit_within_paper.py X J` | `within_paper_*.csv`, `audit_within_paper.out` |
| Tables 8–9, Section 4.9 | Hierarchical model | `hierarchical_comparison.py X J 11` | `hierarchical_*.csv` |
| Table S1 | Selection counts | `build_intermediates.py`, `reproduce.py` | `build.out`, `pairs_primary.csv` |
| Figure S1 | Bootstrap calibration | `run_all.sh --with-coverage-test` (`coverage_test.py`), `make_supplement_figures.py` | `cov_*.out`, `supplement/FigS1_*` |
| Figure S2 | Bootstrap alternatives | `coverage_conservatism.py X J 1 11`, `make_supplement_figures.py` | `cc_primary_strict1_audit_bootstrap.csv` |
| Figures S3–S5 | Leave-one-paper-out; weighting; stress ratio and source size | `influence_checks.py cc_primary_strict1_records.csv X J`, `make_supplement_figures.py` | `influence_*.csv`, `cc_*audit_weighting.csv` |
| Tables S2, S3 | Full frontier; matched contrasts | as Figure 4 | `cc_primary_strict1_frontier.csv`, `matched_contrasts.csv` |
| Table S4 | Varying-sample variance | `tau_explained.py X J` | `tau_explained.out` |
| Table S5a, S5b | Exclusions; reading error | as Table 5 | as Table 5 |
| Table S6a, S6b | Hierarchical by alloy; decision rule | `hierarchical_comparison.py`; `decision_policy.py cc_primary_strict1_records.csv cc_primary_strict1_frontier.csv 11` | `hierarchical_by_alloy.csv`, `decision_policy_*.csv` |
| Section S5, data screens | Database-wide screens | `data_screens.py X J` | `data_screens.out` |
| Supplementary Information | All S figures and tables | `make_supplement_figures.py`; `make_supplement_doc.py analysis_table_primary.csv` | `supplement/` |

Figures 3, 4, 5 and 6 are redrawn from saved outputs by `figures/make_fig4_frontier.py` and `figures/make_fig3_fig5_fig6.py`. Figures 1 and 2 are schematics without data values and are not included (`figures/README.md`).
