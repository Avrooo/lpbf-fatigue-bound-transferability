# Results

Outputs written by `run_all.sh` (shipped here so that the supplement can be rebuilt and the tests run without the database). All files are derived summaries: paper DOIs, curve identifiers, counts, bound heights and rates. They contain no raw specimen measurements.

| File(s) | Content |
|---|---|
| `primary.out`, `all.out`, `primary_*.out` | Headline rates for the primary set, all 147 papers and the data-exclusion sensitivity runs |
| `pairs_*.csv` | Pair-curve records of each run: alloy, stress ratio, source paper, source curve, target paper, scored specimens and breach counts per method |
| `alloy_results.csv` | Breach rates by alloy (main Table 3) |
| `cc_primary_strict{0,1}_*.csv`, `cc_strict*.out` | Coverage–conservatism records, 13-point frontier, contrasts and audits. `strict1` (covariates chosen from training papers only) is used in the paper; `strict0` reproduces the headline |
| `matched_contrasts.csv` | Paired matched-bound-life contrasts (main Table 2, Table S3) |
| `decision_policy_*.csv`, `hierarchical_*.csv`, `within_paper_*.csv`, `influence_*.csv`, `heterogeneity_fixed_sample*.csv` | Decision rule, hierarchical model, within-paper control, influence checks and between-paper variance |
| `learning_curves*.out`, `comparability.out`, `tau_explained.out`, `data_screens.out`, `build.out` | Learning curves, matched conditions, varying-sample variance, data screens, data counts |
| `digitization/` | Reading-error runs (20 runs, two modes) and `summary.csv` |
| `cov_*.out` | Bootstrap-coverage simulations behind Figure S1 (slow; archived) |
| `supplement/` | Figures S1–S5 and `Supplementary_Information.docx` |
| `run_manifest.json` | Versions, input checksums and completion status of the last `run_all.sh` |

Intervals are 95% percentile intervals from the independent-role paper bootstrap (2,000 resamples) unless stated. Point estimates are deterministic.
