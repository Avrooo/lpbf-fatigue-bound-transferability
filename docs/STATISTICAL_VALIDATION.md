# Statistical validation

## Bootstrap-interval calibration (`coverage_test.py`, Figure S1)
A data-generating process is calibrated to one real alloy and stress-ratio group (curve slopes, between-paper, between-curve and within-curve scatter) using the real stress designs. In each replication K papers are drawn, new lives are simulated, the same estimator and weights as in the main analysis are applied, and the 95% percentile bootstrap interval is compared with the truth at the same K. Five settings (`results/cov_*.out`; configuration on the command lines in `run_all.sh`, seeds 1–5): 300 simulated data sets each (200 for K = 120, 500 for 316L), 500 bootstrap resamples per interval, and the truth at each K estimated from at least 1,500 independent replications. The coverage figures below therefore carry Monte Carlo error of about ±0.3 to ±1.2 points (see the `MC s.e.` column of the `.out` files).

| Setting | Single-study: independent-role / same-paper-set coverage | Pooled failures-only: independent-role / same-paper-set coverage |
|---|---|---|
| Ti-6Al-4V, R = −1, K = 20 | 99.7% / 94.7% | 99.0% / 98.3% |
| Ti-6Al-4V, R = −1, K = 60 | 100.0% / 93.7% | 99.7% / 99.3% |
| Ti-6Al-4V, R = −1, K = 120 | 100.0% / 92.0% | 99.5% / 99.5% |
| Ti-6Al-4V, R = 0.1, K = 30 | 99.7% / 93.7% | 95.7% / 95.3% |
| 316L, R = −1, K = 22 | 99.2% / 93.2% | 94.8% / 93.8% |

**What the simulation covers:** the 5% single-study prediction bound and the pooled failures-only 5% bound, one group at a time, no covariates, no runouts, Gaussian errors.
**What it does not cover:** the 95/95 tolerance bound, the censored pooled model, the hierarchical model, the coverage–conservatism frontier and the matched-bound contrasts. Their intervals use the same bootstrap but were not tested by simulation.

## Alternative resampling (`results/cc_primary_strict1_audit_bootstrap.csv`, Figure S2)
Independent roles (primary), same paper set for both roles, source papers only, target papers only. The independent-role interval is the widest for single-study rates and is used throughout.

## Other checks
Weighting alternatives, leave-one-source/target-paper-out (`influence_jackknife.csv`), stress-ratio and source-size splits (`influence_by_size_and_R.csv`), ambiguous-runout audit, strict-versus-loose covariate selection (identical rates), flagged-curve exclusions, reading-error runs in two eligibility modes (`results/digitization/`), within-paper audit (`audit_within_paper.out`).

## Automated tests (`tests/`)
They check counts, source–target separation, breach-definition invariants of the saved records, the weighting hierarchy (re-derived from `pairs_primary.csv`), the 13-point frontier and the matched-height contrast (re-derived from the records), and the published numbers frozen in `tests/expected_results.json`. They test the saved outputs and the logic that can be separated from the scripts; the scripts themselves are monolithic, so there are no unit tests of their internals.
