# Methods (implementation notes)

Full descriptions are in Section 3 of the manuscript and Section S3 of the supplement. This file records the choices that matter for reproducing the code.

**S–N model.** Within an alloy and stress-ratio group, y = log₁₀ N = β₀ + β₁ x + γᵀz + ε, x = log₁₀(S / S_ref), ε ~ N(0, σ²). z holds recorded covariates (surface condition, HIP, heat treatment ≥ 700 °C) in the pooled and hierarchical models.

**Single-study bounds** (`reproduce.py`). One eligible curve is fitted by least squares to its failures. The 5% prediction bound is the one-sided 95% prediction limit; the approximate 95/95 tolerance bound uses a non-central t factor with n_eff = 1/h. Both are used only inside the stress range of the source curve.

**Pooled censored bound.** Leave-one-paper-out fit per alloy and stress-ratio group, runouts right-censored (EM). Bound at probability p: μ̂ + z_p σ̂. A failures-only variant is reported for comparison.

**Hierarchical censored bound** (`hierarchical_comparison.py`, `coverage_conservatism.py`). Paper-level random intercept, 21-node Gauss–Hermite quadrature, two starting values for τ. Predictive SD for an unseen paper: (σ̂² + τ̂²)^½. All 128 leave-one-paper-out fits converged.

**Confirmed breach.** A specimen that failed with a life below the bound at its applied stress. Runouts stay in the denominator and are never breaches; a runout stopped below the bound is "ambiguous" (11 of 57,483 scored specimens, single-study only) and is not counted.

**Weighting.** Breach rates are averaged by curve within a source paper, then over source papers within a target paper, then over target papers (both stress ratios averaged first when a pair appears in both). Alternatives (source-paper, pair, specimen weights) are in `results/cc_primary_strict1_audit_weighting.csv`.

**Bootstrap.** Two-way cluster bootstrap: source and target papers resampled independently, 2,000 resamples, percentile intervals, models not refitted. Differences between methods are taken within resamples.

**Coverage–conservatism** (`coverage_conservatism.py`, `matched_contrasts.py`). The bound probability p runs over the grid 10⁻⁵, 10⁻⁴, 3·10⁻⁴, 10⁻³, 2.5·10⁻³, 5·10⁻³, 0.01, 0.025, 0.05, 0.1, 0.2, 0.3, 0.5. For each p and method the weighted breach rate and the weighted mean of log₁₀ N of the bound over the scored specimens are recorded. Methods are compared at the mean bound life of the pooled bound at p = 1%, 2.5%, 5%, 10%, 20% by linear interpolation inside each frontier's sampled range; the reference height is recomputed in every bootstrap resample so the interval of a difference is paired. STRICT=1 chooses covariate columns from training papers only.

**Within-paper control** (`within_paper_control.py`). 213 curves from 89 papers, 100 random splits, 30% of failures held out; the bound is fitted on the remaining 70%.

**Stress equivalent.** For a bound lower by ΔL log₁₀ cycles, the stress at the same life is lower by the factor 10^(ΔL/b), b the S–N slope (mean pooled slope −2.79). An illustration, not a design factor.
