# Known limitations

* **Not allowables.** Confirmed-breach rates describe early failures across published studies. They are not component failure probabilities, and no pooled tolerance bound or qualified design allowable is estimated.
* **Runouts.** A runout is never a breach; the breach rate is therefore not an exact failure probability. Eleven runouts stopped below a single-study bound are ambiguous and not counted (effect on the specimen-weighted rate: 19.006% → 19.025%).
* **Dependence between papers.** Source and target papers are different publications, but they can share authors, laboratories, machines or test campaigns. Results are cross-publication evaluations, not independent between-laboratory validation.
* **Secondary data.** Many values were digitized from figures and the database does not mark which. Random reading error was simulated (Table S5b); systematic extraction errors cannot be excluded.
* **Post hoc life window.** 10³–2 × 10⁷ cycles was chosen after looking at the data; all 147 papers are reported as a sensitivity set.
* **Scope.** Four alloys, uniaxial loading, R = −1 and 0.1, one database, a lognormal S–N model. No variable-amplitude, multiaxial or very-high-cycle data.
* **Intervals** come from resampling papers with fitted models fixed; they exclude refitting uncertainty, digitization error and dependence between related papers. The simulation check covers only part of the methods (see `STATISTICAL_VALIDATION.md`).
* **Frontier.** Mean bound life can hide differences in shape across stress. The single-study frontier is sparse below about 3.1 log₁₀ cycles, so its 5% crossing (2.86, 95% CI 2.38–3.40) is indicative only.
* **Hierarchical bound** is a plug-in value without parameter uncertainty.
* **Machine-specific details.** AM machine identity is closely tied to paper identity, so its step in the variance table is not interpretable alone. The variance comparisons are descriptive.
* **Figure versions.** Figures 3, 5 and 6 are redrawn from saved outputs, but the versions printed in the submitted manuscript were drawn earlier with a different font and layout; the plotted values are the same. Figures 1 and 2 are schematics and are not in the repository.
