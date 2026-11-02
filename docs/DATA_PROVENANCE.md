# Data provenance

## Source
FatigueData-AM2022 (Zhang and Xu), figshare DOI 10.6084/m9.figshare.22337629. The S–N part holds 15,146 observations in 1,610 datasets from 459 articles; many values were read from published figures, and the database does not say which. The analysed files are pinned by SHA-256 in `data/SHA256SUMS`. The files are not redistributed here (their licence has not been verified for this repository). `scripts/fetch_data.py` prints the licence reported by the figshare API when it downloads them.

## Filtering (`build_intermediates.py`, `reproduce.py`)
| Step | Papers | Curves | Specimens | Failures | Runouts |
|---|---|---|---|---|---|
| L-PBF, uniaxial, Ti-6Al-4V / AlSi10Mg / IN718 / 316L, R = −1 or 0.1 (the "all" set) | 147 | 492 | 4,609 | 4,083 | 526 |
| Primary set: papers with any selected observation outside 10³–2 × 10⁷ cycles are removed as a whole | 123 | 424 | 3,924 | 3,464 | 460 |
| Eligible source curves (≥ 6 failures at ≥ 3 stress levels) | 99 | 296 | | | |

The life window was chosen after looking at the data. All 147 papers are analysed as a sensitivity set (`results/all.out`).

## Pairs
A source curve is compared only with target papers of the same alloy and stress ratio. The primary analysis has **1,562 distinct source–target paper pairs**; a pair compared under both stress ratios appears in two alloy/stress-ratio groups, so the group-level count is **1,564**. 295 of the 296 eligible source curves appear in at least one pair.

## Derived data
`build_intermediates.py` writes `analysis_table_{primary,all}.csv` (one row per specimen: curve id, paper DOI, stress, life, runout flag, alloy, process, test type, stress ratio, documentation rating, surface, HIP, heat-treatment class). These are subsets of the source database and are not shipped (`.gitignore`). The files in `results/` are summaries (counts, rates, bound heights, paper DOIs).

## Data-quality screens (`data_screens.py`)
Over the whole database: no non-positive or missing values; 38 duplicate points within curves; 23 rows in 12 curves with stress above the recorded ultimate strength; 24 curves with a slope above −1 among curves with at least six failures. The paper flags curves by slope (≥ 0 by default; > −1 as a sensitivity) and removes duplicates; see Table 5 and Supplement S5.
