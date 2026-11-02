# Cross-study transferability of lower-tail S–N bounds for laser powder bed fusion alloys

Reproducibility repository for the manuscript *Cross-study transferability of lower-tail S–N bounds for laser powder bed fusion alloys: coverage and conservatism*, prepared for submission to *Archives of Civil and Mechanical Engineering*.

## Overview
This repository holds the analysis code and the saved results of a cross-publication assessment of lower-tail stress–life (S–N) bounds in laser powder bed fusion (L-PBF) alloys. It asks whether a bound fitted to one published study keeps its lower-tail coverage on independent publications of the same alloy and stress ratio, and whether pooling evidence from several publications improves external coverage beyond simply producing more conservative bounds.

## Data
The only external input is the public **FatigueData-AM2022** database (https://doi.org/10.6084/m9.figshare.22337629). It is not redistributed here; `scripts/fetch_data.py` downloads it and `data/SHA256SUMS` pins the analysed files. Primary analysis set:

* 123 publications, 424 S–N curves, 3,924 specimens (3,464 failures, 460 right-censored runouts)
* 296 eligible source curves from 99 publications
* 1,562 distinct source–target paper pairs (1,564 alloy/stress-ratio group occurrences)
* four alloys (Ti-6Al-4V, AlSi10Mg, IN718, 316L), uniaxial loading, R = −1 and R = 0.1

The full filtered database has 147 publications and 4,609 specimens; see `docs/DATA_PROVENANCE.md`.

## Main results
* Within their source publications, nominal 5% prediction bounds had a confirmed-breach rate of 4.12%. Under cross-publication evaluation of the same eligible source curves the rate rose to 26.87%, an empirical transfer gap of 22.75 percentage points. Runouts stay in the denominator and are never counted as breaches.
* Primary analysis: 27.2% (single-study) against 6.5% (pooled censored).
* At a matched mean lower-bound life of 4.20 log₁₀ cycles the confirmed-breach rates are 15.6% (single-study) and 6.5% (pooled), a paired difference of 9.1 percentage points (95% CI 6.5–11.8). The pooled advantage is therefore not only a result of more conservative bounds. A hierarchical model lies on the pooled frontier.
* The bounds are not qualified design allowables (`docs/KNOWN_LIMITATIONS.md`).

## Quick start
```
git clone https://github.com/USERNAME/lpbf-fatigue-bound-transferability.git
cd lpbf-fatigue-bound-transferability
python3 -m venv .venv && source .venv/bin/activate        # Windows: use Git Bash or WSL for the shell scripts
python -m pip install -r requirements-lock.txt -r requirements-dev.txt      # pins the tested environment; requirements.txt lists the direct dependencies
python scripts/fetch_data.py data                          # or place the two files in data/ by hand
sh run_all.sh                                              # main pipeline
python -m pytest tests/ -v                                 # or: python tests/run_tests.py
```
Optional slow steps: `sh run_all.sh --with-digitization` (reading-error runs, Table S5b) and `sh run_all.sh --with-coverage-test` (bootstrap-coverage simulations, Figure S1). Their archived outputs are shipped in `results/`.

`run_all.sh` stops at the first failing step and writes `results/run_manifest.json` (versions, input checksums, git commit, completion status). It overwrites the files in `results/`; commit or copy them first if you want to compare with the shipped outputs. Details, timings and expected outputs: `docs/REPRODUCIBILITY.md`. Which command makes which table or figure: `docs/FIGURE_TABLE_MAP.md`.

## Repository layout
| Path | Content |
|---|---|
| `reproduce.py` and the other `*.py` scripts at the root | Analysis (see `docs/FIGURE_TABLE_MAP.md`); they are run from the root and read `reproduce.py` by relative path |
| `scripts/` | Data download, checksum check, run manifest |
| `results/` | Saved outputs, supplement figures and `Supplementary_Information.docx` |
| `figures/` | Main Figures 3-6 regenerated from saved outputs (Figures 1-2 are schematics; see `figures/README.md`) |
| `docs/` | Provenance, methods, statistical validation, limitations, reproducibility, figure/table map |
| `tests/` | Automated checks of counts, source–target separation, breach definition, weighting, frontier and reported numbers |
| `paper/` | Reading copies of the manuscript and supplement (added when final) |

## Verification status
See `docs/REPRODUCIBILITY.md` (section "Verification record"). The record states exactly what was run, where and with what outcome.

## Limitations
Confirmed-breach rates are not component-level failure probabilities. Source and target papers are separate publications, but several may come from related laboratories or test campaigns. Results apply to the evaluated alloys, stress ratios, database and statistical assumptions. The bootstrap-coverage simulation covers only part of the methods (`docs/STATISTICAL_VALIDATION.md`).

## Citation
Cite the article and the archived software release. See `CITATION.cff`. The publication reference and the Zenodo DOI will be added when available.

## License
Code and derived result files: MIT (`LICENSE`). The FatigueData-AM2022 database keeps its own terms.
