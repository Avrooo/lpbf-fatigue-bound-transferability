# Reproducibility

## Requirements
* Linux or macOS with `sh`, Python 3.13 (tested with 3.13.16). On Windows use Git Bash or WSL for the shell scripts.
* `python -m pip install -r requirements-lock.txt -r requirements-dev.txt`. `requirements.txt` lists the direct dependencies (numpy, scipy, pandas, openpyxl, matplotlib, python-docx).
* Internet access only for `scripts/fetch_data.py`; the analysis itself is offline. If the download fails, place `FatigueData-AM2022.xlsx` and `.json` from https://doi.org/10.6084/m9.figshare.22337629 in `data/`.
* About 2 CPU cores and a few GB of RAM were available in the verification run; peak memory was not measured.

## Commands
```
python scripts/fetch_data.py data           # download the pinned input files
python scripts/check_data.py data           # compare with data/SHA256SUMS (run_all.sh does this too)
sh run_all.sh                               # main pipeline
sh run_all.sh --with-digitization           # adds the 20 reading-error runs (Table S5b)
sh run_all.sh --with-coverage-test          # adds the bootstrap-coverage simulations (Figure S1)
python -m pytest tests/ -v                  # or python tests/run_tests.py when pytest is not installed
```
All random steps have fixed seeds on the command lines in `run_all.sh`. Point estimates are deterministic. Re-running overwrites `results/`; copy it first to compare. `run_all.sh` stops at the first failing step and records the outcome in `results/run_manifest.json`.

## Expected outputs
Headline check (`results/primary.out`): 123 papers, 424 curves, 3,924 specimens; single-study 5% prediction bound 27.21%, pooled censored 6.50%, paired difference +20.71 points; 1,562 pairs. `results/matched_contrasts.csv`: at the pooled 5% height (4.20 log₁₀ cycles) single-study 15.57%, pooled 6.50%, difference 9.07 points (95% CI 6.45–11.77). `results/within_paper_summary.csv`: within-paper 4.12%, cross-paper 26.87%.

## Measured run times (Linux container, 2 CPU cores, Python 3.13.16)
| Step | Wall time |
|---|---|
| `sh run_all.sh` (main pipeline, including supplement figures and document) | about 12 minutes |
| `sh digitization_sensitivity.sh both` (20 runs) | about 16 minutes |
| `--with-coverage-test` | not re-run in the verification below; run times were not recorded (the original run took hours) |

## Verification record
* **What was run.** A fresh copy of the repository with `results/` emptied (except the archived outputs of the two optional steps) and the two database files placed in `data/`; then `sh run_all.sh`, then `sh digitization_sensitivity.sh both` after removing the archived reading-error outputs, then the tests.
* **Outcome.** `run_all.sh` completed (`results/run_manifest.json`: status "complete"; inputs matched `data/SHA256SUMS`). Of the 86 regenerated text and CSV output files compared with the shipped ones, 83 were identical (numeric CSV values to 1e-6) and the other 3 (`cc_strict0.out`, `cc_strict1.out`, `hierarchical.out`) differ only in their "done N s" timing lines. The 20 regenerated reading-error runs gave a `summary.csv` identical to the shipped one. The tests: 27 passed, 0 skipped, 0 failed on the regenerated outputs.
* **Tests without the database.** In a fresh checkout (no `analysis_table_primary.csv`) 24 tests pass and 3 are skipped; the 3 need the derived table, which `run_all.sh` writes.
* **Not verified.** The bootstrap-coverage simulations (`--with-coverage-test`) were not re-run; `results/cov_*.out` are the archived outputs of the original run. `fetch_data.py` could not be tested from the verification environment (no access to figshare from there); the files were supplied by hand and their checksums checked. The tests were run with `tests/run_tests.py` because pytest could not be installed in that environment; they were written as plain pytest-style functions. No run was made on Windows or macOS, and no run used a second, independent machine.
* **Environment limits.** The lock file lists the packages that were present in the verification environment, not a fresh resolution from PyPI.
