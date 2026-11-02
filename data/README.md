# Input data

The only external scientific input is **FatigueData-AM2022** (Zhang and Xu), a public database of fatigue results from published papers on additively manufactured metals.

* Source: https://doi.org/10.6084/m9.figshare.22337629
* Files used: `FatigueData-AM2022.xlsx` and `FatigueData-AM2022.json`
* Checksums of the analysed files: `SHA256SUMS` (computed from the files used for the reported results)

The files are **not included** in this repository, because their redistribution terms have not been confirmed here. Get them with

```
python3 scripts/fetch_data.py data
```

or download them by hand from the DOI above into this folder. `python3 scripts/check_data.py data` compares them with `SHA256SUMS`, and `run_all.sh` stops if they differ. A different figshare version can give different counts and results, so the pinned version is not replaced silently. To inspect another version on purpose, set `ALLOW_DATA_MISMATCH=1` and do not compare the output with the paper.

See `docs/DATA_PROVENANCE.md` for the filtering from the full database to the analysis set.
