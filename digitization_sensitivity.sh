#!/bin/sh
# Digitization-error sensitivity. Needs the data in ./data. About 16 minutes on two CPU cores (20 runs).
# Every stress and life value is multiplied by random reading error (worst case: all points treated as digitized).
#   moderate: stress sd 1%, log10(life) sd 0.02 (about 5%);  severe: stress sd 3%, log10(life) sd 0.05 (about 12%); seeds 1-5.
# dynamic: eligibility and stress-range masks are recomputed from the perturbed values (comparison sets can change).
# fixed:   same source curves and the same target specimens as the unperturbed analysis; only the measured values change.
# Usage: sh digitization_sensitivity.sh [dynamic|fixed|both]   (default both)
set -e
MODE=${1:-both}; X=data/FatigueData-AM2022.xlsx; J=data/FatigueData-AM2022.json
for m in dynamic fixed; do
  [ "$MODE" = both ] || [ "$MODE" = $m ] || continue
  mkdir -p results/digitization/$m
  for sd in 1 2 3 4 5; do
    for sc in "0.01,0.02" "0.03,0.05"; do
      tag=$(echo $sc | tr ',' '_')_$sd
      if [ $m = fixed ]; then DIGIT_FIXED_ELIGIBILITY=1 DIGIT_NOISE="$sc,$sd" python3 -I reproduce.py $X $J primary 11 > results/digitization/$m/out_$tag.txt
      else DIGIT_NOISE="$sc,$sd" python3 -I reproduce.py $X $J primary 11 > results/digitization/$m/out_$tag.txt; fi
    done
  done
done
rm -f pairs_primary_digit_*.csv
python3 -I summarize_digitization.py
