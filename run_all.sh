#!/bin/sh
# Full reproduction of the manuscript and supplement results. Run from the repository root, with the dependencies of requirements.txt installed.
#   sh run_all.sh                          main pipeline (about 12 minutes on two CPU cores)
#   sh run_all.sh --with-digitization      also the 20 reading-error runs of Table S5b (about 16 more minutes)
#   sh run_all.sh --with-coverage-test     also the bootstrap-coverage simulations of Figure S1 (hours)
# The pipeline stops at the first failing step; results/run_manifest.json then says "failed". Archived outputs of the two optional
# steps are shipped in results/, so the supplement can be rebuilt without them.
set -eu
mkdir -p results
trap 'python3 -I scripts/write_manifest.py failed "$@" || true' EXIT
[ -f data/FatigueData-AM2022.xlsx ] && [ -f data/FatigueData-AM2022.json ] || python3 -I scripts/fetch_data.py data
python3 -I scripts/check_data.py data
X=data/FatigueData-AM2022.xlsx; J=data/FatigueData-AM2022.json
python3 -I build_intermediates.py $X $J > results/build.out
python3 -I data_screens.py $X $J > results/data_screens.out
for w in primary all primary_noflag primary_nopaperflag primary_highrated primary_lowrated; do python3 -I reproduce.py $X $J $w 11 > results/$w.out; done
python3 -I reproduce.py $X $J primary_noflag 11 -1 > results/primary_noflag_slope_gt_minus1.out
python3 -I reproduce.py $X $J primary_nopaperflag 11 -1 > results/primary_nopaperflag_slope_gt_minus1.out
python3 -I learning_curves.py $X $J 51 > results/learning_curves.out
python3 -I learning_curves.py $X $J 52 > results/learning_curves_seed52.out
python3 -I within_paper_control.py $X $J 100 5 > results/within_paper.out
python3 -I audit_within_paper.py $X $J > results/audit_within_paper.out
python3 -I hierarchical_comparison.py $X $J 11 > results/hierarchical.out
python3 -I heterogeneity_fixed_sample.py $X $J > results/heterogeneity_fixed_sample.out
python3 -I heterogeneity_fixed_sample.py $X $J 'AlSi10Mg|0.1' > results/heterogeneity_fixed_sample_excl_AlSi10Mg_R0.1.out
python3 -I coverage_conservatism.py $X $J 0 11 > results/cc_strict0.out     # STRICT=0 reproduces the headline rates
python3 -I coverage_conservatism.py $X $J 1 11 > results/cc_strict1.out     # STRICT=1: covariates chosen from training papers only (used in the paper)
python3 -I decision_policy.py cc_primary_strict1_records.csv cc_primary_strict1_frontier.csv 11 > results/decision_policy.out
python3 -W ignore -I influence_checks.py cc_primary_strict1_records.csv $X $J > results/influence_checks.out
python3 -I matched_contrasts.py cc_primary_strict1_records.csv 11 > results/matched_contrasts.out
python3 -I comparability.py $X $J > results/comparability.out
python3 -I tau_explained.py $X $J > results/tau_explained.out
mv -f matched_contrasts.csv cc_*.csv decision_policy_*.csv influence_*.csv results/
mv -f within_paper_*.csv hierarchical_*.csv heterogeneity_fixed_sample*.csv pairs_*.csv results/
python3 -I alloy_results.py 11 > results/alloy_results.out
for a in "$@"; do
  [ "$a" = "--with-digitization" ] && sh digitization_sensitivity.sh both
  if [ "$a" = "--with-coverage-test" ]; then   # slow
    python3 -I coverage_test.py $X $J "Ti-6Al-4V|-1" 20 300 500 1 > results/cov_Ti_m1_K20.out
    python3 -I coverage_test.py $X $J "Ti-6Al-4V|-1" 60 300 500 2 > results/cov_Ti_m1_K60.out
    python3 -I coverage_test.py $X $J "Ti-6Al-4V|-1" 120 200 500 5 > results/cov_Ti_m1_K120.out
    python3 -I coverage_test.py $X $J "Ti-6Al-4V|0.1" 30 300 500 3 > results/cov_Ti_01_K30.out
    python3 -I coverage_test.py $X $J "316L|-1" 22 500 500 4 > results/cov_316L_m1_K22.out
  fi
done
python3 -I figures/make_fig4_frontier.py
python3 -I figures/make_fig3_fig5_fig6.py
python3 -I make_supplement_figures.py
python3 -I make_supplement_doc.py analysis_table_primary.csv
trap - EXIT
python3 -I scripts/write_manifest.py complete "$@"
