#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "$0")" && pwd)

bash "$ROOT/01_windowed_fst.sh" "$ROOT/windowed_fst_jobs.tsv"
bash "$ROOT/02_windowed_pi.sh" "$ROOT/windowed_pi_jobs.tsv"
python3 "$ROOT/03_calc_pi_ratio.py" \
  --pairs "$ROOT/pi_ratio_pairs.tsv" \
  --out-dir results/pi_ratio
bash "$ROOT/04_run_selscan_xpehh.sh" "$ROOT/xpehh_jobs.tsv"
python3 "$ROOT/05_summarize_xpehh_windows.py" \
  --job-table "$ROOT/xpehh_window_jobs.tsv"
python3 "$ROOT/06_make_candidate_regions.py" \
  --match "$ROOT/selective_signals_file_match.tsv" \
  --outdir results/selective_signal
python3 "$ROOT/08_make_fst_deciles.py" \
  --job-table "$ROOT/windowed_fst_jobs.tsv" \
  --outdir results/fst_deciles
python3 "$ROOT/07_make_candidate_supplement.py" \
  --input-dir results/selective_signal/04.candidate_merged \
  --output results/selective_signal/04.candidate_merged/Supp.Table.S6.Candidate_selective_sweep_regions_identified_based_on_the_reference_population.tsv
