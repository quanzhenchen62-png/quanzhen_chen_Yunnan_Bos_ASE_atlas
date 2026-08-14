#!/bin/bash

set -euo pipefail

# Step 1: build individual-level D-statistic tests for reference-panel screening
python3 01_build_reference_dstat_tests.py \
  --candidate-table data/reference_candidates.tsv \
  --donor-file taurine=data/taurine_reference_entities.txt \
  --donor-file yak=data/yak_reference_entities.txt \
  --donor-file gayal=data/gayal_reference_entities.txt \
  --outgroup-id swamp_buffalo \
  --output results/reference_panel_dstat_tests.tsv

# Step 2: run qpDstat for reference-panel screening
python3 03_run_qpdstat.py \
  --qpdstat-bin tools/AdmixTools/bin/qpDstat \
  --geno data/all_samples.geno \
  --snp data/all_samples.snp \
  --ind data/all_samples.ind \
  --tests results/reference_panel_dstat_tests.tsv \
  --out-prefix results/reference_panel_dstat

# Step 3: retain putatively unadmixed reference individuals
python3 04_filter_reference_panel.py \
  --results results/reference_panel_dstat.results.tsv \
  --z-threshold 3 \
  --audit-output results/reference_panel_retention_audit.tsv \
  --retained-output results/reference_panel_retained.tsv

# Step 4: build target-breed D-statistic tests
python3 02_build_target_dstat_tests.py \
  --source-id retained_SAI_panel \
  --target-file data/target_breeds.txt \
  --donor-file gayal=data/gayal_reference_panel.txt \
  --donor-file yak=data/yak_reference_panel.txt \
  --outgroup-id swamp_buffalo \
  --output results/target_breed_dstat_tests.tsv

# Step 5: run qpDstat for target breeds
python3 03_run_qpdstat.py \
  --qpdstat-bin tools/AdmixTools/bin/qpDstat \
  --geno data/all_samples.geno \
  --snp data/all_samples.snp \
  --ind data/all_samples.ind \
  --tests results/target_breed_dstat_tests.tsv \
  --out-prefix results/target_breed_dstat

# Step 6: calculate U20 in non-overlapping 50-kb windows
python3 05_scan_u20_windows.py \
  --vcf data/filtered_biallelic_snps.vcf.gz \
  --sample-table data/sample_population_table.tsv \
  --source-pop SAI \
  --target-pop target_breed \
  --donor-pop gayal \
  --window-size 50000 \
  --source-max-af 0.01 \
  --target-min-af 0.20 \
  --donor-fixed-af 1.00 \
  --site-output results/u20_sites.tsv \
  --window-output results/u20_windows.tsv \
  --top-output results/u20_top1pct_windows.tsv
