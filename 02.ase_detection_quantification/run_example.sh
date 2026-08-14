#!/bin/bash

set -euo pipefail

sample_id="sample_001"

bash 01_run_star_wasp.sh \
  "${sample_id}" \
  data/${sample_id}_R1.fastq.gz \
  data/${sample_id}_R2.fastq.gz \
  data/${sample_id}.het.vcf.gz \
  refs/star_index \
  refs/annotation.gtf \
  12 \
  results/star_wasp/${sample_id} \
  results/star_wasp/${sample_id}.stats.tsv

bash 02_run_ase_read_counter.sh \
  "${sample_id}" \
  results/star_wasp/${sample_id}/${sample_id}.WASP_filtered.bam \
  data/${sample_id}.het.vcf.gz \
  refs/reference.fa \
  tools/GenomeAnalysisTK.jar \
  results/ase_read_counter/${sample_id}.ASEReadCounter.tsv \
  results/tmp_rg

python3 03_filter_test_ase_sites.py \
  --input results/ase_read_counter/${sample_id}.ASEReadCounter.tsv \
  --filtered-output results/ase_sites/${sample_id}.filtered.tsv \
  --significant-output results/ase_sites/${sample_id}.significant.tsv

python3 04_call_brd_ase.py \
  --metadata data/brd_ase_metadata.tsv \
  --output-prefix results/brd_ase/brd_ase

python3 05_call_pop_ase.py \
  --metadata data/pop_ase_metadata.tsv \
  --output-prefix results/pop_ase/pop_ase
