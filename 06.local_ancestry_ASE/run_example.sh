#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RESULTS="results/example"

mkdir -p "$RESULTS/rfmix"

bash "$ROOT/00_run_rfmix.sh" \
  --reference-bcf data/rfmix_inputs/chr1.reference.clean.bcf \
  --query-bcf data/rfmix_inputs/chr1.query.clean.bcf \
  --sample-map data/rfmix_inputs/reference_sample_map.tsv \
  --genetic-map data/rfmix_inputs/chr1.uniform_fullchr_1Mb_1cM.map.tsv \
  --chromosome 1 \
  --output-prefix "$RESULTS/rfmix/chr1.rfmix.fullchr"

python3 "$ROOT/01_merge_filter_tracts.py" \
  --msp "$RESULTS/rfmix/chr1.rfmix.fullchr.msp.tsv" \
  --fb "$RESULTS/rfmix/chr1.rfmix.fullchr.fb.tsv" \
  --output "$RESULTS/tracts.all.tsv" \
  --strict-output "$RESULTS/tracts.strict.tsv"

python3 "$ROOT/02_make_gene_intervals.py" \
  --tracts "$RESULTS/tracts.strict.tsv" \
  --gtf data/annotation.gtf \
  --output "$RESULTS/gene_associated_intervals.tsv"

python3 "$ROOT/03_map_phaser_blocks.py" \
  --phaser-blocks data/phaser/sample.haplotypic_counts.tsv \
  --phased-genotypes data/phased/sample.vcf.gz \
  --gene-intervals "$RESULTS/gene_associated_intervals.tsv" \
  --rna-sample RNA_001 \
  --wgs-individual WGS_001 \
  --tissue Lung \
  --audit-output "$RESULTS/RNA_001.blocks.audit.tsv" \
  --strict-output "$RESULTS/RNA_001.blocks.strict.tsv"

python3 "$ROOT/04_test_interval_ASE.py" \
  --strict-blocks "$RESULTS/RNA_001.blocks.strict.tsv" \
  --all-output "$RESULTS/interval_ASE.all.tsv" \
  --tested-output "$RESULTS/interval_ASE.tested.tsv" \
  --significant-output "$RESULTS/interval_ASE.significant.tsv"

python3 "$ROOT/05_check_site_ASE.py" \
  --site-ase data/site_ASE/significant.tsv \
  --interval-ase "$RESULTS/interval_ASE.significant.tsv" \
  --phased-genotypes data/phased/sample.vcf.gz \
  --output "$RESULTS/site_ASE.interval_annotated.tsv"
