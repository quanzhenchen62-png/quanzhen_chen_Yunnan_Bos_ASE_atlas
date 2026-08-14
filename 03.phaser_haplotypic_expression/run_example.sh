#!/bin/bash

set -euo pipefail

sample_id="sample_001"

bash 01_run_phaser.sh \
  "${sample_id}" \
  data/${sample_id}.het.vcf.gz \
  data/${sample_id}.WASP_filtered.bam \
  tools/phaser/phaser.py \
  refs/phaser_blacklist.bed \
  12 \
  1 \
  results/phaser_haplo_count/${sample_id}/${sample_id}

python3 02_make_gene_feature_bed.py \
  --gtf refs/annotation.gtf \
  --output refs/genes.genebody.bed

bash 03_run_phaser_gene_ae.sh \
  results/phaser_haplo_count/${sample_id}/${sample_id}.haplotypic_counts.txt \
  refs/genes.genebody.bed \
  tools/phaser/phaser_gene_ae.py \
  results/phaser_gene_ae/${sample_id}_phaser_gene_ae.txt

bash 04_run_phaser_expr_matrix.sh \
  results/phaser_gene_ae \
  refs/genes.genebody.bed \
  tools/phaser/phaser_expr_matrix.py \
  4 \
  results/phaser_expr_matrix/tissue_A.haplotype_count_matrix

python3 05_make_phaser_cis_map.py \
  --vcf data/tissue_A.significant_ase_pairs.vcf.gz \
  --bed data/tissue_A.gene_ae.bed.gz \
  --output results/phaser_cis_var/tissue_A.phaser_map.txt

bash 06_run_phaser_cis_var.sh \
  data/tissue_A.significant_ase_variant_gene_pairs.bed \
  data/tissue_A.gene_ae.bed.gz \
  data/tissue_A.significant_ase_pairs.vcf.gz \
  results/phaser_cis_var/tissue_A.phaser_map.txt \
  tools/phaser/phaser_cis_var.py \
  results/phaser_cis_var/tissue_A.phaser_cis_var.txt
