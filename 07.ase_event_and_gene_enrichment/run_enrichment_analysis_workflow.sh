#!/usr/bin/env bash
set -euo pipefail

python3 01_prepare_input_tables_for_ASE_event_enrichment.py   --sample-meta metadata/sample_meta.tsv   --target-manifest metadata/selection_target_manifest.tsv   --outdir results/01_ASE_event_input_tables

python3 02_bin_background_SNPs_for_ASE_event_enrichment.py   --sample-meta results/01_ASE_event_input_tables/sample_meta.with_prepared_snp_tables.tsv   --outdir results/02_binned_ASE_event_background

python3 03_test_enrichment_of_ASE_events_in_selection_regions.py   --sample-meta results/02_binned_ASE_event_background/sample_meta.with_clean_binned.tsv   --outdir results/03_ASE_event_enrichment_in_selection_regions

python3 04_prepare_gene_overlap_table_for_enrichment.py   --background-genes metadata/background_genes.bed   --target-gene-classes metadata/target_gene_classes.bed   --target-manifest metadata/selection_target_manifest.tsv   --outdir results/04_gene_overlap_tables

python3 05_test_enrichment_of_gene_sets_in_selection_regions.py   --gene-table results/04_gene_overlap_tables/background_gene_selection_table.tsv   --output results/05_gene_set_enrichment/gene_set_enrichment.tsv

python3 06_test_enrichment_of_ASE_events_in_local_ancestry_segments.py   --meta metadata/local_ancestry_sample_meta.tsv   --prepared-dir results/01_ASE_event_input_tables/prepared_snp_tables   --tracts metadata/high_confidence_local_ancestry_tracts.tsv.gz   --outdir results/06_ASE_event_enrichment_in_local_ancestry_segments
