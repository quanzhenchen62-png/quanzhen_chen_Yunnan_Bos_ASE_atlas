### `01.expression_phylogeny`

Expression-based phylogeny analysis and matched-background testing.

`01_define_targets_and_fixed_bins.R`
  defines target gene sets and fixed background bins
`02_match_backgrounds_and_build_trees.R`
  matches background genes and builds expression trees
`03_empirical_tests_and_violin_plots.R`
  performs empirical comparisons between observed and matched-background tree statistics and generates violin plots

### `02.ase_detection_quantification`

Site-level ASE detection and quantification, followed by breed-level and population-level ASE summarization.

`01_run_star_wasp.sh`
  STAR alignment with WASP filtering
`02_run_ase_read_counter.sh`
  GATK ASEReadCounter wrapper
`03_filter_test_ase_sites.py`
  filters sites, runs per-sample binomial tests, applies BH correction, and calculates individual `aFC`
`04_call_brd_ase.py`
  generates `brd-ASE-1`, `brd-ASE-2`, and `brd-ASE-3`
`05_call_pop_ase.py`
  generates tissue-level `pop-ASE`
`run_example.sh`
  minimal example workflow

### `03.phaser_haplotypic_expression`

phASER-based linking of site-level ASE to gene-level haplotypic expression.

`01_run_phaser.sh`
  runs `phaser.py` per sample
`02_make_gene_feature_bed.py`
  creates gene-body BED features from GTF
`03_run_phaser_gene_ae.sh`
  runs `phaser_gene_ae.py`
`04_run_phaser_expr_matrix.sh`
  runs `phaser_expr_matrix.py`
`05_make_phaser_cis_map.py`
  builds the VCF-to-BED sample map for cis-variant analysis
`06_run_phaser_cis_var.sh`
  runs `phaser_cis_var.py`
`run_example.sh`
  minimal example workflow

### `04.selective_sweep_scan`

Selective sweep analyses based on windowed population-genetic statistics.

`01_windowed_fst.sh`
  runs windowed `FST`
`02_windowed_pi.sh`
  runs windowed `theta-pi`
`03_calc_pi_ratio.py`
  calculates `pi_ratio`
`04_run_selscan_xpehh.sh`
  runs chromosome-level `XP-EHH`
`05_summarize_xpehh_windows.py`
  summarizes normalized `XP-EHH` into sliding windows
`06_make_candidate_regions.py`
  creates candidate selective-sweep regions
`07_make_fst_deciles.py`
  assigns `FST` windows to deciles
`run_workflow.sh`
  example workflow order

### `05.introgression_analysis`

Introgression analyses based on D-statistics and U20.

`01_build_reference_dstat_tests.py`
  builds individual-level D-statistic tests for reference-panel screening
`02_build_target_dstat_tests.py`
  builds target-breed D-statistic tests
`03_run_qpdstat.py`
  runs ADMIXTOOLS `qpDstat`
`04_filter_reference_panel.py`
  retains putatively unadmixed reference individuals
`05_scan_u20_windows.py`
  calculates U20 in non-overlapping 50-kb windows and identifies top windows
`run_example.sh`
  minimal example workflow

### `06.local_ancestry_ASE`

Local-ancestry-resolved ASE workflow using RFMix and phASER outputs.

`00_run_rfmix.sh`
  runs chromosome-level RFMix local ancestry inference using prepared reference/query BCFs and genetic maps
`01_merge_filter_tracts.py`
  merges and filters RFMix tracts
`02_make_gene_intervals.py`
  intersects tracts with genes
`03_map_phaser_blocks.py`
  maps phASER blocks to ancestry-resolved haplotypes
`04_test_interval_ASE.py`
  tests ancestry-resolved interval ASE
`05_check_site_ASE.py`
  annotates interval ASE with site-level ASE evidence
`config.example.yaml`
  example configuration
`run_example.sh`
  minimal example workflow

### `07.ase_event_and_gene_enrichment`

Downstream enrichment analyses for ASE events and ASE-related gene sets.

`01_prepare_input_tables_for_ASE_event_enrichment.py`
  prepares SNP-level input tables
`02_bin_background_SNPs_for_ASE_event_enrichment.py`
  bins matched background SNPs
`03_test_enrichment_of_ASE_events_in_selection_regions.py`
  tests ASE-event enrichment in selection regions and FST windows
`04_prepare_gene_overlap_table_for_enrichment.py`
  prepares gene-overlap tables
`05_test_enrichment_of_gene_sets_in_selection_regions.py`
  tests gene-set enrichment in selection regions and FST windows
`06_test_enrichment_of_ASE_events_in_local_ancestry_segments.py`
  tests ASE-event enrichment in local-ancestry segments
`run_enrichment_analysis_workflow.sh`
  example workflow order

### `08.popdiff_ASE_variants_annotation_enrichment`

Functional annotation and enrichment analysis for population-differentiation-associated ASE variants.

`ASE_four_class_annotation_enrichment.py`
  analyzes four ASE-variant classes across chromatin states, QTL annotations, and genomic functional categories

### `09.phenotype_altitude_environment_association`

Association analyses linking ASE variants to phenotypes and environmental variables.

`phenotype_environment_association.R`
  phenotype-environment association
`phenotype_genotype_association.py`
  phenotype-associated ASE variant identification with linear mixed models
`genotype_environment_association.R`
  genotype-environment association for phenotype-associated ASE variants
`run_example.sh`
  minimal example workflow
