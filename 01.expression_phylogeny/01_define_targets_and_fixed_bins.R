#!/usr/bin/env Rscript

# Minimal input files expected under data/:
#
# 1) analysis_manifest.tsv
#    analysis_id, analysis_level, tissue, total_group_n
#    One row per analysis level x tissue. Exclude unwanted tissues before this step.
#
# 2) support_calls.tsv
#    analysis_id, target_type, gene_id, support_n
#    target_type is tissue_specific or broadly_expressed.
#
# 3) g0_features.tsv
#    analysis_id, gene_id, mean_expression, gene_length, gc_content, exon_number
#    G0 contains genes that passed expression filtering, have complete expression in all
#    participating groups, and have complete matching features. For cross-species analyses,
#    gene_id is a strict one-to-one ortholog ID and structural features are species medians.
#
# This script does not run differential expression, TAU, count normalization, orthology
# inference, or gene-feature calculation. Those inputs are assumed to have been prepared.

suppressPackageStartupMessages(library(data.table))

dir.create("results", showWarnings = FALSE, recursive = TRUE)

n_bins <- 10L

support_proportions <- c(0.20)

manifest <- fread("data/analysis_manifest.tsv")
support_calls <- fread("data/support_calls.tsv")
g0 <- fread("data/g0_features.tsv")

stopifnot(
  !anyDuplicated(manifest$analysis_id),
  !anyDuplicated(g0[, paste(analysis_id, gene_id)]),
  all(support_calls$target_type %in% c("tissue_specific", "broadly_expressed"))
)

thresholds <- data.table(
  support_proportion = support_proportions,
  support_threshold = as.integer(100 * support_proportions)
)

design <- merge(manifest, thresholds, by = NULL, allow.cartesian = TRUE)
design[, required_support_n := ceiling(support_proportion * total_group_n)]

targets <- merge(
  support_calls,
  design[, .(analysis_id, analysis_level, tissue, total_group_n,
             support_threshold, support_proportion, required_support_n)],
  by = "analysis_id",
  allow.cartesian = TRUE
)
targets <- targets[support_n >= required_support_n]

# Restrict all thresholds to the same fixed, fully annotated and expressed G0 universe.
targets <- merge(
  targets,
  unique(g0[, .(analysis_id, gene_id)]),
  by = c("analysis_id", "gene_id")
)

rank_bin <- function(values, n_bins = 10L) {
  ranks <- frank(values, ties.method = "average") / length(values)
  pmin(n_bins, as.integer(ceiling(n_bins * ranks)))
}

# Decile bin boundaries are defined once from G0 and are not recomputed for different target thresholds.
g0[, expression_bin := rank_bin(mean_expression, n_bins), by = analysis_id]
g0[, length_bin := rank_bin(log10(gene_length), n_bins), by = analysis_id]
g0[, gc_bin := rank_bin(gc_content, n_bins), by = analysis_id]
g0[, exon_bin := rank_bin(log10(exon_number), n_bins), by = analysis_id]
g0[, joint_bin := paste(expression_bin, length_bin, gc_bin, exon_bin, sep = "_")]

fwrite(targets, "results/target_sets.tsv", sep = "\t")
fwrite(g0, "results/g0_features_with_fixed_bins.tsv", sep = "\t")
fwrite(
  targets[, .(
    n_target_genes = uniqueN(gene_id)
  ), by = .(support_threshold, required_support_n, analysis_id,
            analysis_level, tissue, target_type)],
  "results/target_gene_counts.tsv",
  sep = "\t"
)
