#!/usr/bin/env Rscript
# Group-exclusion sensitivity analysis for the expression phylogenies
# (e.g. excluding a group quantified on a non-species-matched reference).
#
# Reads the same minimal inputs as scripts 01-02 under data/ and writes
# transformed copies under data_excluded/; then re-run scripts 01-03 with
# data/ replaced by data_excluded/ (or copied over) to obtain the
# sensitivity-analysis trees and tests with all procedures otherwise identical.
#
# Additional expected column:
#   support_calls.tsv: supporting_groups (semicolon-separated group names),
#   so that support counts can be recomputed without the excluded group.
#
# Note: the upstream group expression statistics, G0 expression filtering, and
# mean_expression matching covariate should be recomputed from the sample set
# without the excluded group before this step, as in the formal pipeline.
#
# Usage: Rscript 04_group_exclusion_sensitivity.R <group_name_to_exclude>

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 04_group_exclusion_sensitivity.R <group_name_to_exclude>")
exclude <- args[1L]

dir.create("data_excluded", showWarnings = FALSE, recursive = TRUE)

manifest <- fread("data/analysis_manifest.tsv")
support <- fread("data/support_calls.tsv")
stats <- fread("data/group_expression_stats.tsv")
g0 <- fread("data/g0_features.tsv")

dropped_groups <- unique(stats[group == exclude, unique(analysis_id)])
manifest[, total_group_n := total_group_n - 1L]
fwrite(manifest, "data_excluded/analysis_manifest.tsv", sep = "\t")

if (!"supporting_groups" %in% names(support)) {
  stop("support_calls.tsv must contain a supporting_groups column for group exclusion")
}
support[, supporting_groups := vapply(strsplit(supporting_groups, ";", fixed = TRUE),
                                      function(x) paste(setdiff(x, exclude), collapse = ";"), "")]
support[, support_n := vapply(strsplit(supporting_groups, ";", fixed = TRUE),
                              function(x) sum(nzchar(x)), 0L)]
fwrite(support, "data_excluded/support_calls.tsv", sep = "\t")

stats <- stats[group != exclude]
fwrite(stats, "data_excluded/group_expression_stats.tsv", sep = "\t")

fwrite(g0, "data_excluded/g0_features.tsv", sep = "\t")

cat("Excluded group:", exclude, "\n")
cat("Analyses affected:", paste(dropped_groups, collapse = ", "), "\n")
cat("Re-run scripts 01-03 against data_excluded/ to obtain sensitivity results.\n")
