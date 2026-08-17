#!/usr/bin/env Rscript

# Additional minimal input expected under data/:
#
# group_expression_stats.tsv
# analysis_id, gene_id, group, median_expression, replicate_variance, replicate_n
#
# median_expression: median log2(TMM-normalized CPM + 1) across biological replicates.
# replicate_variance: variance of that gene across biological replicates within the group.
# replicate_n: number of biological replicates in that group.
#
# results/target_sets.tsv and results/g0_features_with_fixed_bins.tsv are produced by script 01.
# The same G0, fixed 10-bin assignments, samples, expression values and denominator are used when
# comparing each observed target set with its matched background gene sets.

suppressPackageStartupMessages({
  library(data.table)
  library(ape)
})

dir.create("results", showWarnings = FALSE, recursive = TRUE)

n_bins <- 10L
background_replicates <- 1000L
base_seed <- 20260728L

targets <- fread("results/target_sets.tsv")
g0 <- fread("results/g0_features_with_fixed_bins.tsv")
group_stats <- fread("data/group_expression_stats.tsv")

bin_columns <- c("expression_bin", "length_bin", "gc_bin", "exon_bin")

compute_fixed_denominator <- function(current_stats, g0_genes) {
  x <- current_stats[gene_id %in% g0_genes]
  components <- x[, {
    n_values <- unique(replicate_n)
    if (length(n_values) != 1L || n_values < 2L) stop("Invalid replicate count")
    theta_g0 <- mean(replicate_variance)
    list(
      replicate_n = n_values,
      theta_g0 = theta_g0,
      group_component = var(median_expression) - theta_g0 / n_values
    )
  }, by = group]
  denominator <- mean(components$group_component)
  if (!is.finite(denominator) || denominator <= 0) stop("Invalid fixed G0 denominator")
  denominator
}

wang_distance <- function(genes, current_stats, denominator) {
  x <- current_stats[gene_id %in% genes]
  median_wide <- dcast(x, gene_id ~ group, value.var = "median_expression")
  variance_wide <- dcast(x, gene_id ~ group, value.var = "replicate_variance")
  groups <- setdiff(names(median_wide), "gene_id")
  if (length(groups) < 3L || nrow(median_wide) < 2L) stop("Insufficient tree input")
  setorder(median_wide, gene_id)
  setorder(variance_wide, gene_id)
  if (!identical(median_wide$gene_id, variance_wide$gene_id)) stop("Gene order mismatch")

  medians <- as.matrix(median_wide[, ..groups])
  replicate_variances <- as.matrix(variance_wide[, ..groups])
  theta <- colMeans(replicate_variances)
  replicate_n <- x[, unique(replicate_n), by = group]
  n_by_group <- setNames(replicate_n$V1, replicate_n$group)[groups]

  distance_matrix <- matrix(0, length(groups), length(groups), dimnames = list(groups, groups))
  pair_qc <- vector("list", choose(length(groups), 2L))
  row_index <- 0L
  for (first in seq_len(length(groups) - 1L)) {
    for (second in seq.int(first + 1L, length(groups))) {
      row_index <- row_index + 1L
      variance_before <- var(medians[, first] - medians[, second])
      correction_first <- theta[first] / n_by_group[first]
      correction_second <- theta[second] / n_by_group[second]
      corrected_raw <- variance_before - correction_first - correction_second
      nonnegative_raw <- max(0, corrected_raw)
      normalized <- nonnegative_raw / denominator
      distance_matrix[first, second] <- normalized
      distance_matrix[second, first] <- normalized
      pair_qc[[row_index]] <- data.table(
        group_i = groups[first],
        group_j = groups[second],
        variance_before_correction = variance_before,
        theta_i_divided_by_n_i = correction_first,
        theta_j_divided_by_n_j = correction_second,
        corrected_raw_distance = corrected_raw,
        corrected_distance_is_negative = corrected_raw < 0,
        nonnegative_raw_distance = nonnegative_raw,
        fixed_G0_denominator = denominator,
        normalized_NJ_distance = normalized
      )
    }
  }
  list(matrix = distance_matrix, pair_qc = rbindlist(pair_qc))
}

draw_without_replacement <- function(indices, number, used) {
  available <- indices[!used[indices]]
  if (!length(available) || number <= 0L) return(integer())
  take <- min(number, length(available))
  available[sample.int(length(available), take, replace = FALSE)]
}

match_one_background <- function(target_table, candidate_table, seed) {
  set.seed(seed)
  candidate_table <- copy(candidate_table)
  candidate_table[, candidate_index := .I]
  used <- rep(FALSE, nrow(candidate_table))
  selected <- integer()

  strata <- target_table[, .N, by = bin_columns]
  deficits <- vector("list", nrow(strata))

  # Phase 1: fill every exact joint stratum before borrowing from neighboring strata.
  for (row_index in seq_len(nrow(strata))) {
    target_bin <- as.integer(strata[row_index, ..bin_columns])
    exact <- which(Reduce(`&`, Map(`==`, candidate_table[, ..bin_columns], target_bin)))
    chosen <- draw_without_replacement(exact, strata$N[row_index], used)
    used[chosen] <- TRUE
    selected <- c(selected, chosen)
    deficits[[row_index]] <- strata$N[row_index] - length(chosen)
  }

  # Phase 2: exact shortages use one-coordinate/one-bin neighbors first, then nearest bins.
  for (row_index in sample(seq_len(nrow(strata)))) {
    deficit <- deficits[[row_index]]
    if (deficit == 0L) next
    target_bin <- as.integer(strata[row_index, ..bin_columns])
    available <- which(!used)
    candidate_bins <- as.matrix(candidate_table[available, ..bin_columns])
    differences <- abs(sweep(candidate_bins, 2L, target_bin, `-`))
    one_bin_neighbor <- rowSums(differences == 1L) == 1L & rowSums(differences == 0L) == 3L
    standardized_distance <- sqrt(rowSums((differences / (n_bins - 1L))^2))
    random_tie_break <- runif(length(available))
    ranked <- available[order(!one_bin_neighbor, standardized_distance, random_tie_break)]
    if (length(ranked) < deficit) stop("Candidate pool is smaller than target set")
    chosen <- ranked[seq_len(deficit)]
    used[chosen] <- TRUE
    selected <- c(selected, chosen)
  }

  if (length(selected) != nrow(target_table) || anyDuplicated(selected)) {
    stop("Background size or uniqueness check failed")
  }
  candidate_table$gene_id[selected]
}

configurations <- unique(targets[, .(
  support_threshold, required_support_n, analysis_id,
  analysis_level, tissue, target_type
)])

tree_length_rows <- list()
background_qc_rows <- list()
newick_rows <- list()
result_index <- 0L

for (configuration_index in seq_len(nrow(configurations))) {
  config <- configurations[configuration_index]
  current_g0 <- g0[analysis_id == config$analysis_id]
  current_stats <- group_stats[analysis_id == config$analysis_id]
  target_genes <- targets[
    analysis_id == config$analysis_id &
      target_type == config$target_type &
      support_threshold == config$support_threshold,
    unique(gene_id)
  ]
  candidate_genes <- setdiff(current_g0$gene_id, target_genes)

  if (length(target_genes) < 2L || length(candidate_genes) < length(target_genes)) next

  denominator <- compute_fixed_denominator(current_stats, current_g0$gene_id)
  target_table <- current_g0[gene_id %in% target_genes]
  candidate_table <- current_g0[gene_id %in% candidate_genes]

  observed_distance <- wang_distance(target_genes, current_stats, denominator)
  observed_tree <- nj(as.dist(observed_distance$matrix))
  result_index <- result_index + 1L
  tree_length_rows[[result_index]] <- data.table(
    config,
    tree_source = "observed",
    background_replicate = 0L,
    n_genes = length(target_genes),
    total_tree_length = sum(observed_tree$edge.length),
    negative_branch_n = sum(observed_tree$edge.length < 0)
  )
  newick_rows[[result_index]] <- data.table(config, tree_source = "observed",
                                             background_replicate = 0L,
                                             newick = write.tree(observed_tree))

  analysis_number <- match(config$analysis_id, unique(configurations$analysis_id))
  type_offset <- if (config$target_type == "tissue_specific") 0L else 5000L

  for (replicate in seq_len(background_replicates)) {
    # The seed depends on the analysis configuration so that target and background
    # sampling are reproducible for each run.
    seed <- base_seed + analysis_number * 1000000L +
      config$required_support_n * 10000L + type_offset + replicate
    background_genes <- match_one_background(target_table, candidate_table, seed)
    background_distance <- wang_distance(background_genes, current_stats, denominator)
    background_tree <- nj(as.dist(background_distance$matrix))

    result_index <- result_index + 1L
    tree_length_rows[[result_index]] <- data.table(
      config,
      tree_source = "background",
      background_replicate = replicate,
      n_genes = length(background_genes),
      total_tree_length = sum(background_tree$edge.length),
      negative_branch_n = sum(background_tree$edge.length < 0)
    )
    newick_rows[[result_index]] <- data.table(config, tree_source = "background",
                                               background_replicate = replicate,
                                               newick = write.tree(background_tree))
    background_qc_rows[[length(background_qc_rows) + 1L]] <- data.table(
      config,
      background_replicate = replicate,
      random_seed = seed,
      expected_gene_n = length(target_genes),
      observed_gene_n = length(background_genes),
      unique_gene_n = uniqueN(background_genes),
      target_overlap_n = sum(background_genes %in% target_genes)
    )
  }
}

fwrite(rbindlist(tree_length_rows), "results/tree_lengths_long.tsv", sep = "\t")
fwrite(rbindlist(newick_rows), "results/trees_newick.tsv", sep = "\t")
fwrite(rbindlist(background_qc_rows), "results/background_set_qc.tsv", sep = "\t")
