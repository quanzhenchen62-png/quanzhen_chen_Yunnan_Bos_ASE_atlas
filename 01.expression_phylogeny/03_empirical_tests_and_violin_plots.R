#!/usr/bin/env Rscript

# Input:
# results/tree_lengths_long.tsv from script 02.
# Each valid analysis level x tissue x target type must contain one observed tree
# and 1,000 matched background trees.
#
# One-sided empirical P values test longer tissue-specific trees and shorter
# broadly expressed trees. Empirical P values are calculated as (1 + r) / 1001,
# where r is the number of matched background trees at least as extreme as the
# observed tree in the expected direction. P values are adjusted together using
# the Benjamini-Hochberg method.
#
# For visualization, relative tree length is defined as log2(L / M), where L is
# an observed or randomized tree length and M is the median of its 1,000 matched
# background tree lengths. Zero therefore denotes the matched-background median.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

dir.create("results/plots", showWarnings = FALSE, recursive = TRUE)

tree_lengths <- fread("results/tree_lengths_long.tsv")
group_columns <- c(
  "support_threshold", "required_support_n", "analysis_id",
  "analysis_level", "tissue", "target_type"
)

test_rows <- list()
configurations <- unique(tree_lengths[, ..group_columns])

for (row_index in seq_len(nrow(configurations))) {
  config <- configurations[row_index]
  current <- tree_lengths[
    support_threshold == config$support_threshold &
      required_support_n == config$required_support_n &
      analysis_id == config$analysis_id &
      analysis_level == config$analysis_level &
      tissue == config$tissue &
      target_type == config$target_type
  ]

  observed <- current[tree_source == "observed", total_tree_length]
  background <- current[tree_source == "background", total_tree_length]
  if (length(observed) != 1L || length(background) != 1000L) next

  r_value <- if (config$target_type == "tissue_specific") {
    sum(background >= observed)
  } else {
    sum(background <= observed)
  }
  empirical_p <- (1 + r_value) / 1001
  background_median <- median(background)

  test_rows[[row_index]] <- data.table(
    config,
    r = r_value,
    observed_total_tree_length = observed,
    background_mean = mean(background),
    background_median = background_median,
    background_2.5_percentile = quantile(background, 0.025, names = FALSE),
    background_97.5_percentile = quantile(background, 0.975, names = FALSE),
    observed_minus_background_median = observed - background_median,
    observed_divided_by_background_median = observed / background_median,
    log2_observed_divided_by_background_median = log2(observed / background_median),
    empirical_P = empirical_p,
    valid_background_n = length(background)
  )
}

tests <- rbindlist(test_rows, fill = TRUE)
tests[, FDR_BH := NA_real_]
if (nrow(tests[is.finite(empirical_P)]) > 0L) {
  tests[is.finite(empirical_P), FDR_BH := p.adjust(empirical_P, method = "BH")]
}
fwrite(tests, "results/empirical_tests.tsv", sep = "	")

background_medians <- tree_lengths[
  tree_source == "background",
  .(background_median = median(total_tree_length)),
  by = group_columns
]

plot_data <- merge(tree_lengths, background_medians, by = group_columns, all.x = TRUE)
plot_data[, relative_tree_length := fifelse(
  is.finite(background_median) & background_median > 0,
  log2(total_tree_length / background_median),
  NA_real_
)]

fwrite(plot_data, "results/tree_lengths_relative_long.tsv", sep = "	")

# Example visualization for the formal 20% threshold.
plot_data <- plot_data[support_threshold == 20 & is.finite(relative_tree_length)]
plot_data[, tissue := factor(
  tissue,
  levels = c("Heart", "Kidney", "Liver", "Lung", "Muscle", "Spleen")
)]
plot_data[, target_type := factor(
  target_type,
  levels = c("tissue_specific", "broadly_expressed")
)]

fill_colors <- c(tissue_specific = "#bc8149", broadly_expressed = "#9e4951")
point_colors <- c(tissue_specific = "#bc8149", broadly_expressed = "#9e4951")

for (current_level in unique(plot_data$analysis_level)) {
  background_data <- plot_data[analysis_level == current_level & tree_source == "background"]
  observed_data <- plot_data[analysis_level == current_level & tree_source == "observed"]

  figure <- ggplot() +
    geom_violin(
      data = background_data,
      aes(
        x = tissue,
        y = relative_tree_length,
        fill = target_type,
        group = interaction(tissue, target_type)
      ),
      position = position_dodge(width = 0.82),
      width = 0.78,
      scale = "width",
      alpha = 0.30,
      color = NA
    ) +
    geom_point(
      data = observed_data,
      aes(
        x = tissue,
        y = relative_tree_length,
        color = target_type,
        group = target_type
      ),
      position = position_dodge(width = 0.82),
      size = 3.0
    ) +
    scale_fill_manual(values = fill_colors, name = NULL) +
    scale_color_manual(values = point_colors, name = NULL) +
    labs(x = NULL, y = "Relative tree length, log2(L/M)", title = current_level) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 35, hjust = 1),
      legend.position = "top"
    )

  ggsave(
    filename = file.path("results/plots", paste0(current_level, "_tree_length_violin.pdf")),
    plot = figure,
    width = 8.5,
    height = 5.5
  )
}
