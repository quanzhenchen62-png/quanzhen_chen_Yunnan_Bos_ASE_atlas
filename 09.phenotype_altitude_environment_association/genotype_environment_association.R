#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(LEA))

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  defaults <- list(
    candidate_lfmm = "data/GENO_candidate.lfmm",
    candidate_snp_info = "data/GENO_candidate.snp.info",
    candidate_ind_order = "data/GENO_candidate.ind.order",
    pruned_lfmm = "data/Bos_pruned.lfmm",
    pruned_ind_order = "data/Bos_pruned.ind.order",
    env_table = "data/Bos_ind_env.info",
    outdir = "results/phenotype_associated_ASE_variant_environment_association",
    individual_col = "individual",
    altitude_col = "altitude",
    temperature_col = "temperature_annual_mean",
    precipitation_col = "precipitation_annual",
    uvb_col = "UVB_annual_mean",
    longitude_col = "longitude",
    latitude_col = "latitude",
    k = 6
  )

  if (length(args) == 0) {
    return(defaults)
  }

  i <- 1
  while (i <= length(args)) {
    key <- args[i]
    if (!startsWith(key, "--")) {
      stop("Unknown argument: ", key)
    }
    if (i == length(args)) {
      stop("Missing value for ", key)
    }

    value <- args[i + 1]
    name <- sub("^--", "", key)
    name <- gsub("-", "_", name)

    if (!(name %in% names(defaults))) {
      stop("Unknown option: ", key)
    }

    defaults[[name]] <- value
    i <- i + 2
  }

  defaults$k <- as.integer(defaults$k)
  defaults
}

ARGS <- parse_args()

dir.create(ARGS$outdir, showWarnings = FALSE, recursive = TRUE)

require_file <- function(path) {
  if (!file.exists(path)) {
    stop("Missing file: ", path)
  }
}

for (path in c(
  ARGS$candidate_lfmm,
  ARGS$candidate_snp_info,
  ARGS$candidate_ind_order,
  ARGS$pruned_lfmm,
  ARGS$pruned_ind_order,
  ARGS$env_table
)) {
  require_file(path)
}

read_ind_order <- function(path) {
  x <- scan(path, what = "character", quiet = TRUE)
  if (length(x) == 0) {
    stop("Empty individual-order file: ", path)
  }
  x
}

read_snp_info <- function(path) {
  x <- read.table(
    path,
    header = FALSE,
    sep = "\t",
    stringsAsFactors = FALSE,
    fill = TRUE,
    quote = "",
    comment.char = ""
  )
  snps <- trimws(x[[1]])
  snps <- snps[nzchar(snps)]
  if (length(snps) == 0) {
    stop("Empty SNP-info file: ", path)
  }
  if (snps[1] %in% c("SNP", "variant_id", "candidate_snp")) {
    snps <- snps[-1]
  }
  snps
}

coerce_to_snp_by_env <- function(x, n_snp, n_env, what = "matrix") {
  if (is.null(dim(x))) {
    x <- as.numeric(x)
    if (length(x) == n_snp * n_env) {
      return(matrix(x, nrow = n_snp, ncol = n_env))
    }
    if (length(x) == n_snp && n_env == 1) {
      return(matrix(x, nrow = n_snp, ncol = 1))
    }
    stop(sprintf(
      "%s length %d does not match n_snp=%d and n_env=%d",
      what, length(x), n_snp, n_env
    ))
  }

  m <- as.matrix(x)
  if (nrow(m) == n_snp && ncol(m) == n_env) return(m)
  if (nrow(m) == n_env && ncol(m) == n_snp) return(t(m))

  stop(sprintf(
    "%s dimensions %d x %d do not match n_snp=%d and n_env=%d",
    what, nrow(m), ncol(m), n_snp, n_env
  ))
}

candidate_ind <- read_ind_order(ARGS$candidate_ind_order)
pruned_ind <- read_ind_order(ARGS$pruned_ind_order)

if (!identical(candidate_ind, pruned_ind)) {
  stop("Candidate and pruned genotype matrices do not have identical individual order.")
}

candidate_snp <- read_snp_info(ARGS$candidate_snp_info)

env_raw <- read.table(
  ARGS$env_table,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (!(ARGS$individual_col %in% colnames(env_raw))) {
  stop("Environmental table is missing individual column: ", ARGS$individual_col)
}

env_cols <- c(
  ARGS$altitude_col,
  ARGS$temperature_col,
  ARGS$precipitation_col,
  ARGS$uvb_col,
  ARGS$longitude_col,
  ARGS$latitude_col
)

missing_cols <- setdiff(c(ARGS$individual_col, env_cols), colnames(env_raw))
if (length(missing_cols) > 0) {
  stop("Environmental table is missing columns: ", paste(missing_cols, collapse = ", "))
}

idx <- match(pruned_ind, env_raw[[ARGS$individual_col]])
if (any(is.na(idx))) {
  stop("Some individuals from genotype matrices are missing in the environmental table.")
}

env_ordered <- env_raw[idx, c(ARGS$individual_col, env_cols), drop = FALSE]

if (!identical(pruned_ind, env_ordered[[ARGS$individual_col]])) {
  stop("Environmental table could not be reordered to match genotype individual order.")
}

for (col in env_cols) {
  env_ordered[[col]] <- as.numeric(env_ordered[[col]])
  if (any(is.na(env_ordered[[col]]))) {
    stop("Environmental predictor contains NA after numeric coercion: ", col)
  }
}

env_mat <- scale(as.matrix(env_ordered[, env_cols, drop = FALSE]))
colnames(env_mat) <- env_cols

lfmm_model <- lfmm2(
  input = ARGS$pruned_lfmm,
  env = env_mat,
  K = ARGS$k
)

saveRDS(lfmm_model, file.path(ARGS$outdir, "lfmm2_model_pruned_genomewide.rds"))

lfmm_test <- lfmm2.test(
  object = lfmm_model,
  input = ARGS$candidate_lfmm,
  env = env_mat,
  linear = TRUE,
  genomic.control = FALSE,
  full = FALSE
)

p_mat <- coerce_to_snp_by_env(
  lfmm_test$pvalues,
  n_snp = length(candidate_snp),
  n_env = length(env_cols),
  what = "pvalue matrix"
)
z_mat <- coerce_to_snp_by_env(
  lfmm_test$zscores,
  n_snp = length(candidate_snp),
  n_env = length(env_cols),
  what = "zscore matrix"
)

colnames(p_mat) <- env_cols
rownames(p_mat) <- candidate_snp
colnames(z_mat) <- env_cols
rownames(z_mat) <- candidate_snp

result <- data.frame(
  variant_id = rep(candidate_snp, times = ncol(p_mat)),
  environmental_variable = rep(colnames(p_mat), each = nrow(p_mat)),
  p_value = as.vector(p_mat),
  z_score = as.vector(z_mat),
  stringsAsFactors = FALSE
)

result$FDR <- ave(
  result$p_value,
  result$environmental_variable,
  FUN = function(x) p.adjust(x, method = "BH")
)
result$significant_FDR_0_1 <- ifelse(result$FDR < 0.1, 1, 0)

write.table(
  result,
  file.path(ARGS$outdir, "phenotype_associated_ASE_variant_environment_association.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

summary_df <- aggregate(
  significant_FDR_0_1 ~ environmental_variable,
  data = result,
  FUN = sum
)
summary_df$tested_variant_n <- aggregate(
  variant_id ~ environmental_variable,
  data = result,
  FUN = length
)$variant_id
summary_df <- summary_df[, c("environmental_variable", "tested_variant_n", "significant_FDR_0_1")]
colnames(summary_df)[3] <- "significant_variant_n"

write.table(
  summary_df,
  file.path(ARGS$outdir, "phenotype_associated_ASE_variant_environment_association.summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
