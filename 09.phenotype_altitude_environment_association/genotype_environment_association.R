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
    units_file = "data/panel_units.tsv",
    outdir = "results/phenotype_associated_ASE_variant_environment_association",
    individual_col = "individual",
    altitude_col = "altitude",
    temperature_col = "temperature_annual_mean",
    precipitation_col = "precipitation_annual",
    uvb_col = "UVB_annual_mean",
    longitude_col = "longitude",
    latitude_col = "latitude",
    k = 6,
    n_perm = 1000,
    seed = 20260924
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
  defaults$n_perm <- as.integer(defaults$n_perm)
  defaults$seed <- as.integer(defaults$seed)
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
  ARGS$env_table,
  ARGS$units_file
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

env6 <- c(
  ARGS$altitude_col,
  ARGS$uvb_col,
  ARGS$temperature_col,
  ARGS$precipitation_col,
  ARGS$longitude_col,
  ARGS$latitude_col
)

missing_cols <- setdiff(c(ARGS$individual_col, env6), colnames(env_raw))
if (length(missing_cols) > 0) {
  stop("Environmental table is missing columns: ", paste(missing_cols, collapse = ", "))
}

idx <- match(pruned_ind, env_raw[[ARGS$individual_col]])
if (any(is.na(idx))) {
  stop("Some individuals from genotype matrices are missing in the environmental table.")
}

env_ordered <- env_raw[idx, c(ARGS$individual_col, env6), drop = FALSE]

if (!identical(pruned_ind, env_ordered[[ARGS$individual_col]])) {
  stop("Environmental table could not be reordered to match genotype individual order.")
}

for (col in env6) {
  env_ordered[[col]] <- as.numeric(env_ordered[[col]])
  if (any(is.na(env_ordered[[col]]))) {
    stop("Environmental predictor contains NA after numeric coercion: ", col)
  }
}

# ---------------------------------------------------------------------
# Restricted permutation units. The units file assigns each individual
# to a breed-locality group and to an independent sampling locality;
# co-sampled breeds share a locality. Environmental values are permuted
# as intact packages across units at each level, with the LFMM2 latent
# factors held at their observed estimates.
# ---------------------------------------------------------------------
units_raw <- read.table(
  ARGS$units_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)
missing_unit_cols <- setdiff(
  c(ARGS$individual_col, "group_breed_locality", "locality"),
  colnames(units_raw)
)
if (length(missing_unit_cols) > 0) {
  stop("Units file is missing columns: ", paste(missing_unit_cols, collapse = ", "))
}

u_idx <- match(pruned_ind, units_raw[[ARGS$individual_col]])
if (any(is.na(u_idx))) {
  stop("Some individuals are missing in the units file.")
}
units_ordered <- units_raw[u_idx, , drop = FALSE]
if (!identical(pruned_ind, units_ordered[[ARGS$individual_col]])) {
  stop("Units file could not be reordered to match genotype individual order.")
}

# ---------------------------------------------------------------------
# Altitude and annual mean temperature are strongly correlated and are
# therefore modeled separately. The two models share the remaining four
# environmental variables so that the altitude and temperature
# coefficients are directly comparable.
# ---------------------------------------------------------------------
models <- list(
  altitude = setdiff(env6, ARGS$temperature_col),
  temperature = setdiff(env6, ARGS$altitude_col)
)

per_env_pieces <- list()

for (focal in names(models)) {
  env_cols <- models[[focal]]
  env_mat <- scale(as.matrix(env_ordered[, env_cols, drop = FALSE]))
  colnames(env_mat) <- env_cols

  lfmm_model <- lfmm2(
    input = ARGS$pruned_lfmm,
    env = env_mat,
    K = ARGS$k
  )
  saveRDS(
    lfmm_model,
    file.path(ARGS$outdir, paste0("lfmm2_model_", focal, "_K", ARGS$k, ".rds"))
  )

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
    focal_environmental_variable = focal,
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

  # -----------------------------------------------------------------
  # Restricted permutation tests at the breed-locality and locality
  # levels. Statistics are two-sided absolute z scores; empirical
  # P values follow (1 + #{T_perm >= T_obs}) / (B + 1).
  # -----------------------------------------------------------------
  for (level in c("group_breed_locality", "locality")) {
    unit_var <- units_ordered[[level]]
    uids <- unique(unit_var)
    ui <- match(unit_var, uids)
    pkg <- as.matrix(env_ordered[match(uids, unit_var), env_cols, drop = FALSE])
    if (anyNA(pkg)) {
      stop("Permutation units contain missing environmental values at level: ", level)
    }
    # each unit must share one environmental package within itself
    for (u in seq_along(uids)) {
      rows <- which(ui == u)
      if (nrow(unique(as.matrix(env_ordered[rows, env_cols, drop = FALSE]))) != 1) {
        stop("Environmental values are not constant within unit: ", uids[u])
      }
    }

    set.seed(ARGS$seed + match(focal, names(models)) * 101L +
             match(level, c("group_breed_locality", "locality")) * 7919L)

    z_perm <- array(NA_real_, dim = c(ARGS$n_perm, length(candidate_snp), length(env_cols)))
    for (b in seq_len(ARGS$n_perm)) {
      perm_map <- sample.int(length(uids))
      env_b <- pkg[perm_map[ui], , drop = FALSE]
      env_perm <- as.matrix(scale(env_b))
      colnames(env_perm) <- env_cols
      pv_perm <- lfmm2.test(
        object = lfmm_model,
        input = ARGS$candidate_lfmm,
        env = env_perm,
        linear = TRUE,
        genomic.control = FALSE,
        full = FALSE
      )
      z_perm[b, , ] <- coerce_to_snp_by_env(
        pv_perm$zscores,
        n_snp = length(candidate_snp),
        n_env = length(env_cols),
        what = "permutation zscore matrix"
      )
    }

    p_perm <- rep(NA_real_, nrow(result))
    env_pos <- match(result$environmental_variable, env_cols)
    snp_pos <- match(result$variant_id, candidate_snp)
    for (i in seq_len(nrow(result))) {
      zp <- abs(z_perm[, snp_pos[i], env_pos[i]])
      p_perm[i] <- (1 + sum(zp >= abs(result$z_score[i]))) / (ARGS$n_perm + 1)
    }
    result[[paste0("P_perm_", level)]] <- p_perm
  }

  per_env_pieces[[focal]] <- result
}

per_env <- do.call(rbind, per_env_pieces)

write.table(
  per_env,
  file.path(ARGS$outdir, "phenotype_associated_ASE_variant_environment_association.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

summary_df <- aggregate(
  significant_FDR_0_1 ~ focal_environmental_variable + environmental_variable,
  data = per_env,
  FUN = sum
)
summary_df$tested_variant_n <- aggregate(
  variant_id ~ focal_environmental_variable + environmental_variable,
  data = per_env,
  FUN = length
)$variant_id
summary_df <- summary_df[, c(
  "focal_environmental_variable",
  "environmental_variable",
  "tested_variant_n",
  "significant_FDR_0_1"
)]
colnames(summary_df)[4] <- "significant_variant_n"

write.table(
  summary_df,
  file.path(ARGS$outdir, "phenotype_associated_ASE_variant_environment_association.summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
