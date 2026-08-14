#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(sommer))

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  defaults <- list(
    lmm_info = "data/LMM.info",
    phe_info = "data/phe.info",
    outdir = "results/phenotype_altitude_environment_association",
    gcta_bin = "gcta64",
    cor_threshold = 0.8,
    n_pca_compute = 20,
    autosome_num = 29
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

  defaults$cor_threshold <- as.numeric(defaults$cor_threshold)
  defaults$n_pca_compute <- as.integer(defaults$n_pca_compute)
  defaults$autosome_num <- as.integer(defaults$autosome_num)

  defaults
}

ARGS <- parse_args()

BUFFALO_BREEDS <- c("BLJSN", "DDNSN", "DHSN", "YJSN")
YAK_BREED <- "ZDMN"
DL_BREED <- "DL"
MIN_SAMPLE_SIZE <- 10

PCA_CACHE <- new.env(parent = emptyenv())

dir.create(ARGS$outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(ARGS$outdir, "diagnostics"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(ARGS$outdir, "genoPC_tmp"), showWarnings = FALSE, recursive = TRUE)

VAR_PRIORITY <- c(
  altitude_km = 1,
  temperature_annual_mean = 1,
  precipitation_annual = 2,
  UVB_annual_mean = 2,
  longitude = 3,
  latitude = 3,
  Age = 4,
  geno_PC1 = 5,
  geno_PC2 = 5,
  geno_PC3 = 5,
  geno_PC4 = 5,
  geno_PC5 = 5,
  geno_PC6 = 5,
  geno_PC7 = 5,
  geno_PC8 = 5,
  geno_PC9 = 5,
  geno_PC10 = 5
)

stop_if_missing_columns <- function(df, cols, label) {
  missing_cols <- setdiff(cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop(label, " is missing columns: ", paste(missing_cols, collapse = ", "))
  }
}

decide_n_geno_pc <- function(class_str) {
  if (grepl("^(HM|JF)", class_str, ignore.case = TRUE)) return(1)
  if (grepl("all", class_str, ignore.case = TRUE)) return(6)
  if (grepl("Bos", class_str, ignore.case = TRUE)) return(3)
  if (grepl("Cattle", class_str, ignore.case = TRUE)) return(2)
  6
}

parse_class <- function(class_str) {
  without_dl <- grepl("without_DL", class_str, ignore.case = TRUE)

  if (grepl("Cattle", class_str, ignore.case = TRUE)) {
    base_group <- "Cattle"
  } else if (grepl("Bos", class_str, ignore.case = TRUE)) {
    base_group <- "Bos"
  } else {
    base_group <- "all"
  }

  list(base_group = base_group, without_dl = without_dl)
}

filter_by_breed <- function(df, base_group, without_dl) {
  out <- df

  if (base_group %in% c("Bos", "Cattle")) {
    out <- out[!(out$Breed %in% BUFFALO_BREEDS), , drop = FALSE]
  }

  if (base_group == "Cattle") {
    out <- out[out$Breed != YAK_BREED, , drop = FALSE]
  }

  if (without_dl) {
    out <- out[out$Breed != DL_BREED, , drop = FALSE]
  }

  out
}

read_gcta_grm <- function(grm_prefix) {
  id_file <- paste0(grm_prefix, ".grm.id")
  bin_file <- paste0(grm_prefix, ".grm.bin")

  ids <- read.table(id_file, stringsAsFactors = FALSE)
  n <- nrow(ids)

  con <- file(bin_file, "rb")
  vals <- readBin(con, what = "numeric", n = n * (n + 1) / 2, size = 4)
  close(con)

  G <- matrix(0, n, n)
  idx <- 1
  for (i in seq_len(n)) {
    for (j in seq_len(i)) {
      G[i, j] <- vals[idx]
      G[j, i] <- vals[idx]
      idx <- idx + 1
    }
  }

  rownames(G) <- ids[, 2]
  colnames(G) <- ids[, 2]
  G
}

guess_grm_prefix <- function(grm_dir) {
  id_files <- list.files(grm_dir, pattern = "\\.grm\\.id$", full.names = TRUE)
  if (length(id_files) == 0) {
    stop("No .grm.id file found in: ", grm_dir)
  }
  sub("\\.grm\\.id$", "", id_files[1])
}

guess_bed_prefix <- function(grm_prefix) {
  pruned <- paste0(grm_prefix, "_pruned")
  if (file.exists(paste0(pruned, ".bed"))) return(pruned)
  if (file.exists(paste0(grm_prefix, ".bed"))) return(grm_prefix)

  beds <- list.files(dirname(grm_prefix), pattern = "_pruned\\.bed$", full.names = TRUE)
  if (length(beds) > 0) {
    return(sub("\\.bed$", "", beds[1]))
  }

  stop("No matching BED file found for GRM prefix: ", grm_prefix)
}

make_pca_key <- function(bed_prefix, keep_ids) {
  paste0(bed_prefix, "||", paste(sort(unique(keep_ids)), collapse = ";"))
}

compute_geno_pc <- function(bed_prefix, keep_ids, tag) {
  key <- make_pca_key(bed_prefix, keep_ids)
  if (exists(key, envir = PCA_CACHE, inherits = FALSE)) {
    return(get(key, envir = PCA_CACHE, inherits = FALSE))
  }

  safe_tag <- gsub("[^A-Za-z0-9_]+", "_", tag)
  tmp_dir <- file.path(ARGS$outdir, "genoPC_tmp")
  keep_file <- file.path(tmp_dir, paste0(safe_tag, ".keep.txt"))
  grm_out <- file.path(tmp_dir, paste0(safe_tag, ".subgrm"))
  pca_out <- file.path(tmp_dir, paste0(safe_tag, ".pca"))

  fam <- read.table(paste0(bed_prefix, ".fam"), stringsAsFactors = FALSE)
  keep_fam <- fam[fam[, 2] %in% keep_ids, c(1, 2), drop = FALSE]

  if (nrow(keep_fam) < MIN_SAMPLE_SIZE) {
    stop("Too few matched individuals in FAM: ", nrow(keep_fam))
  }

  write.table(
    keep_fam,
    keep_file,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )

  n_compute <- min(ARGS$n_pca_compute, nrow(keep_fam) - 1)

  cmd1 <- sprintf(
    "%s --bfile %s --keep %s --autosome --autosome-num %d --make-grm --out %s --thread-num 4",
    shQuote(ARGS$gcta_bin),
    shQuote(bed_prefix),
    shQuote(keep_file),
    ARGS$autosome_num,
    shQuote(grm_out)
  )

  cmd2 <- sprintf(
    "%s --grm %s --pca %d --out %s --thread-num 4",
    shQuote(ARGS$gcta_bin),
    shQuote(grm_out),
    n_compute,
    shQuote(pca_out)
  )

  if (system(cmd1, ignore.stdout = TRUE, ignore.stderr = FALSE) != 0) {
    stop("GCTA --make-grm failed")
  }
  if (system(cmd2, ignore.stdout = TRUE, ignore.stderr = FALSE) != 0) {
    stop("GCTA --pca failed")
  }

  ev <- read.table(paste0(pca_out, ".eigenvec"), stringsAsFactors = FALSE)
  pc_n <- ncol(ev) - 2
  colnames(ev) <- c("FID", "individual", paste0("geno_PC", seq_len(pc_n)))
  ev$FID <- NULL

  result <- list(pc_df = ev, n_pc = pc_n)
  assign(key, result, envir = PCA_CACHE)
  result
}

stabilize_grm <- function(G, eig_floor = 1e-6) {
  G <- (G + t(G)) / 2
  eig <- eigen(G, symmetric = TRUE)
  if (min(eig$values) < eig_floor) {
    eig$values[eig$values < eig_floor] <- eig_floor
    G <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
    G <- (G + t(G)) / 2
  }
  G
}

resolve_collinearity <- function(df, candidate_vars, focal_var, threshold = 0.8) {
  vars <- candidate_vars[candidate_vars %in% colnames(df)]

  if (!(focal_var %in% vars)) {
    stop("Focal variable not found: ", focal_var)
  }

  if (length(vars) < 2) return(vars)

  complete_df <- df[, vars, drop = FALSE]
  if (sum(complete.cases(complete_df)) < 5) return(vars)

  cor_mat <- cor(complete_df, use = "complete.obs")

  focal_cor <- abs(cor_mat[focal_var, ])
  focal_cor[focal_var] <- 0
  removed <- names(focal_cor)[focal_cor >= threshold]
  remaining <- setdiff(vars, c(focal_var, removed))

  if (length(remaining) >= 2) {
    repeat {
      if (length(remaining) < 2) break

      sub_cor <- cor_mat[remaining, remaining, drop = FALSE]
      pair <- NULL
      max_r <- 0

      for (i in seq_len(length(remaining) - 1)) {
        for (j in (i + 1):length(remaining)) {
          r <- abs(sub_cor[i, j])
          if (!is.na(r) && r >= threshold && r > max_r) {
            max_r <- r
            pair <- c(remaining[i], remaining[j])
          }
        }
      }

      if (is.null(pair)) break

      p1 <- VAR_PRIORITY[pair[1]]
      p2 <- VAR_PRIORITY[pair[2]]
      if (is.na(p1)) p1 <- 99
      if (is.na(p2)) p2 <- 99

      drop_var <- if (p1 > p2) {
        pair[1]
      } else if (p2 > p1) {
        pair[2]
      } else {
        pair[2]
      }

      remaining <- setdiff(remaining, drop_var)
    }
  }

  c(focal_var, remaining)
}

count_unique_non_na <- function(x) {
  length(unique(x[!is.na(x)]))
}

drop_non_estimable_terms <- function(df, terms) {
  keep <- c()
  dropped <- c()

  for (term in terms) {
    x <- df[[term]]
    n_unique <- count_unique_non_na(x)
    if (n_unique < 2) {
      dropped <- c(dropped, term)
    } else {
      keep <- c(keep, term)
    }
  }

  list(keep = keep, dropped = dropped)
}

fit_model <- function(df, G_sub, phe_col, fixed_terms, common_ids) {
  fixed_formula <- as.formula(
    paste0(phe_col, " ~ ", paste(fixed_terms, collapse = " + "))
  )

  df$individual <- factor(df$individual, levels = common_ids)

  model <- mmer(
    fixed = fixed_formula,
    random = ~ vsr(individual, Gu = G_sub),
    rcov = ~ units,
    data = df,
    verbose = FALSE
  )

  list(
    model = model,
    beta_table = summary(model)$betas,
    formula = fixed_formula
  )
}

extract_wald <- function(beta_table, cls, phe_col, model_name, focal_var,
                         n_individuals, terms_used, dropped_terms) {
  out <- beta_table
  if (!("Effect" %in% colnames(out))) {
    out$Effect <- rownames(out)
  }

  out$wald_chisq <- (out$Estimate / out$Std.Error) ^ 2
  out$p_value <- pchisq(out$wald_chisq, df = 1, lower.tail = FALSE)
  out$ci_lower <- out$Estimate - 1.96 * out$Std.Error
  out$ci_upper <- out$Estimate + 1.96 * out$Std.Error
  out$class <- cls
  out$phe <- phe_col
  out$model <- model_name
  out$focal_variable <- focal_var
  out$n_individuals <- n_individuals
  out$terms_used <- paste(terms_used, collapse = ",")
  out$dropped_terms <- paste(dropped_terms, collapse = ",")

  out <- out[, c(
    "class", "phe", "model", "focal_variable", "Effect",
    "Estimate", "Std.Error", "wald_chisq", "p_value",
    "ci_lower", "ci_upper", "n_individuals", "terms_used", "dropped_terms"
  )]

  rownames(out) <- NULL
  out
}

plot_residual_diagnostics <- function(model, out_prefix, title_text) {
  res <- tryCatch(as.numeric(unlist(residuals(model))), error = function(e) NULL)
  if (is.null(res)) {
    res <- as.numeric(unlist(model$residuals))
  }
  res <- res[is.finite(res)]

  if (length(res) < 5) {
    return(NULL)
  }

  write.table(
    data.frame(residual = res),
    paste0(out_prefix, ".residuals.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  png(paste0(out_prefix, ".png"), width = 1800, height = 900, res = 180)
  par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3.5, 1))

  hist(
    res,
    breaks = 20,
    probability = TRUE,
    main = paste0(title_text, "\nResidual histogram"),
    xlab = "Residuals",
    col = "grey80",
    border = "grey40"
  )

  curve(
    dnorm(x, mean = mean(res), sd = sd(res)),
    add = TRUE,
    lwd = 2,
    col = "#2C74B3"
  )

  qqnorm(res, main = paste0(title_text, "\nNormal Q-Q plot"), pch = 16, cex = 0.8)
  qqline(res, col = "#2C74B3", lwd = 2)
  dev.off()
}

required_phe_cols <- c(
  "individual", "Breed", "Sex", "Age", "altitude",
  "temperature_annual_mean", "precipitation_annual", "UVB_annual_mean",
  "longitude", "latitude"
)

phe_info <- read.table(ARGS$phe_info, header = TRUE, sep = "\t", check.names = FALSE)
lmm_info <- read.table(ARGS$lmm_info, header = TRUE, sep = "\t", check.names = FALSE)

stop_if_missing_columns(phe_info, required_phe_cols, "phe.info")
stop_if_missing_columns(lmm_info, c("class", "phe", "GRM"), "LMM.info")

phe_info$individual <- trimws(as.character(phe_info$individual))
phe_info$Breed <- trimws(as.character(phe_info$Breed))
phe_info$Sex <- factor(trimws(as.character(phe_info$Sex)))

for (v in c(
  "Age", "altitude", "temperature_annual_mean", "precipitation_annual",
  "UVB_annual_mean", "longitude", "latitude"
)) {
  phe_info[[v]] <- as.numeric(trimws(as.character(phe_info[[v]])))
}

phe_info$altitude_km <- phe_info$altitude / 1000

focal_models <- list(
  altitude_model = list(
    focal = "altitude_km",
    env_covariates = c("precipitation_annual", "UVB_annual_mean")
  ),
  temperature_model = list(
    focal = "temperature_annual_mean",
    env_covariates = c("precipitation_annual", "UVB_annual_mean", "altitude_km")
  )
)

all_results <- list()
focal_results <- list()
meta_results <- list()
run_logs <- list()

for (i in seq_len(nrow(lmm_info))) {
  cls <- lmm_info$class[i]
  phe_col <- lmm_info$phe[i]
  grm_dir <- lmm_info$GRM[i]

  cat("\nProcessing:", cls, "| phenotype:", phe_col, "\n")

  tryCatch({
    parsed <- parse_class(cls)
    df <- filter_by_breed(phe_info, parsed$base_group, parsed$without_dl)

    if (!(phe_col %in% colnames(df))) {
      stop("Phenotype column not found: ", phe_col)
    }

    df[[phe_col]] <- as.numeric(trimws(as.character(df[[phe_col]])))
    df <- df[!is.na(df[[phe_col]]), , drop = FALSE]

    if (nrow(df) < MIN_SAMPLE_SIZE) {
      stop("Too few non-missing phenotype records: ", nrow(df))
    }

    grm_prefix <- guess_grm_prefix(grm_dir)
    G <- read_gcta_grm(grm_prefix)

    common_ids <- intersect(df$individual, rownames(G))
    if (length(common_ids) < MIN_SAMPLE_SIZE) {
      stop("Too few individuals shared by phenotype and GRM: ", length(common_ids))
    }

    df <- df[df$individual %in% common_ids, , drop = FALSE]
    df <- df[match(common_ids, df$individual), , drop = FALSE]

    bed_prefix <- guess_bed_prefix(grm_prefix)
    pc_info <- compute_geno_pc(
      bed_prefix = bed_prefix,
      keep_ids = common_ids,
      tag = paste(cls, phe_col, sep = "__")
    )

    n_pc_use <- min(decide_n_geno_pc(cls), pc_info$n_pc)
    geno_pc_names <- paste0("geno_PC", seq_len(n_pc_use))
    pc_df <- pc_info$pc_df[, c("individual", geno_pc_names), drop = FALSE]

    df <- merge(df, pc_df, by = "individual", all.x = TRUE, sort = FALSE)
    df <- df[match(common_ids, df$individual), , drop = FALSE]

    for (model_name in names(focal_models)) {
      focal_var <- focal_models[[model_name]]$focal

      numeric_candidates <- c(
        focal_var,
        "Age",
        "longitude",
        "latitude",
        focal_models[[model_name]]$env_covariates,
        geno_pc_names
      )
      numeric_candidates <- unique(numeric_candidates)
      numeric_candidates <- numeric_candidates[numeric_candidates %in% colnames(df)]

      keep_numeric <- resolve_collinearity(
        df = df,
        candidate_vars = numeric_candidates,
        focal_var = focal_var,
        threshold = ARGS$cor_threshold
      )

      base_terms <- c("Sex", keep_numeric)
      vars_needed <- c(phe_col, base_terms)
      df_model <- df[complete.cases(df[, vars_needed, drop = FALSE]), , drop = FALSE]

      if (nrow(df_model) < MIN_SAMPLE_SIZE) {
        stop(paste(model_name, "has too few complete cases:", nrow(df_model)))
      }

      term_info <- drop_non_estimable_terms(df_model, base_terms)
      if (!(focal_var %in% term_info$keep)) {
        stop(paste("Focal predictor is not estimable in", model_name))
      }

      common_ids_model <- intersect(df_model$individual, rownames(G))
      df_model <- df_model[match(common_ids_model, df_model$individual), , drop = FALSE]
      G_sub <- stabilize_grm(G[common_ids_model, common_ids_model, drop = FALSE])

      fit <- fit_model(
        df = df_model,
        G_sub = G_sub,
        phe_col = phe_col,
        fixed_terms = term_info$keep,
        common_ids = common_ids_model
      )

      wald_df <- extract_wald(
        beta_table = fit$beta_table,
        cls = cls,
        phe_col = phe_col,
        model_name = model_name,
        focal_var = focal_var,
        n_individuals = nrow(df_model),
        terms_used = term_info$keep,
        dropped_terms = term_info$dropped
      )

      focal_row <- wald_df[wald_df$Effect == focal_var, , drop = FALSE]

      diag_prefix <- file.path(
        ARGS$outdir, "diagnostics",
        paste0(gsub("[^A-Za-z0-9_]+", "_", cls), "__",
               gsub("[^A-Za-z0-9_]+", "_", phe_col), "__", model_name)
      )
      plot_residual_diagnostics(fit$model, diag_prefix, paste(cls, phe_col, model_name, sep = " | "))

      all_results[[length(all_results) + 1]] <- wald_df
      focal_results[[length(focal_results) + 1]] <- focal_row
      meta_results[[length(meta_results) + 1]] <- data.frame(
        class = cls,
        phe = phe_col,
        model = model_name,
        focal_variable = focal_var,
        n_individuals = nrow(df_model),
        n_geno_pc = n_pc_use,
        fixed_terms = paste(term_info$keep, collapse = ","),
        dropped_terms = paste(term_info$dropped, collapse = ","),
        stringsAsFactors = FALSE
      )
      run_logs[[length(run_logs) + 1]] <- data.frame(
        class = cls,
        phe = phe_col,
        model = model_name,
        status = "OK",
        note = "",
        stringsAsFactors = FALSE
      )
    }
  }, error = function(e) {
    msg <- conditionMessage(e)
    cat("[ERROR]", msg, "\n")

    for (model_name in names(focal_models)) {
      run_logs[[length(run_logs) + 1]] <<- data.frame(
        class = cls,
        phe = phe_col,
        model = model_name,
        status = "FAILED",
        note = msg,
        stringsAsFactors = FALSE
      )
    }
  })
}

all_results_df <- if (length(all_results) > 0) do.call(rbind, all_results) else data.frame()
focal_results_df <- if (length(focal_results) > 0) do.call(rbind, focal_results) else data.frame()
meta_results_df <- if (length(meta_results) > 0) do.call(rbind, meta_results) else data.frame()
run_logs_df <- if (length(run_logs) > 0) do.call(rbind, run_logs) else data.frame()

write.table(
  all_results_df,
  file.path(ARGS$outdir, "wald_results_all_fixed_terms.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  focal_results_df,
  file.path(ARGS$outdir, "wald_results_focal_predictors.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  meta_results_df,
  file.path(ARGS$outdir, "model_metadata.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  run_logs_df,
  file.path(ARGS$outdir, "run_log.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat("Done.\n")
