#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(sommer))

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  get_val <- function(key, default = NULL) {
    idx <- match(key, args)
    if (is.na(idx)) return(default)
    args[idx + 1]
  }
  list(
    phe_info = get_val("--phe-info", "data/phe.info"),
    lmm_info = get_val("--lmm-info", "data/LMM.info"),
    outdir   = get_val("--outdir", "results/phenotype_environment_association"),
    gcta_bin = get_val("--gcta-bin", "gcta64"),
    cor_threshold = as.numeric(get_val("--cor-threshold", "0.8")),
    n_pca_compute = as.integer(get_val("--n-pca-compute", "20")),
    autosome_num  = as.integer(get_val("--autosome-num", "29")),
    B   = as.integer(get_val("--b", "1000")),
    seed = as.integer(get_val("--seed", "1"))
  )
}

ARGS <- parse_args()
set.seed(ARGS$seed)
dir.create(ARGS$outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(ARGS$outdir, "genoPC_tmp"), recursive = TRUE, showWarnings = FALSE)

BUFFALO_BREEDS <- c("BLJSN", "DDNSN", "DHSN", "YJSN")
YAK_BREED <- "ZDMN"
MIN_SAMPLE_SIZE <- 10

ENV_GEO <- c("altitude_km", "temperature_annual_mean", "precipitation_annual",
             "UVB_annual_mean", "longitude", "latitude")
COV_PRIORITY <- c(precipitation_annual = 1, UVB_annual_mean = 1,
                  longitude = 2, latitude = 2, Age = 3)

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
    out <- out[!grepl("^DL\\.", out$individual), , drop = FALSE]
  }
  out
}

read_gcta_grm <- function(grm_prefix) {
  id_files <- list.files(dirname(grm_prefix), pattern = "\\.grm\\.id$", full.names = TRUE)
  if (length(id_files) == 0) stop("No .grm.id found in: ", dirname(grm_prefix))
  prefix <- sub("\\.grm\\.id$", "", id_files[1])
  ids <- read.table(paste0(prefix, ".grm.id"), stringsAsFactors = FALSE)
  n <- nrow(ids)
  con <- file(paste0(prefix, ".grm.bin"), "rb")
  vals <- readBin(con, what = "numeric", n = n * (n + 1) / 2, size = 4)
  close(con)
  G <- matrix(0, n, n)
  k <- 1
  for (i in seq_len(n)) for (j in seq_len(i)) {
    G[i, j] <- vals[k]; G[j, i] <- vals[k]; k <- k + 1
  }
  dimnames(G) <- list(ids[, 2], ids[, 2])
  G
}

compute_geno_pc <- function(grm_prefix, keep_ids, tag) {
  grm <- read_gcta_grm(grm_prefix)
  ids <- intersect(keep_ids, rownames(grm))
  G <- grm[ids, ids, drop = FALSE]
  G <- (G + t(G)) / 2
  eig <- eigen(G, symmetric = TRUE, only.values = TRUE)
  if (min(eig$values) < 1e-6) {
    eig <- eigen(G, symmetric = TRUE)
    eig$values[eig$values < 1e-6] <- 1e-6
    G <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
    G <- (G + t(G)) / 2
    dimnames(G) <- list(ids, ids)
  }
  ev <- eigen(G, symmetric = TRUE)$vectors[, 1:min(ARGS$n_pca_compute, length(ids) - 1),
                                           drop = FALSE]
  pc <- as.data.frame(ev * sqrt(nrow(G)))
  names(pc) <- paste0("geno_PC", seq_len(ncol(pc)))
  cbind(data.frame(individual = ids, stringsAsFactors = FALSE), pc)
}

filter_collinear <- function(df, candidates, focal, threshold) {
  keep <- candidates
  cc <- complete.cases(df[, c(focal, candidates)])
  r_f <- sapply(candidates, function(v)
    abs(cor(df[[v]][cc], df[[focal]][cc], use = "complete.obs")))
  keep <- keep[is.na(r_f) | r_f < threshold]
  repeat {
    if (length(keep) < 2) break
    cc2 <- complete.cases(df[, keep])
    cm <- cor(df[cc2, keep, drop = FALSE], use = "complete.obs")
    worst <- 0; pair <- NULL
    for (i in seq_len(length(keep) - 1)) for (j in (i + 1):length(keep)) {
      r0 <- abs(cm[i, j])
      if (!is.na(r0) && r0 >= threshold && r0 > worst) {
        worst <- r0; pair <- c(keep[i], keep[j])
      }
    }
    if (is.null(pair)) break
    p1 <- COV_PRIORITY[[pair[1]]]; p2 <- COV_PRIORITY[[pair[2]]]
    if (is.na(p1)) p1 <- 99
    if (is.na(p2)) p2 <- 99
    drop_v <- if (p1 > p2) pair[1] else if (p2 > p1) pair[2] else pair[2]
    keep <- setdiff(keep, drop_v)
  }
  keep
}

phe <- read.table(ARGS$phe_info, header = TRUE, sep = "\t",
                  stringsAsFactors = FALSE, check.names = FALSE,
                  quote = "", comment.char = "")
phe$individual <- trimws(as.character(phe$individual))
phe$Breed <- trimws(as.character(phe$Breed))
phe$site <- paste(phe$longitude, phe$latitude, sep = "_")
phe$breed_locality <- paste(phe$Breed, phe$site, sep = "::")
for (v in c("Age", "altitude", "temperature_annual_mean", "precipitation_annual",
            "UVB_annual_mean", "longitude", "latitude")) {
  phe[[v]] <- suppressWarnings(as.numeric(trimws(as.character(phe[[v]]))))
}
phe$altitude_km <- phe$altitude / 1000
phe$Sex <- factor(trimws(as.character(phe$Sex)))

lmm_info <- read.table(ARGS$lmm_info, header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE, check.names = FALSE)

focal_models <- list(
  altitude = list(focal = "altitude_km"),
  temperature = list(focal = "temperature_annual_mean")
)

all_rows <- list(); focal_rows <- list(); meta_rows <- list(); log_rows <- list()

for (i in seq_len(nrow(lmm_info))) {
  cls <- lmm_info$class[i]
  phe_col <- lmm_info$phe[i]
  grm_dir <- lmm_info$GRM[i]
  parsed <- parse_class(cls)
  df0 <- filter_by_breed(phe, parsed$base_group, parsed$without_dl)
  df0 <- df0[!is.na(df0[[phe_col]]) & trimws(as.character(df0[[phe_col]])) != "", , drop = FALSE]
  df0[[phe_col]] <- as.numeric(trimws(as.character(df0[[phe_col]])))

  G_full <- read_gcta_grm(file.path(grm_dir, "grm"))
  common <- intersect(df0$individual, rownames(G_full))
  if (length(common) < MIN_SAMPLE_SIZE) {
    log_rows[[length(log_rows) + 1]] <- data.frame(class = cls, phe = phe_col,
      model = NA, status = "FAILED", note = "too few individuals")
    next
  }
  df0 <- df0[match(common, df0$individual), , drop = FALSE]
  G <- G_full[common, common, drop = FALSE]
  G <- (G + t(G)) / 2
  eig <- eigen(G, symmetric = TRUE)
  if (min(eig$values) < 1e-6) {
    eig$values[eig$values < 1e-6] <- 1e-6
    G <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
    G <- (G + t(G)) / 2
    dimnames(G) <- list(common, common)
  }

  n_pc <- decide_n_geno_pc(cls)
  pc_tab <- compute_geno_pc(file.path(grm_dir, "grm"), common, cls)
  gpc <- paste0("geno_PC", seq_len(min(n_pc, ncol(pc_tab) - 1)))

  df <- data.frame(individual = df0$individual, Breed = df0$Breed,
                   site = df0$site, breed_locality = df0$breed_locality,
                   phe = df0[[phe_col]], Sex = df0$Sex, Age = df0$Age,
                   stringsAsFactors = FALSE)
  for (v in ENV_GEO) df[[v]] <- df0[[v]]
  df <- cbind(df, pc_tab[match(df$individual, pc_tab$individual),
                         gpc, drop = FALSE])

  for (model_name in names(focal_models)) {
    focal <- focal_models[[model_name]]$focal
    cand <- c("Age", "longitude", "latitude", "precipitation_annual",
              "UVB_annual_mean", gpc)
    cand <- cand[sapply(cand, function(t) length(unique(df[[t]])) >= 2)]
    cand <- filter_collinear(df, cand, focal, ARGS$cor_threshold)
    fixed_terms <- c(focal, "Sex", cand)

    fit_once <- function(d) {
      tryCatch({
        d$individual <- factor(d$individual, levels = rownames(G))
        m <- mmer(as.formula(paste("phe ~", paste(fixed_terms, collapse = " + "))),
                  random = ~ vs(individual, Gu = G), rcov = ~ units,
                  data = d, verbose = FALSE)
        bt <- summary(m)$betas
        rn <- if ("Effect" %in% colnames(bt)) bt$Effect else rownames(bt)
        k <- which(rn == focal)
        if (!length(k)) return(c(NA_real_, NA_real_))
        c(bt$Estimate[k], bt$Std.Error[k])
      }, error = function(e) c(NA_real_, NA_real_))
    }

    obs <- fit_once(df)
    if (any(is.na(obs))) {
      log_rows[[length(log_rows) + 1]] <- data.frame(class = cls, phe = phe_col,
        model = model_name, status = "FAILED", note = "model fit failed")
      next
    }
    t_obs <- (obs[1] / obs[2])^2
    p_wald <- pchisq(t_obs, df = 1, lower.tail = FALSE)

    g19_tab <- unique(df[, c("breed_locality", ENV_GEO)])
    s17_tab <- unique(df[, c("site", ENV_GEO)])
    p_perm <- c(g19 = NA_real_, s17 = NA_real_); n_ok <- c(g19 = 0L, s17 = 0L)

    for (lvl in c("g19", "s17")) {
      chis <- rep(NA_real_, ARGS$B)
      for (b in seq_len(ARGS$B)) {
        if (lvl == "g19") {
          lab <- sample(g19_tab$breed_locality)
          pt <- g19_tab
          pt[, ENV_GEO] <- g19_tab[match(lab, g19_tab$breed_locality), ENV_GEO, drop = FALSE]
          dp <- df
          dp[, ENV_GEO] <- pt[match(df$breed_locality, pt$breed_locality), ENV_GEO, drop = FALSE]
        } else {
          lab <- sample(s17_tab$site)
          pt <- s17_tab
          pt[, ENV_GEO] <- s17_tab[match(lab, s17_tab$site), ENV_GEO, drop = FALSE]
          dp <- df
          dp[, ENV_GEO] <- pt[match(df$site, pt$site), ENV_GEO, drop = FALSE]
        }
        r <- fit_once(dp)
        if (!any(is.na(r))) chis[b] <- (r[1] / r[2])^2
      }
      n_ok[lvl] <- sum(!is.na(chis))
      p_perm[lvl] <- (1 + sum(chis >= t_obs, na.rm = TRUE)) / (1 + n_ok[lvl])
    }

    meta_rows[[length(meta_rows) + 1]] <- data.frame(
      class = cls, phe = phe_col, model = model_name,
      n_individuals = nrow(df),
      n_breed_locality_groups = nrow(g19_tab),
      n_unique_localities = nrow(s17_tab),
      n_geno_pc = length(gpc),
      fixed_terms = paste(fixed_terms, collapse = ","),
      B = ARGS$B, seed = ARGS$seed)
    focal_rows[[length(focal_rows) + 1]] <- data.frame(
      class = cls, phe = phe_col, model = model_name, focal = focal,
      Estimate = obs[1], SE = obs[2],
      Wald_chisq = t_obs, Wald_P = p_wald,
      P_perm_breed_locality = p_perm["g19"],
      P_perm_locality = p_perm["s17"],
      B_valid_g19 = n_ok["g19"], B_valid_s17 = n_ok["s17"])
    log_rows[[length(log_rows) + 1]] <- data.frame(class = cls, phe = phe_col,
      model = model_name, status = "OK", note = "")
  }
}

wr <- function(x, f) if (length(x)) write.table(do.call(rbind, x), f,
  sep = "\t", row.names = FALSE, quote = FALSE)
wr(focal_rows, file.path(ARGS$outdir, "wald_results_focal_predictors.tsv"))
wr(meta_rows, file.path(ARGS$outdir, "model_metadata.tsv"))
wr(log_rows, file.path(ARGS$outdir, "run_log.tsv"))
cat("Done.\n")
