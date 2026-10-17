#!/usr/bin/env Rscript
## Simulate paired-end reads with Polyester simulate_experiment_countmat().
## Usage: Rscript 05_simulate_reads.R <species> <haplotype_fasta> <readmat_dir>
##   <error_model_dir> <error_model_prefix> <outroot> <fraglen> <fragsd> <readlen> <seed_base>
suppressMessages(library(polyester))
source(file.path(Sys.getenv("PIPE_DIR"), "04_patch_polyester_150.R"))

a <- commandArgs(trailingOnly = TRUE)
stopifnot(length(a) == 10)
sp <- a[1]; hapfa <- a[2]; rmdir <- a[3]; emdir <- a[4]; empre <- a[5]
outroot <- a[6]; fraglen <- as.numeric(a[7]); fragsd <- as.numeric(a[8])
readlen <- as.integer(a[9]); seed0 <- as.integer(a[10])

rmfiles <- list.files(rmdir, pattern = sprintf("^readmat_.*\\.%s\\.tsv$", sp), full.names = TRUE)
for (rmf in rmfiles) {
  ratio <- sub(sprintf("^readmat_(.*)\\.%s\\.tsv$", sp), "\\1", basename(rmf))
  rmdf <- read.delim(rmf, row.names = 1, check.names = FALSE)
  readmat <- as.matrix(rmdf); storage.mode(readmat) <- "integer"
  for (k in seq_len(ncol(readmat))) {
    rep <- colnames(readmat)[k]
    seed <- seed0 + k
    repdir <- file.path(outroot, sp, paste0("ratio", ratio), rep)
    dir.create(repdir, recursive = TRUE, showWarnings = FALSE)
    t0 <- Sys.time()
    simulate_experiment_countmat(
      fasta = hapfa, readmat = readmat[, k, drop = FALSE],
      outdir = repdir, paired = TRUE, readlen = readlen,
      fraglen = fraglen, fragsd = fragsd, distr = "normal",
      error_model = "custom", model_path = emdir, model_prefix = empre,
      bias = "none", gzip = TRUE, seed = seed)
    base <- file.path(repdir, sprintf("%s_%s_%s", sp, ratio, rep))
    file.rename(file.path(repdir, "sample_01_1.fasta.gz"), paste0(base, "_1.fasta.gz"))
    file.rename(file.path(repdir, "sample_01_2.fasta.gz"), paste0(base, "_2.fasta.gz"))
    for (mate in c(1, 2)) {
      plain <- file.path(repdir, sprintf("sample_01_%d.fasta", mate))
      if (file.exists(plain)) {
        system(sprintf("mv '%s' '%s_%d.fasta' && gzip -9 '%s_%d.fasta'", plain, base, mate, base, mate))
      }
    }
    cat(sprintf("[%s %s %s seed=%d] %.1f sec\n", sp, ratio, rep, seed,
                as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }
}
