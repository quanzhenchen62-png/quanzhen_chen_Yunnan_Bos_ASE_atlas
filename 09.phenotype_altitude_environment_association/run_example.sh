#!/bin/bash

set -euo pipefail

# Example 1: phenotype-environment association (LMM + cluster-level permutation)
Rscript phenotype_environment_association.R \
  --phe-info data/phe.info \
  --lmm-info data/LMM.info \
  --outdir results/phenotype_environment_association \
  --gcta-bin gcta64 \
  --cor-threshold 0.8 \
  --b 1000 \
  --seed 1

# Example 2: identify phenotype-associated ASE variants
python3 phenotype_genotype_association.py \
  --ld-ase-table data/ld_brd_ASE3_candidates.tsv \
  --hm-ase-table data/hm_brd_ASE3_candidates.tsv \
  --ld-bfile data/LD_Bos \
  --hm-bfile data/HM_cattle \
  --ld-diameter-grm data/LD_diameter_GRM \
  --ld-density-grm data/LD_density_GRM \
  --hm-relative-weight-grm data/HM_relative_hump_weight_GRM \
  --ld-diameter-pheno data/pheno_LD_diameter.txt \
  --ld-density-pheno data/pheno_LD_density.txt \
  --hm-relative-weight-pheno data/pheno_relative_hump_weight.txt \
  --ld-covar data/covar_LD.txt \
  --ld-qcovar data/qcovar_LD.txt \
  --hm-covar data/covar_HM.txt \
  --hm-qcovar data/qcovar_HM.txt \
  --outdir results/phenotype_associated_ASE_variants \
  --gcta-bin gcta64 \
  --plink-bin plink

# Example 3: environmental association of phenotype-associated ASE variants
# Altitude and annual mean temperature are modeled separately because they
# are strongly correlated; both focal models include the remaining four
# environmental variables. Significance is assessed by BH-adjusted FDR and
# by restricted permutation tests at the breed-locality and locality levels
# (environmental packages permuted across units, latent factors held fixed).
Rscript genotype_environment_association.R \
  --candidate-lfmm data/GENO_candidate.lfmm \
  --candidate-snp-info data/GENO_candidate.snp.info \
  --candidate-ind-order data/GENO_candidate.ind.order \
  --pruned-lfmm data/Bos_pruned.lfmm \
  --pruned-ind-order data/Bos_pruned.ind.order \
  --env-table data/Bos_ind_env.info \
  --units-file data/panel_units.tsv \
  --outdir results/phenotype_associated_ASE_variant_environment_association \
  --individual-col individual \
  --altitude-col altitude \
  --temperature-col temperature_annual_mean \
  --precipitation-col precipitation_annual \
  --uvb-col UVB_annual_mean \
  --longitude-col longitude \
  --latitude-col latitude \
  --k 6 \
  --n-perm 1000 \
  --seed 20260924
