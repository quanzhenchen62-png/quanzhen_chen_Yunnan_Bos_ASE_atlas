#!/bin/bash
# Example workflow order (paths are placeholders).
set -eo pipefail

# 1. empirical per-cycle error model from a representative BAM (PE150, MD tags)
python 01_build_error_model.py sample.bam -o error_model -s 0.03 --readlen 150

# 2. per species: target het-SNP loci + Ref/Alt haplotype FASTA
#    (cigar map from splice-aligning species transcripts to the reference genome)
python 02_select_targets_build_haplotypes.py --species cattle \
  --transcripts cattle.transcripts.fa --cigar-map cattle.cigar_map.json.gz \
  --vcf sample.het.snp.vcf.gz --genome ref.fa --outdir 02.haplotypes \
  --n-targets 120 --require-native-match
python 02_select_targets_build_haplotypes.py --species yak    --transcripts yak.transcripts.fa    --cigar-map yak.cigar_map.json.gz    --vcf yak.het.snp.vcf.gz    --genome ref.fa --outdir 02.haplotypes --n-targets 120
python 02_select_targets_build_haplotypes.py --species gayal  --transcripts gayal.transcripts.fa  --cigar-map gayal.cigar_map.json.gz  --vcf gayal.het.snp.vcf.gz  --genome ref.fa --outdir 02.haplotypes --n-targets 120

# 3. count matrices (fixed fragments per haplotype per replicate)
python 03_build_readmat.py --haplotypes 02.haplotypes/cattle.haplotypes.fa --species cattle --species-code 1 --outdir 04.truth_tables
python 03_build_readmat.py --haplotypes 02.haplotypes/yak.haplotypes.fa    --species yak    --species-code 2 --outdir 04.truth_tables
python 03_build_readmat.py --haplotypes 02.haplotypes/gayal.haplotypes.fa  --species gayal  --species-code 3 --outdir 04.truth_tables

# 4-5. simulate paired-end 150-bp reads (fragment lengths from real libraries)
PIPE_DIR=. Rscript 05_simulate_reads.R cattle 02.haplotypes/cattle.haplotypes.fa 04.truth_tables error_model empirical150 03.simulated_reads 238 30 150 1001
PIPE_DIR=. Rscript 05_simulate_reads.R yak    02.haplotypes/yak.haplotypes.fa    04.truth_tables error_model empirical150 03.simulated_reads 263 30 150 2001
PIPE_DIR=. Rscript 05_simulate_reads.R gayal  02.haplotypes/gayal.haplotypes.fa  04.truth_tables error_model empirical150 03.simulated_reads 269 30 150 3001

# 6. truth verification
python 06_verify_truth.py --simroot 03.simulated_reads \
  --target-loci 02.haplotypes/cattle.target_loci.tsv 02.haplotypes/yak.target_loci.tsv 02.haplotypes/gayal.target_loci.tsv \
  --seeds 04.truth_tables/seeds.cattle.tsv 04.truth_tables/seeds.yak.tsv 04.truth_tables/seeds.gayal.tsv \
  --outdir 04.truth_tables

# 7. FASTA -> FASTQ (uniform Q40), then analyze with the identical ASE pipeline
python 07_fasta_to_fastq.py --simroot 03.simulated_reads --outdir 01.fastq
