# Polyester simulation for ASE reference-bias benchmarking

Paired-end 150-bp reads are simulated from cattle/yak/gayal Ref/Alt haplotype
transcripts at known allelic ratios and processed through the production ASE
pipeline (STAR + WASP + ASEReadCounter) to measure residual reference bias
and detection power.

Requires Python 3, R with polyester (tested 1.38.0), samtools, bcftools,
and minimap2. Run per species:

    # 1. per-cycle error model from a representative BAM
    python 01_build_error_model.py sample.bam -o error_model -s 0.03 --readlen 150

    # 2. pick target heterozygous SNPs and build Ref/Alt haplotype FASTA;
    #    the cigar map comes from a splice alignment of the species
    #    transcripts, e.g.
    #    minimap2 -ax splice -uf -N 10 --secondary=no ref.fa transcripts.fa
    python 02_select_targets_build_haplotypes.py --species yak \
        --transcripts yak.transcripts.fa --cigar-map yak.cigar_map.json.gz \
        --vcf sample.het.snp.vcf.gz --genome ref.fa --outdir 02.haplotypes \
        --n-targets 120 [--require-native-match]   # cattle only

    # 3. fragment count matrices per haplotype and replicate
    python 03_build_readmat.py --haplotypes 02.haplotypes/yak.haplotypes.fa \
        --species yak --species-code 2 --outdir 04.truth_tables

    # 4-5. simulate reads (source the patch from this directory first)
    PIPE_DIR=. Rscript 05_simulate_reads.R yak 02.haplotypes/yak.haplotypes.fa \
        04.truth_tables error_model empirical150 03.simulated_reads 263 30 150 2001

    # 6. verify truth: fragment counts, SNP-covering reads, pairing
    python 06_verify_truth.py --simroot 03.simulated_reads \
        --target-loci 02.haplotypes/cattle.target_loci.tsv \
        --seeds 04.truth_tables/seeds.cattle.tsv --outdir 04.truth_tables

    # 7. FASTA to FASTQ (Q40 placeholder; the errors live in the sequences)
    python 07_fasta_to_fastq.py --simroot 03.simulated_reads --outdir 01.fastq

04_patch_polyester_150.R lifts polyester's 101-bp read-length guard so that
custom error models work at 150 bp; nothing else in the function is touched.
The 70:30 ratio is fixed per haplotype in the count matrix; the realized
SNP-covering ratios reported by step 6 are the truth values for any
downstream benchmark.
