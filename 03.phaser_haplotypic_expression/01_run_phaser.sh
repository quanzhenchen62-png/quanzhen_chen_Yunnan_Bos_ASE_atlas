#!/bin/bash

set -euo pipefail

if [[ $# -ne 8 ]]; then
  echo "Usage: $0 <sample_id> <het.vcf.gz> <wasp_filtered.bam> <phaser.py> <blacklist.bed> <threads> <paired_end:0|1> <output_prefix>" >&2
  exit 1
fi

sample_id="$1"
het_vcf="$2"
bam="$3"
phaser_py="$4"
blacklist_bed="$5"
threads="$6"
paired_end="$7"
output_prefix="$8"

if [[ ! -s "${bam}.bai" ]]; then
  samtools index -@ "${threads}" "${bam}"
fi

python "${phaser_py}" \
  --vcf "${het_vcf}" \
  --bam "${bam}" \
  --paired_end "${paired_end}" \
  --mapq 255 \
  --baseq 10 \
  --sample "${sample_id}" \
  --haplo_count_blacklist "${blacklist_bed}" \
  --threads "${threads}" \
  --o "${output_prefix}"
