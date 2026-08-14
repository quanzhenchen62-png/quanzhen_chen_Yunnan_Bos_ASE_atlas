#!/bin/bash

set -euo pipefail

if [[ $# -ne 7 ]]; then
  echo "Usage: $0 <sample_id> <wasp_filtered.bam> <het.vcf.gz> <reference.fa> <gatk_jar> <output_table.tsv> <temp_rg_dir>" >&2
  exit 1
fi

sample_id="$1"
bam="$2"
het_vcf="$3"
reference_fasta="$4"
gatk_jar="$5"
output_table="$6"
temp_rg_dir="$7"

mkdir -p "$(dirname "${output_table}")" "${temp_rg_dir}"

input_bam="${bam}"
rg_bam="${temp_rg_dir}/${sample_id}.WASP_filtered.RG.bam"

if [[ ! -s "${bam}.bai" ]]; then
  samtools index "${bam}"
fi

if ! samtools view -H "${bam}" | grep -q '^@RG'; then
  samtools addreplacerg \
    -O BAM \
    -r "ID:${sample_id}" \
    -r "SM:${sample_id}" \
    -r "LB:${sample_id}" \
    -r "PL:ILLUMINA" \
    -r "PU:${sample_id}" \
    -o "${rg_bam}" \
    "${bam}"
  samtools index "${rg_bam}"
  input_bam="${rg_bam}"
fi

java -Xmx30G -jar "${gatk_jar}" \
  -T ASEReadCounter \
  -R "${reference_fasta}" \
  -I "${input_bam}" \
  -sites "${het_vcf}" \
  -o "${output_table}" \
  -U ALLOW_N_CIGAR_READS \
  -minDepth 1 \
  -mmq 255 \
  -mbq 10

if [[ "${input_bam}" == "${rg_bam}" ]]; then
  rm -f "${rg_bam}" "${rg_bam}.bai"
fi
