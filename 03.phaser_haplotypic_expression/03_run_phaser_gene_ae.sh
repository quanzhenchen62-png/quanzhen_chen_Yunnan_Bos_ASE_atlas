#!/bin/bash

set -euo pipefail

if [[ $# -ne 4 ]]; then
  echo "Usage: $0 <haplotypic_counts.txt> <gene_features.bed> <phaser_gene_ae.py> <output.txt>" >&2
  exit 1
fi

haplotypic_counts="$1"
gene_features="$2"
phaser_gene_ae_py="$3"
output_txt="$4"

tmp_clean="${output_txt}.haplotypic_counts.clean.tmp"
tmp_bad="${output_txt}.haplotypic_counts.bad_lines.tmp"

awk -F'\t' 'NF==18 {print}' "${haplotypic_counts}" > "${tmp_clean}"
awk -F'\t' 'NF!=18 {print}' "${haplotypic_counts}" > "${tmp_bad}" || true

python "${phaser_gene_ae_py}" \
  --haplotypic_counts "${tmp_clean}" \
  --features "${gene_features}" \
  --o "${output_txt}"

rm -f "${tmp_clean}" "${tmp_bad}"
