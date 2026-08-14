#!/bin/bash

set -euo pipefail

if [[ $# -ne 5 ]]; then
  echo "Usage: $0 <gene_ae_dir> <gene_features.bed> <phaser_expr_matrix.py> <threads> <output_prefix>" >&2
  exit 1
fi

gene_ae_dir="$1"
gene_features="$2"
phaser_expr_matrix_py="$3"
threads="$4"
output_prefix="$5"

tmp_dir="${output_prefix}.tmp_gene_ae_dir"
rm -rf "${tmp_dir}"
mkdir -p "${tmp_dir}"

find "${gene_ae_dir}" -maxdepth 1 -type f -name '*_phaser_gene_ae.txt' -exec cp -f {} "${tmp_dir}"/ \;

for file in "${tmp_dir}"/*_phaser_gene_ae.txt; do
  header=$(grep -m1 -P '^contig\tstart\tstop\tname\t' "${file}" || true)
  grep -v -P '^contig\tstart\tstop\tname\t' "${file}" | sort -k1,1V -k2,2n > "${file}.body"
  {
    printf '%s\n' "${header}"
    cat "${file}.body"
  } > "${file}.fixed"
  mv "${file}.fixed" "${file}"
  rm -f "${file}.body"
done

python "${phaser_expr_matrix_py}" \
  --gene_ae_dir "${tmp_dir}" \
  --features "${gene_features}" \
  --t "${threads}" \
  --o "${output_prefix}"

rm -rf "${tmp_dir}"
