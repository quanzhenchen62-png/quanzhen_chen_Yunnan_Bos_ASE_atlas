#!/bin/bash

set -euo pipefail

if [[ $# -ne 6 ]]; then
  echo "Usage: $0 <pair_file> <bed.gz> <vcf.gz> <map_file> <phaser_cis_var.py> <output.txt>" >&2
  exit 1
fi

pair_file="$1"
bed_gz="$2"
vcf_gz="$3"
map_file="$4"
phaser_cis_var_py="$5"
output_txt="$6"

if [[ ! -f "${vcf_gz}.csi" && ! -f "${vcf_gz}.tbi" ]]; then
  tabix -p vcf "${vcf_gz}"
fi

python "${phaser_cis_var_py}" \
  --bed "${bed_gz}" \
  --vcf "${vcf_gz}" \
  --pair "${pair_file}" \
  --map "${map_file}" \
  --o "${output_txt}" \
  --ignore_v 1 \
  --bs 1000 \
  --t 4
