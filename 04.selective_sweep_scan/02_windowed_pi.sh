#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
    echo "Usage: bash 02_windowed_pi.sh <job_table.tsv> [vcftools_bin]" >&2
    exit 1
fi

JOB_TABLE=$1
VCFTOOLS_BIN=${2:-vcftools}

awk -v bin="$VCFTOOLS_BIN" 'BEGIN{FS=OFS="\t"} NR==1{next} NF>0{
    population=$1
    keep=$2
    win=$3
    step=$4
    vcf=$5
    out=$6

    cmd=sprintf("%s --gzvcf %s --keep %s --window-pi %s --window-pi-step %s --out %s",
        bin, vcf, keep, win, step, out)
    print cmd
}' "$JOB_TABLE" | while IFS= read -r cmd
do
    eval "$cmd"
done
