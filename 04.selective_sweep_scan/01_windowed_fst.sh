#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
    echo "Usage: bash 01_windowed_fst.sh <job_table.tsv> [vcftools_bin]" >&2
    exit 1
fi

JOB_TABLE=$1
VCFTOOLS_BIN=${2:-vcftools}

awk -v bin="$VCFTOOLS_BIN" 'BEGIN{FS=OFS="\t"} NR==1{next} NF>0{
    comparison=$1
    pop1=$2
    pop2=$3
    win=$4
    step=$5
    vcf=$6
    out=$7

    cmd=sprintf("%s --gzvcf %s --weir-fst-pop %s --weir-fst-pop %s --fst-window-size %s --fst-window-step %s --out %s",
        bin, vcf, pop1, pop2, win, step, out)
    print cmd
}' "$JOB_TABLE" | while IFS= read -r cmd
do
    eval "$cmd"
done
