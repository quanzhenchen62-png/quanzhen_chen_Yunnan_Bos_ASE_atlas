#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 || $# -gt 3 ]]; then
    echo "Usage: bash 04_run_selscan_xpehh.sh <job_table.tsv> [selscan_bin] [norm_bin]" >&2
    exit 1
fi

JOB_TABLE=$1
SELSCAN_BIN=${2:-selscan}
NORM_BIN=${3:-norm}

awk -v selscan_bin="$SELSCAN_BIN" -v norm_bin="$NORM_BIN" '
BEGIN{FS=OFS="\t"}
NR==1{next}
NF>0{
    comparison=$1
    chrom_list=$2
    query_pattern=$3
    ref_pattern=$4
    map_pattern=$5
    out_pattern=$6

    n=split(chrom_list, chroms, ",")
    for(i=1; i<=n; i++){
        chrom=chroms[i]
        gsub(/^[[:space:]]+|[[:space:]]+$/, "", chrom)
        if(chrom=="") continue

        query=query_pattern
        ref=ref_pattern
        mapf=map_pattern
        out=out_pattern

        gsub(/\{chrom\}/, chrom, query)
        gsub(/\{chrom\}/, chrom, ref)
        gsub(/\{chrom\}/, chrom, mapf)
        gsub(/\{chrom\}/, chrom, out)
        gsub(/\{comparison\}/, comparison, out)

        print sprintf("%s --xpehh --hap %s --ref %s --map %s --out %s", selscan_bin, query, ref, mapf, out)
        print sprintf("%s --xpehh --files %s.xpehh.out --out %s", norm_bin, out, out)
    }
}' "$JOB_TABLE" | while IFS= read -r cmd
do
    eval "$cmd"
done
