#!/usr/bin/env python3

import argparse
import gzip
import subprocess
from collections import defaultdict


KNOWN_SUFFIXES = [
    "_WASP_filtered",
    "_Aligned.sortedByCoord.out",
]


def parse_args():
    parser = argparse.ArgumentParser(description="Match VCF sample IDs to phASER BED sample IDs.")
    parser.add_argument("--vcf", required=True, help="VCF used for phaser_cis_var.")
    parser.add_argument("--bed", required=True, help="phASER bed.gz file.")
    parser.add_argument("--output", required=True, help="Output map file.")
    return parser.parse_args()


def clean(value):
    return value.strip().replace("\r", "")


def core_name(value):
    value = clean(value)
    for suffix in KNOWN_SUFFIXES:
        if value.endswith(suffix):
            return value[: -len(suffix)]
    return value


def main():
    args = parse_args()
    vcf_samples = subprocess.check_output(["bcftools", "query", "-l", args.vcf], text=True).splitlines()
    vcf_samples = [clean(x) for x in vcf_samples if clean(x)]

    with gzip.open(args.bed, "rt", encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").rstrip("\r").split("\t")
    bed_samples = [clean(x) for x in header[4:] if clean(x)]

    bed_exact = set(bed_samples)
    bed_by_core = defaultdict(list)
    for sample in bed_samples:
        bed_by_core[core_name(sample)].append(sample)

    used_bed = set()
    pairs = []
    for vcf_sample in vcf_samples:
        if vcf_sample in bed_exact:
            pairs.append((vcf_sample, vcf_sample))
            used_bed.add(vcf_sample)
            continue
        candidates = [x for x in bed_by_core.get(core_name(vcf_sample), []) if x not in used_bed]
        if len(candidates) == 1:
            pairs.append((vcf_sample, candidates[0]))
            used_bed.add(candidates[0])

    if not pairs:
        raise SystemExit("no matched samples between VCF and BED")

    with open(args.output, "w", encoding="utf-8") as handle:
        handle.write("vcf_sample\tbed_sample\n")
        for vcf_sample, bed_sample in pairs:
            handle.write(f"{vcf_sample}\t{bed_sample}\n")


if __name__ == "__main__":
    main()
