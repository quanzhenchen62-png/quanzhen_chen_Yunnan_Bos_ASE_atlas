#!/usr/bin/env python3

import argparse
import re


def parse_args():
    parser = argparse.ArgumentParser(description="Convert gene-body features from GTF to BED for phASER.")
    parser.add_argument("--gtf", required=True, help="Gene annotation GTF.")
    parser.add_argument("--output", required=True, help="Output BED file.")
    return parser.parse_args()


def extract_gene_id(attributes):
    patterns = [
        r'gene_id "([^"]+)"',
        r'ID=gene:([^;]+)',
        r'ID=([^;]+)',
    ]
    for pattern in patterns:
        match = re.search(pattern, attributes)
        if match:
            return match.group(1)
    return None


def main():
    args = parse_args()
    records = []
    with open(args.gtf, "r", encoding="utf-8") as handle:
      for line in handle:
        if not line or line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9 or fields[2] != "gene":
            continue
        gene_id = extract_gene_id(fields[8])
        if gene_id is None:
            continue
        chrom = fields[0]
        start = int(fields[3]) - 1
        end = int(fields[4])
        records.append((chrom, start, end, gene_id))

    records.sort(key=lambda x: (x[0], x[1], x[2], x[3]))
    with open(args.output, "w", encoding="utf-8") as handle:
        for chrom, start, end, gene_id in records:
            handle.write(f"{chrom}\t{start}\t{end}\t{gene_id}\n")


if __name__ == "__main__":
    main()
