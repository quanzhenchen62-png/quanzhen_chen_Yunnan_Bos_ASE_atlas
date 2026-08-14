#!/usr/bin/env python3

import argparse
import csv
import math


def parse_args():
    parser = argparse.ArgumentParser(
        description="Filter ASEReadCounter sites and test allelic imbalance per sample."
    )
    parser.add_argument("--input", required=True, help="ASEReadCounter output table.")
    parser.add_argument("--filtered-output", required=True, help="Filtered site table.")
    parser.add_argument("--significant-output", required=True, help="Significant ASE site table.")
    parser.add_argument("--min-allele-count", type=int, default=3, help="Minimum REF and ALT count. Default: 3.")
    parser.add_argument("--min-minor-ratio", type=float, default=0.01, help="Minimum minor-allele ratio. Default: 0.01.")
    parser.add_argument("--fdr-threshold", type=float, default=0.05, help="BH FDR threshold. Default: 0.05.")
    return parser.parse_args()


def binom_two_sided_p(k, n):
    if n <= 0:
        return float("nan")
    m = min(k, n - k)
    log2 = math.log(2.0)
    logps = []
    for i in range(m + 1):
        logp = math.lgamma(n + 1) - math.lgamma(i + 1) - math.lgamma(n - i + 1) - n * log2
        logps.append(logp)
    maxlog = max(logps)
    cdf = math.exp(maxlog) * sum(math.exp(x - maxlog) for x in logps)
    return min(1.0, 2.0 * cdf)


def bh_fdr(pvals):
    indexed = sorted(enumerate(pvals), key=lambda x: x[1])
    n = len(indexed)
    qvals = [1.0] * n
    running = 1.0
    for rank in range(n, 0, -1):
        idx, p = indexed[rank - 1]
        q = p * n / rank
        if q < running:
            running = q
        qvals[idx] = min(running, 1.0)
    return qvals


def main():
    args = parse_args()
    with open(args.input, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"contig", "position", "variantID", "refAllele", "altAllele", "refCount", "altCount", "totalCount"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError("ASEReadCounter table is missing required columns")
        rows = list(reader)

    filtered = []
    for row in rows:
        ref_count = int(float(row["refCount"]))
        alt_count = int(float(row["altCount"]))
        total_count = int(float(row["totalCount"]))
        if total_count <= 0:
            continue
        minor = min(ref_count, alt_count)
        minor_ratio = minor / total_count
        if ref_count > args.min_allele_count and alt_count > args.min_allele_count and minor_ratio > args.min_minor_ratio:
            out = dict(row)
            out["minorAlleleCount"] = str(minor)
            out["minorAlleleRatio"] = f"{minor_ratio:.6f}"
            out["pval"] = f"{binom_two_sided_p(alt_count, total_count):.12g}"
            alt_over_ref = alt_count / ref_count
            out["aFC"] = f"{math.log2(alt_over_ref):.12g}"
            filtered.append(out)

    if filtered:
        qvals = bh_fdr([float(row["pval"]) for row in filtered])
        for row, q in zip(filtered, qvals):
            row["FDR"] = f"{q:.12g}"
    else:
        qvals = []

    fieldnames = [
        "contig",
        "position",
        "variantID",
        "refAllele",
        "altAllele",
        "refCount",
        "altCount",
        "totalCount",
        "minorAlleleCount",
        "minorAlleleRatio",
        "pval",
        "FDR",
        "aFC",
    ]

    with open(args.filtered_output, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(filtered)

    significant = [row for row in filtered if float(row["FDR"]) < args.fdr_threshold]
    with open(args.significant_output, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(significant)


if __name__ == "__main__":
    main()
