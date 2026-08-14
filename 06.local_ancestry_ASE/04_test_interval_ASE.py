#!/usr/bin/env python3

import argparse
from collections import defaultdict

from lib.common import benjamini_hochberg, exact_binomial_two_sided, chromosome_sort_key, read_tsv, write_tsv


REQUIRED_COLUMNS = [
    "RNA_sample",
    "WGS_individual",
    "tissue",
    "gene_associated_interval_id",
    "gene_associated_interval",
    "gene_id",
    "gene_symbol",
    "tract_id",
    "block_id",
    "chromosome",
    "block_start",
    "block_end",
    "n_informative_variants",
    "introgressed_ancestry",
    "background_ancestry",
    "introgressed_rfmix_haplotype",
    "background_rfmix_haplotype",
    "rfmix_hap1_ancestry",
    "rfmix_hap2_ancestry",
    "introgressed_haplotype_count",
    "background_haplotype_count",
    "totalCount",
    "strict_main_pass",
]

OUTPUT_COLUMNS = [
    "RNA_sample",
    "WGS_individual",
    "tissue",
    "gene_associated_interval_id",
    "gene_associated_interval",
    "gene_id",
    "gene_symbol",
    "tract_id",
    "chromosome",
    "interval_start",
    "interval_end",
    "introgressed_ancestry",
    "background_ancestry",
    "rfmix_hap1_ancestry",
    "rfmix_hap2_ancestry",
    "introgressed_rfmix_haplotype",
    "background_rfmix_haplotype",
    "introgressed_haplotype_count",
    "background_haplotype_count",
    "total_haplotype_count",
    "n_blocks",
    "block_ids",
    "n_informative_variants",
    "count_conservation_status",
    "minimum_total_count",
    "test_status",
    "binomial_p",
    "BH_FDR",
    "ASE_direction",
    "significant",
]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Aggregate strict phASER blocks and test ancestry-resolved interval ASE."
    )
    parser.add_argument("--strict-blocks", nargs="+", required=True)
    parser.add_argument("--all-output", required=True)
    parser.add_argument("--tested-output", required=True)
    parser.add_argument("--significant-output", required=True)
    parser.add_argument("--minimum-total-count", type=int, default=10)
    return parser.parse_args()


def parse_interval(interval):
    coordinate = interval.split(":", 1)[1]
    start, end = coordinate.split("-", 1)
    return int(start), int(end)


def main():
    args = parse_args()
    if args.minimum_total_count < 0:
        raise ValueError("--minimum-total-count must be non-negative")
    grouped = {}
    order = []
    seen_blocks = set()
    for path in args.strict_blocks:
        rows, _ = read_tsv(path, REQUIRED_COLUMNS, "strict block table")
        for row in rows:
            if row["strict_main_pass"] != "TRUE":
                raise ValueError(f"non-strict block supplied to aggregation: {row['block_id']}")
            block_record_id = (row["RNA_sample"], row["WGS_individual"], row["block_id"])
            if block_record_id in seen_blocks:
                raise ValueError(f"duplicate phASER block record: {block_record_id}")
            seen_blocks.add(block_record_id)
            key = (
                row["RNA_sample"],
                row["WGS_individual"],
                row["tissue"],
                row["gene_id"],
                row["gene_associated_interval_id"],
                row["tract_id"],
                row["introgressed_ancestry"],
                row["background_ancestry"],
                row["introgressed_rfmix_haplotype"],
                row["background_rfmix_haplotype"],
            )
            if key not in grouped:
                grouped[key] = {"rows": [], "first": row}
                order.append(key)
            grouped[key]["rows"].append(row)

    output = []
    for key in order:
        rows = grouped[key]["rows"]
        first = grouped[key]["first"]
        intro_values = [int(row["introgressed_haplotype_count"]) for row in rows]
        background_values = [int(row["background_haplotype_count"]) for row in rows]
        if any(value < 0 for value in intro_values + background_values):
            raise ValueError(f"negative haplotype count in interval {key}")
        intro_count = sum(intro_values)
        background_count = sum(background_values)
        total_count = intro_count + background_count
        expected_total = sum(int(row["totalCount"]) for row in rows)
        if total_count != expected_total:
            raise ValueError(f"block count conservation failed in interval {key}")
        if intro_count > background_count:
            direction = "introgressed_higher"
        elif intro_count < background_count:
            direction = "introgressed_lower"
        else:
            direction = "equal"
        if total_count < args.minimum_total_count:
            test_status, p_value = "LOW_TOTAL_COUNT", None
        else:
            test_status, p_value = "TESTED", exact_binomial_two_sided(intro_count, total_count)
        interval_start, interval_end = parse_interval(first["gene_associated_interval"])
        output.append(
            {
                "RNA_sample": first["RNA_sample"],
                "WGS_individual": first["WGS_individual"],
                "tissue": first["tissue"],
                "gene_associated_interval_id": first["gene_associated_interval_id"],
                "gene_associated_interval": first["gene_associated_interval"],
                "gene_id": first["gene_id"],
                "gene_symbol": first["gene_symbol"],
                "tract_id": first["tract_id"],
                "chromosome": first["chromosome"],
                "interval_start": interval_start,
                "interval_end": interval_end,
                "introgressed_ancestry": first["introgressed_ancestry"],
                "background_ancestry": first["background_ancestry"],
                "rfmix_hap1_ancestry": first["rfmix_hap1_ancestry"],
                "rfmix_hap2_ancestry": first["rfmix_hap2_ancestry"],
                "introgressed_rfmix_haplotype": first["introgressed_rfmix_haplotype"],
                "background_rfmix_haplotype": first["background_rfmix_haplotype"],
                "introgressed_haplotype_count": intro_count,
                "background_haplotype_count": background_count,
                "total_haplotype_count": total_count,
                "n_blocks": len(rows),
                "block_ids": ",".join(sorted(row["block_id"] for row in rows)),
                "n_informative_variants": sum(int(row["n_informative_variants"]) for row in rows),
                "count_conservation_status": "PASS",
                "minimum_total_count": args.minimum_total_count,
                "test_status": test_status,
                "binomial_p": f"{p_value:.17g}" if p_value is not None else "NA",
                "BH_FDR": "NA",
                "ASE_direction": direction,
                "significant": "FALSE",
            }
        )

    tested_by_sample = defaultdict(list)
    for index, row in enumerate(output):
        if row["test_status"] == "TESTED":
            tested_by_sample[row["RNA_sample"]].append(index)
    for sample, indexes in tested_by_sample.items():
        p_values = [float(output[index]["binomial_p"]) for index in indexes]
        for index, fdr in zip(indexes, benjamini_hochberg(p_values)):
            output[index]["BH_FDR"] = f"{fdr:.17g}"
            output[index]["significant"] = "TRUE" if fdr < 0.05 else "FALSE"

    output.sort(
        key=lambda row: (
            row["RNA_sample"],
            chromosome_sort_key(row["chromosome"]),
            int(row["interval_start"]),
            row["gene_id"],
            row["tract_id"],
        )
    )
    tested = [row for row in output if row["test_status"] == "TESTED"]
    significant = [row for row in tested if row["significant"] == "TRUE"]
    write_tsv(args.all_output, OUTPUT_COLUMNS, output)
    write_tsv(args.tested_output, OUTPUT_COLUMNS, tested)
    write_tsv(args.significant_output, OUTPUT_COLUMNS, significant)
    print(
        f"PASS interval_records={len(output)} tested={len(tested)} significant={len(significant)}"
    )


if __name__ == "__main__":
    main()
