#!/usr/bin/env python3

import argparse
import csv
import math
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(
        description="Call tissue-level population ASE sites from per-sample ASE outputs."
    )
    parser.add_argument(
        "--metadata",
        required=True,
        help=(
            "TSV with columns: sample_id, breed, tissue, filtered_ase_path, significant_ase_path."
        ),
    )
    parser.add_argument(
        "--output-prefix",
        required=True,
        help="Prefix for pop-ASE outputs.",
    )
    parser.add_argument(
        "--min-heterozygous-samples",
        type=int,
        default=5,
        help="Minimum number of heterozygous samples passing ASE filters. Default: 5.",
    )
    parser.add_argument(
        "--min-significant-fraction",
        type=float,
        default=0.60,
        help="Minimum fraction of significant samples among filtered heterozygous samples. Default: 0.60.",
    )
    parser.add_argument(
        "--min-breeds",
        type=int,
        default=2,
        help="Minimum number of breeds among significant samples. Default: 2.",
    )
    return parser.parse_args()


def safe_float(value):
    try:
        return float(value)
    except Exception:
        return float("nan")


def median(values):
    ordered = sorted(values)
    n = len(ordered)
    if n % 2 == 1:
        return ordered[n // 2]
    return (ordered[n // 2 - 1] + ordered[n // 2]) / 2.0


def signed_group_afc(afc_values):
    if len(afc_values) == 1:
        return afc_values[0], "single"
    if all(x > 0 for x in afc_values) or all(x < 0 for x in afc_values):
        median_fold = median([2 ** x for x in afc_values])
        return math.log2(median_fold), "consistent"
    return median([abs(x) for x in afc_values]), "inconsistent"


def load_tsv(path, required):
    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"file is missing required columns: {path}")
        return list(reader)


def main():
    args = parse_args()
    metadata = load_tsv(
        args.metadata,
        {"sample_id", "breed", "tissue", "filtered_ase_path", "significant_ase_path"},
    )

    grouped = defaultdict(list)
    for row in metadata:
        grouped[row["tissue"]].append(row)

    fieldnames = [
        "tissue",
        "contig",
        "position",
        "variantID",
        "refAllele",
        "altAllele",
        "ori_count",
        "sig_count",
        "sig_fraction",
        "sig_breed_count",
        "pop_aFC",
        "direction_consistency",
        "ori_samples",
        "sig_samples",
    ]

    with open(f"{args.output_prefix}.all.tsv", "w", encoding="utf-8", newline="") as all_handle:
        all_writer = csv.DictWriter(all_handle, fieldnames=fieldnames, delimiter="\t")
        all_writer.writeheader()

        with open(f"{args.output_prefix}.significant.tsv", "w", encoding="utf-8", newline="") as sig_handle:
            sig_writer = csv.DictWriter(sig_handle, fieldnames=fieldnames, delimiter="\t")
            sig_writer.writeheader()

            for tissue, sample_rows in sorted(grouped.items()):
                filtered_sites = defaultdict(dict)
                significant_sites = defaultdict(dict)
                sample_to_breed = {}

                for sample_row in sample_rows:
                    sample_id = sample_row["sample_id"]
                    sample_to_breed[sample_id] = sample_row["breed"]

                    filtered = load_tsv(
                        sample_row["filtered_ase_path"],
                        {"contig", "position", "variantID", "refAllele", "altAllele", "aFC"},
                    )
                    for row in filtered:
                        filtered_sites[row["variantID"]][sample_id] = row

                    significant = load_tsv(
                        sample_row["significant_ase_path"],
                        {"contig", "position", "variantID", "refAllele", "altAllele", "aFC"},
                    )
                    for row in significant:
                        significant_sites[row["variantID"]][sample_id] = row

                for variant_id in sorted(filtered_sites):
                    ori_samples = sorted(filtered_sites[variant_id])
                    sig_samples = sorted(significant_sites.get(variant_id, {}))
                    ori_count = len(ori_samples)
                    sig_count = len(sig_samples)
                    sig_fraction = sig_count / ori_count if ori_count > 0 else 0.0
                    sig_breeds = sorted({sample_to_breed[s] for s in sig_samples if s in sample_to_breed})
                    afc_values = [
                        safe_float(significant_sites[variant_id][sample_id]["aFC"])
                        for sample_id in sig_samples
                        if not math.isnan(safe_float(significant_sites[variant_id][sample_id]["aFC"]))
                    ]

                    if afc_values:
                        pop_afc, direction = signed_group_afc(afc_values)
                        pop_afc_text = f"{pop_afc:.12g}"
                    else:
                        pop_afc_text = ""
                        direction = ""

                    first = filtered_sites[variant_id][ori_samples[0]]
                    out_row = {
                        "tissue": tissue,
                        "contig": first["contig"],
                        "position": first["position"],
                        "variantID": variant_id,
                        "refAllele": first["refAllele"],
                        "altAllele": first["altAllele"],
                        "ori_count": ori_count,
                        "sig_count": sig_count,
                        "sig_fraction": f"{sig_fraction:.6f}",
                        "sig_breed_count": len(sig_breeds),
                        "pop_aFC": pop_afc_text,
                        "direction_consistency": direction,
                        "ori_samples": ",".join(ori_samples),
                        "sig_samples": ",".join(sig_samples),
                    }
                    all_writer.writerow(out_row)

                    if (
                        ori_count >= args.min_heterozygous_samples
                        and sig_fraction >= args.min_significant_fraction
                        and len(sig_breeds) >= args.min_breeds
                    ):
                        sig_writer.writerow(out_row)


if __name__ == "__main__":
    main()
