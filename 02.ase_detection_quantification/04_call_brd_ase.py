#!/usr/bin/env python3

import argparse
import csv
import math
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(
        description="Call breed-level ASE sites within breed-tissue groups."
    )
    parser.add_argument(
        "--metadata",
        required=True,
        help=(
            "TSV with columns: sample_id, breed, tissue, significant_ase_path. "
            "Each significant_ase_path should point to the per-sample ind-ASE table."
        ),
    )
    parser.add_argument(
        "--output-prefix",
        required=True,
        help="Prefix for brd-ASE outputs.",
    )
    return parser.parse_args()


def safe_float(value):
    try:
        return float(value)
    except Exception:
        return float("nan")


def signed_group_afc(afc_values):
    if len(afc_values) == 1:
        return afc_values[0], "single"
    if all(x > 0 for x in afc_values) or all(x < 0 for x in afc_values):
        median_fold = median([2 ** x for x in afc_values])
        return math.log2(median_fold), "consistent"
    return median([abs(x) for x in afc_values]), "inconsistent"


def median(values):
    ordered = sorted(values)
    n = len(ordered)
    if n % 2 == 1:
        return ordered[n // 2]
    return (ordered[n // 2 - 1] + ordered[n // 2]) / 2.0


def load_metadata(path):
    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"sample_id", "breed", "tissue", "significant_ase_path"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError("metadata is missing required columns")
        return list(reader)


def load_significant_sites(path):
    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "contig",
            "position",
            "variantID",
            "refAllele",
            "altAllele",
            "aFC",
        }
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"significant ASE file is missing required columns: {path}")
        rows = []
        for row in reader:
            row["aFC"] = safe_float(row["aFC"])
            rows.append(row)
        return rows


def main():
    args = parse_args()
    metadata = load_metadata(args.metadata)

    grouped = defaultdict(list)
    for row in metadata:
        grouped[(row["breed"], row["tissue"])].append(row)

    all_rows = []
    threshold_rows = {1: [], 2: [], 3: []}

    for (breed, tissue), sample_rows in sorted(grouped.items()):
        variant_records = defaultdict(list)
        for sample_row in sample_rows:
            sample_id = sample_row["sample_id"]
            sig_rows = load_significant_sites(sample_row["significant_ase_path"])
            seen = set()
            for row in sig_rows:
                variant_id = row["variantID"]
                if variant_id in seen:
                    continue
                seen.add(variant_id)
                variant_records[variant_id].append(
                    {
                        "sample_id": sample_id,
                        "aFC": row["aFC"],
                        "contig": row["contig"],
                        "position": row["position"],
                        "refAllele": row["refAllele"],
                        "altAllele": row["altAllele"],
                    }
                )

        for variant_id, records in sorted(variant_records.items()):
            afc_values = [record["aFC"] for record in records if not math.isnan(record["aFC"])]
            if not afc_values:
                continue
            pop_afc, direction = signed_group_afc(afc_values)
            sig_samples = sorted({record["sample_id"] for record in records})
            base = {
                "breed": breed,
                "tissue": tissue,
                "contig": records[0]["contig"],
                "position": records[0]["position"],
                "variantID": variant_id,
                "refAllele": records[0]["refAllele"],
                "altAllele": records[0]["altAllele"],
                "sig_sample_count": len(sig_samples),
                "pop_aFC": f"{pop_afc:.12g}",
                "direction_consistency": direction,
                "sig_samples": ",".join(sig_samples),
            }
            all_rows.append(base)
            for threshold in (1, 2, 3):
                if len(sig_samples) >= threshold:
                    threshold_rows[threshold].append(dict(base))

    fieldnames = [
        "breed",
        "tissue",
        "contig",
        "position",
        "variantID",
        "refAllele",
        "altAllele",
        "sig_sample_count",
        "pop_aFC",
        "direction_consistency",
        "sig_samples",
    ]

    with open(f"{args.output_prefix}.all.tsv", "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(all_rows)

    for threshold in (1, 2, 3):
        with open(
            f"{args.output_prefix}.brd_ASE_{threshold}.tsv",
            "w",
            encoding="utf-8",
            newline="",
        ) as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(threshold_rows[threshold])


if __name__ == "__main__":
    main()
