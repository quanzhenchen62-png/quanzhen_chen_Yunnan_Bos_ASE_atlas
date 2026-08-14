#!/usr/bin/env python3

import argparse
import csv
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Retain putatively unadmixed reference individuals that satisfy "
            "|Z-score| < threshold across all individual-level D-statistic tests."
        )
    )
    parser.add_argument(
        "--results",
        required=True,
        help=(
            "TSV produced by 03_run_qpdstat.py from reference-panel tests. "
            "Must contain focal_population, focal_entity, and Zscore."
        ),
    )
    parser.add_argument(
        "--z-threshold",
        type=float,
        default=3.0,
        help="Absolute Z-score threshold. Default: 3.0.",
    )
    parser.add_argument("--audit-output", required=True, help="Full audit TSV.")
    parser.add_argument("--retained-output", required=True, help="Retained individuals TSV.")
    return parser.parse_args()


def main():
    args = parse_args()
    with open(args.results, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"focal_population", "focal_entity", "Zscore"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(
                "results file must contain columns: focal_population, focal_entity, Zscore"
            )
        grouped = defaultdict(list)
        for row in reader:
            focal_population = row["focal_population"].strip()
            focal_entity = row["focal_entity"].strip()
            z_raw = row["Zscore"].strip()
            if not focal_population or not focal_entity or not z_raw:
                continue
            grouped[(focal_population, focal_entity)].append(abs(float(z_raw)))

    audit_rows = []
    retained_rows = []
    for (population, entity_id), abs_zscores in sorted(grouped.items()):
        max_abs_z = max(abs_zscores)
        n_tests = len(abs_zscores)
        keep = max_abs_z < args.z_threshold
        row = {
            "population": population,
            "entity_id": entity_id,
            "n_tests": n_tests,
            "max_abs_zscore": f"{max_abs_z:.6f}",
            "z_threshold": f"{args.z_threshold:.6f}",
            "retain_reference_individual": "TRUE" if keep else "FALSE",
        }
        audit_rows.append(row)
        if keep:
            retained_rows.append({"population": population, "entity_id": entity_id})

    with open(args.audit_output, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "population",
                "entity_id",
                "n_tests",
                "max_abs_zscore",
                "z_threshold",
                "retain_reference_individual",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(audit_rows)

    with open(args.retained_output, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["population", "entity_id"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(retained_rows)


if __name__ == "__main__":
    main()
