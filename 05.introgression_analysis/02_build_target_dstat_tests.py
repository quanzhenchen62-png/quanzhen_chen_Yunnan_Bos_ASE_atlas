#!/usr/bin/env python3

import argparse
import csv


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Build D-statistic test tuples for target breeds under "
            "D(source, target; donor, outgroup)."
        )
    )
    parser.add_argument(
        "--source-id",
        required=True,
        help="Source label to use as W, for example the retained SAI reference panel label.",
    )
    parser.add_argument(
        "--target-file",
        required=True,
        help="Text file with one target label per line.",
    )
    parser.add_argument(
        "--donor-file",
        required=True,
        action="append",
        help="Donor specification in the form donor_population=path/to/entities.txt.",
    )
    parser.add_argument(
        "--outgroup-id",
        required=True,
        help="Outgroup label to use as Z in D(W, X; Y, Z).",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output TSV of D-statistic tests.",
    )
    return parser.parse_args()


def load_list(path):
    values = []
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            value = line.strip()
            if value and not value.startswith("#"):
                values.append(value)
    return values


def parse_donor_specs(specs):
    donors = []
    for spec in specs:
        if "=" not in spec:
            raise ValueError(f"invalid donor specification: {spec}")
        donor_population, path = spec.split("=", 1)
        donor_population = donor_population.strip()
        path = path.strip()
        if not donor_population or not path:
            raise ValueError(f"invalid donor specification: {spec}")
        for entity_id in load_list(path):
            donors.append((donor_population, entity_id))
    return donors


def main():
    args = parse_args()
    targets = load_list(args.target_file)
    donors = parse_donor_specs(args.donor_file)

    fieldnames = [
        "test_scope",
        "source_entity",
        "target_entity",
        "donor_population",
        "donor_entity",
        "outgroup_entity",
        "W",
        "X",
        "Y",
        "Z",
    ]

    with open(args.output, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for target_entity in targets:
            for donor_population, donor_entity in donors:
                writer.writerow(
                    {
                        "test_scope": "target_breed_introgression",
                        "source_entity": args.source_id,
                        "target_entity": target_entity,
                        "donor_population": donor_population,
                        "donor_entity": donor_entity,
                        "outgroup_entity": args.outgroup_id,
                        "W": args.source_id,
                        "X": target_entity,
                        "Y": donor_entity,
                        "Z": args.outgroup_id,
                    }
                )


if __name__ == "__main__":
    main()
