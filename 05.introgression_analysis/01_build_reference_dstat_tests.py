#!/usr/bin/env python3

import argparse
import csv
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Build individual-level D-statistic test tuples for screening putatively "
            "unadmixed reference individuals."
        )
    )
    parser.add_argument(
        "--candidate-table",
        required=True,
        help="TSV with columns: population, entity_id.",
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


def load_candidates(path):
    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"population", "entity_id"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError("candidate table must contain columns: population, entity_id")
        grouped = defaultdict(list)
        for row in reader:
            population = row["population"].strip()
            entity_id = row["entity_id"].strip()
            if population and entity_id:
                grouped[population].append(entity_id)
    return grouped


def load_entity_list(path):
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
        for entity_id in load_entity_list(path):
            donors.append((donor_population, entity_id))
    return donors


def main():
    args = parse_args()
    candidates = load_candidates(args.candidate_table)
    donors = parse_donor_specs(args.donor_file)

    rows = []
    for focal_population in sorted(candidates):
        peers = sorted(set(candidates[focal_population]))
        for focal_entity in peers:
            for peer_entity in peers:
                if peer_entity == focal_entity:
                    continue
                for donor_population, donor_entity in donors:
                    rows.append(
                        {
                            "test_scope": "reference_panel_screen",
                            "focal_population": focal_population,
                            "focal_entity": focal_entity,
                            "peer_entity": peer_entity,
                            "donor_population": donor_population,
                            "donor_entity": donor_entity,
                            "outgroup_entity": args.outgroup_id,
                            "W": focal_entity,
                            "X": peer_entity,
                            "Y": donor_entity,
                            "Z": args.outgroup_id,
                        }
                    )

    fieldnames = [
        "test_scope",
        "focal_population",
        "focal_entity",
        "peer_entity",
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
        writer.writerows(rows)


if __name__ == "__main__":
    main()
