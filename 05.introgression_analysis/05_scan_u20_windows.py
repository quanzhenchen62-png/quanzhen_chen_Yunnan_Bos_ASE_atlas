#!/usr/bin/env python3

import argparse
import csv
import gzip
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Calculate U20-style donor-allele counts in non-overlapping windows and "
            "identify windows in the top 1%% of the window-level distribution."
        )
    )
    parser.add_argument("--vcf", required=True, help="Input VCF or VCF.GZ.")
    parser.add_argument(
        "--sample-table",
        required=True,
        help="TSV with columns sample and population.",
    )
    parser.add_argument("--source-pop", required=True, help="Source population label.")
    parser.add_argument("--target-pop", required=True, help="Target population label.")
    parser.add_argument("--donor-pop", required=True, help="Donor population label.")
    parser.add_argument(
        "--window-size",
        type=int,
        default=50000,
        help="Non-overlapping window size in bp. Default: 50000.",
    )
    parser.add_argument(
        "--source-max-af",
        type=float,
        default=0.01,
        help="Maximum donor-allele frequency allowed in the source population. Default: 0.01.",
    )
    parser.add_argument(
        "--target-min-af",
        type=float,
        default=0.20,
        help="Minimum donor-allele frequency required in the target population. Default: 0.20.",
    )
    parser.add_argument(
        "--donor-fixed-af",
        type=float,
        default=1.00,
        help="Required donor-allele frequency in the donor population. Default: 1.00.",
    )
    parser.add_argument(
        "--site-output",
        required=True,
        help="Per-site audit TSV.",
    )
    parser.add_argument(
        "--window-output",
        required=True,
        help="Window-level U20 summary TSV.",
    )
    parser.add_argument(
        "--top-output",
        required=True,
        help="Top-1%% putative introgression windows TSV.",
    )
    return parser.parse_args()


def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def load_sample_table(path):
    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"sample", "population"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError("sample table must contain columns: sample, population")
        sample_to_population = {}
        for row in reader:
            sample = row["sample"].strip()
            population = row["population"].strip()
            if sample:
                sample_to_population[sample] = population
    return sample_to_population


def parse_gt(field):
    gt = field.split(":", 1)[0]
    if gt in {".", "./.", ".|."}:
        return None
    sep = "|" if "|" in gt else "/"
    alleles = gt.split(sep)
    if len(alleles) != 2 or "." in alleles:
        return None
    try:
        return [int(alleles[0]), int(alleles[1])]
    except ValueError:
        return None


def alt_allele_frequency(sample_fields, sample_indices):
    alt_count = 0
    total_alleles = 0
    for idx in sample_indices:
        gt = parse_gt(sample_fields[idx])
        if gt is None:
            continue
        alt_count += gt[0] + gt[1]
        total_alleles += 2
    if total_alleles == 0:
        return None
    return alt_count / total_alleles


def donor_allele_state(donor_alt_af, donor_fixed_af):
    if donor_alt_af is None:
        return None
    if donor_alt_af >= donor_fixed_af:
        return "ALT"
    if (1.0 - donor_alt_af) >= donor_fixed_af:
        return "REF"
    return None


def donor_allele_frequency(alt_af, donor_allele):
    if alt_af is None:
        return None
    return alt_af if donor_allele == "ALT" else 1.0 - alt_af


def percentile_threshold(values, percentile):
    ordered = sorted(values)
    if not ordered:
        return 0
    rank = int(percentile * len(ordered))
    if percentile * len(ordered) > rank:
        rank += 1
    rank = min(max(rank, 1), len(ordered))
    return ordered[rank - 1]


def top_n_count(total_count, fraction):
    if total_count <= 0:
        return 0
    rank = int(total_count * fraction)
    if total_count * fraction > rank:
        rank += 1
    return max(rank, 1)


def main():
    args = parse_args()
    sample_to_population = load_sample_table(args.sample_table)

    site_rows = []
    windows = defaultdict(lambda: {"u20_count": 0, "informative_sites": 0, "qualifying_sites": 0})

    with open_text(args.vcf) as handle:
        header_samples = None
        population_indices = {}

        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header_samples = line.rstrip("\n").split("\t")[9:]
                pop_to_indices = defaultdict(list)
                for i, sample in enumerate(header_samples):
                    population = sample_to_population.get(sample)
                    if population is not None:
                        pop_to_indices[population].append(i)
                for population in [args.source_pop, args.target_pop, args.donor_pop]:
                    population_indices[population] = pop_to_indices.get(population, [])
                    if not population_indices[population]:
                        raise ValueError(f"no VCF samples found for population: {population}")
                continue

            fields = line.rstrip("\n").split("\t")
            chrom, pos, variant_id, ref, alt = fields[:5]
            if "," in alt:
                continue
            sample_fields = fields[9:]
            donor_alt_af = alt_allele_frequency(sample_fields, population_indices[args.donor_pop])
            source_alt_af = alt_allele_frequency(sample_fields, population_indices[args.source_pop])
            target_alt_af = alt_allele_frequency(sample_fields, population_indices[args.target_pop])

            donor_allele = donor_allele_state(donor_alt_af, args.donor_fixed_af)
            if donor_allele is None:
                continue

            source_donor_af = donor_allele_frequency(source_alt_af, donor_allele)
            target_donor_af = donor_allele_frequency(target_alt_af, donor_allele)
            informative = source_donor_af is not None and target_donor_af is not None
            qualifies = (
                informative
                and source_donor_af < args.source_max_af
                and target_donor_af > args.target_min_af
            )

            pos_int = int(pos)
            window_start = ((pos_int - 1) // args.window_size) * args.window_size + 1
            window_end = window_start + args.window_size - 1
            window_key = (chrom, window_start, window_end)

            if informative:
                windows[window_key]["informative_sites"] += 1
            if donor_allele is not None:
                windows[window_key]["qualifying_sites"] += 1
            if qualifies:
                windows[window_key]["u20_count"] += 1

            site_rows.append(
                {
                    "chromosome": chrom,
                    "position": pos_int,
                    "variant_id": variant_id,
                    "REF": ref,
                    "ALT": alt,
                    "source_population": args.source_pop,
                    "target_population": args.target_pop,
                    "donor_population": args.donor_pop,
                    "donor_allele": donor_allele,
                    "source_donor_allele_frequency": "" if source_donor_af is None else f"{source_donor_af:.6f}",
                    "target_donor_allele_frequency": "" if target_donor_af is None else f"{target_donor_af:.6f}",
                    "donor_alt_allele_frequency": "" if donor_alt_af is None else f"{donor_alt_af:.6f}",
                    "passes_u20_definition": "TRUE" if qualifies else "FALSE",
                    "window_start": window_start,
                    "window_end": window_end,
                }
            )

    u20_values = [payload["u20_count"] for payload in windows.values()]
    threshold = percentile_threshold(u20_values, 0.99)

    with open(args.site_output, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "chromosome",
                "position",
                "variant_id",
                "REF",
                "ALT",
                "source_population",
                "target_population",
                "donor_population",
                "donor_allele",
                "source_donor_allele_frequency",
                "target_donor_allele_frequency",
                "donor_alt_allele_frequency",
                "passes_u20_definition",
                "window_start",
                "window_end",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(site_rows)

    sorted_window_keys = sorted(windows)
    ranked_window_keys = sorted(
        sorted_window_keys,
        key=lambda key: (
            -windows[key]["u20_count"],
            key[0],
            key[1],
            key[2],
        ),
    )
    n_top_windows = top_n_count(len(ranked_window_keys), 0.01)
    top_window_keys = set(ranked_window_keys[:n_top_windows])

    window_rows = []
    for chrom, start, end in sorted_window_keys:
        payload = windows[(chrom, start, end)]
        row = {
            "chromosome": chrom,
            "window_start": start,
            "window_end": end,
            "source_population": args.source_pop,
            "target_population": args.target_pop,
            "donor_population": args.donor_pop,
            "window_size_bp": args.window_size,
            "u20_count": payload["u20_count"],
            "informative_sites": payload["informative_sites"],
            "donor_fixed_sites": payload["qualifying_sites"],
            "u20_top_1pct_threshold": threshold,
            "is_top_1pct_window": (
                "TRUE" if (chrom, start, end) in top_window_keys else "FALSE"
            ),
        }
        window_rows.append(row)

    with open(args.window_output, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "chromosome",
                "window_start",
                "window_end",
                "source_population",
                "target_population",
                "donor_population",
                "window_size_bp",
                "u20_count",
                "informative_sites",
                "donor_fixed_sites",
                "u20_top_1pct_threshold",
                "is_top_1pct_window",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(window_rows)

    top_rows = [row for row in window_rows if row["is_top_1pct_window"] == "TRUE"]
    with open(args.top_output, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "chromosome",
                "window_start",
                "window_end",
                "source_population",
                "target_population",
                "donor_population",
                "window_size_bp",
                "u20_count",
                "informative_sites",
                "donor_fixed_sites",
                "u20_top_1pct_threshold",
                "is_top_1pct_window",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(top_rows)


if __name__ == "__main__":
    main()
