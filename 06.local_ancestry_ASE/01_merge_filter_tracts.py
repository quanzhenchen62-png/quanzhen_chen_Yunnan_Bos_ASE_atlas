#!/usr/bin/env python3

import argparse
import csv
import math
import re
from collections import defaultdict
from pathlib import Path

from lib.common import (
    NA,
    bool_text,
    chromosome_sort_key,
    normalize_chromosome,
    open_text,
    split_fields,
    write_tsv,
)


OUTPUT_COLUMNS = [
    "WGS_individual",
    "chromosome",
    "tract_id",
    "tract_start",
    "tract_end",
    "tract_length_bp",
    "tract_n_snps",
    "tract_n_msp_intervals",
    "first_msp_row",
    "last_msp_row",
    "rfmix_hap1_ancestry",
    "rfmix_hap2_ancestry",
    "ordered_ancestry_pair",
    "n_fb_markers",
    "rfmix_hap1_assigned_posterior_mean",
    "rfmix_hap2_assigned_posterior_mean",
    "posterior_threshold",
    "posterior_pass",
    "minimum_marker_count",
    "marker_count_pass",
    "tract_high_confidence",
    "putatively_introgressed",
    "strict_ancestry_heterozygous_target",
    "introgressed_ancestry",
    "background_ancestry",
    "introgressed_rfmix_haplotype",
    "background_rfmix_haplotype",
    "tract_filter_reason",
]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Merge phased RFMix MSP intervals and filter high-confidence ancestry tracts."
    )
    parser.add_argument("--msp", nargs="+", required=True, help="RFMix MSP TSV files")
    parser.add_argument("--fb", nargs="+", required=True, help="Matching RFMix FB TSV files")
    parser.add_argument("--output", required=True, help="All merged tracts TSV")
    parser.add_argument("--strict-output", required=True, help="Strict donor/background tracts TSV")
    parser.add_argument("--posterior-threshold", type=float, default=0.80)
    parser.add_argument("--minimum-markers", type=int, default=30)
    parser.add_argument("--donor-ancestries", default="Yak,Gayal")
    parser.add_argument("--background-ancestries", default="SAI,Taurine")
    parser.add_argument("--posterior-sum-tolerance", type=float, default=1e-4)
    return parser.parse_args()


def parse_msp(path):
    ancestry_codes = {}
    header = None
    raw_rows = []
    with open_text(path) as handle:
        for line in handle:
            if not line.strip():
                continue
            if line.startswith("#Subpopulation order/codes:"):
                for ancestry, code in re.findall(r"([^\s=]+)=(\d+)", line):
                    ancestry_codes[int(code)] = ancestry
                continue
            if line.startswith("#chm"):
                header = split_fields(line)
                header[0] = header[0].lstrip("#")
                continue
            if line.startswith("#"):
                continue
            if header is None:
                raise ValueError(f"MSP header beginning with #chm was not found: {path}")
            fields = split_fields(line)
            if len(fields) != len(header):
                raise ValueError(f"MSP field-count mismatch in {path}")
            raw_rows.append(fields)
    if not ancestry_codes:
        raise ValueError(f"MSP ancestry code declaration was not found: {path}")
    if not raw_rows:
        raise ValueError(f"MSP contains no data rows: {path}")

    sample_columns = defaultdict(dict)
    for index, name in enumerate(header[6:], start=6):
        match = re.fullmatch(r"(.+)\.([01])", name)
        if not match:
            raise ValueError(f"invalid MSP haplotype column: {name}")
        sample_columns[match.group(1)][int(match.group(2))] = index
    for sample, columns in sample_columns.items():
        if set(columns) != {0, 1}:
            raise ValueError(f"MSP sample lacks .0 or .1 haplotype column: {sample}")

    rows = []
    for row_number, fields in enumerate(raw_rows, start=1):
        chromosome = normalize_chromosome(fields[0])
        start, end = int(float(fields[1])), int(float(fields[2]))
        if start < 1 or end < start:
            raise ValueError(f"invalid MSP coordinates at data row {row_number}")
        n_snps = int(float(fields[5]))
        if n_snps < 0:
            raise ValueError(f"negative MSP marker count at data row {row_number}")
        ancestry_by_sample = {}
        for sample, columns in sample_columns.items():
            hap1_code, hap2_code = int(fields[columns[0]]), int(fields[columns[1]])
            if hap1_code not in ancestry_codes or hap2_code not in ancestry_codes:
                raise ValueError(f"unknown ancestry code at MSP data row {row_number}")
            ancestry_by_sample[sample] = (
                ancestry_codes[hap1_code],
                ancestry_codes[hap2_code],
            )
        rows.append(
            {
                "row_number": row_number,
                "chromosome": chromosome,
                "start": start,
                "end": end,
                "n_snps": n_snps,
                "ancestry_by_sample": ancestry_by_sample,
            }
        )
    for previous, current in zip(rows, rows[1:]):
        if current["chromosome"] != previous["chromosome"]:
            raise ValueError("one MSP input must contain exactly one chromosome")
        if current["start"] < previous["start"]:
            raise ValueError("MSP rows are not coordinate sorted")
    return rows, sorted(sample_columns), ancestry_codes


def build_tracts(msp_rows, samples):
    tracts = defaultdict(list)
    row_to_tract = {sample: [] for sample in samples}
    for sample in samples:
        current = None
        for row in msp_rows:
            ancestry_pair = row["ancestry_by_sample"][sample]
            contiguous = current is not None and row["start"] <= current["tract_end"] + 1
            if current is None or not contiguous or current["_pair"] != ancestry_pair:
                current = {
                    "WGS_individual": sample,
                    "chromosome": row["chromosome"],
                    "tract_start": row["start"],
                    "tract_end": row["end"],
                    "tract_n_snps": row["n_snps"],
                    "tract_n_msp_intervals": 1,
                    "first_msp_row": row["row_number"],
                    "last_msp_row": row["row_number"],
                    "rfmix_hap1_ancestry": ancestry_pair[0],
                    "rfmix_hap2_ancestry": ancestry_pair[1],
                    "_pair": ancestry_pair,
                    "_hap1_posteriors": [],
                    "_hap2_posteriors": [],
                }
                tracts[sample].append(current)
            else:
                current["tract_end"] = row["end"]
                current["tract_n_snps"] += row["n_snps"]
                current["tract_n_msp_intervals"] += 1
                current["last_msp_row"] = row["row_number"]
            row_to_tract[sample].append(len(tracts[sample]) - 1)
        for tract_number, tract in enumerate(tracts[sample], start=1):
            tract["tract_id"] = (
                f"{sample}_chr{tract['chromosome']}_tract{tract_number:05d}"
            )
    return tracts, row_to_tract


def parse_fb_header(path):
    populations = []
    with open_text(path) as handle:
        for line in handle:
            if not line.strip():
                continue
            if line.startswith("#reference_panel_population:") or line.startswith(
                "#reference population:"
            ):
                populations = re.split(r"\s+", line.split(":", 1)[1].strip())
                continue
            if line.startswith("#"):
                continue
            return split_fields(line), populations
    raise ValueError(f"FB header was not found: {path}")


def add_fb_posteriors(path, msp_rows, samples, ancestry_codes, tracts, row_to_tract, tolerance):
    header, declared_populations = parse_fb_header(path)
    populations = declared_populations or [ancestry_codes[index] for index in sorted(ancestry_codes)]
    if set(populations) != set(ancestry_codes.values()):
        raise ValueError(f"FB/MSP ancestry populations differ in {path}")
    column_index = {name: index for index, name in enumerate(header)}
    if "physical_position" not in column_index:
        raise ValueError(f"FB lacks physical_position: {path}")
    if "chromosome" not in column_index:
        raise ValueError(f"FB lacks chromosome: {path}")
    for sample in samples:
        for haplotype in ("hap1", "hap2"):
            for ancestry in populations:
                key = f"{sample}:::{haplotype}:::{ancestry}"
                if key not in column_index:
                    raise ValueError(f"FB lacks posterior column: {key}")

    msp_index = 0
    previous_position = None
    posterior_sum_failures = 0
    with open_text(path) as handle:
        seen_header = False
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = split_fields(line)
            if not seen_header:
                if fields != header:
                    raise ValueError(f"FB header changed while reading {path}")
                seen_header = True
                continue
            if len(fields) != len(header):
                raise ValueError(f"FB field-count mismatch in {path}")
            chromosome = normalize_chromosome(fields[column_index["chromosome"]])
            position = int(float(fields[column_index["physical_position"]]))
            if chromosome != msp_rows[0]["chromosome"]:
                raise ValueError(f"FB chromosome does not match MSP input: {chromosome}")
            if previous_position is not None and position < previous_position:
                raise ValueError(f"FB markers are not coordinate sorted in {path}")
            previous_position = position
            while msp_index + 1 < len(msp_rows) and position > msp_rows[msp_index]["end"]:
                msp_index += 1
            msp_row = msp_rows[msp_index]
            if not (msp_row["start"] <= position <= msp_row["end"]):
                raise ValueError(f"FB marker {position} does not match an MSP interval")
            for sample in samples:
                for haplotype in ("hap1", "hap2"):
                    posterior_values = [
                        float(fields[column_index[f"{sample}:::{haplotype}:::{ancestry}"]])
                        for ancestry in populations
                    ]
                    if any(
                        not math.isfinite(value) or value < 0 or value > 1
                        for value in posterior_values
                    ):
                        raise ValueError(
                            f"invalid FB posterior for {sample} {haplotype} at {position}"
                        )
                    posterior_sum = math.fsum(posterior_values)
                    if abs(posterior_sum - 1.0) > tolerance:
                        posterior_sum_failures += 1
                ancestry1, ancestry2 = msp_row["ancestry_by_sample"][sample]
                tract = tracts[sample][row_to_tract[sample][msp_index]]
                tract["_hap1_posteriors"].append(
                    float(fields[column_index[f"{sample}:::hap1:::{ancestry1}"]])
                )
                tract["_hap2_posteriors"].append(
                    float(fields[column_index[f"{sample}:::hap2:::{ancestry2}"]])
                )
    if posterior_sum_failures:
        raise ValueError(
            f"FB posterior sums differ from 1 beyond tolerance in {posterior_sum_failures} sample-haplotype marker rows"
        )


def finalize_tracts(tracts, posterior_threshold, minimum_markers, donors, backgrounds):
    output = []
    for sample in sorted(tracts):
        for tract in tracts[sample]:
            hap1_values = tract.pop("_hap1_posteriors")
            hap2_values = tract.pop("_hap2_posteriors")
            tract.pop("_pair")
            if len(hap1_values) != len(hap2_values):
                raise ValueError(f"FB marker count differs between haplotypes for {tract['tract_id']}")
            hap1_mean = math.fsum(hap1_values) / len(hap1_values) if hap1_values else None
            hap2_mean = math.fsum(hap2_values) / len(hap2_values) if hap2_values else None
            posterior_pass = (
                hap1_mean is not None
                and hap2_mean is not None
                and hap1_mean >= posterior_threshold
                and hap2_mean >= posterior_threshold
            )
            marker_pass = tract["tract_n_snps"] >= minimum_markers
            high_confidence = posterior_pass and marker_pass
            ancestry1, ancestry2 = tract["rfmix_hap1_ancestry"], tract["rfmix_hap2_ancestry"]
            donor_haplotypes = [
                haplotype
                for haplotype, ancestry in (("hap1", ancestry1), ("hap2", ancestry2))
                if ancestry in donors
            ]
            background_haplotypes = [
                haplotype
                for haplotype, ancestry in (("hap1", ancestry1), ("hap2", ancestry2))
                if ancestry in backgrounds
            ]
            strict_target = high_confidence and len(donor_haplotypes) == 1 and len(background_haplotypes) == 1
            if not hap1_values:
                reason = "NO_FB_MARKERS"
            elif not posterior_pass:
                reason = "POSTERIOR_BELOW_THRESHOLD"
            elif not marker_pass:
                reason = "MARKER_COUNT_BELOW_MINIMUM"
            elif strict_target:
                reason = "PASS_STRICT_TARGET"
            elif donor_haplotypes:
                reason = "PASS_PUTATIVELY_INTROGRESSED_NON_TARGET_CONFIGURATION"
            else:
                reason = "PASS_NON_DONOR_TRACT"
            introgressed_haplotype = donor_haplotypes[0] if strict_target else NA
            background_haplotype = background_haplotypes[0] if strict_target else NA
            introgressed_ancestry = (
                ancestry1 if introgressed_haplotype == "hap1" else ancestry2
            ) if strict_target else NA
            background_ancestry = (
                ancestry1 if background_haplotype == "hap1" else ancestry2
            ) if strict_target else NA
            tract.update(
                {
                    "tract_length_bp": tract["tract_end"] - tract["tract_start"] + 1,
                    "ordered_ancestry_pair": f"{ancestry1}|{ancestry2}",
                    "n_fb_markers": len(hap1_values),
                    "rfmix_hap1_assigned_posterior_mean": f"{hap1_mean:.10g}" if hap1_mean is not None else NA,
                    "rfmix_hap2_assigned_posterior_mean": f"{hap2_mean:.10g}" if hap2_mean is not None else NA,
                    "posterior_threshold": f"{posterior_threshold:.10g}",
                    "posterior_pass": bool_text(posterior_pass),
                    "minimum_marker_count": minimum_markers,
                    "marker_count_pass": bool_text(marker_pass),
                    "tract_high_confidence": bool_text(high_confidence),
                    "putatively_introgressed": bool_text(high_confidence and bool(donor_haplotypes)),
                    "strict_ancestry_heterozygous_target": bool_text(strict_target),
                    "introgressed_ancestry": introgressed_ancestry,
                    "background_ancestry": background_ancestry,
                    "introgressed_rfmix_haplotype": introgressed_haplotype,
                    "background_rfmix_haplotype": background_haplotype,
                    "tract_filter_reason": reason,
                }
            )
            output.append(tract)
    return output


def main():
    args = parse_args()
    if len(args.msp) != len(args.fb):
        raise ValueError("--msp and --fb must contain the same number of files")
    if not 0 <= args.posterior_threshold <= 1:
        raise ValueError("--posterior-threshold must be between 0 and 1")
    if args.minimum_markers < 1:
        raise ValueError("--minimum-markers must be positive")
    donors = {value.strip() for value in args.donor_ancestries.split(",") if value.strip()}
    backgrounds = {value.strip() for value in args.background_ancestries.split(",") if value.strip()}
    if donors & backgrounds:
        raise ValueError("donor and background ancestry sets must be disjoint")

    all_tracts = []
    observed_sample_chromosomes = set()
    for msp_path, fb_path in zip(args.msp, args.fb):
        msp_rows, samples, ancestry_codes = parse_msp(msp_path)
        chromosome = msp_rows[0]["chromosome"]
        duplicates = {(sample, chromosome) for sample in samples} & observed_sample_chromosomes
        if duplicates:
            raise ValueError(f"duplicate sample-chromosome inputs: {sorted(duplicates)}")
        observed_sample_chromosomes.update((sample, chromosome) for sample in samples)
        tracts, row_to_tract = build_tracts(msp_rows, samples)
        add_fb_posteriors(
            fb_path,
            msp_rows,
            samples,
            ancestry_codes,
            tracts,
            row_to_tract,
            args.posterior_sum_tolerance,
        )
        all_tracts.extend(
            finalize_tracts(
                tracts,
                args.posterior_threshold,
                args.minimum_markers,
                donors,
                backgrounds,
            )
        )

    all_tracts.sort(
        key=lambda row: (
            row["WGS_individual"],
            chromosome_sort_key(row["chromosome"]),
            row["tract_start"],
            row["tract_end"],
        )
    )
    strict = [row for row in all_tracts if row["strict_ancestry_heterozygous_target"] == "TRUE"]
    write_tsv(args.output, OUTPUT_COLUMNS, all_tracts)
    write_tsv(args.strict_output, OUTPUT_COLUMNS, strict)
    print(f"PASS merged_tracts={len(all_tracts)} strict_target_tracts={len(strict)}")


if __name__ == "__main__":
    main()
