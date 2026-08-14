#!/usr/bin/env python3

import argparse
from collections import defaultdict

from lib.common import (
    NA,
    chromosome_sort_key,
    load_phased_genotypes,
    normalize_chromosome,
    parse_variant_token,
    phased_alleles,
    read_tsv,
    require_columns,
    write_tsv,
)


INTERVAL_COLUMNS = [
    "gene_associated_interval_id",
    "WGS_individual",
    "chromosome",
    "gene_associated_start",
    "gene_associated_end",
    "gene_associated_interval",
    "gene_id",
    "gene_symbol",
    "gene_biotype",
    "tract_id",
    "tract_start",
    "tract_end",
    "rfmix_hap1_ancestry",
    "rfmix_hap2_ancestry",
    "ordered_ancestry_pair",
    "introgressed_ancestry",
    "background_ancestry",
    "introgressed_rfmix_haplotype",
    "background_rfmix_haplotype",
]

OUTPUT_COLUMNS = [
    "RNA_sample",
    "WGS_individual",
    "tissue",
    "block_id",
    "chromosome",
    "block_start",
    "block_end",
    "variants",
    "variantCount",
    "haplotypeA",
    "haplotypeB",
    "aCount",
    "bCount",
    "totalCount",
    "blockGWPhase",
    "n_genotypes_found",
    "n_informative_variants",
    "n_A_to_G0",
    "n_A_to_G1",
    "n_missing_genotypes",
    "n_allele_mismatches",
    "n_non_phased_heterozygous",
    "informative_variant_ids",
    "RNA_haplotype_A_genomic_haplotype",
    "RNA_haplotype_B_genomic_haplotype",
    "RNA_haplotype_A_rfmix_haplotype",
    "RNA_haplotype_B_rfmix_haplotype",
    "RNA_haplotype_A_ancestry",
    "RNA_haplotype_B_ancestry",
    "phase_mapping_status",
    "tract_mapping_status",
    "gene_interval_mapping_status",
    "candidate_gene_interval_ids",
    "gene_associated_interval_id",
    "gene_associated_interval",
    "gene_id",
    "gene_symbol",
    "gene_biotype",
    "tract_id",
    "tract_start",
    "tract_end",
    "rfmix_hap1_ancestry",
    "rfmix_hap2_ancestry",
    "ordered_ancestry_pair",
    "introgressed_ancestry",
    "background_ancestry",
    "introgressed_rfmix_haplotype",
    "background_rfmix_haplotype",
    "rfmix_hap1_count",
    "rfmix_hap2_count",
    "introgressed_haplotype_count",
    "background_haplotype_count",
    "count_conservation_status",
    "strict_main_pass",
    "block_analysis_status",
]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Map phASER RNA haplotypes to phased WGS haplotypes and local ancestry."
    )
    parser.add_argument("--phaser-blocks", required=True)
    parser.add_argument("--phased-genotypes", required=True, help="Phased VCF/VCF.GZ or normalized genotype TSV")
    parser.add_argument("--gene-intervals", required=True)
    parser.add_argument("--rna-sample", required=True)
    parser.add_argument("--wgs-individual", required=True)
    parser.add_argument("--tissue", required=True)
    parser.add_argument("--audit-output", required=True)
    parser.add_argument("--strict-output", required=True)
    return parser.parse_args()


def resolve_block_schema(fieldnames):
    aliases = {
        "chromosome": ("chromosome", "contig"),
        "start": ("block_start", "start"),
        "end": ("block_end", "stop"),
    }
    resolved = {}
    observed = set(fieldnames or [])
    for label, choices in aliases.items():
        matches = [choice for choice in choices if choice in observed]
        if not matches:
            raise ValueError(f"phASER block table lacks a {label} column; accepted names: {choices}")
        resolved[label] = matches[0]
    require_columns(
        fieldnames,
        ["variants", "variantCount", "haplotypeA", "haplotypeB", "aCount", "bCount", "totalCount"],
        "phASER block table",
    )
    return resolved


def load_intervals(path, individual):
    rows, _ = read_tsv(path, INTERVAL_COLUMNS, "gene-associated interval table")
    intervals = [row for row in rows if row["WGS_individual"] == individual]
    if not intervals:
        raise ValueError(f"no gene-associated intervals found for {individual}")
    for row in intervals:
        row["_chromosome"] = normalize_chromosome(row["chromosome"])
        row["_interval_start"] = int(row["gene_associated_start"])
        row["_interval_end"] = int(row["gene_associated_end"])
        row["_tract_start"] = int(row["tract_start"])
        row["_tract_end"] = int(row["tract_end"])
    intervals.sort(
        key=lambda row: (
            chromosome_sort_key(row["_chromosome"]),
            row["_tract_start"],
            row["_tract_end"],
            row["gene_id"],
        )
    )
    return intervals


def map_phase(block, schema, genotypes, individual):
    variant_tokens = [token for token in block["variants"].split(",") if token]
    haplotype_a = [allele.upper() for allele in block["haplotypeA"].split(",") if allele]
    haplotype_b = [allele.upper() for allele in block["haplotypeB"].split(",") if allele]
    declared_count = int(block["variantCount"])
    if not (len(variant_tokens) == len(haplotype_a) == len(haplotype_b) == declared_count):
        return None, "VARIANT_LIST_LENGTH_MISMATCH", []

    found = informative = a_to_g0 = a_to_g1 = missing = mismatches = non_phased = 0
    informative_ids = []
    parsed_variants = []
    for token, allele_a, allele_b in zip(variant_tokens, haplotype_a, haplotype_b):
        chromosome, position, ref, alt = parse_variant_token(token)
        parsed_variants.append((chromosome, position, ref, alt))
        genotype = genotypes.get((individual, chromosome, position))
        if genotype is None:
            missing += 1
            continue
        found += 1
        if genotype["REF"] != ref or genotype["ALT"] != alt:
            mismatches += 1
            continue
        alleles = phased_alleles(ref, alt, genotype["GT"])
        if alleles is None:
            non_phased += 1
            continue
        g0, g1 = alleles
        if allele_a == g0 and allele_b == g1:
            informative += 1
            a_to_g0 += 1
            informative_ids.append(token)
        elif allele_a == g1 and allele_b == g0:
            informative += 1
            a_to_g1 += 1
            informative_ids.append(token)
        else:
            mismatches += 1

    if mismatches:
        status = "ALLELE_MISMATCH"
        orientation = None
    elif non_phased:
        status = "NON_PHASED_HETEROZYGOUS_VARIANT"
        orientation = None
    elif a_to_g0 and a_to_g1:
        status = "PHASE_CONFLICT"
        orientation = None
    elif informative == 0:
        status = "NO_INFORMATIVE_PHASED_VARIANT"
        orientation = None
    elif informative == 1:
        status = "PASS_SINGLE_SNP"
        orientation = "A_TO_G0" if a_to_g0 else "A_TO_G1"
    else:
        status = "PASS_MULTI_SNP"
        orientation = "A_TO_G0" if a_to_g0 else "A_TO_G1"
    return {
        "n_genotypes_found": found,
        "n_informative_variants": informative,
        "n_A_to_G0": a_to_g0,
        "n_A_to_G1": a_to_g1,
        "n_missing_genotypes": missing,
        "n_allele_mismatches": mismatches,
        "n_non_phased_heterozygous": non_phased,
        "informative_variant_ids": ",".join(informative_ids) if informative_ids else NA,
        "orientation": orientation,
        "parsed_variants": parsed_variants,
    }, status, parsed_variants


def map_tract_and_gene(chromosome, start, end, parsed_variants, intervals):
    chromosome_intervals = [row for row in intervals if row["_chromosome"] == chromosome]
    tract_records = {}
    for row in chromosome_intervals:
        tract_records[row["tract_id"]] = row
    overlapping_tracts = [
        row
        for row in tract_records.values()
        if start <= row["_tract_end"] and row["_tract_start"] <= end
    ]
    containing = [
        row for row in overlapping_tracts if row["_tract_start"] <= start and end <= row["_tract_end"]
    ]
    if len(containing) == 1:
        tract = containing[0]
        tract_status = "UNIQUE_FULLY_CONTAINED"
    elif len(containing) > 1:
        tract = None
        tract_status = "BOUNDARY_AMBIGUOUS"
    elif len(overlapping_tracts) > 1:
        tract = None
        tract_status = "CROSSOVER_BLOCK"
    elif len(overlapping_tracts) == 1:
        tract = None
        tract_status = "BOUNDARY_AMBIGUOUS"
    else:
        tract = None
        tract_status = "NO_TARGET_TRACT"

    if tract is None:
        return None, tract_status, [], "NOT_EVALUATED"
    candidates = []
    for row in chromosome_intervals:
        if row["tract_id"] != tract["tract_id"]:
            continue
        if not (start <= row["_interval_end"] and row["_interval_start"] <= end):
            continue
        if not any(
            variant_chromosome == chromosome
            and row["_interval_start"] <= position <= row["_interval_end"]
            for variant_chromosome, position, _, _ in parsed_variants
        ):
            continue
        candidates.append(row)
    if len(candidates) == 1:
        return tract, tract_status, candidates, "UNIQUE_GENE_INTERVAL"
    if not candidates:
        return tract, tract_status, candidates, "NO_GENE_INTERVAL"
    return tract, tract_status, candidates, "MULTIPLE_GENE_INTERVALS"


def build_row(block, schema, phase_result, phase_status, tract, tract_status, candidates, gene_status, args, index):
    chromosome = normalize_chromosome(block[schema["chromosome"]])
    start, end = int(block[schema["start"]]), int(block[schema["end"]])
    block_id = block.get("block_id") or f"{args.rna_sample}:chr{chromosome}:{start}-{end}:{index}"
    a_count, b_count, total_count = int(block["aCount"]), int(block["bCount"]), int(block["totalCount"])
    if min(a_count, b_count, total_count) < 0:
        raise ValueError(f"negative phASER count for {block_id}")
    if a_count + b_count != total_count:
        raise ValueError(f"phASER count conservation failed for {block_id}")
    row = {field: NA for field in OUTPUT_COLUMNS}
    row.update(
        {
            "RNA_sample": args.rna_sample,
            "WGS_individual": args.wgs_individual,
            "tissue": args.tissue,
            "block_id": block_id,
            "chromosome": chromosome,
            "block_start": start,
            "block_end": end,
            "variants": block["variants"],
            "variantCount": block["variantCount"],
            "haplotypeA": block["haplotypeA"],
            "haplotypeB": block["haplotypeB"],
            "aCount": block["aCount"],
            "bCount": block["bCount"],
            "totalCount": block["totalCount"],
            "blockGWPhase": block.get("blockGWPhase") or NA,
            "phase_mapping_status": phase_status,
            "tract_mapping_status": tract_status,
            "gene_interval_mapping_status": gene_status,
            "candidate_gene_interval_ids": ",".join(
                candidate["gene_associated_interval_id"] for candidate in candidates
            ) if candidates else NA,
        }
    )
    if phase_result is not None:
        for field in (
            "n_genotypes_found",
            "n_informative_variants",
            "n_A_to_G0",
            "n_A_to_G1",
            "n_missing_genotypes",
            "n_allele_mismatches",
            "n_non_phased_heterozygous",
            "informative_variant_ids",
        ):
            row[field] = phase_result[field]
    interval = candidates[0] if len(candidates) == 1 else tract
    if interval is not None:
        for field in (
            "gene_associated_interval_id",
            "gene_associated_interval",
            "gene_id",
            "gene_symbol",
            "gene_biotype",
            "tract_id",
            "tract_start",
            "tract_end",
            "rfmix_hap1_ancestry",
            "rfmix_hap2_ancestry",
            "ordered_ancestry_pair",
            "introgressed_ancestry",
            "background_ancestry",
            "introgressed_rfmix_haplotype",
            "background_rfmix_haplotype",
        ):
            row[field] = interval.get(field, NA)

    strict = phase_status == "PASS_MULTI_SNP" and tract_status == "UNIQUE_FULLY_CONTAINED" and gene_status == "UNIQUE_GENE_INTERVAL"
    if phase_result is not None and phase_result["orientation"] is not None:
        if phase_result["orientation"] == "A_TO_G0":
            a_genomic, b_genomic, a_rfmix, b_rfmix = "G0", "G1", "hap1", "hap2"
        else:
            a_genomic, b_genomic, a_rfmix, b_rfmix = "G1", "G0", "hap2", "hap1"
        row.update(
            {
                "RNA_haplotype_A_genomic_haplotype": a_genomic,
                "RNA_haplotype_B_genomic_haplotype": b_genomic,
                "RNA_haplotype_A_rfmix_haplotype": a_rfmix,
                "RNA_haplotype_B_rfmix_haplotype": b_rfmix,
            }
        )
        if interval is not None:
            ancestry = {
                "hap1": interval["rfmix_hap1_ancestry"],
                "hap2": interval["rfmix_hap2_ancestry"],
            }
            row["RNA_haplotype_A_ancestry"] = ancestry[a_rfmix]
            row["RNA_haplotype_B_ancestry"] = ancestry[b_rfmix]
            counts = {a_rfmix: a_count, b_rfmix: b_count}
            row["rfmix_hap1_count"] = counts["hap1"]
            row["rfmix_hap2_count"] = counts["hap2"]
            row["introgressed_haplotype_count"] = counts[interval["introgressed_rfmix_haplotype"]]
            row["background_haplotype_count"] = counts[interval["background_rfmix_haplotype"]]
            row["count_conservation_status"] = "PASS" if sum(counts.values()) == total_count else "FAIL"
    row["strict_main_pass"] = "TRUE" if strict else "FALSE"
    if strict:
        row["block_analysis_status"] = "PASS_STRICT_MULTI_SNP"
    elif phase_status == "PASS_SINGLE_SNP":
        row["block_analysis_status"] = "AUDIT_ONLY_SINGLE_SNP"
    elif phase_status != "PASS_MULTI_SNP":
        row["block_analysis_status"] = phase_status
    elif tract_status != "UNIQUE_FULLY_CONTAINED":
        row["block_analysis_status"] = tract_status
    else:
        row["block_analysis_status"] = gene_status
    return row


def main():
    args = parse_args()
    block_rows, block_header = read_tsv(args.phaser_blocks, label="phASER block table")
    schema = resolve_block_schema(block_header)
    intervals = load_intervals(args.gene_intervals, args.wgs_individual)
    requested = set()
    for block in block_rows:
        for token in block["variants"].split(","):
            if token:
                chromosome, position, _, _ = parse_variant_token(token)
                requested.add((args.wgs_individual, chromosome, position))
    genotypes = load_phased_genotypes(
        args.phased_genotypes,
        expected_individual=args.wgs_individual,
        requested_positions=requested,
    )

    output = []
    observed_ids = set()
    for index, block in enumerate(block_rows, start=1):
        chromosome = normalize_chromosome(block[schema["chromosome"]])
        start, end = int(block[schema["start"]]), int(block[schema["end"]])
        phase_result, phase_status, parsed_variants = map_phase(
            block, schema, genotypes, args.wgs_individual
        )
        tract, tract_status, candidates, gene_status = map_tract_and_gene(
            chromosome, start, end, parsed_variants, intervals
        )
        row = build_row(
            block,
            schema,
            phase_result,
            phase_status,
            tract,
            tract_status,
            candidates,
            gene_status,
            args,
            index,
        )
        if row["block_id"] in observed_ids:
            raise ValueError(f"duplicate block_id: {row['block_id']}")
        observed_ids.add(row["block_id"])
        output.append(row)
    output.sort(
        key=lambda row: (
            chromosome_sort_key(row["chromosome"]),
            int(row["block_start"]),
            int(row["block_end"]),
            row["block_id"],
        )
    )
    strict = [row for row in output if row["strict_main_pass"] == "TRUE"]
    write_tsv(args.audit_output, OUTPUT_COLUMNS, output)
    write_tsv(args.strict_output, OUTPUT_COLUMNS, strict)
    print(f"PASS audited_blocks={len(output)} strict_multi_snp_blocks={len(strict)}")


if __name__ == "__main__":
    main()
