#!/usr/bin/env python3

import argparse
from collections import defaultdict

from lib.common import load_phased_genotypes, normalize_chromosome, phased_alleles, read_tsv, write_tsv


INTERVAL_REQUIRED_COLUMNS = [
    "RNA_sample",
    "WGS_individual",
    "tissue",
    "gene_id",
    "gene_symbol",
    "gene_associated_interval_id",
    "gene_associated_interval",
    "tract_id",
    "chromosome",
    "introgressed_ancestry",
    "background_ancestry",
    "rfmix_hap1_ancestry",
    "rfmix_hap2_ancestry",
    "introgressed_rfmix_haplotype",
    "background_rfmix_haplotype",
    "introgressed_haplotype_count",
    "background_haplotype_count",
    "total_haplotype_count",
    "binomial_p",
    "BH_FDR",
    "ASE_direction",
    "significant",
]

OUTPUT_APPEND_COLUMNS = [
    "has_matched_significant_ASE_site",
    "n_matched_significant_ASE_sites",
    "n_site_direction_concordant",
    "n_site_direction_discordant",
    "site_row_type",
    "site_significance_status",
    "site_mapping_status",
    "site_id",
    "site_chromosome",
    "site_position",
    "REF_allele",
    "ALT_allele",
    "REF_count",
    "ALT_count",
    "site_FDR",
    "GT",
    "G0_allele",
    "G1_allele",
    "REF_allele_ancestry",
    "ALT_allele_ancestry",
    "introgressed_allele",
    "background_allele",
    "introgressed_allele_count",
    "background_allele_count",
    "site_direction",
    "gene_interval_direction",
    "direction_concordance",
]

SITE_NA_COLUMNS = [
    "site_significance_status",
    "site_mapping_status",
    "site_id",
    "site_chromosome",
    "site_position",
    "REF_allele",
    "ALT_allele",
    "REF_count",
    "ALT_count",
    "site_FDR",
    "GT",
    "G0_allele",
    "G1_allele",
    "REF_allele_ancestry",
    "ALT_allele_ancestry",
    "introgressed_allele",
    "background_allele",
    "introgressed_allele_count",
    "background_allele_count",
    "site_direction",
    "gene_interval_direction",
    "direction_concordance",
]


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Annotate significant interval-level local-ancestry ASE events with site-level ASE evidence. "
            "Site-level concordance is reported as an audit field and is not used to filter interval events."
        )
    )
    parser.add_argument("--site-ase", required=True)
    parser.add_argument("--interval-ase", required=True, help="Significant interval ASE table from 04_test_interval_ASE.py")
    parser.add_argument("--phased-genotypes", required=True, help="Phased VCF/VCF.GZ or normalized genotype TSV")
    parser.add_argument("--output", required=True)
    parser.add_argument("--site-fdr-threshold", type=float, default=0.05)
    return parser.parse_args()


def pick(row, names, label, required=True):
    for name in names:
        if name in row and row[name] not in ("", "NA"):
            return row[name]
    if required:
        raise ValueError(f"site ASE input lacks usable {label}; accepted columns: {names}")
    return "NA"


def normalized_site(row):
    chromosome = normalize_chromosome(pick(row, ["chromosome", "contig"], "chromosome"))
    position = int(pick(row, ["position", "POS"], "position"))
    ref = pick(row, ["REF", "REF_allele", "refAllele"], "REF allele").upper()
    alt = pick(row, ["ALT", "ALT_allele", "altAllele"], "ALT allele").upper()
    return {
        "RNA_sample": pick(row, ["RNA_sample"], "RNA_sample"),
        "WGS_individual": pick(row, ["WGS_individual", "WGS_sample"], "WGS individual"),
        "tissue": pick(row, ["tissue"], "tissue"),
        "gene_id": pick(row, ["gene_id"], "gene_id"),
        "tract_id": pick(row, ["tract_id"], "tract_id", required=False),
        "chromosome": chromosome,
        "position": position,
        "REF": ref,
        "ALT": alt,
        "ref_count": int(pick(row, ["refCount", "ref_count", "REF_count"], "REF count")),
        "alt_count": int(pick(row, ["altCount", "alt_count", "ALT_count"], "ALT count")),
        "fdr": float(pick(row, ["site_FDR", "ASE_site_BH_FDR", "FDR"], "site FDR")),
        "site_id": pick(row, ["ASE_site_id", "variantID", "site_id"], "site ID", required=False)
        if any(name in row for name in ["ASE_site_id", "variantID", "site_id"])
        else f"{chromosome}_{position}_{ref}_{alt}",
    }


def parse_interval(interval):
    chromosome, coordinate = interval.split(":", 1)
    start, end = coordinate.split("-", 1)
    return normalize_chromosome(chromosome), int(start), int(end)


def interval_key(row):
    return (
        row["RNA_sample"],
        row["WGS_individual"],
        row["gene_id"],
        row["gene_associated_interval_id"],
        row["tract_id"],
        row["introgressed_ancestry"],
        row["background_ancestry"],
    )


def candidate_key(site):
    return (site["RNA_sample"], site["WGS_individual"], site["tissue"], site["gene_id"], site["chromosome"])


def interval_candidate_key(row):
    return (
        row["RNA_sample"],
        row["WGS_individual"],
        row["tissue"],
        row["gene_id"],
        normalize_chromosome(row["chromosome"]),
    )


def site_na_payload(row_type, mapping_status="NA"):
    payload = {column: "NA" for column in SITE_NA_COLUMNS}
    payload["site_row_type"] = row_type
    payload["site_mapping_status"] = mapping_status
    payload["site_significance_status"] = "NO_MATCHED_SIGNIFICANT_ASE_SITE"
    return payload


def map_site_to_interval(site, interval, genotype):
    result = {
        "site_row_type": "MATCHED_SIGNIFICANT_ASE_SITE",
        "site_significance_status": "SIGNIFICANT_INPUT",
        "site_id": site["site_id"],
        "site_chromosome": site["chromosome"],
        "site_position": site["position"],
        "REF_allele": site["REF"],
        "ALT_allele": site["ALT"],
        "REF_count": site["ref_count"],
        "ALT_count": site["alt_count"],
        "site_FDR": f"{site['fdr']:.17g}",
        "gene_interval_direction": interval["ASE_direction"],
    }
    if genotype is None:
        result.update({column: "NA" for column in SITE_NA_COLUMNS if column not in result})
        result["site_mapping_status"] = "PHASED_GENOTYPE_NOT_FOUND"
        return result
    result["GT"] = genotype["GT"]
    if genotype["REF"] != site["REF"] or genotype["ALT"] != site["ALT"]:
        result.update({column: "NA" for column in SITE_NA_COLUMNS if column not in result})
        result["site_mapping_status"] = "ALLELE_MISMATCH"
        return result
    alleles = phased_alleles(site["REF"], site["ALT"], genotype["GT"])
    if alleles is None:
        result.update({column: "NA" for column in SITE_NA_COLUMNS if column not in result})
        result["site_mapping_status"] = "NON_PHASED_HETEROZYGOUS_GT"
        return result

    g0, g1 = alleles
    ancestry_by_haplotype = {
        "hap1": interval["rfmix_hap1_ancestry"],
        "hap2": interval["rfmix_hap2_ancestry"],
    }
    intro_haplotype = interval["introgressed_rfmix_haplotype"]
    background_haplotype = interval["background_rfmix_haplotype"]
    if intro_haplotype not in ancestry_by_haplotype or background_haplotype not in ancestry_by_haplotype:
        result.update({column: "NA" for column in SITE_NA_COLUMNS if column not in result})
        result["site_mapping_status"] = "INTERVAL_ANCESTRY_DIRECTION_UNAVAILABLE"
        return result

    allele_by_haplotype = {"hap1": g0, "hap2": g1}
    count_by_allele = {site["REF"]: site["ref_count"], site["ALT"]: site["alt_count"]}
    intro_allele = allele_by_haplotype[intro_haplotype]
    background_allele = allele_by_haplotype[background_haplotype]
    intro_count = count_by_allele[intro_allele]
    background_count = count_by_allele[background_allele]
    if intro_count > background_count:
        site_direction = "introgressed_higher"
    elif intro_count < background_count:
        site_direction = "introgressed_lower"
    else:
        site_direction = "equal"
    concordance = site_direction == interval["ASE_direction"] and site_direction != "equal"
    result.update(
        {
            "site_mapping_status": "PASS",
            "G0_allele": g0,
            "G1_allele": g1,
            "REF_allele_ancestry": ancestry_by_haplotype["hap1"]
            if g0 == site["REF"]
            else ancestry_by_haplotype["hap2"],
            "ALT_allele_ancestry": ancestry_by_haplotype["hap1"]
            if g0 == site["ALT"]
            else ancestry_by_haplotype["hap2"],
            "introgressed_allele": intro_allele,
            "background_allele": background_allele,
            "introgressed_allele_count": intro_count,
            "background_allele_count": background_count,
            "site_direction": site_direction,
            "direction_concordance": "TRUE" if concordance else "FALSE",
        }
    )
    return result


def main():
    args = parse_args()
    if not 0 <= args.site_fdr_threshold <= 1:
        raise ValueError("--site-fdr-threshold must be between 0 and 1")
    site_rows, _ = read_tsv(args.site_ase, label="site ASE table")
    interval_rows, interval_header = read_tsv(
        args.interval_ase, INTERVAL_REQUIRED_COLUMNS, "significant interval ASE table"
    )

    intervals = []
    intervals_by_key = defaultdict(list)
    for row in interval_rows:
        if row["significant"] != "TRUE":
            continue
        chromosome, start, end = parse_interval(row["gene_associated_interval"])
        row["_chromosome"], row["_start"], row["_end"] = chromosome, start, end
        intervals.append(row)
        intervals_by_key[interval_candidate_key(row)].append(row)

    normalized_sites = []
    requested = set()
    for row in site_rows:
        site = normalized_site(row)
        if min(site["ref_count"], site["alt_count"]) < 0:
            raise ValueError("site allele counts must be non-negative")
        if not 0 <= site["fdr"] <= 1:
            raise ValueError("site FDR must be between 0 and 1")
        if site["fdr"] >= args.site_fdr_threshold:
            continue
        normalized_sites.append(site)
        requested.add((site["WGS_individual"], site["chromosome"], site["position"]))
    genotypes = load_phased_genotypes(args.phased_genotypes, requested_positions=requested)

    output_by_interval = defaultdict(list)
    for site in normalized_sites:
        candidates = [
            row
            for row in intervals_by_key.get(candidate_key(site), [])
            if row["_start"] <= site["position"] <= row["_end"]
            and (site["tract_id"] == "NA" or row["tract_id"] == site["tract_id"])
        ]
        if len(candidates) != 1:
            continue
        interval = candidates[0]
        genotype = genotypes.get((site["WGS_individual"], site["chromosome"], site["position"]))
        mapped = map_site_to_interval(site, interval, genotype)
        output_by_interval[interval_key(interval)].append(mapped)

    output = []
    for interval in intervals:
        key = interval_key(interval)
        mapped_sites = output_by_interval.get(key, [])
        n_sites = len(mapped_sites)
        n_concordant = sum(site.get("direction_concordance") == "TRUE" for site in mapped_sites)
        n_discordant = sum(site.get("direction_concordance") == "FALSE" for site in mapped_sites)
        shared = {
            "has_matched_significant_ASE_site": "TRUE" if n_sites else "FALSE",
            "n_matched_significant_ASE_sites": str(n_sites),
            "n_site_direction_concordant": str(n_concordant),
            "n_site_direction_discordant": str(n_discordant),
        }
        if not mapped_sites:
            output.append({**interval, **shared, **site_na_payload("NO_MATCHED_SIGNIFICANT_ASE_SITE")})
            continue
        for site in mapped_sites:
            output.append({**interval, **shared, **site})

    output_header = interval_header + [column for column in OUTPUT_APPEND_COLUMNS if column not in interval_header]
    write_tsv(args.output, output_header, output)
    print(
        "PASS "
        f"significant_intervals={len(intervals)} "
        f"output_rows={len(output)} "
        f"matched_site_rows={sum(row['site_row_type'] == 'MATCHED_SIGNIFICANT_ASE_SITE' for row in output)} "
        f"no_site_placeholder_rows={sum(row['site_row_type'] == 'NO_MATCHED_SIGNIFICANT_ASE_SITE' for row in output)}"
    )


if __name__ == "__main__":
    main()
