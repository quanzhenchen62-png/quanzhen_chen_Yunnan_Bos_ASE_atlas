#!/usr/bin/env python3

import argparse
import re
from bisect import bisect_left
from collections import defaultdict

from lib.common import chromosome_sort_key, normalize_chromosome, open_text, read_tsv, write_tsv


TRACT_COLUMNS = [
    "WGS_individual",
    "chromosome",
    "tract_id",
    "tract_start",
    "tract_end",
    "tract_n_snps",
    "rfmix_hap1_ancestry",
    "rfmix_hap2_ancestry",
    "ordered_ancestry_pair",
    "introgressed_ancestry",
    "background_ancestry",
    "introgressed_rfmix_haplotype",
    "background_rfmix_haplotype",
    "tract_high_confidence",
    "strict_ancestry_heterozygous_target",
]

OUTPUT_COLUMNS = [
    "gene_associated_interval_id",
    "WGS_individual",
    "chromosome",
    "gene_associated_start",
    "gene_associated_end",
    "gene_associated_interval",
    "gene_id",
    "gene_symbol",
    "gene_biotype",
    "gene_start",
    "gene_end",
    "gene_strand",
    "tract_id",
    "tract_start",
    "tract_end",
    "tract_n_snps",
    "rfmix_hap1_ancestry",
    "rfmix_hap2_ancestry",
    "ordered_ancestry_pair",
    "introgressed_ancestry",
    "background_ancestry",
    "introgressed_rfmix_haplotype",
    "background_rfmix_haplotype",
]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Intersect high-confidence ancestry-heterozygous tracts with protein-coding gene bodies."
    )
    parser.add_argument("--tracts", required=True, help="Strict tract TSV from step 1")
    parser.add_argument("--gtf", required=True, help="ARS-UCD1.2-compatible GTF or GTF.GZ")
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def parse_attributes(text):
    return dict(re.findall(r'([A-Za-z0-9_]+)\s+"([^"]*)"', text))


def read_protein_coding_genes(path):
    genes = defaultdict(list)
    gene_ids = set()
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"GTF line {line_number} does not have nine columns")
            if fields[2] != "gene":
                continue
            attributes = parse_attributes(fields[8])
            biotype = attributes.get("gene_biotype") or attributes.get("gene_type")
            if biotype != "protein_coding":
                continue
            gene_id = attributes.get("gene_id")
            if not gene_id:
                raise ValueError(f"protein-coding gene lacks gene_id at GTF line {line_number}")
            if gene_id in gene_ids:
                raise ValueError(f"duplicate protein-coding gene_id in GTF: {gene_id}")
            gene_ids.add(gene_id)
            chromosome = normalize_chromosome(fields[0])
            start, end = int(fields[3]), int(fields[4])
            if start < 1 or end < start:
                raise ValueError(f"invalid GTF coordinates at line {line_number}")
            genes[chromosome].append(
                {
                    "chromosome": chromosome,
                    "gene_start": start,
                    "gene_end": end,
                    "gene_id": gene_id,
                    "gene_symbol": attributes.get("gene_name") or "NA",
                    "gene_biotype": biotype,
                    "gene_strand": fields[6],
                }
            )
    if not gene_ids:
        raise ValueError("no protein-coding gene features were found in the GTF")
    starts = {}
    for chromosome, rows in genes.items():
        rows.sort(key=lambda row: (row["gene_start"], row["gene_end"], row["gene_id"]))
        starts[chromosome] = [row["gene_start"] for row in rows]
    return genes, starts


def main():
    args = parse_args()
    tract_rows, _ = read_tsv(args.tracts, TRACT_COLUMNS, "strict tract table")
    genes, gene_starts = read_protein_coding_genes(args.gtf)
    output = []
    observed_ids = set()
    for tract in tract_rows:
        if tract["tract_high_confidence"] != "TRUE" or tract["strict_ancestry_heterozygous_target"] != "TRUE":
            raise ValueError(f"non-strict tract found in strict input: {tract['tract_id']}")
        chromosome = normalize_chromosome(tract["chromosome"])
        tract_start, tract_end = int(tract["tract_start"]), int(tract["tract_end"])
        chromosome_genes = genes.get(chromosome, [])
        limit = bisect_left(gene_starts.get(chromosome, []), tract_end + 1)
        for gene in chromosome_genes[:limit]:
            if gene["gene_end"] < tract_start:
                continue
            overlap_start = max(gene["gene_start"], tract_start)
            overlap_end = min(gene["gene_end"], tract_end)
            if overlap_start > overlap_end:
                continue
            interval_id = f"{tract['WGS_individual']}|{gene['gene_id']}|{tract['tract_id']}"
            if interval_id in observed_ids:
                raise ValueError(f"duplicate gene-associated interval: {interval_id}")
            observed_ids.add(interval_id)
            output.append(
                {
                    "gene_associated_interval_id": interval_id,
                    "WGS_individual": tract["WGS_individual"],
                    "chromosome": chromosome,
                    "gene_associated_start": overlap_start,
                    "gene_associated_end": overlap_end,
                    "gene_associated_interval": f"{chromosome}:{overlap_start}-{overlap_end}",
                    "gene_id": gene["gene_id"],
                    "gene_symbol": gene["gene_symbol"],
                    "gene_biotype": gene["gene_biotype"],
                    "gene_start": gene["gene_start"],
                    "gene_end": gene["gene_end"],
                    "gene_strand": gene["gene_strand"],
                    "tract_id": tract["tract_id"],
                    "tract_start": tract_start,
                    "tract_end": tract_end,
                    "tract_n_snps": tract["tract_n_snps"],
                    "rfmix_hap1_ancestry": tract["rfmix_hap1_ancestry"],
                    "rfmix_hap2_ancestry": tract["rfmix_hap2_ancestry"],
                    "ordered_ancestry_pair": tract["ordered_ancestry_pair"],
                    "introgressed_ancestry": tract["introgressed_ancestry"],
                    "background_ancestry": tract["background_ancestry"],
                    "introgressed_rfmix_haplotype": tract["introgressed_rfmix_haplotype"],
                    "background_rfmix_haplotype": tract["background_rfmix_haplotype"],
                }
            )
    output.sort(
        key=lambda row: (
            row["WGS_individual"],
            chromosome_sort_key(row["chromosome"]),
            row["gene_associated_start"],
            row["gene_associated_end"],
            row["gene_id"],
            row["tract_id"],
        )
    )
    write_tsv(args.output, OUTPUT_COLUMNS, output)
    print(f"PASS protein_coding_gene_associated_intervals={len(output)}")


if __name__ == "__main__":
    main()
