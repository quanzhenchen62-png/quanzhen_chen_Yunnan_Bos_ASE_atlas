#!/usr/bin/env python3

"""
Identify phenotype-associated ASE variants with GCTA MLMA.

This script implements the analysis described in the manuscript section
"Identification of phenotype-associated ASE variants".

Implemented analyses
1. LD muscle fiber diameter
2. LD muscle fiber density
3. Relative hump weight

Method-constrained filtering
1. LD traits:
   concordant brd-ASE-3 variants located in altitude-associated
   expression-cluster genes that are also candidate introgressed genes.
2. Relative hump weight:
   concordant HM brd-ASE-3 variants located in HM-specific genes.

Assumptions
1. Input ASE candidate tables are already prepared.
2. Genotype BED/BIM/FAM, GRM, phenotype, covariate, and qcovariate files
   are already prepared.
3. LD phenotype/covariate files already correspond to Bos individuals.
4. Relative hump weight phenotype/covariate files already correspond to
   cattle with HM phenotypes.
5. qcovar files already encode the manuscript-specified PCs:
   LD traits = age + PC1-3
   Relative hump weight = age + PC1
"""

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Sequence


TRUTHY = {"1", "true", "t", "yes", "y"}
CONSISTENT = TRUTHY | {"consistent", "concordant"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run GCTA MLMA for phenotype-associated ASE variants."
    )
    parser.add_argument("--ld-ase-table", required=True,
                        help="TSV of LD brd-ASE-3 candidates with header.")
    parser.add_argument("--hm-ase-table", required=True,
                        help="TSV of HM brd-ASE-3 candidates with header.")
    parser.add_argument("--ld-bfile", required=True,
                        help="PLINK prefix for LD Bos genotypes.")
    parser.add_argument("--hm-bfile", required=True,
                        help="PLINK prefix for HM cattle genotypes.")
    parser.add_argument("--ld-diameter-grm", required=True,
                        help="GRM prefix for LD muscle fiber diameter.")
    parser.add_argument("--ld-density-grm", required=True,
                        help="GRM prefix for LD muscle fiber density.")
    parser.add_argument("--hm-relative-weight-grm", required=True,
                        help="GRM prefix for relative hump weight.")
    parser.add_argument("--ld-diameter-pheno", required=True,
                        help="GCTA phenotype file for LD muscle fiber diameter.")
    parser.add_argument("--ld-density-pheno", required=True,
                        help="GCTA phenotype file for LD muscle fiber density.")
    parser.add_argument("--hm-relative-weight-pheno", required=True,
                        help="GCTA phenotype file for relative hump weight.")
    parser.add_argument("--ld-covar", required=True,
                        help="GCTA covariate file for LD traits.")
    parser.add_argument("--ld-qcovar", required=True,
                        help="GCTA qcovariate file for LD traits.")
    parser.add_argument("--hm-covar", required=True,
                        help="GCTA covariate file for relative hump weight.")
    parser.add_argument("--hm-qcovar", required=True,
                        help="GCTA qcovariate file for relative hump weight.")
    parser.add_argument("--outdir", required=True,
                        help="Output directory.")
    parser.add_argument("--plink-bin", default="plink",
                        help="PLINK executable. Default: plink")
    parser.add_argument("--gcta-bin", default="gcta64",
                        help="GCTA executable. Default: gcta64")
    parser.add_argument("--threads", type=int, default=4,
                        help="Threads for PLINK/GCTA. Default: 4")
    parser.add_argument("--variant-id-col", default="variant_id",
                        help="Column name containing ASE variant ID.")
    parser.add_argument("--ld-consistency-col", default="direction_consistency",
                        help="LD concordance column.")
    parser.add_argument("--ld-altitude-gene-col", default="in_altitude_cluster_gene",
                        help="LD altitude-associated expression-cluster gene flag.")
    parser.add_argument("--ld-introgressed-gene-col", default="in_introgressed_gene",
                        help="LD candidate introgressed gene flag.")
    parser.add_argument("--hm-consistency-col", default="direction_consistency",
                        help="HM concordance column.")
    parser.add_argument("--hm-specific-gene-col", default="in_hm_specific_gene",
                        help="HM-specific gene flag.")
    return parser.parse_args()


def prefix_file(prefix: Path, suffix: str) -> Path:
    return Path(f"{prefix}{suffix}")


def must_exist(path_str: str, label: str) -> Path:
    path = Path(path_str)
    if not path.exists():
        raise FileNotFoundError(f"{label} not found: {path}")
    return path


def must_exist_prefix(prefix_str: str, label: str) -> Path:
    prefix = Path(prefix_str)
    required = [prefix_file(prefix, ".bed"), prefix_file(prefix, ".bim"), prefix_file(prefix, ".fam")]
    missing = [str(p) for p in required if not p.exists()]
    if missing:
        raise FileNotFoundError(f"{label} missing required PLINK files: {', '.join(missing)}")
    return prefix


def must_exist_grm_prefix(prefix_str: str, label: str) -> Path:
    prefix = Path(prefix_str)
    required = [prefix_file(prefix, ".grm.bin"), prefix_file(prefix, ".grm.N.bin"), prefix_file(prefix, ".grm.id")]
    missing = [str(p) for p in required if not p.exists()]
    if missing:
        raise FileNotFoundError(f"{label} missing required GRM files: {', '.join(missing)}")
    return prefix


def normalize_flag(value: str) -> str:
    return value.strip().lower()


def is_true(value: str) -> bool:
    return normalize_flag(value) in TRUTHY


def is_consistent(value: str) -> bool:
    return normalize_flag(value) in CONSISTENT


def read_tsv(path: Path) -> List[Dict[str, str]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"No header found in {path}")
        return list(reader)


def ensure_columns(rows: Sequence[Dict[str, str]], required: Sequence[str], label: str) -> None:
    if not rows:
        raise ValueError(f"{label} is empty")
    missing = [col for col in required if col not in rows[0]]
    if missing:
        raise ValueError(f"{label} missing columns: {', '.join(missing)}")


def unique_preserve_order(values: Iterable[str]) -> List[str]:
    seen = set()
    ordered: List[str] = []
    for value in values:
        if value not in seen:
            seen.add(value)
            ordered.append(value)
    return ordered


def filter_ld_variants(args: argparse.Namespace, rows: Sequence[Dict[str, str]]) -> List[Dict[str, str]]:
    required = [
        args.variant_id_col,
        args.ld_consistency_col,
        args.ld_altitude_gene_col,
        args.ld_introgressed_gene_col,
    ]
    ensure_columns(rows, required, "LD ASE table")
    return [
        row for row in rows
        if row[args.variant_id_col].strip()
        and is_consistent(row[args.ld_consistency_col])
        and is_true(row[args.ld_altitude_gene_col])
        and is_true(row[args.ld_introgressed_gene_col])
    ]


def filter_hm_variants(args: argparse.Namespace, rows: Sequence[Dict[str, str]]) -> List[Dict[str, str]]:
    required = [
        args.variant_id_col,
        args.hm_consistency_col,
        args.hm_specific_gene_col,
    ]
    ensure_columns(rows, required, "HM ASE table")
    return [
        row for row in rows
        if row[args.variant_id_col].strip()
        and is_consistent(row[args.hm_consistency_col])
        and is_true(row[args.hm_specific_gene_col])
    ]


def write_variant_table(rows: Sequence[Dict[str, str]], out_path: Path) -> None:
    if not rows:
        raise ValueError(f"No rows to write for {out_path}")
    fieldnames = list(rows[0].keys())
    with out_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_variant_ids(rows: Sequence[Dict[str, str]], variant_col: str, out_path: Path) -> List[str]:
    variant_ids = unique_preserve_order(row[variant_col].strip() for row in rows if row[variant_col].strip())
    with out_path.open("w") as handle:
        for variant_id in variant_ids:
            handle.write(f"{variant_id}\n")
    return variant_ids


def run_command(command: Sequence[str]) -> None:
    subprocess.run(command, check=True)


def extract_genotypes(plink_bin: str, bfile_prefix: Path, extract_file: Path, out_prefix: Path, threads: int) -> None:
    command = [
        plink_bin,
        "--bfile", str(bfile_prefix),
        "--extract", str(extract_file),
        "--make-bed",
        "--out", str(out_prefix),
        "--threads", str(threads),
    ]
    run_command(command)


def read_bim_variant_ids(bim_path: Path) -> List[str]:
    if not bim_path.exists():
        return []
    variant_ids: List[str] = []
    with bim_path.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                fields = line.rstrip("\n").split()
            if len(fields) >= 2:
                variant_ids.append(fields[1])
    return variant_ids


def write_missing_variants(requested: Sequence[str], observed: Sequence[str], out_path: Path) -> None:
    observed_set = set(observed)
    with out_path.open("w") as handle:
        handle.write("variant_id\n")
        for variant_id in requested:
            if variant_id not in observed_set:
                handle.write(f"{variant_id}\n")


def run_mlma(
    gcta_bin: str,
    bfile_prefix: Path,
    grm_prefix: Path,
    pheno: Path,
    covar: Path,
    qcovar: Path,
    out_prefix: Path,
    threads: int,
) -> None:
    command = [
        gcta_bin,
        "--mlma",
        "--bfile", str(bfile_prefix),
        "--grm", str(grm_prefix),
        "--pheno", str(pheno),
        "--covar", str(covar),
        "--qcovar", str(qcovar),
        "--out", str(out_prefix),
        "--thread-num", str(threads),
    ]
    run_command(command)


def bh_fdr(p_values: Sequence[float]) -> List[float]:
    n = len(p_values)
    order = sorted(range(n), key=lambda i: p_values[i])
    adjusted = [1.0] * n
    prev = 1.0
    for rank in range(n, 0, -1):
        idx = order[rank - 1]
        q = min(1.0, p_values[idx] * n / rank)
        prev = min(prev, q)
        adjusted[idx] = prev
    return adjusted


def add_fdr_to_mlma(mlma_path: Path, out_path: Path, fdr_threshold: float = 0.1) -> int:
    with mlma_path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"No header in {mlma_path}")
        rows = list(reader)
        fieldnames = list(reader.fieldnames)

    p_col = None
    for candidate in ("p", "P"):
        if candidate in fieldnames:
            p_col = candidate
            break
    if p_col is None:
        raise ValueError(f"No P-value column found in {mlma_path}")

    p_values = [float(row[p_col]) for row in rows]
    fdr_values = bh_fdr(p_values)

    out_fields = fieldnames + ["FDR", "phenotype_associated_ASE_variant"]
    significant = 0
    with out_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=out_fields, delimiter="\t")
        writer.writeheader()
        for row, fdr in zip(rows, fdr_values):
            row = dict(row)
            row["FDR"] = f"{fdr:.10g}"
            row["phenotype_associated_ASE_variant"] = "1" if fdr < fdr_threshold else "0"
            if fdr < fdr_threshold:
                significant += 1
            writer.writerow(row)
    return significant


def write_summary(summary_rows: Sequence[Dict[str, str]], out_path: Path) -> None:
    fieldnames = [
        "trait",
        "tested_variant_n",
        "genotyped_variant_n",
        "missing_from_genotype_n",
        "fdr_significant_variant_n",
    ]
    with out_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(summary_rows)


def require_executable(program: str) -> None:
    if shutil.which(program) is None and not Path(program).exists():
        raise FileNotFoundError(f"Executable not found: {program}")


def main() -> int:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    require_executable(args.plink_bin)
    require_executable(args.gcta_bin)

    ld_ase_table = must_exist(args.ld_ase_table, "LD ASE table")
    hm_ase_table = must_exist(args.hm_ase_table, "HM ASE table")
    ld_bfile = must_exist_prefix(args.ld_bfile, "LD genotype prefix")
    hm_bfile = must_exist_prefix(args.hm_bfile, "HM genotype prefix")
    ld_diameter_grm = must_exist_grm_prefix(args.ld_diameter_grm, "LD diameter GRM prefix")
    ld_density_grm = must_exist_grm_prefix(args.ld_density_grm, "LD density GRM prefix")
    hm_relative_weight_grm = must_exist_grm_prefix(args.hm_relative_weight_grm, "HM relative weight GRM prefix")
    ld_diameter_pheno = must_exist(args.ld_diameter_pheno, "LD diameter phenotype")
    ld_density_pheno = must_exist(args.ld_density_pheno, "LD density phenotype")
    hm_relative_weight_pheno = must_exist(args.hm_relative_weight_pheno, "HM relative weight phenotype")
    ld_covar = must_exist(args.ld_covar, "LD covariate file")
    ld_qcovar = must_exist(args.ld_qcovar, "LD qcovariate file")
    hm_covar = must_exist(args.hm_covar, "HM covariate file")
    hm_qcovar = must_exist(args.hm_qcovar, "HM qcovariate file")

    ld_rows = read_tsv(ld_ase_table)
    hm_rows = read_tsv(hm_ase_table)

    ld_filtered = filter_ld_variants(args, ld_rows)
    hm_filtered = filter_hm_variants(args, hm_rows)
    if not ld_filtered:
        raise ValueError("No LD variants passed manuscript-defined filters.")
    if not hm_filtered:
        raise ValueError("No HM variants passed manuscript-defined filters.")

    ld_dir = outdir / "LD_muscle"
    hm_dir = outdir / "HM_relative_hump_weight"
    ld_dir.mkdir(exist_ok=True)
    hm_dir.mkdir(exist_ok=True)

    write_variant_table(ld_filtered, ld_dir / "tested_brd_ASE3_variants.tsv")
    write_variant_table(hm_filtered, hm_dir / "tested_brd_ASE3_variants.tsv")

    ld_variant_ids = write_variant_ids(
        ld_filtered,
        args.variant_id_col,
        ld_dir / "tested_brd_ASE3_variant_ids.txt",
    )
    hm_variant_ids = write_variant_ids(
        hm_filtered,
        args.variant_id_col,
        hm_dir / "tested_brd_ASE3_variant_ids.txt",
    )

    ld_subset_prefix = ld_dir / "tested_brd_ASE3_genotypes"
    hm_subset_prefix = hm_dir / "tested_brd_ASE3_genotypes"
    extract_genotypes(args.plink_bin, ld_bfile, ld_dir / "tested_brd_ASE3_variant_ids.txt", ld_subset_prefix, args.threads)
    extract_genotypes(args.plink_bin, hm_bfile, hm_dir / "tested_brd_ASE3_variant_ids.txt", hm_subset_prefix, args.threads)

    ld_observed = read_bim_variant_ids(prefix_file(ld_subset_prefix, ".bim"))
    hm_observed = read_bim_variant_ids(prefix_file(hm_subset_prefix, ".bim"))
    write_missing_variants(ld_variant_ids, ld_observed, ld_dir / "variants_missing_from_genotype.tsv")
    write_missing_variants(hm_variant_ids, hm_observed, hm_dir / "variants_missing_from_genotype.tsv")

    if not ld_observed:
        raise ValueError("No LD variants were extracted from the genotype dataset.")
    if not hm_observed:
        raise ValueError("No HM variants were extracted from the genotype dataset.")

    ld_diameter_prefix = ld_dir / "LD_diameter"
    ld_density_prefix = ld_dir / "LD_density"
    hm_relative_weight_prefix = hm_dir / "HM_relative_hump_weight"

    run_mlma(
        args.gcta_bin,
        ld_subset_prefix,
        ld_diameter_grm,
        ld_diameter_pheno,
        ld_covar,
        ld_qcovar,
        ld_diameter_prefix,
        args.threads,
    )
    run_mlma(
        args.gcta_bin,
        ld_subset_prefix,
        ld_density_grm,
        ld_density_pheno,
        ld_covar,
        ld_qcovar,
        ld_density_prefix,
        args.threads,
    )
    run_mlma(
        args.gcta_bin,
        hm_subset_prefix,
        hm_relative_weight_grm,
        hm_relative_weight_pheno,
        hm_covar,
        hm_qcovar,
        hm_relative_weight_prefix,
        args.threads,
    )

    ld_diameter_sig = add_fdr_to_mlma(
        prefix_file(ld_diameter_prefix, ".mlma"),
        prefix_file(ld_diameter_prefix, ".FDR.tsv"),
    )
    ld_density_sig = add_fdr_to_mlma(
        prefix_file(ld_density_prefix, ".mlma"),
        prefix_file(ld_density_prefix, ".FDR.tsv"),
    )
    hm_relative_weight_sig = add_fdr_to_mlma(
        prefix_file(hm_relative_weight_prefix, ".mlma"),
        prefix_file(hm_relative_weight_prefix, ".FDR.tsv"),
    )

    summary_rows = [
        {
            "trait": "LD_muscle_fiber_diameter",
            "tested_variant_n": str(len(ld_variant_ids)),
            "genotyped_variant_n": str(len(ld_observed)),
            "missing_from_genotype_n": str(len(ld_variant_ids) - len(ld_observed)),
            "fdr_significant_variant_n": str(ld_diameter_sig),
        },
        {
            "trait": "LD_muscle_fiber_density",
            "tested_variant_n": str(len(ld_variant_ids)),
            "genotyped_variant_n": str(len(ld_observed)),
            "missing_from_genotype_n": str(len(ld_variant_ids) - len(ld_observed)),
            "fdr_significant_variant_n": str(ld_density_sig),
        },
        {
            "trait": "relative_hump_weight",
            "tested_variant_n": str(len(hm_variant_ids)),
            "genotyped_variant_n": str(len(hm_observed)),
            "missing_from_genotype_n": str(len(hm_variant_ids) - len(hm_observed)),
            "fdr_significant_variant_n": str(hm_relative_weight_sig),
        },
    ]
    write_summary(summary_rows, outdir / "phenotype_associated_ASE_variant_summary.tsv")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
