#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Prepare background-gene selection-overlap tables for tissue-specific and broadly expressed gene enrichment.")
    parser.add_argument("--background-genes", required=True, help="BED/TSV with gene coordinates after low-expression filtering.")
    parser.add_argument("--target-gene-classes", required=True, help="BED/TSV or gene-class TSV for tissue-specific and broadly expressed genes.")
    parser.add_argument("--target-manifest", required=True, help="TSV describing candidate-region and FST-decile files.")
    parser.add_argument("--outdir", required=True, help="Output directory.")
    parser.add_argument("--min-overlap-frac", type=float, default=0.20, help="Minimum gene-body overlap fraction.")
    return parser.parse_args()


def normalize_chr(series: pd.Series) -> pd.Series:
    return series.astype(str).str.replace(r"^chr", "", regex=True)


def normalize_comparison(x: str) -> str:
    return str(x).replace("vs.taurus", "vs.taurine")


def read_table(path: str | Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", compression="infer", dtype=str)


def load_manifest(path: str | Path) -> pd.DataFrame:
    df = read_table(path)
    required = {"target_kind", "comparison", "window_kb", "step_kb", "path"}
    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError(f"Target manifest missing columns: {', '.join(missing)}")
    df = df.copy()
    df["comparison"] = df["comparison"].map(normalize_comparison)
    df["window_kb"] = pd.to_numeric(df["window_kb"], errors="coerce").astype("Int64")
    df["step_kb"] = pd.to_numeric(df["step_kb"], errors="coerce").astype("Int64")
    if "threshold_label" not in df.columns:
        df["threshold_label"] = pd.NA
    return df


def load_background_genes(path: str | Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", compression="infer", header=None)
    if df.shape[1] < 4:
        raise ValueError("Background gene file must contain at least 4 columns: chr, start0, end0, gene")
    out = pd.DataFrame(
        {
            "chr": normalize_chr(df.iloc[:, 0]),
            "start0": pd.to_numeric(df.iloc[:, 1], errors="coerce"),
            "end0": pd.to_numeric(df.iloc[:, 2], errors="coerce"),
            "gene": df.iloc[:, 3].astype(str).str.replace(r"\..*$", "", regex=True),
        }
    )
    if df.shape[1] >= 5:
        out["length"] = pd.to_numeric(df.iloc[:, 4], errors="coerce")
    else:
        out["length"] = out["end0"] - out["start0"]
    out = out.dropna(subset=["chr", "start0", "end0", "gene"]).copy()
    out["start0"] = out["start0"].astype(int)
    out["end0"] = out["end0"].astype(int)
    out["length"] = pd.to_numeric(out["length"], errors="coerce").fillna(out["end0"] - out["start0"]).astype(int)
    out = out.sort_values(["chr", "start0", "end0", "gene"]).drop_duplicates(subset=["gene"], keep="first").reset_index(drop=True)
    return out


def load_target_gene_classes(path: str | Path) -> tuple[set[str], set[str]]:
    df = pd.read_csv(path, sep="\t", compression="infer", header=None)
    if df.shape[1] >= 5:
        cls = df.iloc[:, 3].astype(str)
        gene = df.iloc[:, 4].astype(str)
    elif df.shape[1] >= 2:
        gene = df.iloc[:, 0].astype(str)
        cls = df.iloc[:, 1].astype(str)
    else:
        raise ValueError("Target gene class file must have either 2 columns (gene, class) or 5 BED-like columns.")
    gene = gene.str.replace(r"\..*$", "", regex=True)
    cls = cls.str.strip()
    tissue_specific = set(gene[cls.eq("tissue_specific")])
    tissue_broadly = set(gene[cls.eq("tissue_broadly")])
    return tissue_specific, tissue_broadly


def read_candidate_bed(path: str | Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", header=None, compression="infer")
    out = df.iloc[:, :3].copy()
    out.columns = ["chr", "start0", "end0"]
    out["chr"] = normalize_chr(out["chr"])
    out["start0"] = pd.to_numeric(out["start0"], errors="coerce").astype("Int64")
    out["end0"] = pd.to_numeric(out["end0"], errors="coerce").astype("Int64")
    out = out.dropna(subset=["chr", "start0", "end0"]).copy()
    out["start0"] = out["start0"].astype(int)
    out["end0"] = out["end0"].astype(int)
    return out.sort_values(["chr", "start0", "end0"]).reset_index(drop=True)


def read_fst_deciles(path: str | Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", header=None, compression="infer")
    if df.shape[1] < 7:
        raise ValueError("FST-decile file must contain at least 7 columns.")
    out = pd.DataFrame(
        {
            "chr": normalize_chr(df.iloc[:, 0]),
            "start0": pd.to_numeric(df.iloc[:, 1], errors="coerce") - 1,
            "end0": pd.to_numeric(df.iloc[:, 2], errors="coerce"),
            "weighted_fst": pd.to_numeric(df.iloc[:, 5], errors="coerce"),
            "decile": pd.to_numeric(df.iloc[:, 6], errors="coerce"),
        }
    )
    out = out.dropna(subset=["chr", "start0", "end0", "weighted_fst", "decile"]).copy()
    out["start0"] = out["start0"].astype(int)
    out["end0"] = out["end0"].astype(int)
    out["decile"] = out["decile"].astype(int)
    return out.sort_values(["chr", "start0", "end0", "weighted_fst"], ascending=[True, True, True, False]).reset_index(drop=True)


def overlap_bp(start1: int, end1: int, start2: int, end2: int) -> int:
    return max(0, min(end1, end2) - max(start1, start2))


def assign_candidate_genes(genes: pd.DataFrame, intervals: pd.DataFrame, min_frac: float) -> pd.Series:
    out = pd.Series(0, index=genes.index, dtype="int64")
    for chrom, gsub in genes.groupby("chr", sort=False):
        isub = intervals.loc[intervals["chr"] == chrom]
        if isub.empty:
            continue
        starts = isub["start0"].to_numpy()
        ends = isub["end0"].to_numpy()
        for idx, row in gsub.iterrows():
            left = np.searchsorted(starts, row["end0"], side="right")
            assigned = 0
            j = left - 1
            while j >= 0 and ends[j] > row["start0"]:
                ov = overlap_bp(row["start0"], row["end0"], int(starts[j]), int(ends[j]))
                if ov / max(int(row["length"]), 1) >= min_frac:
                    assigned = 1
                    break
                j -= 1
            out.at[idx] = assigned
    return out


def assign_fst_gene_deciles(genes: pd.DataFrame, windows: pd.DataFrame, min_frac: float, window_bp: int) -> pd.Series:
    out = pd.Series(pd.NA, index=genes.index, dtype="Int64")
    for chrom, gsub in genes.groupby("chr", sort=False):
        wsub = windows.loc[windows["chr"] == chrom]
        if wsub.empty:
            continue
        starts = wsub["start0"].to_numpy()
        ends = wsub["end0"].to_numpy()
        weights = wsub["weighted_fst"].to_numpy(dtype=float)
        deciles = wsub["decile"].to_numpy(dtype=int)
        for idx, row in gsub.iterrows():
            left = np.searchsorted(starts, int(row["end0"]), side="right")
            cand = []
            j = left - 1
            while j >= 0 and starts[j] < int(row["end0"]):
                ov = overlap_bp(int(row["start0"]), int(row["end0"]), int(starts[j]), int(ends[j]))
                if ov / max(int(row["length"]), 1) >= min_frac:
                    cand.append(j)
                if int(starts[j]) < int(row["start0"]) - window_bp and int(ends[j]) <= int(row["start0"]):
                    break
                j -= 1
            if not cand:
                continue
            cand = np.array(cand, dtype=int)
            best = cand[np.argmax(weights[cand])]
            out.at[idx] = int(deciles[best])
    return out


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    genes = load_background_genes(args.background_genes)
    tissue_specific, tissue_broadly = load_target_gene_classes(args.target_gene_classes)
    manifest = load_manifest(args.target_manifest)

    genes["is_tissue_specific"] = genes["gene"].isin(tissue_specific).astype(int)
    genes["is_tissue_broadly"] = genes["gene"].isin(tissue_broadly).astype(int)

    for row in manifest.itertuples(index=False):
        comp = normalize_comparison(row.comparison)
        prefix = f"{comp}_{int(row.window_kb)}_{int(row.step_kb)}"
        path = Path(str(row.path))
        if row.target_kind == "candidate_region":
            label = str(row.threshold_label)
            if label in {"", "nan", "<NA>"}:
                raise ValueError(f"threshold_label is required for candidate-region rows: {path}")
            col = f"in_candidate_{prefix}_{label}"
            genes[col] = assign_candidate_genes(genes, read_candidate_bed(path), args.min_overlap_frac)
        elif row.target_kind == "fst_decile":
            deciles = assign_fst_gene_deciles(genes, read_fst_deciles(path), args.min_overlap_frac, int(row.window_kb) * 1000)
            for d in range(1, 11):
                genes[f"in_fst_decile_{prefix}_d{d}"] = deciles.eq(d).astype(int)
        else:
            raise ValueError(f"Unsupported target_kind: {row.target_kind}")

    genes.to_csv(outdir / "background_gene_selection_table.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
