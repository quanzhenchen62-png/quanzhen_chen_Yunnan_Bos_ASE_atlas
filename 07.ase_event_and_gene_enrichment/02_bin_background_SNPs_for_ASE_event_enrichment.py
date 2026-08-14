#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


REQUIRED_COLS = {"sample", "chr", "pos", "SNP_ID", "gene", "ASE_status", "MAF", "heterozygosity", "TPM", "totalCount"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Filter and decile-bin per-sample SNP tables for ind-ASE enrichment.")
    parser.add_argument("--sample-meta", required=True, help="TSV with Sample and prepared_file columns.")
    parser.add_argument("--outdir", required=True, help="Output directory.")
    parser.add_argument("--n-bins", type=int, default=10, help="Number of rank bins per covariate.")
    return parser.parse_args()


def rank_bin(values: pd.Series, n_bins: int) -> pd.Series:
    x = pd.to_numeric(values, errors="coerce")
    out = pd.Series(pd.NA, index=values.index, dtype="Int64")
    ok = x.notna() & np.isfinite(x)
    if ok.sum() == 0:
        return out
    ranks = x.loc[ok].rank(method="average")
    bins = np.ceil(ranks / ranks.max() * n_bins).astype(int)
    bins[bins < 1] = 1
    bins[bins > n_bins] = n_bins
    out.loc[ok] = bins.astype("int64")
    return out


def clean_one(path: str | Path, n_bins: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    df = pd.read_csv(path, sep="\t", compression="infer")
    missing = sorted(REQUIRED_COLS - set(df.columns))
    if missing:
        raise ValueError(f"Input SNP table missing columns: {', '.join(missing)}")

    df = df.copy()
    df["MAF"] = pd.to_numeric(df["MAF"], errors="coerce")
    df["heterozygosity"] = pd.to_numeric(df["heterozygosity"], errors="coerce")
    df["TPM"] = pd.to_numeric(df["TPM"], errors="coerce")
    df["totalCount"] = pd.to_numeric(df["totalCount"], errors="coerce")
    df["ASE_status"] = pd.to_numeric(df["ASE_status"], errors="coerce").fillna(0).astype(int)

    keep = (
        df["gene"].notna()
        & df["MAF"].between(0, 0.5, inclusive="both")
        & df["heterozygosity"].between(0, 1, inclusive="both")
        & df["TPM"].notna()
        & (df["TPM"] >= 0)
        & df["totalCount"].notna()
        & (df["totalCount"] > 0)
    )
    clean = df.loc[keep].copy()
    clean["log1p_TPM"] = np.log1p(clean["TPM"])
    clean["log1p_totalCount"] = np.log1p(clean["totalCount"])
    clean["MAF_bin"] = rank_bin(clean["MAF"], n_bins)
    clean["heterozygosity_bin"] = rank_bin(clean["heterozygosity"], n_bins)
    clean["TPM_bin"] = rank_bin(clean["log1p_TPM"], n_bins)
    clean["totalCount_bin"] = rank_bin(clean["log1p_totalCount"], n_bins)
    clean["match_bin"] = (
        "MAF" + clean["MAF_bin"].astype(str)
        + "_Ho" + clean["heterozygosity_bin"].astype(str)
        + "_TPM" + clean["TPM_bin"].astype(str)
        + "_Depth" + clean["totalCount_bin"].astype(str)
    )

    target_cols = sorted([c for c in clean.columns if c.startswith("in_candidate_") or c.startswith("in_fst_decile_")])
    keep_cols = [
        "sample", "chr", "pos", "SNP_ID", "gene", "ASE_status", *target_cols,
        "MAF", "heterozygosity", "TPM", "totalCount",
        "log1p_TPM", "log1p_totalCount",
        "MAF_bin", "heterozygosity_bin", "TPM_bin", "totalCount_bin", "match_bin",
    ]
    clean = clean[keep_cols].copy()

    summary = pd.DataFrame(
        [
            {
                "sample": clean["sample"].iloc[0] if not clean.empty else df["sample"].astype(str).iloc[0],
                "n_raw": int(df.shape[0]),
                "n_clean": int(clean.shape[0]),
                "n_removed": int(df.shape[0] - clean.shape[0]),
                "n_ASE_clean": int(clean["ASE_status"].sum()) if not clean.empty else 0,
                "n_match_bins": int(clean["match_bin"].nunique()) if not clean.empty else 0,
                "n_bins_setting": n_bins,
            }
        ]
    )
    return clean, summary


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    clean_dir = outdir / "clean_binned"
    summary_dir = outdir / "summary"
    clean_dir.mkdir(parents=True, exist_ok=True)
    summary_dir.mkdir(parents=True, exist_ok=True)

    meta = pd.read_csv(args.sample_meta, sep="\t", compression="infer", dtype=str)
    required = {"Sample", "prepared_file"}
    missing = sorted(required - set(meta.columns))
    if missing:
        raise ValueError(f"Sample metadata missing columns: {', '.join(missing)}")

    status_rows: list[dict[str, object]] = []
    for row in meta.itertuples(index=False):
        sample = str(row.Sample)
        infile = Path(str(row.prepared_file))
        clean_out = clean_dir / f"{sample}.clean_binned.tsv"
        try:
            clean, summary = clean_one(infile, args.n_bins)
            clean.to_csv(clean_out, sep="\t", index=False)
            summary.to_csv(summary_dir / f"{sample}.summary.tsv", sep="\t", index=False)
            status_rows.append(
                {
                    "Sample": sample,
                    "clean_binned_file": str(clean_out),
                    "n_clean": int(clean.shape[0]),
                    "n_ASE_clean": int(clean["ASE_status"].sum()) if not clean.empty else 0,
                    "status": "OK",
                    "message": "",
                }
            )
        except Exception as exc:
            status_rows.append(
                {
                    "Sample": sample,
                    "clean_binned_file": str(clean_out),
                    "n_clean": pd.NA,
                    "n_ASE_clean": pd.NA,
                    "status": "FAILED",
                    "message": str(exc),
                }
            )

    status = pd.DataFrame(status_rows)
    status.to_csv(outdir / "clean_binning_status.tsv", sep="\t", index=False)
    meta_out = meta.merge(status, on="Sample", how="left")
    meta_out.to_csv(outdir / "sample_meta.with_clean_binned.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
