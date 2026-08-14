#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prepare per-sample SNP tables for ind-ASE enrichment in candidate selection regions and FST-decile windows."
    )
    parser.add_argument("--sample-meta", required=True, help="TSV with Sample, breed, tissue, detectable_sites, sig_ase_sites, maf_heterozygosity_file and optionally tpm_file.")
    parser.add_argument("--target-manifest", required=True, help="TSV describing candidate-region and FST-decile interval files.")
    parser.add_argument("--outdir", required=True, help="Output directory.")
    return parser.parse_args()


def normalize_chr(series: pd.Series) -> pd.Series:
    return series.astype(str).str.replace(r"^chr", "", regex=True)


def normalize_comparison(x: str) -> str:
    return str(x).replace("vs.taurus", "vs.taurine")


def first_present(df: pd.DataFrame, names: Iterable[str], label: str) -> str:
    for name in names:
        if name in df.columns:
            return name
    raise ValueError(f"Missing required column for {label}: tried {', '.join(names)}")


def read_table(path: str | Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", compression="infer", dtype=str)


def load_detectable_sites(path: str | Path) -> pd.DataFrame:
    df = read_table(path)
    chrom_col = first_present(df, ["chr", "chrom", "CHROM", "contig"], "detectable chromosome")
    pos_col = first_present(df, ["pos", "position", "POS"], "detectable position")
    snp_col = first_present(df, ["SNP_ID", "variantID", "variant_id", "ID"], "detectable SNP ID")
    gene_col = first_present(df, ["gene", "Gene"], "detectable gene")
    depth_col = first_present(df, ["totalCount", "total_count", "depth"], "detectable total count")
    out = df[[chrom_col, pos_col, snp_col, gene_col, depth_col]].copy()
    out.columns = ["chr", "pos", "SNP_ID", "gene", "totalCount"]
    out["chr"] = normalize_chr(out["chr"])
    out["pos"] = pd.to_numeric(out["pos"], errors="coerce").astype("Int64")
    out["totalCount"] = pd.to_numeric(out["totalCount"], errors="coerce")
    out["gene"] = out["gene"].astype(str).str.replace(r"\..*$", "", regex=True)
    out = out.dropna(subset=["chr", "pos", "SNP_ID"]).copy()
    out["pos"] = out["pos"].astype(int)
    out["point0"] = out["pos"] - 1
    out = out.drop_duplicates(subset=["chr", "pos", "SNP_ID"], keep="first")
    return out


def load_sig_ase_ids(path: str | Path) -> set[str]:
    df = read_table(path)
    snp_col = first_present(df, ["SNP_ID", "variantID", "variant_id", "ID"], "significant ASE SNP ID")
    return set(df[snp_col].dropna().astype(str))


def load_maf_heterozygosity(path: str | Path) -> pd.DataFrame:
    df = read_table(path)
    chrom_col = first_present(df, ["chr", "chrom", "CHROM"], "MAF chromosome")
    pos_col = first_present(df, ["pos", "position", "POS"], "MAF position")
    maf_col = first_present(df, ["MAF", "maf"], "MAF")
    het_col = first_present(df, ["heterozygosity", "Ho", "ho"], "heterozygosity")
    out = df[[chrom_col, pos_col, maf_col, het_col]].copy()
    out.columns = ["chr", "pos", "MAF", "heterozygosity"]
    out["chr"] = normalize_chr(out["chr"])
    out["pos"] = pd.to_numeric(out["pos"], errors="coerce").astype("Int64")
    out["MAF"] = pd.to_numeric(out["MAF"], errors="coerce")
    out["heterozygosity"] = pd.to_numeric(out["heterozygosity"], errors="coerce")
    out = out.dropna(subset=["chr", "pos"]).copy()
    out["pos"] = out["pos"].astype(int)
    out = out.drop_duplicates(subset=["chr", "pos"], keep="first")
    return out


def load_tpm(path: str | Path) -> pd.DataFrame:
    df = read_table(path)
    gene_col = first_present(df, ["gene", "Gene"], "TPM gene")
    tpm_col = first_present(df, ["TPM", "tpm"], "TPM")
    out = df[[gene_col, tpm_col]].copy()
    out.columns = ["gene", "TPM"]
    out["gene"] = out["gene"].astype(str).str.replace(r"\..*$", "", regex=True)
    out["TPM"] = pd.to_numeric(out["TPM"], errors="coerce")
    out = out.drop_duplicates(subset=["gene"], keep="first")
    return out


def load_target_manifest(path: str | Path) -> pd.DataFrame:
    df = read_table(path)
    required = {"target_kind", "comparison", "window_kb", "step_kb", "path"}
    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError(f"Target manifest missing columns: {', '.join(missing)}")
    df = df.copy()
    df["target_kind"] = df["target_kind"].astype(str)
    df["comparison"] = df["comparison"].map(normalize_comparison)
    df["window_kb"] = pd.to_numeric(df["window_kb"], errors="coerce").astype("Int64")
    df["step_kb"] = pd.to_numeric(df["step_kb"], errors="coerce").astype("Int64")
    if "threshold_label" not in df.columns:
        df["threshold_label"] = pd.NA
    return df


def read_candidate_bed(path: str | Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", header=None, comment="#", compression="infer")
    if df.shape[1] < 3:
        raise ValueError(f"Candidate region file must have at least 3 columns: {path}")
    out = df.iloc[:, :3].copy()
    out.columns = ["chr", "start0", "end0"]
    out["chr"] = normalize_chr(out["chr"])
    out["start0"] = pd.to_numeric(out["start0"], errors="coerce").astype("Int64")
    out["end0"] = pd.to_numeric(out["end0"], errors="coerce").astype("Int64")
    out = out.dropna(subset=["chr", "start0", "end0"]).copy()
    out["start0"] = out["start0"].astype(int)
    out["end0"] = out["end0"].astype(int)
    out = out.sort_values(["chr", "start0", "end0"]).drop_duplicates().reset_index(drop=True)
    return out


def read_fst_decile_windows(path: str | Path) -> pd.DataFrame:
    header = pd.read_csv(path, sep="\t", nrows=5, compression="infer")
    lower = {str(c).lower(): c for c in header.columns}
    if {"chrom", "start_0based", "end_0based", "weighted_fst", "decile"} <= set(lower):
        df = pd.read_csv(path, sep="\t", compression="infer")
        out = pd.DataFrame(
            {
                "chr": normalize_chr(df[lower["chrom"]]),
                "start0": pd.to_numeric(df[lower["start_0based"]], errors="coerce"),
                "end0": pd.to_numeric(df[lower["end_0based"]], errors="coerce"),
                "weighted_fst": pd.to_numeric(df[lower["weighted_fst"]], errors="coerce"),
                "decile": pd.to_numeric(df[lower["decile"]], errors="coerce"),
            }
        )
    else:
        df = pd.read_csv(path, sep="\t", header=None, compression="infer")
        if df.shape[1] < 7:
            raise ValueError(f"FST-decile file must have >=7 columns or a header: {path}")
        out = pd.DataFrame(
            {
                "chr": normalize_chr(df.iloc[:, 0]),
                "start1": pd.to_numeric(df.iloc[:, 1], errors="coerce"),
                "end1": pd.to_numeric(df.iloc[:, 2], errors="coerce"),
                "weighted_fst": pd.to_numeric(df.iloc[:, 5], errors="coerce"),
                "decile": pd.to_numeric(df.iloc[:, 6], errors="coerce"),
            }
        )
        out["start0"] = out["start1"] - 1
        out["end0"] = out["end1"]
        out = out.drop(columns=["start1", "end1"])
    out = out.dropna(subset=["chr", "start0", "end0", "weighted_fst", "decile"]).copy()
    out["start0"] = out["start0"].astype(int)
    out["end0"] = out["end0"].astype(int)
    out["decile"] = out["decile"].astype(int)
    out = out.sort_values(["chr", "start0", "end0", "weighted_fst"], ascending=[True, True, True, False]).reset_index(drop=True)
    return out


def assign_candidate_points(point_df: pd.DataFrame, interval_df: pd.DataFrame) -> pd.Series:
    result = pd.Series(0, index=point_df.index, dtype="int64")
    if point_df.empty or interval_df.empty:
        return result
    for chrom, pts in point_df.groupby("chr", sort=False):
        sub = interval_df.loc[interval_df["chr"] == chrom].copy()
        if sub.empty:
            continue
        starts = sub["start0"].to_numpy()
        ends = sub["end0"].to_numpy()
        point0 = pts["point0"].to_numpy()
        idx = np.searchsorted(starts, point0, side="right") - 1
        ok = (idx >= 0) & (point0 < ends[np.clip(idx, 0, len(ends) - 1)])
        if ok.any():
            result.loc[pts.index[ok]] = 1
    return result


def assign_fst_deciles(point_df: pd.DataFrame, interval_df: pd.DataFrame, window_bp: int) -> pd.Series:
    result = pd.Series(pd.NA, index=point_df.index, dtype="Int64")
    if point_df.empty or interval_df.empty:
        return result
    for chrom, pts in point_df.groupby("chr", sort=False):
        sub = interval_df.loc[interval_df["chr"] == chrom].copy()
        if sub.empty:
            continue
        starts = sub["start0"].to_numpy()
        ends = sub["end0"].to_numpy()
        weights = sub["weighted_fst"].to_numpy(dtype=float)
        deciles = sub["decile"].to_numpy(dtype=int)
        for idx_row, point0 in zip(pts.index, pts["point0"].to_numpy()):
            left = np.searchsorted(starts, max(point0 - window_bp + 1, 0), side="left")
            right = np.searchsorted(starts, point0, side="right")
            if right <= left:
                continue
            cand = np.arange(left, right)
            cand = cand[point0 < ends[cand]]
            if len(cand) == 0:
                continue
            best = cand[np.argmax(weights[cand])]
            result.at[idx_row] = int(deciles[best])
    return result


def build_output(sample_row: pd.Series, manifest: pd.DataFrame) -> pd.DataFrame:
    detectable = load_detectable_sites(sample_row["detectable_sites"])
    sig_ids = load_sig_ase_ids(sample_row["sig_ase_sites"])
    maf = load_maf_heterozygosity(sample_row["maf_heterozygosity_file"])
    tpm = load_tpm(sample_row["tpm_file"])

    out = detectable.merge(maf, on=["chr", "pos"], how="left").merge(tpm, on="gene", how="left")
    out.insert(0, "sample", str(sample_row["Sample"]))
    out["ASE_status"] = out["SNP_ID"].astype(str).isin(sig_ids).astype(int)

    for row in manifest.itertuples(index=False):
        comp = normalize_comparison(row.comparison)
        prefix = f"{comp}_{int(row.window_kb)}_{int(row.step_kb)}"
        target_path = Path(str(row.path))
        kind = str(row.target_kind)
        if kind == "candidate_region":
            label = str(row.threshold_label)
            if label in {"", "nan", "<NA>"}:
                raise ValueError(f"threshold_label is required for candidate_region rows: {target_path}")
            col = f"in_candidate_{prefix}_{label}"
            intervals = read_candidate_bed(target_path)
            out[col] = assign_candidate_points(out[["chr", "point0"]], intervals)
        elif kind == "fst_decile":
            windows = read_fst_decile_windows(target_path)
            deciles = assign_fst_deciles(out[["chr", "point0"]], windows, int(row.window_kb) * 1000)
            for decile in range(1, 11):
                out[f"in_fst_decile_{prefix}_d{decile}"] = deciles.eq(decile).astype(int)
        else:
            raise ValueError(f"Unsupported target_kind: {kind}")

    keep_cols = ["sample", "chr", "pos", "SNP_ID", "gene", "ASE_status", "MAF", "heterozygosity", "TPM", "totalCount"]
    target_cols = [c for c in out.columns if c.startswith("in_candidate_") or c.startswith("in_fst_decile_")]
    keep_cols.extend(sorted(target_cols))
    return out[keep_cols].copy()


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    meta = pd.read_csv(args.sample_meta, sep="\t", compression="infer", dtype=str)
    required = {"Sample", "breed", "tissue", "detectable_sites", "sig_ase_sites", "maf_heterozygosity_file", "tpm_file"}
    missing = sorted(required - set(meta.columns))
    if missing:
        raise ValueError(f"Sample metadata missing columns: {', '.join(missing)}")

    manifest = load_target_manifest(args.target_manifest)
    status_rows: list[dict[str, object]] = []

    prepared_dir = outdir / "prepared_snp_tables"
    prepared_dir.mkdir(parents=True, exist_ok=True)

    for row in meta.itertuples(index=False):
        sample = str(row.Sample)
        outfile = prepared_dir / f"{sample}.selection_input.tsv"
        try:
            out = build_output(pd.Series(row._asdict()), manifest)
            out.to_csv(outfile, sep="\t", index=False)
            status_rows.append(
                {
                    "Sample": sample,
                    "breed": row.breed,
                    "tissue": row.tissue,
                    "prepared_file": str(outfile),
                    "n_detectable_sites": int(out.shape[0]),
                    "n_sig_ase_sites": int(out["ASE_status"].sum()),
                    "status": "OK",
                    "message": "",
                }
            )
        except Exception as exc:
            status_rows.append(
                {
                    "Sample": sample,
                    "breed": row.breed,
                    "tissue": row.tissue,
                    "prepared_file": str(outfile),
                    "n_detectable_sites": pd.NA,
                    "n_sig_ase_sites": pd.NA,
                    "status": "FAILED",
                    "message": str(exc),
                }
            )

    status = pd.DataFrame(status_rows)
    status.to_csv(outdir / "snp_preparation_status.tsv", sep="\t", index=False)

    meta_out = meta.merge(status[["Sample", "prepared_file", "status"]], on="Sample", how="left")
    meta_out.to_csv(outdir / "sample_meta.with_prepared_snp_tables.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
