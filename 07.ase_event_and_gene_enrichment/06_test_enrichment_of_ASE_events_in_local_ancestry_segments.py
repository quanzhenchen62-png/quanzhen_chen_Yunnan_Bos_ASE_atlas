#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


TEST_CLASSES = ("introgressed", "non_introgressed")
REQUIRED_META_COLS = {"breed", "tissue", "Sample", "Individual"}
REQUIRED_SITE_COLS = {"chr", "pos", "SNP_ID", "ASE_status", "MAF", "heterozygosity", "TPM", "totalCount"}
REQUIRED_TRACT_COLS = {
    "sample_id", "chromosome", "tract_start", "tract_end", "tract_n_snps",
    "rfmix_hap1_ancestry", "rfmix_hap2_ancestry", "ordered_ancestry_pair",
    "unordered_ancestry_pair", "yak_dosage", "gayal_dosage", "tract_high_confidence",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run matched ind-ASE enrichment in high-confidence local-ancestry segments.")
    parser.add_argument("--meta", required=True, help="TSV with breed, tissue, Sample and Individual.")
    parser.add_argument("--prepared-dir", required=True, help="Directory containing <Sample>.selection_input.tsv or equivalent per-sample SNP tables.")
    parser.add_argument("--tracts", required=True, help="High-confidence diploid tract table.")
    parser.add_argument("--outdir", required=True, help="Output directory.")
    parser.add_argument("--n-iter", type=int, default=10000, help="Number of null iterations.")
    parser.add_argument("--n-bins", type=int, default=10, help="Number of rank bins per covariate.")
    parser.add_argument("--min-tract-snps", type=int, default=30, help="Minimum marker count per eligible tract.")
    parser.add_argument("--seed", type=int, default=20260803, help="Random seed.")
    return parser.parse_args()


def mkdir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def validate_columns(df: pd.DataFrame, required: set[str], label: str) -> None:
    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError(f"{label} missing columns: {', '.join(missing)}")


def normalize_chr(series: pd.Series) -> pd.Series:
    return series.astype(str).str.replace(r"^chr", "", regex=True)


def as_bool(series: pd.Series) -> pd.Series:
    return series.astype(str).str.upper().isin({"TRUE", "T", "1", "YES"})


def bh_adjust(pvals: pd.Series) -> np.ndarray:
    x = pd.to_numeric(pvals, errors="coerce").to_numpy(dtype=float)
    out = np.full(len(x), np.nan)
    ok = np.isfinite(x)
    if ok.sum() == 0:
        return out
    idx = np.where(ok)[0]
    order = idx[np.argsort(x[idx])]
    ranked = x[order] * ok.sum() / np.arange(1, ok.sum() + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out[order] = np.minimum(ranked, 1.0)
    return out


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


def build_match_bins(df: pd.DataFrame, n_bins: int) -> pd.DataFrame:
    out = df.copy()
    out["MAF"] = pd.to_numeric(out["MAF"], errors="coerce")
    out["heterozygosity"] = pd.to_numeric(out["heterozygosity"], errors="coerce")
    out["TPM"] = pd.to_numeric(out["TPM"], errors="coerce")
    out["totalCount"] = pd.to_numeric(out["totalCount"], errors="coerce")
    keep = (
        out["MAF"].between(0, 0.5, inclusive="both")
        & out["heterozygosity"].between(0, 1, inclusive="both")
        & out["TPM"].notna()
        & (out["TPM"] >= 0)
        & out["totalCount"].notna()
        & (out["totalCount"] > 0)
    )
    out = out.loc[keep].copy()
    out["maf_bin"] = rank_bin(out["MAF"], n_bins)
    out["het_bin"] = rank_bin(out["heterozygosity"], n_bins)
    out["tpm_bin"] = rank_bin(np.log1p(pd.to_numeric(out["TPM"], errors="coerce")), n_bins)
    out["depth_bin"] = rank_bin(np.log1p(pd.to_numeric(out["totalCount"], errors="coerce")), n_bins)
    out["match_bin"] = (
        "MAF" + out["maf_bin"].astype(str)
        + "_Ho" + out["het_bin"].astype(str)
        + "_TPM" + out["tpm_bin"].astype(str)
        + "_Depth" + out["depth_bin"].astype(str)
    )
    return out


def assign_tracts_to_sites(site_df: pd.DataFrame, tract_df: pd.DataFrame, min_tract_snps: int) -> pd.DataFrame:
    out = site_df.copy()
    add_cols = [
        "tract_n_snps", "rfmix_hap1_ancestry", "rfmix_hap2_ancestry",
        "ordered_ancestry_pair", "unordered_ancestry_pair",
        "yak_dosage", "gayal_dosage", "tract_high_confidence",
    ]
    for col in add_cols:
        out[col] = pd.NA
    if tract_df.empty:
        out["ancestry_status"] = "unknown_no_tracts"
        return out

    tract_df = tract_df.copy()
    tract_df["chromosome"] = normalize_chr(tract_df["chromosome"])
    out["chr"] = normalize_chr(out["chr"])

    pieces = []
    for chrom, gsub in out.groupby("chr", sort=False):
        tsub = tract_df.loc[tract_df["chromosome"] == chrom].copy()
        if tsub.empty:
            gsub = gsub.copy()
            gsub["ancestry_status"] = "unknown_no_chr_tract"
            pieces.append(gsub)
            continue
        tsub["tract_start"] = pd.to_numeric(tsub["tract_start"], errors="coerce").astype(int)
        tsub["tract_end"] = pd.to_numeric(tsub["tract_end"], errors="coerce").astype(int)
        starts = tsub["tract_start"].to_numpy()
        ends = tsub["tract_end"].to_numpy()
        rows = []
        for idx, row in gsub.iterrows():
            pos = int(row["pos"])
            center = np.searchsorted(starts, pos, side="right") - 1
            hits = [i for i in (center - 1, center, center + 1) if 0 <= i < len(starts) and starts[i] <= pos <= ends[i]]
            if len(hits) != 1:
                rows.append({"index": idx, "ancestry_status": "boundary_ambiguous" if len(hits) > 1 else "unknown_gap"})
                continue
            hit = tsub.iloc[hits[0]]
            status = "eligible" if as_bool(pd.Series([hit["tract_high_confidence"]])).iloc[0] and int(hit["tract_n_snps"]) >= min_tract_snps else "unknown_low_confidence_or_short"
            rec = {"index": idx, "ancestry_status": status}
            for col in add_cols:
                rec[col] = hit[col]
            rows.append(rec)
        hit_df = pd.DataFrame(rows).set_index("index")
        block = gsub.join(hit_df, how="left", rsuffix="_hit")
        for col in add_cols + ["ancestry_status"]:
            if f"{col}_hit" in block.columns:
                block[col] = block[f"{col}_hit"]
                block = block.drop(columns=[f"{col}_hit"])
        pieces.append(block)
    return pd.concat(pieces, ignore_index=True)


def classify_segments(df: pd.DataFrame) -> pd.Series:
    yak = pd.to_numeric(df["yak_dosage"], errors="coerce").fillna(0).astype(int)
    gayal = pd.to_numeric(df["gayal_dosage"], errors="coerce").fillna(0).astype(int)
    total = yak + gayal
    hap1 = df["rfmix_hap1_ancestry"].astype(str)
    hap2 = df["rfmix_hap2_ancestry"].astype(str)
    non_intro = total.eq(0) & hap1.isin({"SAI", "Taurine"}) & hap2.isin({"SAI", "Taurine"})
    return pd.Series(np.where(total.gt(0), "introgressed", np.where(non_intro, "non_introgressed", "other")), index=df.index)


def empirical_test(df: pd.DataFrame, mask: pd.Series, n_iter: int, rng: np.random.Generator) -> tuple[int, float, float, int]:
    observed = int(((df["ASE_status"] == 1) & mask).sum())
    null = np.zeros(n_iter, dtype=float)
    n_strata = 0
    for (_, _), group in df.groupby(["Sample", "match_bin"], sort=False, dropna=False):
        n_total = int(group.shape[0])
        n_ase = int((group["ASE_status"] == 1).sum())
        n_target = int(mask.loc[group.index].sum())
        if n_total == 0 or n_ase == 0:
            continue
        null += rng.hypergeometric(ngood=n_target, nbad=n_total - n_target, nsample=n_ase, size=n_iter)
        n_strata += 1
    if n_strata == 0:
        return observed, np.nan, np.nan, 0
    mean_null = float(null.mean())
    p_emp = float((1 + np.sum(null >= observed)) / (n_iter + 1))
    return observed, mean_null, p_emp, n_strata


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    mkdir(outdir)

    meta = pd.read_csv(args.meta, sep="\t", compression="infer")
    validate_columns(meta, REQUIRED_META_COLS, "metadata")

    tracts = pd.read_csv(args.tracts, sep="\t", compression="infer")
    validate_columns(tracts, REQUIRED_TRACT_COLS, "tract table")

    result_rows: list[dict[str, object]] = []
    qc_rows: list[dict[str, object]] = []

    for (breed, tissue), group_meta in meta.groupby(["breed", "tissue"], sort=True):
        sample_frames = []
        missing_files = 0
        for row in group_meta.itertuples(index=False):
            site_path = Path(args.prepared_dir) / f"{row.Sample}.selection_input.tsv"
            if not site_path.exists():
                alt = Path(args.prepared_dir) / f"{row.Sample}.enrichment_input.tsv"
                site_path = alt if alt.exists() else site_path
            if not site_path.exists():
                missing_files += 1
                continue
            site_df = pd.read_csv(site_path, sep="\t", compression="infer")
            validate_columns(site_df, REQUIRED_SITE_COLS, f"site table {site_path}")
            site_df = site_df[list(REQUIRED_SITE_COLS)].copy()
            site_df["Sample"] = row.Sample
            site_df["Individual"] = row.Individual
            site_df["breed"] = breed
            site_df["tissue"] = tissue
            tract_df = tracts.loc[tracts["sample_id"].astype(str) == str(row.Individual)].copy()
            sample_frames.append(assign_tracts_to_sites(site_df, tract_df, args.min_tract_snps))
        if not sample_frames:
            qc_rows.append({"breed": breed, "tissue": tissue, "status": "no_usable_samples", "missing_input_files": missing_files})
            continue

        group_df = pd.concat(sample_frames, ignore_index=True)
        group_df["ASE_status"] = pd.to_numeric(group_df["ASE_status"], errors="coerce").fillna(0).astype(int)
        group_df["segment_class"] = classify_segments(group_df)
        eligible = group_df.loc[group_df["ancestry_status"].eq("eligible") & group_df["segment_class"].isin(TEST_CLASSES)].copy()
        qc_rows.append(
            {
                "breed": breed,
                "tissue": tissue,
                "status": "OK" if not eligible.empty else "no_eligible_sites",
                "n_total_sites": int(group_df.shape[0]),
                "n_eligible_sites": int(eligible.shape[0]),
                "n_samples": int(group_df["Sample"].nunique()),
                "missing_input_files": missing_files,
            }
        )
        if eligible.empty:
            continue

        eligible = build_match_bins(eligible, args.n_bins)
        rng = np.random.default_rng(args.seed + len(result_rows) + hash((breed, tissue)) % 100000)
        for segment_class in TEST_CLASSES:
            mask = eligible["segment_class"].eq(segment_class)
            observed, mean_null, p_emp, n_strata = empirical_test(eligible, mask, args.n_iter, rng)
            result_rows.append(
                {
                    "breed": breed,
                    "tissue": tissue,
                    "segment_class": segment_class,
                    "n_samples": int(eligible["Sample"].nunique()),
                    "n_eligible_sites": int(eligible.shape[0]),
                    "n_ASE_sites": int((eligible["ASE_status"] == 1).sum()),
                    "observed_overlap": observed,
                    "mean_random_overlap": mean_null,
                    "fold_enrichment": observed / mean_null if pd.notna(mean_null) and mean_null > 0 else np.nan,
                    "empirical_p": p_emp,
                    "n_sample_strata": n_strata,
                    "n_iter": args.n_iter,
                }
            )

    results = pd.DataFrame(result_rows)
    if not results.empty:
        results["BH_FDR"] = bh_adjust(results["empirical_p"])
    pd.DataFrame(qc_rows).to_csv(outdir / "local_ancestry_group_qc.tsv", sep="\t", index=False)
    results.to_csv(outdir / "local_ancestry_ase_enrichment.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
