#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run breed-by-tissue ind-ASE enrichment for candidate selection regions and FST-decile windows.")
    parser.add_argument("--sample-meta", required=True, help="TSV with Sample, breed, tissue and clean_binned_file.")
    parser.add_argument("--outdir", required=True, help="Output directory.")
    parser.add_argument("--n-iter", type=int, default=10000, help="Number of null iterations.")
    parser.add_argument("--seed", type=int, default=20260803, help="Random seed.")
    return parser.parse_args()


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


def parse_target_column(col: str) -> dict[str, object]:
    if col.startswith("in_candidate_"):
        body = col.removeprefix("in_candidate_")
        comparison, win, step, threshold = body.rsplit("_", 3)
        return {
            "analysis_type": "candidate_region",
            "comparison": comparison,
            "window_kb": int(win),
            "step_kb": int(step),
            "threshold_label": threshold,
            "fst_decile": pd.NA,
        }
    if col.startswith("in_fst_decile_"):
        body = col.removeprefix("in_fst_decile_")
        comparison, win, step, decile = body.rsplit("_", 3)
        return {
            "analysis_type": "fst_decile",
            "comparison": comparison,
            "window_kb": int(win),
            "step_kb": int(step),
            "threshold_label": pd.NA,
            "fst_decile": int(decile.removeprefix("d")),
        }
    raise ValueError(f"Unsupported target column: {col}")


def load_clean_table(path: str | Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", compression="infer")
    required = {"sample", "SNP_ID", "ASE_status", "match_bin"}
    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError(f"Clean table missing columns: {', '.join(missing)}")
    df["ASE_status"] = pd.to_numeric(df["ASE_status"], errors="coerce").fillna(0).astype(int)
    df["match_bin"] = df["match_bin"].astype(str)
    return df


def empirical_test(df: pd.DataFrame, target_col: str, n_iter: int, rng: np.random.Generator) -> tuple[int, float, float, int]:
    observed = int(((df["ASE_status"] == 1) & (df[target_col] == 1)).sum())
    null = np.zeros(n_iter, dtype=float)
    n_strata = 0
    for (_, _), group in df.groupby(["sample", "match_bin"], sort=False, dropna=False):
        n_total = int(group.shape[0])
        n_ase = int((group["ASE_status"] == 1).sum())
        n_target = int((group[target_col] == 1).sum())
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
    outdir.mkdir(parents=True, exist_ok=True)

    meta = pd.read_csv(args.sample_meta, sep="\t", compression="infer", dtype=str)
    required = {"Sample", "breed", "tissue", "clean_binned_file"}
    missing = sorted(required - set(meta.columns))
    if missing:
        raise ValueError(f"Sample metadata missing columns: {', '.join(missing)}")

    result_rows: list[dict[str, object]] = []
    sample_qc_rows: list[dict[str, object]] = []

    group_counter = 0
    for (breed, tissue), group_meta in meta.groupby(["breed", "tissue"], sort=True):
        sample_frames = []
        for row in group_meta.itertuples(index=False):
            path = Path(str(row.clean_binned_file))
            if not path.exists():
                sample_qc_rows.append({"breed": breed, "tissue": tissue, "Sample": row.Sample, "status": "missing_clean_file"})
                continue
            df = load_clean_table(path)
            sample_frames.append(df)
            sample_qc_rows.append({"breed": breed, "tissue": tissue, "Sample": row.Sample, "status": "OK", "n_sites": int(df.shape[0]), "n_ASE": int(df["ASE_status"].sum())})

        if not sample_frames:
            continue

        group_df = pd.concat(sample_frames, ignore_index=True)
        target_cols = sorted([c for c in group_df.columns if c.startswith("in_candidate_") or c.startswith("in_fst_decile_")])
        rng = np.random.default_rng(args.seed + group_counter)
        group_counter += 1

        for col in target_cols:
            parsed = parse_target_column(col)
            observed, mean_null, p_emp, n_strata = empirical_test(group_df, col, args.n_iter, rng)
            result_rows.append(
                {
                    "breed": breed,
                    "tissue": tissue,
                    "analysis_type": parsed["analysis_type"],
                    "comparison": parsed["comparison"],
                    "window_kb": parsed["window_kb"],
                    "step_kb": parsed["step_kb"],
                    "threshold_label": parsed["threshold_label"],
                    "fst_decile": parsed["fst_decile"],
                    "target_column": col,
                    "n_samples": int(group_df["sample"].nunique()),
                    "n_detectable_sites": int(group_df.shape[0]),
                    "n_ASE_sites": int((group_df["ASE_status"] == 1).sum()),
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
        results["BH_FDR"] = np.nan
        for comparison, idx in results.groupby("comparison", sort=False).groups.items():
            results.loc[idx, "BH_FDR"] = bh_adjust(results.loc[idx, "empirical_p"])
        results = results.sort_values(
            ["analysis_type", "comparison", "window_kb", "step_kb", "threshold_label", "fst_decile", "breed", "tissue"],
            na_position="last",
        ).reset_index(drop=True)

    pd.DataFrame(sample_qc_rows).to_csv(outdir / "indase_selection_sample_qc.tsv", sep="\t", index=False)
    results.to_csv(outdir / "indase_selection_enrichment.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
