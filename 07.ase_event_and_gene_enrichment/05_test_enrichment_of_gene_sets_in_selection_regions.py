#!/usr/bin/env python3

from __future__ import annotations

import argparse
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run Fisher-exact enrichment for tissue-specific and broadly expressed genes.")
    parser.add_argument("--gene-table", required=True, help="TSV produced by 04_prepare_gene_selection_table.py.")
    parser.add_argument("--output", required=True, help="Output TSV.")
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


def main() -> None:
    args = parse_args()
    df = pd.read_csv(args.gene_table, sep="\t", compression="infer")
    required = {"gene", "is_tissue_specific", "is_tissue_broadly"}
    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError(f"Gene table missing columns: {', '.join(missing)}")

    target_cols = sorted([c for c in df.columns if c.startswith("in_candidate_") or c.startswith("in_fst_decile_")])
    result_rows: list[dict[str, object]] = []

    for gene_set, flag_col in [("tissue_specific", "is_tissue_specific"), ("tissue_broadly", "is_tissue_broadly")]:
        target_flag = pd.to_numeric(df[flag_col], errors="coerce").fillna(0).astype(int)
        for col in target_cols:
            select_flag = pd.to_numeric(df[col], errors="coerce").fillna(0).astype(int)
            a = int(((target_flag == 1) & (select_flag == 1)).sum())
            b = int(((target_flag == 1) & (select_flag == 0)).sum())
            c = int(((target_flag == 0) & (select_flag == 1)).sum())
            d = int(((target_flag == 0) & (select_flag == 0)).sum())
            odds_ratio, p_value = fisher_exact([[a, b], [c, d]], alternative="two-sided")
            expected = ((a + b) * (a + c) / max(a + b + c + d, 1))
            parsed = parse_target_column(col)
            result_rows.append(
                {
                    "gene_set": gene_set,
                    "analysis_type": parsed["analysis_type"],
                    "comparison": parsed["comparison"],
                    "window_kb": parsed["window_kb"],
                    "step_kb": parsed["step_kb"],
                    "threshold_label": parsed["threshold_label"],
                    "fst_decile": parsed["fst_decile"],
                    "target_column": col,
                    "a_target_and_selected": a,
                    "b_target_only": b,
                    "c_selected_only": c,
                    "d_background_only": d,
                    "observed_overlap": a,
                    "expected_overlap": expected,
                    "fold_enrichment": a / expected if expected > 0 else np.nan,
                    "odds_ratio": odds_ratio,
                    "fisher_p": p_value,
                }
            )

    out = pd.DataFrame(result_rows)
    if not out.empty:
        out["BH_FDR"] = np.nan
        for (comparison, gene_set), idx in out.groupby(["comparison", "gene_set"], sort=False).groups.items():
            out.loc[idx, "BH_FDR"] = bh_adjust(out.loc[idx, "fisher_p"])
        out = out.sort_values(
            ["gene_set", "analysis_type", "comparison", "window_kb", "step_kb", "threshold_label", "fst_decile"],
            na_position="last",
        ).reset_index(drop=True)
    out.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
