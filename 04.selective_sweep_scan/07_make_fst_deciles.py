#!/usr/bin/env python3

import argparse
import math
import os
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description="Assign FST windows to deciles within each comparison and window scheme.")
    parser.add_argument("--job-table", required=True, help="TSV with comparison, window_bp, step_bp, output_prefix.")
    parser.add_argument("--outdir", required=True)
    return parser.parse_args()


def mode_value(values):
    vals = [v for v in values if pd.notna(v)]
    if not vals:
        return None
    return pd.Series(vals).mode().iloc[0]


def detect_coord_type(start, end, expected_bp):
    inclusive = mode_value((end - start + 1).tolist())
    halfopen = mode_value((end - start).tolist())
    if inclusive == expected_bp:
        return "1-based_inclusive"
    if halfopen == expected_bp:
        return "0-based_half-open_BED"
    raise ValueError(
        f"Cannot determine coordinate type for expected_bp={expected_bp}; "
        f"mode inclusive={inclusive}, mode half-open={halfopen}"
    )


def standardize_to_bed(df, expected_bp):
    df = df.copy()
    df["CHROM"] = df["CHROM"].astype(str)
    df["BIN_START"] = pd.to_numeric(df["BIN_START"], errors="coerce")
    df["BIN_END"] = pd.to_numeric(df["BIN_END"], errors="coerce")
    df = df.dropna(subset=["BIN_START", "BIN_END"]).copy()
    df["BIN_START"] = df["BIN_START"].astype(int)
    df["BIN_END"] = df["BIN_END"].astype(int)
    coord_type = detect_coord_type(df["BIN_START"], df["BIN_END"], expected_bp)
    if coord_type == "1-based_inclusive":
        df["BED_START"] = df["BIN_START"] - 1
        df["BED_END"] = df["BIN_END"]
    else:
        df["BED_START"] = df["BIN_START"]
        df["BED_END"] = df["BIN_END"]
        df["BIN_START"] = df["BED_START"] + 1
        df["BIN_END"] = df["BED_END"]
    df["coord_type"] = coord_type
    return df


def chrom_sort_key(chrom):
    c = str(chrom)
    c2 = c[3:] if c.lower().startswith("chr") else c
    return (0, int(c2)) if c2.isdigit() else (1, c2)


def assign_deciles(df):
    x = df.copy()
    x["WEIGHTED_FST"] = pd.to_numeric(x["WEIGHTED_FST"], errors="coerce")
    x = x.dropna(subset=["WEIGHTED_FST"]).copy()
    x["_chrom_key"] = x["CHROM"].map(chrom_sort_key)
    x = x.sort_values(by=["WEIGHTED_FST", "_chrom_key", "BED_START", "BED_END"]).reset_index(drop=True)
    n = len(x)
    if n == 0:
        return x
    x["rank_order"] = np.arange(1, n + 1)
    x["fst_decile"] = x["rank_order"].map(lambda r: int(math.ceil(r * 10.0 / n)))
    x.loc[x["fst_decile"] < 1, "fst_decile"] = 1
    x.loc[x["fst_decile"] > 10, "fst_decile"] = 10
    x["fst_decile_label"] = "decile" + x["fst_decile"].astype(str)
    x = x.drop(columns="_chrom_key")
    return x


def write_decile_beds(df, prefix, outdir):
    bed_dir = Path(outdir) / "bed_by_decile"
    bed_dir.mkdir(parents=True, exist_ok=True)
    for decile in range(1, 11):
        sub = df[df["fst_decile"] == decile].copy()
        out = bed_dir / f"{prefix}.fst_decile{decile}.bed"
        if sub.empty:
            out.write_text("")
            continue
        sub["window_id"] = (
            prefix
            + "|decile"
            + sub["fst_decile"].astype(str)
            + "|"
            + sub["CHROM"].astype(str)
            + ":"
            + sub["BED_START"].astype(str)
            + "-"
            + sub["BED_END"].astype(str)
        )
        sub[["CHROM", "BED_START", "BED_END", "window_id"]].to_csv(out, sep="\t", index=False, header=False)


def main():
    args = parse_args()
    jobs = pd.read_csv(args.job_table, sep="\t")
    required = ["comparison", "window_bp", "step_bp", "output_prefix"]
    missing = [c for c in required if c not in jobs.columns]
    if missing:
        raise ValueError(f"Missing columns in job table: {', '.join(missing)}")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    table_dir = outdir / "tables"
    table_dir.mkdir(parents=True, exist_ok=True)
    summary_rows = []

    for row in jobs.itertuples(index=False):
        comparison = row.comparison
        window_bp = int(row.window_bp)
        step_bp = int(row.step_bp)
        output_prefix = str(row.output_prefix)
        fst_file = output_prefix + ".windowed.weir.fst"
        if not os.path.exists(fst_file):
            raise FileNotFoundError(f"FST file not found: {fst_file}")

        df = pd.read_csv(fst_file, sep="\t", dtype={"CHROM": str})
        required_fst = ["CHROM", "BIN_START", "BIN_END", "WEIGHTED_FST"]
        missing_fst = [c for c in required_fst if c not in df.columns]
        if missing_fst:
            raise ValueError(f"Missing columns in {fst_file}: {', '.join(missing_fst)}")

        df = standardize_to_bed(df, expected_bp=window_bp)
        df["comparison"] = comparison
        df["window_kb"] = window_bp // 1000
        df["step_kb"] = step_bp // 1000
        df = assign_deciles(df)

        prefix = f"{comparison}_{window_bp // 1000}_{step_bp // 1000}"
        out_table = table_dir / f"{prefix}.fst_deciles.tsv"
        df.to_csv(out_table, sep="\t", index=False)
        write_decile_beds(df, prefix, outdir)

        for decile, n in df["fst_decile"].value_counts().sort_index().items():
            summary_rows.append({
                "comparison": comparison,
                "window_kb": window_bp // 1000,
                "step_kb": step_bp // 1000,
                "fst_decile": int(decile),
                "n_windows": int(n),
                "source_fst_file": fst_file,
                "output_table": str(out_table),
            })

    pd.DataFrame(summary_rows).to_csv(outdir / "fst_decile_summary.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
