#!/usr/bin/env python3

import argparse
import os
import re

import numpy as np
import pandas as pd


def normalize_chr(x):
    return re.sub(r"^chr", "", str(x).strip(), flags=re.IGNORECASE)


def parse_args():
    parser = argparse.ArgumentParser(description="Summarize normalized XP-EHH values into sliding windows.")
    parser.add_argument("--job-table", help="TSV with comparison, input_dir, chrom_lengths, window_bp, step_bp, output.")
    parser.add_argument("--comparison")
    parser.add_argument("--input-dir", help="Directory containing chromosome-specific *.norm files.")
    parser.add_argument("--chrom-lengths", help="TSV/TXT with chromosome and length.")
    parser.add_argument("--window-bp", type=int)
    parser.add_argument("--step-bp", type=int)
    parser.add_argument("--output")
    return parser.parse_args()


def read_chrom_lengths(path):
    out = {}
    with open(path) as handle:
        for line in handle:
            if not line.strip():
                continue
            arr = line.strip().split()
            if len(arr) < 2:
                continue
            out[normalize_chr(arr[0])] = int(float(arr[1]))
    if not out:
        raise ValueError("No chromosome lengths were loaded.")
    return out


def find_norm_file(input_dir, chrom, comparison):
    candidates = [
        os.path.join(input_dir, f"Chr{chrom}.{comparison}.norm"),
        os.path.join(input_dir, f"chr{chrom}.{comparison}.norm"),
        os.path.join(input_dir, f"{chrom}.{comparison}.norm"),
    ]
    for path in candidates:
        if os.path.exists(path):
            return path
    return None


def read_norm_file(path, chrom):
    df = pd.read_csv(path, sep=r"\s+", engine="python")
    df.columns = [str(c).strip().lower() for c in df.columns]
    if "pos" not in df.columns or "normxpehh" not in df.columns:
        raise ValueError(f"Missing pos or normxpehh in {path}")
    if "chr" in df.columns:
        target = normalize_chr(chrom)
        df = df[df["chr"].map(lambda x: normalize_chr(x) == target)].copy()
    out = df[["pos", "normxpehh"]].copy()
    out["pos"] = pd.to_numeric(out["pos"], errors="coerce")
    out["normxpehh"] = pd.to_numeric(out["normxpehh"], errors="coerce")
    out = out.dropna(subset=["pos", "normxpehh"])
    out = out.sort_values("pos").drop_duplicates().reset_index(drop=True)
    return out


def sliding_window(df, chrom, chrom_len, window_bp, step_bp):
    if df.empty:
        return pd.DataFrame(columns=["chrom", "start", "end", "count", "mean_normxpehh"])
    pos = df["pos"].to_numpy(dtype=np.int64)
    val = df["normxpehh"].to_numpy(dtype=float)
    starts = np.arange(1, chrom_len + 1, step_bp, dtype=np.int64)
    rows = []
    for start in starts:
        end = min(start + window_bp - 1, chrom_len)
        left = np.searchsorted(pos, start, side="left")
        right = np.searchsorted(pos, end, side="right")
        if right > left:
            x = val[left:right]
            rows.append([chrom, int(start), int(end), int(len(x)), float(np.nanmean(x))])
    return pd.DataFrame(rows, columns=["chrom", "start", "end", "count", "mean_normxpehh"])


def chrom_sort_key(chrom):
    c = normalize_chr(chrom)
    return (0, int(c)) if c.isdigit() else (1, c)


def run_one(comparison, input_dir, chrom_lengths_file, window_bp, step_bp, output):
    chrom_lengths = read_chrom_lengths(chrom_lengths_file)
    all_windows = []
    for chrom in sorted(chrom_lengths, key=chrom_sort_key):
        path = find_norm_file(input_dir, chrom, comparison)
        if path is None:
            continue
        site_df = read_norm_file(path, chrom)
        win_df = sliding_window(site_df, chrom, chrom_lengths[chrom], window_bp, step_bp)
        if not win_df.empty:
            all_windows.append(win_df)
    if all_windows:
        out = pd.concat(all_windows, ignore_index=True)
    else:
        out = pd.DataFrame(columns=["chrom", "start", "end", "count", "mean_normxpehh"])
    os.makedirs(os.path.dirname(output) or ".", exist_ok=True)
    out.to_csv(output, sep="\t", index=False)


def main():
    args = parse_args()
    if args.job_table:
        jobs = pd.read_csv(args.job_table, sep="\t")
        required = ["comparison", "input_dir", "chrom_lengths", "window_bp", "step_bp", "output"]
        missing = [c for c in required if c not in jobs.columns]
        if missing:
            raise ValueError(f"Missing columns in job table: {', '.join(missing)}")
        for row in jobs.itertuples(index=False):
            run_one(
                comparison=row.comparison,
                input_dir=row.input_dir,
                chrom_lengths_file=row.chrom_lengths,
                window_bp=int(row.window_bp),
                step_bp=int(row.step_bp),
                output=row.output,
            )
        return

    required_args = {
        "comparison": args.comparison,
        "input_dir": args.input_dir,
        "chrom_lengths": args.chrom_lengths,
        "window_bp": args.window_bp,
        "step_bp": args.step_bp,
        "output": args.output,
    }
    missing = [k for k, v in required_args.items() if v is None]
    if missing:
        raise ValueError("Missing required arguments without --job-table: " + ", ".join(missing))
    run_one(
        comparison=args.comparison,
        input_dir=args.input_dir,
        chrom_lengths_file=args.chrom_lengths,
        window_bp=args.window_bp,
        step_bp=args.step_bp,
        output=args.output,
    )


if __name__ == "__main__":
    main()
