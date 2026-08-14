#!/usr/bin/env python3

import argparse
import glob
import os

import pandas as pd


THRESHOLD_ORDER = {"top1": 1, "top2.5": 2, "top5": 3, "top10": 4}


def parse_args():
    parser = argparse.ArgumentParser(description="Combine candidate merged-region tables into Supplementary Table S6.")
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def threshold_definition(label):
    mapping = {
        "top1": "top 1% windows within the current comparison and window scheme",
        "top2.5": "top 2.5% windows within the current comparison and window scheme, including top1",
        "top5": "top 5% windows within the current comparison and window scheme, including top1 and top2.5",
        "top10": "top 10% windows within the current comparison and window scheme, including top1, top2.5, and top5",
    }
    return mapping.get(label, label)


def normalize_comparison_label(label):
    x = str(label)
    x = x.replace(".vs.taurus", ".vs.taurine")
    return x


def chrom_sort_key(chrom):
    c = str(chrom)
    c2 = c[3:] if c.lower().startswith("chr") else c
    return (0, int(c2)) if c2.isdigit() else (1, c2)


def main():
    args = parse_args()
    files = sorted(glob.glob(os.path.join(args.input_dir, "*.candidate_merged.with_metrics.tsv")))
    if not files:
        raise FileNotFoundError(f"No candidate_merged.with_metrics.tsv files found in {args.input_dir}")

    frames = []
    for path in files:
        df = pd.read_csv(path, sep="\t", dtype={"CHROM": str})
        if df.empty:
            continue
        frames.append(df)

    if not frames:
        out = pd.DataFrame(columns=[
            "comparison", "window_kb", "step_kb", "cumulative_threshold_label", "cumulative_threshold_definition",
            "chromosome", "start_0based", "end_0based", "start_1based", "end_1based",
            "merged_region_id", "n_windows", "significant_statistics", "max_n_significant_statistics",
            "mean_weighted_fst", "max_weighted_fst", "mean_pi_ratio", "max_pi_ratio",
            "mean_xpehh", "mean_abs_xpehh", "max_abs_xpehh"
        ])
    else:
        out = pd.concat(frames, ignore_index=True)
        out["comparison"] = out["group"].astype(str).str.rsplit("_", n=2).str[0].map(normalize_comparison_label)
        out["window_kb"] = pd.to_numeric(out["win"], errors="coerce")
        out["step_kb"] = pd.to_numeric(out["step"], errors="coerce")
        out["cumulative_threshold_label"] = out["threshold_label"]
        out["cumulative_threshold_definition"] = out["cumulative_threshold_label"].map(threshold_definition)
        out["chromosome"] = out["CHROM"].astype(str)
        out["start_0based"] = pd.to_numeric(out["BED_START"], errors="coerce").astype("Int64")
        out["end_0based"] = pd.to_numeric(out["BED_END"], errors="coerce").astype("Int64")
        out["start_1based"] = pd.to_numeric(out["BIN_START_1based"], errors="coerce").astype("Int64")
        out["end_1based"] = pd.to_numeric(out["BIN_END_1based"], errors="coerce").astype("Int64")
        out["significant_statistics"] = out["sig_methods_merged"]
        out["max_n_significant_statistics"] = pd.to_numeric(out["max_n_sig"], errors="coerce").astype("Int64")
        out = out.rename(columns={
            "mean_WEIGHTED_FST": "mean_weighted_fst",
            "max_WEIGHTED_FST": "max_weighted_fst",
            "mean_mean_normxpehh": "mean_xpehh",
            "mean_abs_mean_normxpehh": "mean_abs_xpehh",
            "max_abs_mean_normxpehh": "max_abs_xpehh",
        })
        out["_chrom_key"] = out["chromosome"].map(chrom_sort_key)
        out["_threshold_order"] = out["cumulative_threshold_label"].map(THRESHOLD_ORDER)
        out = out.sort_values(
            by=["comparison", "window_kb", "step_kb", "_threshold_order", "_chrom_key", "start_0based", "end_0based"]
        )
        out = out[[
            "comparison",
            "window_kb",
            "step_kb",
            "cumulative_threshold_label",
            "cumulative_threshold_definition",
            "chromosome",
            "start_0based",
            "end_0based",
            "start_1based",
            "end_1based",
            "merged_region_id",
            "n_windows",
            "significant_statistics",
            "max_n_significant_statistics",
            "mean_weighted_fst",
            "max_weighted_fst",
            "mean_pi_ratio",
            "max_pi_ratio",
            "mean_xpehh",
            "mean_abs_xpehh",
            "max_abs_xpehh",
        ]]
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    out.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
