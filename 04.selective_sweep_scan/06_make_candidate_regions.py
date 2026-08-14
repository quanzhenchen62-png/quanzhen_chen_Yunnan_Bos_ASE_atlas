#!/usr/bin/env python3

import argparse
import os
from collections import Counter

import numpy as np
import pandas as pd


TOP_INFO = [
    ("top1", 1.0, 0.99),
    ("top2.5", 2.5, 0.975),
    ("top5", 5.0, 0.95),
    ("top10", 10.0, 0.90),
]


def parse_args():
    parser = argparse.ArgumentParser(description="Create candidate selective-sweep regions from FST, pi_ratio, and XP-EHH windows.")
    parser.add_argument("--match", required=True, help="TSV with group, win, step, Fst, XPEHH, Pi_ratio.")
    parser.add_argument("--outdir", required=True)
    return parser.parse_args()


def mode_value(values):
    values = [v for v in values if pd.notna(v)]
    if not values:
        return None
    return Counter(values).most_common(1)[0][0]


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


def standardize_to_bed(df, chr_col, start_col, end_col, expected_bp):
    df = df.copy()
    df["CHROM"] = df[chr_col].astype(str)
    df["START_RAW"] = pd.to_numeric(df[start_col], errors="coerce").astype("Int64")
    df["END_RAW"] = pd.to_numeric(df[end_col], errors="coerce").astype("Int64")
    df = df.dropna(subset=["START_RAW", "END_RAW"]).copy()
    df["START_RAW"] = df["START_RAW"].astype(int)
    df["END_RAW"] = df["END_RAW"].astype(int)
    coord_type = detect_coord_type(df["START_RAW"], df["END_RAW"], expected_bp)
    if coord_type == "1-based_inclusive":
        df["BED_START"] = df["START_RAW"] - 1
        df["BED_END"] = df["END_RAW"]
    else:
        df["BED_START"] = df["START_RAW"]
        df["BED_END"] = df["END_RAW"]
    df["coord_type"] = coord_type
    df["window_bed"] = (
        df["CHROM"].astype(str) + "_" + df["BED_START"].astype(str) + "_" + df["BED_END"].astype(str)
    )
    return df


def read_fst(path, expected_bp):
    df = pd.read_csv(path, sep="\t", dtype={"CHROM": str})
    required = ["CHROM", "BIN_START", "BIN_END", "WEIGHTED_FST"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns in FST file {path}: {', '.join(missing)}")
    df = standardize_to_bed(df, "CHROM", "BIN_START", "BIN_END", expected_bp)
    out = df[["window_bed", "CHROM", "BED_START", "BED_END", "START_RAW", "END_RAW", "coord_type", "WEIGHTED_FST"]].copy()
    out = out.rename(columns={"START_RAW": "FST_START_RAW", "END_RAW": "FST_END_RAW", "coord_type": "FST_coord_type"})
    if "N_VARIANTS" in df.columns:
        out["FST_N_VARIANTS"] = pd.to_numeric(df["N_VARIANTS"], errors="coerce")
    else:
        out["FST_N_VARIANTS"] = np.nan
    if "MEAN_FST" in df.columns:
        out["MEAN_FST"] = pd.to_numeric(df["MEAN_FST"], errors="coerce")
    else:
        out["MEAN_FST"] = np.nan
    out["WEIGHTED_FST"] = pd.to_numeric(out["WEIGHTED_FST"], errors="coerce")
    return out.drop_duplicates(subset=["window_bed"], keep="first")


def read_pi(path, expected_bp):
    df = pd.read_csv(path, sep="\t", dtype={"CHROM": str})
    required = ["CHROM", "BIN_START", "BIN_END", "pi_ratio"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns in pi_ratio file {path}: {', '.join(missing)}")
    df = standardize_to_bed(df, "CHROM", "BIN_START", "BIN_END", expected_bp)
    out = df[["window_bed", "START_RAW", "END_RAW", "coord_type", "pi_ratio"]].copy()
    out = out.rename(columns={"START_RAW": "PI_START_RAW", "END_RAW": "PI_END_RAW", "coord_type": "PI_coord_type"})
    out["pi_ratio"] = pd.to_numeric(out["pi_ratio"], errors="coerce")
    out.loc[~np.isfinite(out["pi_ratio"]), "pi_ratio"] = np.nan
    return out.drop_duplicates(subset=["window_bed"], keep="first")


def read_xpehh(path, expected_bp):
    df = pd.read_csv(path, sep="\t", dtype={"chrom": str})
    required = ["chrom", "start", "end", "mean_normxpehh"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns in XPEHH file {path}: {', '.join(missing)}")
    df = standardize_to_bed(df, "chrom", "start", "end", expected_bp)
    out = df[["window_bed", "START_RAW", "END_RAW", "coord_type", "mean_normxpehh"]].copy()
    out = out.rename(columns={"START_RAW": "XPEHH_START_RAW", "END_RAW": "XPEHH_END_RAW", "coord_type": "XPEHH_coord_type"})
    if "count" in df.columns:
        out["XPEHH_COUNT"] = pd.to_numeric(df["count"], errors="coerce")
    else:
        out["XPEHH_COUNT"] = np.nan
    out["mean_normxpehh"] = pd.to_numeric(out["mean_normxpehh"], errors="coerce")
    out["abs_mean_normxpehh"] = out["mean_normxpehh"].abs()
    return out.drop_duplicates(subset=["window_bed"], keep="first")


def quantile_top(series, prob):
    x = series.replace([np.inf, -np.inf], np.nan).dropna()
    if x.empty:
        return np.nan
    return float(np.quantile(x, prob))


def make_sig_methods(df):
    methods = []
    for fst_sig, pi_sig, xpehh_sig in zip(df["fst_significant"], df["pi_significant"], df["xpehh_significant"]):
        parts = []
        if bool(fst_sig):
            parts.append("FST")
        if bool(pi_sig):
            parts.append("Pi_ratio")
        if bool(xpehh_sig):
            parts.append("abs_XPEHH")
        methods.append(";".join(parts))
    return methods


def collapse_methods(values):
    parts = []
    for value in values:
        if pd.isna(value) or value == "":
            continue
        parts.extend(str(value).split(";"))
    return ";".join(sorted(set(parts)))


def chrom_sort_key(chrom):
    c = str(chrom)
    c2 = c[3:] if c.lower().startswith("chr") else c
    return (0, int(c2)) if c2.isdigit() else (1, c2)


def sort_frame(df):
    if df.empty:
        return df
    tmp = df.copy()
    tmp["_chrom_key"] = tmp["CHROM"].map(chrom_sort_key)
    return tmp.sort_values(by=["_chrom_key", "BED_START", "BED_END"]).drop(columns="_chrom_key")


def merge_regions(unique_dt):
    merged_rows = []
    for threshold_label, x in unique_dt.groupby("threshold_label", sort=False):
        x = sort_frame(x)
        current_chrom = None
        current_start = None
        current_end = None
        bucket = []
        for row in x.itertuples(index=False):
            if current_chrom is None:
                current_chrom = row.CHROM
                current_start = int(row.BED_START)
                current_end = int(row.BED_END)
                bucket = [row]
                continue
            if row.CHROM == current_chrom and int(row.BED_START) <= current_end:
                current_end = max(current_end, int(row.BED_END))
                bucket.append(row)
            else:
                merged_rows.append(summarize_bucket(bucket, current_chrom, current_start, current_end, threshold_label))
                current_chrom = row.CHROM
                current_start = int(row.BED_START)
                current_end = int(row.BED_END)
                bucket = [row]
        if bucket:
            merged_rows.append(summarize_bucket(bucket, current_chrom, current_start, current_end, threshold_label))
    if not merged_rows:
        return pd.DataFrame()
    out = pd.DataFrame(merged_rows)
    out = sort_frame(out)
    return out


def safe_mean(series):
    s = pd.to_numeric(pd.Series(series), errors="coerce")
    return np.nan if s.dropna().empty else float(s.mean())


def safe_max(series):
    s = pd.to_numeric(pd.Series(series), errors="coerce")
    return np.nan if s.dropna().empty else float(s.max())


def summarize_bucket(bucket, chrom, start, end, threshold_label):
    df = pd.DataFrame([r._asdict() for r in bucket])
    row = {
        "group": df["group"].iloc[0],
        "win": int(df["win"].iloc[0]),
        "step": int(df["step"].iloc[0]),
        "threshold_label": threshold_label,
        "top_percent": float(df["top_percent"].iloc[0]),
        "CHROM": chrom,
        "BED_START": int(start),
        "BED_END": int(end),
        "BIN_START_1based": int(start) + 1,
        "BIN_END_1based": int(end),
        "merged_region_id": f"{chrom}:{start}-{end}",
        "n_windows": int(df.shape[0]),
        "mean_WEIGHTED_FST": safe_mean(df["WEIGHTED_FST"]),
        "max_WEIGHTED_FST": safe_max(df["WEIGHTED_FST"]),
        "mean_pi_ratio": safe_mean(df["pi_ratio"]),
        "max_pi_ratio": safe_max(df["pi_ratio"]),
        "mean_mean_normxpehh": safe_mean(df["mean_normxpehh"]),
        "mean_abs_mean_normxpehh": safe_mean(df["abs_mean_normxpehh"]),
        "max_abs_mean_normxpehh": safe_max(df["abs_mean_normxpehh"]),
        "any_fst_significant": bool(df["fst_significant"].fillna(False).any()),
        "any_pi_significant": bool(df["pi_significant"].fillna(False).any()),
        "any_xpehh_significant": bool(df["xpehh_significant"].fillna(False).any()),
        "max_n_sig": int(pd.to_numeric(df["n_sig"], errors="coerce").max()),
        "sig_methods_merged": collapse_methods(df["sig_methods"]),
        "windows_bed": ",".join(df["window_bed"].astype(str)),
    }
    return row


def write_bed(df, path, id_col):
    if df.empty:
        pd.DataFrame(columns=["CHROM", "BED_START", "BED_END", "ID"]).to_csv(path, sep="\t", header=False, index=False)
        return
    bed = df[["CHROM", "BED_START", "BED_END", id_col]].copy()
    bed.to_csv(path, sep="\t", header=False, index=False)


def process_one(row, outdir):
    group = row["group"]
    win = int(row["win"])
    step = int(row["step"])
    expected_bp = win * 1000
    prefix = str(group)

    fst = read_fst(row["Fst"], expected_bp)
    pi = read_pi(row["Pi_ratio"], expected_bp)
    xpehh = read_xpehh(row["XPEHH"], expected_bp)

    merged = fst.merge(pi, on="window_bed", how="inner").merge(xpehh, on="window_bed", how="inner")
    if merged.empty:
        return None

    merged["group"] = group
    merged["win"] = win
    merged["step"] = step
    merged["BIN_START_1based"] = merged["BED_START"] + 1
    merged["BIN_END_1based"] = merged["BED_END"]

    all_list = []
    raw_list = []
    unique_list = []

    for label, pct, qcut in TOP_INFO:
        tmp = merged.copy()
        fst_thr = quantile_top(tmp["WEIGHTED_FST"], qcut)
        pi_thr = quantile_top(tmp["pi_ratio"], qcut)
        xpehh_thr = quantile_top(tmp["abs_mean_normxpehh"], qcut)

        tmp["threshold_label"] = label
        tmp["top_percent"] = pct
        tmp["fst_threshold"] = fst_thr
        tmp["pi_ratio_threshold"] = pi_thr
        tmp["abs_xpehh_threshold"] = xpehh_thr
        tmp["fst_significant"] = tmp["WEIGHTED_FST"].ge(fst_thr) & tmp["WEIGHTED_FST"].notna() if pd.notna(fst_thr) else False
        tmp["pi_significant"] = tmp["pi_ratio"].ge(pi_thr) & tmp["pi_ratio"].notna() if pd.notna(pi_thr) else False
        tmp["xpehh_significant"] = tmp["abs_mean_normxpehh"].ge(xpehh_thr) & tmp["abs_mean_normxpehh"].notna() if pd.notna(xpehh_thr) else False
        tmp["n_sig"] = tmp[["fst_significant", "pi_significant", "xpehh_significant"]].astype(int).sum(axis=1)
        tmp["candidate_signal"] = tmp["n_sig"] >= 2
        tmp["sig_methods"] = make_sig_methods(tmp)
        all_list.append(tmp)

        cand = tmp[tmp["candidate_signal"]].copy()
        raw_list.append(cand)
        unique_list.append(cand.drop_duplicates(subset=["group", "win", "step", "threshold_label", "CHROM", "BED_START", "BED_END"]))

    all_dt = pd.concat(all_list, ignore_index=True) if all_list else pd.DataFrame()
    raw_dt = pd.concat(raw_list, ignore_index=True) if raw_list else pd.DataFrame()
    unique_dt = pd.concat(unique_list, ignore_index=True) if unique_list else pd.DataFrame()

    if not raw_dt.empty:
        raw_dt["raw_id"] = raw_dt["group"].astype(str) + "|" + raw_dt["threshold_label"].astype(str) + "|" + raw_dt["window_bed"].astype(str) + "|" + raw_dt["sig_methods"].astype(str)
    if not unique_dt.empty:
        unique_dt = sort_frame(unique_dt)
        unique_dt["unique_id"] = unique_dt["group"].astype(str) + "|" + unique_dt["threshold_label"].astype(str) + "|" + unique_dt["window_bed"].astype(str) + "|" + unique_dt["sig_methods"].astype(str)

    merged_dt = merge_regions(unique_dt)
    if not merged_dt.empty:
        merged_dt["merged_id"] = (
            merged_dt["group"].astype(str)
            + "|"
            + merged_dt["threshold_label"].astype(str)
            + "|"
            + merged_dt["merged_region_id"].astype(str)
            + "|nwin="
            + merged_dt["n_windows"].astype(str)
            + "|"
            + merged_dt["sig_methods_merged"].astype(str)
        )

    os.makedirs(os.path.join(outdir, "01.all_windows"), exist_ok=True)
    os.makedirs(os.path.join(outdir, "02.candidate_raw"), exist_ok=True)
    os.makedirs(os.path.join(outdir, "03.candidate_unique"), exist_ok=True)
    os.makedirs(os.path.join(outdir, "04.candidate_merged"), exist_ok=True)
    os.makedirs(os.path.join(outdir, "05.summary"), exist_ok=True)

    all_dt.to_csv(os.path.join(outdir, "01.all_windows", f"{prefix}.all_windows.with_thresholds.tsv"), sep="\t", index=False)
    raw_dt.to_csv(os.path.join(outdir, "02.candidate_raw", f"{prefix}.candidate_raw.tsv"), sep="\t", index=False)
    unique_dt.to_csv(os.path.join(outdir, "03.candidate_unique", f"{prefix}.candidate_unique.tsv"), sep="\t", index=False)
    merged_dt.to_csv(os.path.join(outdir, "04.candidate_merged", f"{prefix}.candidate_merged.with_metrics.tsv"), sep="\t", index=False)

    write_bed(raw_dt, os.path.join(outdir, "02.candidate_raw", f"{prefix}.candidate_raw.bed"), "raw_id")
    write_bed(unique_dt, os.path.join(outdir, "03.candidate_unique", f"{prefix}.candidate_unique.bed"), "unique_id")
    write_bed(merged_dt, os.path.join(outdir, "04.candidate_merged", f"{prefix}.candidate_merged.bed"), "merged_id")

    summary_rows = []
    merged_counts = merged_dt.groupby("threshold_label").size().to_dict() if not merged_dt.empty else {}
    raw_counts = raw_dt.groupby("threshold_label").size().to_dict() if not raw_dt.empty else {}
    unique_counts = unique_dt.groupby("threshold_label").size().to_dict() if not unique_dt.empty else {}
    for label, pct, _ in TOP_INFO:
        summary_rows.append({
            "group": group,
            "win": win,
            "step": step,
            "n_fst_windows": int(fst.shape[0]),
            "n_pi_windows": int(pi.shape[0]),
            "n_xpehh_windows": int(xpehh.shape[0]),
            "n_common_windows": int(merged.shape[0]),
            "threshold_label": label,
            "top_percent": pct,
            "n_candidate_raw_windows": int(raw_counts.get(label, 0)),
            "n_candidate_unique_windows": int(unique_counts.get(label, 0)),
            "n_merged_regions": int(merged_counts.get(label, 0)),
        })
    pd.DataFrame(summary_rows).to_csv(os.path.join(outdir, "05.summary", f"{prefix}.summary.tsv"), sep="\t", index=False)
    return pd.DataFrame(summary_rows)


def main():
    args = parse_args()
    match_dt = pd.read_csv(args.match, sep="\t")
    required = ["group", "win", "step", "Fst", "XPEHH", "Pi_ratio"]
    missing = [c for c in required if c not in match_dt.columns]
    if missing:
        raise ValueError(f"Missing columns in match table: {', '.join(missing)}")

    summaries = []
    for _, row in match_dt.iterrows():
        ans = process_one(row, args.outdir)
        if ans is not None:
            summaries.append(ans)
    if summaries:
        pd.concat(summaries, ignore_index=True).to_csv(
            os.path.join(args.outdir, "05.summary", "all_groups.summary.tsv"),
            sep="\t",
            index=False,
        )


if __name__ == "__main__":
    main()
