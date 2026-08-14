#!/usr/bin/env python3

import argparse
import os
import sys

import numpy as np
import pandas as pd


REQUIRED_PI_COLUMNS = ["CHROM", "BIN_START", "BIN_END", "N_VARIANTS", "PI"]


def parse_args():
    parser = argparse.ArgumentParser(description="Calculate pi_ratio from vcftools windowed.pi outputs.")
    parser.add_argument("--pairs", required=True, help="TSV with comparison, numerator_pop, denominator_pop, numerator_file, denominator_file, output_file.")
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def read_pi_file(path, pop):
    df = pd.read_csv(path, sep="\t", dtype={"CHROM": str})
    missing = [c for c in REQUIRED_PI_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns in {path}: {', '.join(missing)}")

    df = df[REQUIRED_PI_COLUMNS].copy()
    df["CHROM"] = df["CHROM"].astype(str).str.strip()
    df["BIN_START"] = pd.to_numeric(df["BIN_START"], errors="coerce")
    df["BIN_END"] = pd.to_numeric(df["BIN_END"], errors="coerce")
    df["N_VARIANTS"] = pd.to_numeric(df["N_VARIANTS"], errors="coerce")
    df["PI"] = pd.to_numeric(df["PI"], errors="coerce")
    df = df.dropna(subset=["CHROM", "BIN_START", "BIN_END"])
    df["BIN_START"] = df["BIN_START"].astype("int64")
    df["BIN_END"] = df["BIN_END"].astype("int64")
    df = df.drop_duplicates(subset=["CHROM", "BIN_START", "BIN_END"], keep="first")
    return df.rename(columns={"N_VARIANTS": f"N_VARIANTS_{pop}", "PI": f"PI_{pop}"})


def calc_one(numerator_pop, denominator_pop, numerator_file, denominator_file, out_file, force):
    if os.path.exists(out_file) and not force:
        return

    num_df = read_pi_file(numerator_file, numerator_pop)
    den_df = read_pi_file(denominator_file, denominator_pop)

    merged = pd.merge(
        num_df,
        den_df,
        on=["CHROM", "BIN_START", "BIN_END"],
        how="inner",
        validate="one_to_one",
    )

    num_col = f"PI_{numerator_pop}"
    den_col = f"PI_{denominator_pop}"
    num_n = f"N_VARIANTS_{numerator_pop}"
    den_n = f"N_VARIANTS_{denominator_pop}"

    merged["pi_ratio"] = np.nan
    ok = merged[num_col].notna() & merged[den_col].notna() & (merged[den_col] != 0)
    merged.loc[ok, "pi_ratio"] = merged.loc[ok, num_col] / merged.loc[ok, den_col]
    merged["log2_pi_ratio"] = np.nan
    ok2 = merged["pi_ratio"].notna() & (merged["pi_ratio"] > 0)
    merged.loc[ok2, "log2_pi_ratio"] = np.log2(merged.loc[ok2, "pi_ratio"])
    merged["abs_log2_pi_ratio"] = merged["log2_pi_ratio"].abs()

    out_cols = [
        "CHROM",
        "BIN_START",
        "BIN_END",
        num_n,
        num_col,
        den_n,
        den_col,
        "pi_ratio",
        "log2_pi_ratio",
        "abs_log2_pi_ratio",
    ]
    merged[out_cols].to_csv(out_file, sep="\t", index=False)


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    pairs = pd.read_csv(args.pairs, sep="\t")
    required = ["comparison", "numerator_pop", "denominator_pop", "numerator_file", "denominator_file", "output_file"]
    missing = [c for c in required if c not in pairs.columns]
    if missing:
        raise ValueError(f"Missing columns in pairs table: {', '.join(missing)}")

    for row in pairs.itertuples(index=False):
        out_file = row.output_file
        if not os.path.isabs(out_file):
            out_file = os.path.join(args.out_dir, out_file)
        os.makedirs(os.path.dirname(out_file) or args.out_dir, exist_ok=True)
        calc_one(
            numerator_pop=row.numerator_pop,
            denominator_pop=row.denominator_pop,
            numerator_file=row.numerator_file,
            denominator_file=row.denominator_file,
            out_file=out_file,
            force=args.force,
        )


if __name__ == "__main__":
    main()
