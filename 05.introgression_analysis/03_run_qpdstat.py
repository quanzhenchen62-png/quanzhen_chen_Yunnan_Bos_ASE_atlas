#!/usr/bin/env python3

import argparse
import csv
import os
import subprocess


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run ADMIXTOOLS qpDstat from a TSV of W/X/Y/Z test tuples."
    )
    parser.add_argument("--qpdstat-bin", required=True, help="Path to qpDstat binary.")
    parser.add_argument("--geno", required=True, help="EIGENSTRAT .geno file.")
    parser.add_argument("--snp", required=True, help="EIGENSTRAT .snp file.")
    parser.add_argument("--ind", required=True, help="EIGENSTRAT .ind file.")
    parser.add_argument("--tests", required=True, help="TSV with columns W, X, Y, Z.")
    parser.add_argument("--out-prefix", required=True, help="Output prefix.")
    parser.add_argument(
        "--blgsize",
        type=float,
        default=0.05,
        help="Block jackknife size in Morgans. Default: 0.05.",
    )
    parser.add_argument(
        "--f4mode",
        choices=["YES", "NO"],
        default="NO",
        help="Use f4mode YES or NO. Default: NO for D-statistics.",
    )
    return parser.parse_args()


def load_tests(path):
    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"W", "X", "Y", "Z"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError("tests file must contain columns: W, X, Y, Z")
        return list(reader)


def write_popfile(rows, path):
    with open(path, "w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(f"{row['W']} {row['X']} {row['Y']} {row['Z']}\n")


def write_parfile(args, popfile, path):
    lines = [
        f"genotypename: {args.geno}",
        f"snpname: {args.snp}",
        f"indivname: {args.ind}",
        "printsd: YES",
        "printname: YES",
        f"f4mode: {args.f4mode}",
        f"blgsize: {args.blgsize}",
        f"popfilename: {popfile}",
    ]
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")


def parse_result_line(line):
    stripped = line.strip()
    if not stripped or stripped.startswith("#"):
        return None
    tokens = stripped.split()
    if tokens[0].lower().startswith("result"):
        tokens = tokens[1:]
    if len(tokens) < 8:
        return None
    try:
        d_value = float(tokens[4])
        z_value = float(tokens[5])
        baba = float(tokens[6])
        abba = float(tokens[7])
        nsnp = int(float(tokens[8])) if len(tokens) > 8 else None
    except ValueError:
        return None
    return {
        "W": tokens[0],
        "X": tokens[1],
        "Y": tokens[2],
        "Z": tokens[3],
        "D": d_value,
        "Zscore": z_value,
        "BABA": baba,
        "ABBA": abba,
        "NSNP": nsnp,
    }


def parse_qpdstat_output(path):
    records = []
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            record = parse_result_line(line)
            if record is not None:
                records.append(record)
    return records


def main():
    args = parse_args()
    tests = load_tests(args.tests)
    if not tests:
        raise SystemExit("tests file contains no D-statistic rows")
    outdir = os.path.dirname(os.path.abspath(args.out_prefix)) or "."
    os.makedirs(outdir, exist_ok=True)

    popfile = args.out_prefix + ".popfilename.txt"
    parfile = args.out_prefix + ".parfile.txt"
    stdout_path = args.out_prefix + ".qpDstat.stdout.txt"
    stderr_path = args.out_prefix + ".qpDstat.stderr.txt"
    results_path = args.out_prefix + ".results.tsv"

    write_popfile(tests, popfile)
    write_parfile(args, popfile, parfile)

    with open(stdout_path, "w", encoding="utf-8") as stdout_handle, open(
        stderr_path, "w", encoding="utf-8"
    ) as stderr_handle:
        completed = subprocess.run(
            [args.qpdstat_bin, "-p", parfile],
            stdout=stdout_handle,
            stderr=stderr_handle,
            check=False,
            text=True,
        )
    if completed.returncode != 0:
        raise SystemExit(f"qpDstat failed with exit code {completed.returncode}")

    parsed = parse_qpdstat_output(stdout_path)
    by_key = {(row["W"], row["X"], row["Y"], row["Z"]): row for row in parsed}

    fieldnames = list(tests[0].keys()) + ["D", "Zscore", "BABA", "ABBA", "NSNP"]
    with open(results_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in tests:
            merged = dict(row)
            parsed_row = by_key.get((row["W"], row["X"], row["Y"], row["Z"]))
            if parsed_row is not None:
                merged.update(parsed_row)
            else:
                merged.update({"D": "", "Zscore": "", "BABA": "", "ABBA": "", "NSNP": ""})
            writer.writerow(merged)


if __name__ == "__main__":
    main()
