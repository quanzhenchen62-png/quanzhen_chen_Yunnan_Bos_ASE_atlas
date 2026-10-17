#!/usr/bin/env python3
"""Concatenate per-chromosome haplotype segment tables, keep high-confidence
segments, and write the genome-wide table with an ancestry code manifest."""
import argparse
import glob
import gzip
import os
import re


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--segments-dir", required=True,
                        help="Directory containing chrN.haplotype_segments.tsv.gz files")
    parser.add_argument("--output", required=True,
                        help="Genome-wide high-confidence segment table (TSV)")
    parser.add_argument("--manifest-output", default=None,
                        help="Ancestry code manifest (TSV)")
    args = parser.parse_args()

    pattern = re.compile(r"^chr(\d+)\.haplotype_segments\.tsv\.gz$")
    files = []
    for path in glob.glob(os.path.join(args.segments_dir, "*.haplotype_segments.tsv.gz")):
        match = pattern.match(os.path.basename(path))
        if match:
            files.append((int(match.group(1)), path))
    files.sort()
    if not files:
        raise SystemExit("no chrN.haplotype_segments.tsv.gz files found in " + args.segments_dir)

    header = None
    n_written = 0
    codes = {}
    with open(args.output, "w") as out:
        for _, path in files:
            with gzip.open(path, "rt") as fh:
                file_header = fh.readline().rstrip("\n").split("\t")
                if header is None:
                    header = file_header
                    out.write("\t".join(header) + "\n")
                elif file_header != header:
                    raise SystemExit("header mismatch in " + path)
                hc_col = header.index("high_confidence_local_ancestry_segment")
                code_col = header.index("ancestry_numeric_code")
                label_col = header.index("ancestry_label")
                for line in fh:
                    fields = line.rstrip("\n").split("\t")
                    if fields[hc_col] == "True":
                        out.write(line if line.endswith("\n") else line + "\n")
                        n_written += 1
                        codes[fields[code_col]] = fields[label_col]

    if args.manifest_output:
        with open(args.manifest_output, "w") as mf:
            mf.write("numeric_code\tancestry_label\n")
            for code in sorted(codes, key=int):
                mf.write("%s\t%s\n" % (code, codes[code]))

    print("wrote %d high-confidence haplotype segments to %s" % (n_written, args.output))


if __name__ == "__main__":
    main()
