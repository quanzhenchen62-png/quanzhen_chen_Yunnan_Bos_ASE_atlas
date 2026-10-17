#!/usr/bin/env python3
"""Verify simulated-read truth: parse Polyester FASTA headers
>read<N>/<hap_id>;mate1:<s>-<e>;mate2:<s>-<e> and count Ref/Alt fragments
and SNP-covering reads per locus, plus pairing and record sanity."""
import os, re, sys, gzip, glob, argparse
from collections import defaultdict

HDR = re.compile(r'^read(\d+)/(.+?);mate1:(-?\d+)-(\d+);mate2:(-?\d+)-(\d+)$')

ap = argparse.ArgumentParser()
ap.add_argument("--simroot", required=True,
                help="root containing <species>/ratio<R>/rep<K>/*_{1,2}.fasta.gz")
ap.add_argument("--target-loci", required=True, nargs="+",
                help="per-species target_loci.tsv from step 02")
ap.add_argument("--seeds", required=True, nargs="+", help="per-species seeds.tsv from step 03")
ap.add_argument("--outdir", required=True)
ap.add_argument("--species", nargs="+", default=["cattle", "yak", "gayal"])
args = ap.parse_args()

by_locus = {}; by_hapref = {}
for p in args.target_loci:
    with open(p) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            d = dict(zip(header, line.rstrip("\n").split("\t")))
            by_locus[d["locus_id"]] = d
            by_hapref[d["haplotype_ref_id"]] = d

seeds = {}
for p in args.seeds:
    with open(p) as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            seeds[(f[0], f[1], f[2])] = f[3]


def read_fasta(path):
    op = gzip.open if path.endswith(".gz") else open
    name = None; seq = []
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name = line[1:].strip(); seq = []
            else:
                seq.append(line.strip())
    if name is not None:
        yield name, "".join(seq)


truth_rows = []; read_rows = []; pair_rows = []; mal_total = 0
for sp in args.species:
    for repdir in sorted(glob.glob(os.path.join(args.simroot, sp, "ratio*", "rep*"))):
        ratio = os.path.basename(os.path.dirname(repdir))[5:]  # strip 'ratio'
        rep = os.path.basename(repdir)
        f1s = glob.glob(repdir + "/*_1.fasta*"); f2s = glob.glob(repdir + "/*_2.fasta*")
        assert len(f1s) == 1 and len(f2s) == 1, (repdir, f1s, f2s)
        reads1 = {}; reads2 = {}
        dup = mal = ncount = negcoord = 0
        for path, store in ((f1s[0], reads1), (f2s[0], reads2)):
            for name, seq in read_fasta(path):
                m = HDR.match(name)
                if not m:
                    mal += 1; continue
                if int(m.group(3)) < 1 or int(m.group(5)) < 1:
                    negcoord += 1
                n = int(m.group(1))
                if n in store: dup += 1
                store[n] = (m.group(2), int(m.group(3)), int(m.group(4)),
                            int(m.group(5)), int(m.group(6)), seq)
                ncount += seq.count("N")
        mal_total += mal
        common = sorted(set(reads1) & set(reads2))
        origin_mismatch = sum(1 for n in common if reads1[n][0] != reads2[n][0])
        agg = defaultdict(lambda: dict(ref_frag=0, alt_frag=0, ref_reads=0, alt_reads=0))
        for n in common:
            hap = reads1[n][0]
            tid, tag = hap.rsplit("__", 1)
            locus = by_hapref.get(tid + "__REF")
            if locus is None: continue
            tp = int(locus["transcript_position"])
            a1, b1 = reads1[n][1], reads1[n][2]
            c1, d1 = reads2[n][1], reads2[n][2]
            fs, fe = min(a1, c1), max(b1, d1)
            if fs <= tp <= fe:
                key = locus["locus_id"]
                agg[key]["ref_frag" if tag == "REF" else "alt_frag"] += 1
                if a1 <= tp <= b1:
                    agg[key]["ref_reads" if tag == "REF" else "alt_reads"] += 1
                if c1 <= tp <= d1:
                    agg[key]["ref_reads" if tag == "REF" else "alt_reads"] += 1
        nom = 0.5 if ratio == "50_50" else 0.7
        seed = seeds.get((sp, ratio, rep), "NA")
        for key in sorted(agg):
            a = agg[key]; locus = by_locus[key]
            tr, ta = a["ref_frag"], a["alt_frag"]
            truth_rows.append([sp, ratio, rep, locus["transcript_id"], key, tr, ta,
                               a["ref_reads"], a["alt_reads"], nom,
                               round(tr / max(1, tr + ta), 4), round(1 - nom, 4),
                               round(ta / max(1, tr + ta), 4), seed])
        read_rows.append([sp, ratio, rep, len(reads1), len(reads2), len(common),
                          origin_mismatch, dup, mal, ncount, negcoord])
        ok = (len(reads1) == len(reads2) == len(common) and origin_mismatch == 0
              and dup == 0 and mal == 0)
        pair_rows.append([sp, ratio, rep, len(reads1), len(reads2), len(common),
                          "PASS" if ok else "CHECK"])

os.makedirs(args.outdir, exist_ok=True)
with open(os.path.join(args.outdir, "replicate_truth_counts.tsv"), "w") as fh:
    fh.write("species\tratio\treplicate\ttranscript_id\tlocus_id\ttotal_ref_fragments"
             "\ttotal_alt_fragments\tref_reads_covering_snp\talt_reads_covering_snp"
             "\tnominal_ref_fraction\trealized_ref_fraction\tnominal_alt_fraction"
             "\trealized_alt_fraction\tseed\n")
    for r in truth_rows: fh.write("\t".join(map(str, r)) + "\n")
with open(os.path.join(args.outdir, "replicate_read_summary.tsv"), "w") as fh:
    fh.write("species\tratio\treplicate\tn_reads_mate1\tn_reads_mate2\tn_pairs"
             "\thap_origin_mismatch_pairs\tduplicate_read_names\tmalformed_records"
             "\tN_bases\tnegative_coord_records\n")
    for r in read_rows: fh.write("\t".join(map(str, r)) + "\n")
with open(os.path.join(args.outdir, "paired_end_validation.tsv"), "w") as fh:
    fh.write("species\tratio\treplicate\tn_mate1\tn_mate2\tn_common_pairs\tstatus\n")
    for r in pair_rows: fh.write("\t".join(map(str, r)) + "\n")
print("verified %d replicates; malformed_total=%d" % (len(read_rows), mal_total))
