#!/usr/bin/env python3
"""Build Polyester count matrices: fixed fragments per haplotype per replicate.
50:50 -> REF=1000, ALT=1000; 70:30 -> REF=700, ALT=300 (alt:ref ratio 0.3)."""
import os, argparse

ap = argparse.ArgumentParser()
ap.add_argument("--haplotypes", required=True, help="haplotype FASTA from step 02")
ap.add_argument("--species", required=True, help="species label used in file names")
ap.add_argument("--species-code", type=int, required=True, help="1/2/3... for seed spacing")
ap.add_argument("--outdir", required=True)
ap.add_argument("--n-rep", type=int, default=10)
ap.add_argument("--ref-alt-counts", default="1000,1000;700,300",
                help="'<ref>,<alt> per ratio, semicolon-separated, matching --ratios'")
ap.add_argument("--ratios", default="50_50,70_30")
args = ap.parse_args()

ids = [line[1:].split()[0] for line in open(args.haplotypes) if line.startswith(">")]
n_ref = sum(1 for i in ids if i.endswith("__REF"))
n_alt = sum(1 for i in ids if i.endswith("__ALT"))
assert n_ref == n_alt and len(ids) == n_ref + n_alt, "REF/ALT entries not paired"

os.makedirs(args.outdir, exist_ok=True)
counts = [tuple(map(int, c.split(","))) for c in args.ref_alt_counts.split(";")]
ratios = args.ratios.split(",")
assert len(counts) == len(ratios)

seed_rows = ["species\tratio\treplicate\tseed"]
for ratio, (rc, ac) in zip(ratios, counts):
    out = os.path.join(args.outdir, "readmat_%s.%s.tsv" % (ratio, args.species))
    with open(out, "w") as fh:
        fh.write("hap_id\t" + "\t".join("rep%02d" % k for k in range(1, args.n_rep + 1)) + "\n")
        for i in ids:
            v = rc if i.endswith("__REF") else ac
            fh.write(i + "\t" + "\t".join([str(v)] * args.n_rep) + "\n")
    for k in range(1, args.n_rep + 1):
        seed_rows.append("%s\t%s\trep%02d\t%d" % (args.species, ratio, k,
                                                   args.species_code * 1000 + k))

with open(os.path.join(args.outdir, "seeds.%s.tsv" % args.species), "w") as fh:
    fh.write("\n".join(seed_rows) + "\n")
print("%s: %d haplotype entries; ratios=%s counts=%s; %d replicates"
      % (args.species, len(ids), ratios, counts, args.n_rep))
