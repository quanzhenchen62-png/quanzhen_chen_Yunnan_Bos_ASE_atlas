#!/usr/bin/env python3
"""Convert simulated FASTA (Polyester output) to gzipped FASTQ with uniform Q40.

Polyester emits FASTA; sequencing errors from the empirical error model are
already embedded in the sequences, so the fixed quality string is a formal
requirement of downstream aligners/callers, not real quality values."""
import os, gzip, glob, argparse

ap = argparse.ArgumentParser()
ap.add_argument("--simroot", required=True,
                help="root containing <species>/ratio<R>/rep<K>/*_{1,2}.fasta.gz")
ap.add_argument("--outdir", required=True)
ap.add_argument("--manifest", default=None, help="output sample manifest path")
args = ap.parse_args()

os.makedirs(args.outdir, exist_ok=True)
manifest = args.manifest or os.path.join(args.outdir, "manifest_sim.tsv")
rows = ["sample_id\tspecies\tratio\treplicate\tR1\tR2"]
n = 0
for repdir in sorted(glob.glob(os.path.join(args.simroot, "*", "ratio*", "rep*"))):
    sp = os.path.basename(os.path.dirname(os.path.dirname(repdir)))
    ratio = os.path.basename(os.path.dirname(repdir))[5:]
    rep = os.path.basename(repdir)
    sid = "%s_%s_%s" % (sp, ratio, rep)
    outs = []
    for mate in (1, 2):
        src = glob.glob(os.path.join(repdir, "*_%d.fasta.gz" % mate))
        assert len(src) == 1, (repdir, src)
        dst = os.path.join(args.outdir, sid + "_%d.fastq.gz" % mate)
        with gzip.open(src[0], "rt") as fin, gzip.open(dst, "wt", compresslevel=6) as fout:
            name = None; seq = []
            for line in fin:
                if line.startswith(">"):
                    if name is not None:
                        fout.write("@%s\n%s\n+\n%s\n" % (name, "".join(seq), "I" * len(seq)))
                    name = line[1:].strip(); seq = []
                else:
                    seq.append(line.strip())
            if name is not None:
                fout.write("@%s\n%s\n+\n%s\n" % (name, "".join(seq), "I" * len(seq)))
        outs.append(dst)
    n += 1
    rows.append("\t".join([sid, sp, ratio, rep] + outs))
with open(manifest, "w") as fh:
    fh.write("\n".join(rows) + "\n")
print("converted %d samples -> %s" % (n, args.outdir))
