#!/usr/bin/env python3
"""Build a per-cycle, per-reference-base empirical error model for Polyester
from a representative aligned BAM (paired-end, fixed read length, MD tags)."""
import subprocess, re, os, argparse
from collections import defaultdict

ap = argparse.ArgumentParser()
ap.add_argument("bam", help="representative BAM aligned to the reference genome")
ap.add_argument("-o", "--outdir", required=True)
ap.add_argument("-s", "--subsample", default="0.03", help="samtools view -s fraction")
ap.add_argument("--samtools", default="samtools")
ap.add_argument("--readlen", type=int, default=150)
ap.add_argument("--prefix", default="empirical150")
args = ap.parse_args()

READLEN = args.readlen
BASES = "ACGT"
NT = ["A", "T", "G", "C", "N"]  # Polyester error-matrix row order
os.makedirs(args.outdir, exist_ok=True)

cnt = {1: [{b: defaultdict(int) for b in BASES} for _ in range(READLEN + 1)],
       2: [{b: defaultdict(int) for b in BASES} for _ in range(READLEN + 1)]}
MD_TOKEN = re.compile(r'(\d+)([A-Z]|\^[A-Z]+|$)')


def md_walk(cigar, md, seq, mate):
    qpos = 0
    events = []
    for num, tail in MD_TOKEN.findall(md):
        for _ in range(int(num)):
            events.append(('M', None))
        if tail.startswith('^'):
            for _ in tail[1:]:
                events.append(('D', None))
        elif tail != '':
            events.append(('X', tail))
    ei = 0
    for n, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar):
        L = int(n)
        if op in 'M=X':
            for k in range(L):
                if ei >= len(events):
                    return
                typ, rb = events[ei]
                ei += 1
                cyc = qpos + k + 1
                if cyc <= READLEN:
                    if typ == 'X':
                        cnt[mate][cyc][rb][seq[qpos + k]] += 1
                    else:
                        cnt[mate][cyc][seq[qpos + k]][seq[qpos + k]] += 1
            qpos += L
        elif op in 'IS':
            qpos += L
        elif op in 'DN':
            for k in range(L):
                if ei >= len(events):
                    return
                if events[ei][0] == 'D':
                    ei += 1
                else:
                    return
    return


proc = subprocess.Popen([args.samtools, "view", "-s", args.subsample, args.bam],
                        stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
nrec = nused = 0
for line in proc.stdout:
    f = line.rstrip("\n").split("\t")
    if len(f) < 11:
        continue
    flag = int(f[1])
    nrec += 1
    if flag & (0x4 | 0x100 | 0x800):
        continue
    if flag & 0x40:
        mate = 1
    elif flag & 0x80:
        mate = 2
    else:
        continue
    seq = f[9].upper()
    if len(seq) != READLEN:
        continue
    md = next((t[5:] for t in f[11:] if t.startswith("MD:Z:")), None)
    if md is None:
        continue
    nused += 1
    md_walk(f[5], md, seq, mate)
proc.stdout.close()
proc.wait()


def write_model(mate, path):
    lines = ["refbase readA readT readG readC readN pos"]
    for cyc in range(1, READLEN + 1):
        for rb in "ATGCN":
            if rb in BASES:
                tot = sum(cnt[mate][cyc][rb].values())
                if tot > 0:
                    probs = [(cnt[mate][cyc][rb].get(ob, 0) + 1e-7) / (tot + 5e-7) for ob in NT]
                else:
                    probs = [1.0 if ob == rb else 0.0 for ob in NT]
                s = sum(probs)
                probs[-1] += (1.0 - s)
            else:
                probs = [1.0 if ob == rb else 0.0 for ob in NT]
            lines.append("ref%s %s %d" % (rb, " ".join("%.8f" % p for p in probs), cyc))
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


for mate, name in ((1, "mate1"), (2, "mate2")):
    write_model(mate, os.path.join(args.outdir, "%s_%s" % (args.prefix, name)))

for mate in (1, 2):
    mm = tot = 0
    for cyc in range(1, READLEN + 1):
        for rb in BASES:
            d = cnt[mate][cyc][rb]
            tot += sum(d.values())
            mm += sum(v for k, v in d.items() if k != rb)
print("records=%d used=%d mismatch_rate=%.5f" % (nrec, nused, mm / max(1, tot)))
