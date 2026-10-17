#!/usr/bin/env python3
"""Select target heterozygous SNP loci and build Ref/Alt haplotype transcript FASTA.

Each transcript yields <tid>__REF and <tid>__ALT: the species transcript sequence
with only the target base set to the reference or alternative allele of a real
heterozygous SNP (reference-genome coordinates). Background sequence divergence
of the species is preserved.
"""
import subprocess, sys, os, re, gzip, random, bisect, json, argparse

CIGAR_RE = re.compile(r'(\d+)([MIDNSHP=X])')
ACGT = set("ACGT")

ap = argparse.ArgumentParser()
ap.add_argument("--species", required=True, help="species label, e.g. cattle")
ap.add_argument("--transcripts", required=True, help="species transcript FASTA")
ap.add_argument("--cigar-map", required=True,
                help="json.gz: {tid: [contig, pos0, cigar, strand, aln_len]} from a splice "
                     "alignment of the transcripts to the reference genome (e.g. minimap2 -ax splice)")
ap.add_argument("--vcf", required=True, help="heterozygous SNP VCF(bcf) in reference-genome coordinates")
ap.add_argument("--genome", required=True, help="reference genome FASTA")
ap.add_argument("--outdir", required=True)
ap.add_argument("--n-targets", type=int, default=120)
ap.add_argument("--n-min", type=int, default=100)
ap.add_argument("--seed", type=int, default=42)
ap.add_argument("--require-native-match", action="store_true",
                help="require transcript base == genome base at the target (use for the reference species)")
ap.add_argument("--bcftools", default="bcftools")
args = ap.parse_args()


def load_fasta(path):
    seqs = {}; tid = None; parts = []
    for line in open(path):
        if line.startswith(">"):
            if tid: seqs[tid] = "".join(parts)
            tid = line[1:].split()[0]; parts = []
        else:
            parts.append(line.strip())
    seqs[tid] = "".join(parts)
    return seqs


def load_het_snps(vcf):
    by_chrom = {}
    p = subprocess.Popen([args.bcftools, "view", "-H", "-v", "snps", "-m2", "-M2", vcf],
                         stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    for line in p.stdout:
        f = line.split("\t", 7)
        if len(f) < 7: continue
        chrom, pos, ref, alt, qual, filt = f[0], int(f[1]), f[3], f[4], f[5], f[6]
        if filt not in ("PASS", ".", ""): continue
        if qual not in (".", ""):
            try:
                if float(qual) < 30: continue
            except ValueError:
                continue
        if len(ref) == 1 and len(alt) == 1 and ref in ACGT and alt in ACGT:
            by_chrom.setdefault(chrom, []).append((pos, ref, alt))
    p.stdout.close(); p.wait()
    for c in by_chrom: by_chrom[c].sort()
    return by_chrom


def cigar_tp_map(pos, cigar):
    """Alignment start -> transcript coordinate map, M blocks, indel/clip events."""
    rp = 0; tp = 0
    mblocks = []; rp2tp = {}; events = []; ref_span = 0
    for n, op in CIGAR_RE.findall(cigar):
        n = int(n)
        if op in "M=X":
            mblocks.append((tp + 1, tp + n, rp + 1, rp + n))
            for k in range(n):
                rp2tp[rp + k] = tp + k + 1
            tp += n; rp += n; ref_span += n
        elif op in "IS":
            events.append((tp + 1, op)); tp += n
        elif op in "DN":
            events.append((tp + 1, op)); rp += n; ref_span += n
    return rp2tp, mblocks, events, ref_span


def low_complexity(seq, tp):
    w = seq[max(0, tp - 16):tp + 15]
    if len(w) < 10: return True
    for b in "ACGT":
        if b * 5 in w: return True
        if w.count(b) / len(w) > 0.7: return True
    return False


def flank_ok(seq, tp, gchr, apos, flank=10, max_mm=2):
    if gchr is None or apos - flank - 1 < 0 or apos + flank > len(gchr):
        return False, None
    mmL = sum(a != b for a, b in zip(seq[tp - 1 - flank:tp - 1], gchr[apos - 1 - flank:apos - 1]))
    mmR = sum(a != b for a, b in zip(seq[tp:tp + flank], gchr[apos:apos + flank]))
    gbase = gchr[apos - 1]
    ok = mmL <= max_mm and mmR <= max_mm and gbase in ACGT
    if args.require_native_match:
        ok = ok and (seq[tp - 1] == gbase)
    return ok, gbase


os.makedirs(args.outdir, exist_ok=True)
seqs = load_fasta(args.transcripts)
with gzip.open(args.cigar_map, "rt") as fh:
    aclass = json.load(fh)
snps = load_het_snps(args.vcf)
genome = load_fasta(args.genome)

cand = []
for tid, (contig, pos, cigar, strand, L) in aclass.items():
    seq = seqs.get(tid)
    if seq is None or L < 600 or "N" in seq: continue
    rp2tp, mblocks, events, ref_span = cigar_tp_map(pos, cigar)
    lst = snps.get(contig)
    if not lst: continue
    lo = bisect.bisect_left(lst, (pos, "", ""))
    hi = bisect.bisect_right(lst, (pos + ref_span + 1, "", ""))
    lo_t = max(200, int(0.25 * L)); hi_t = min(L - 200, int(0.75 * L))
    for (p, ref, alt) in lst[lo:hi]:
        rel = p - pos
        tp = rp2tp.get(rel)
        if tp is None or not (lo_t <= tp <= hi_t): continue
        ok = False
        for (t1, t2, r1, r2) in mblocks:
            if t1 <= tp <= t2:
                ok = (tp - t1 >= 10) and (t2 - tp >= 10)
                break
        if not ok: continue
        if any(abs(e_tp - tp) <= 10 for e_tp, _op in events): continue
        if low_complexity(seq, tp): continue
        native = seq[tp - 1]
        if native not in ACGT: continue
        ok, gbase = flank_ok(seq, tp, genome.get(contig), p)
        if not ok: continue
        cand.append(dict(tid=tid, tp=tp, contig=contig, pos=p, ref=ref, alt=alt,
                         native=native, L=L, strand=strand))
        break  # one locus per transcript

if len(cand) < args.n_min:
    sys.exit("fewer than %d target loci available" % args.n_min)
random.Random(args.seed).shuffle(cand)
chosen = cand[:args.n_targets]

fa_out = []; loci_rows = []
for locus in chosen:
    tid, tp = locus["tid"], locus["tp"]
    seq = seqs[tid]
    fa_out.append((tid + "__REF", seq[:tp - 1] + locus["ref"] + seq[tp:]))
    fa_out.append((tid + "__ALT", seq[:tp - 1] + locus["alt"] + seq[tp:]))
    loci_rows.append([args.species, tid, tid + "__REF", tid + "__ALT", tp,
                      locus["contig"], locus["pos"], locus["ref"], locus["alt"],
                      locus["native"], "yes" if locus["native"] != locus["ref"] else "no",
                      locus["L"], locus["strand"],
                      "%s_%s_%s_%d" % (args.species, tid, locus["contig"], locus["pos"])])

fa_path = os.path.join(args.outdir, args.species + ".haplotypes.fa")
with open(fa_path, "w") as fh:
    for name, s in fa_out:
        fh.write(">%s\n" % name)
        for i in range(0, len(s), 60):
            fh.write(s[i:i + 60] + "\n")

tsv_path = os.path.join(args.outdir, args.species + ".target_loci.tsv")
with open(tsv_path, "w") as fh:
    fh.write("species\ttranscript_id\thaplotype_ref_id\thaplotype_alt_id\ttranscript_position"
             "\tcontig\tposition\tref\talt\tnative_base\tref_hap_modified\ttranscript_length"
             "\tstrand\tlocus_id\n")
    for r in loci_rows:
        fh.write("\t".join(map(str, r)) + "\n")
print("%s: selected %d loci -> %s" % (args.species, len(chosen), fa_path))
