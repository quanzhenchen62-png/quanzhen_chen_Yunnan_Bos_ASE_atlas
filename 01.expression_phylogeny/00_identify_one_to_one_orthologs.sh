#!/bin/bash
# Identify strict one-to-one orthologs among three species from protein sets.
# Usage:
#   00_identify_one_to_one_orthologs.sh <species1.pep.fa> <species2.pep.fa> <species3.pep.fa> \
#     <species1.counts.tsv> <species2.counts.tsv> <species3.counts.tsv> \
#     <species1.tpm.tsv> <species2.tpm.tsv> <species3.tpm.tsv> <output_dir> [species3.gtf]

set -euo pipefail

if [ "$#" -lt 10 ] || [ "$#" -gt 11 ]; then
  echo "Usage: $0 <species1.pep.fa> <species2.pep.fa> <species3.pep.fa> \\" >&2
  echo "         <species1.counts.tsv> <species2.counts.tsv> <species3.counts.tsv> \\" >&2
  echo "         <species1.tpm.tsv> <species2.tpm.tsv> <species3.tpm.tsv> <output_dir> [species3.gtf]" >&2
  exit 1
fi

PEP1=$1; PEP2=$2; PEP3=$3
CNT1=$4; CNT2=$5; CNT3=$6
TPM1=$7; TPM2=$8; TPM3=$9
OUT=${10}
GTF3=${11:-}
mkdir -p "$OUT/proteins"

# ---- step 1: longest protein per gene ----
extract_longest() {
  # $1 = protein FASTA, $2 = output FASTA, $3 = optional GTF (protein->gene map)
  python3 - "$1" "$2" "$3" <<'PYEOF'
import gzip, re, sys
fa, out_path, gtf = sys.argv[1], sys.argv[2], sys.argv[3]
prot2gene = {}
if gtf:
    opener = gzip.open if gtf.endswith(".gz") else open
    for line in opener(gtf, "rt"):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) != 9 or f[2] != "CDS":
            continue
        gm = re.search(r'(?:^|; )gene_id "([^"]+)"', f[8])
        pm = re.search(r'(?:^|; )protein_id "([^"]+)"', f[8])
        if not gm or not pm:
            continue
        if pm.group(1) in prot2gene and prot2gene[pm.group(1)] != gm.group(1):
            raise ValueError("protein maps to multiple genes: %s" % pm.group(1))
        prot2gene[pm.group(1)] = gm.group(1)
best = {}          # gene -> (protein_id, sequence)
unmapped = 0
def store(header, seq):
    global unmapped
    s = "".join(seq).replace("*", "").upper()
    if not s:
        return
    m = re.search(r"(?:^| )gene:([^ ]+)", header)
    tokens = header.split()
    if m:
        # Ensembl path: require protein_coding gene_biotype, as in the original analysis
        bt = re.search(r"(?:^| )gene_biotype:([^ ]+)", header)
        if not bt or bt.group(1) != "protein_coding":
            return
        gene = re.sub(r"\.\d+$", "", m.group(1))
        rank = tokens[0] if tokens else ""
    else:
        tok = tokens[0] if tokens else None
        if gtf:
            if tok not in prot2gene:
                unmapped += 1
                return
            gene, rank = prot2gene[tok], tok
        elif tok:
            gene, rank = tok, tok
        else:
            return
    cur = best.get(gene)
    if cur is None or len(s) > len(cur[1]) or (len(s) == len(cur[1]) and rank < cur[0]):
        best[gene] = (rank, s)
name, seq = None, []
opener = gzip.open if fa.endswith(".gz") else open
for line in opener(fa, "rt"):
    if line.startswith(">"):
        if name is not None:
            store(name, seq)
        name, seq = line[1:].rstrip("\n"), []
    else:
        seq.append(line.strip())
if name is not None:
    store(name, seq)
with open(out_path, "w") as out:
    for g, (rank, s) in best.items():
        out.write(">%s\n%s\n" % (g, s))
print("genes retained: %d%s" % (len(best), (" (unmapped proteins skipped: %d)" % unmapped) if gtf else ""))
PYEOF
}
extract_longest "$PEP1" "$OUT/proteins/species1.fa" ""
extract_longest "$PEP2" "$OUT/proteins/species2.fa" ""
extract_longest "$PEP3" "$OUT/proteins/species3.fa" "$GTF3"

# ---- step 2: OrthoFinder ----
OF_OUT="$OUT/orthofinder_results"
rm -rf "$OF_OUT"
orthofinder -f "$OUT/proteins" -o "$OF_OUT" -t 8

# ---- step 3: keep orthogroups with exactly one gene per species, present in all matrices ----
OG=$(find "$OF_OUT" -type f -name Orthogroups.tsv -print -quit)
if [ -z "$OG" ]; then
  echo "Orthogroups.tsv not found under $OF_OUT" >&2
  exit 1
fi
python3 - "$OG" "$CNT1" "$CNT2" "$CNT3" "$TPM1" "$TPM2" "$TPM3" "$OUT/one_to_one_orthologs.tsv" <<'PYEOF'
import csv, gzip, sys
og_path = sys.argv[1]
def ids(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as h:
        r = csv.reader(h, delimiter="\t")
        return {row[0] for row in r if row and row[0]}
cnt = [ids(p) for p in sys.argv[2:5]]
tpm = [ids(p) for p in sys.argv[5:8]]
with open(og_path) as h:
    header = h.readline().rstrip("\n").split("\t")
    species_cols = header[1:]
    rows = []
    for line in h:
        f = line.rstrip("\n").split("\t")
        genes = [set(x.split(", ")) if x.strip() else set() for x in f[1:]]
        if all(len(g) == 1 for g in genes):
            g = [list(x)[0] for x in genes]
            if all(g[i] in cnt[i] and g[i] in tpm[i] for i in range(3)):
                rows.append([f[0]] + g)
with open(sys.argv[8], "w", newline="") as out:
    w = csv.writer(out, delimiter="\t", lineterminator="\n")
    w.writerow(["ortholog_id"] + [c.replace(" ", "_") for c in species_cols])
    w.writerows(rows)
print("strict one-to-one orthologs retained:", len(rows))
PYEOF
