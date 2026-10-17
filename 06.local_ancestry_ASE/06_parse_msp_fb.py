#!/usr/bin/env python3
"""Parse RFMix MSP/FB output for one chromosome into diploid tracts and
haplotype-level segments with per-ancestry posterior summaries.

A tract is high-confidence when both haplotypes have mean assigned-ancestry
posterior >= 0.80 (with at least one FB marker) and the tract contains
>= 30 RFMix SNPs.
"""
import argparse
import bisect
import gzip
import statistics

BREED_PREFIX = [("HHHNY.", "Honghe"), ("HHHN.", "Honghe"), ("JCHN.", "Jiangcheng"),
                ("DC.", "Dengchuan"), ("DQ.", "Diqing"), ("DZ.", "Dianzhong"),
                ("WS.", "Wenshan"), ("ZT.", "Zhaotong")]


def breed_of(sample):
    for prefix, breed in BREED_PREFIX:
        if sample.startswith(prefix):
            return breed
    return "NA"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chromosome", required=True)
    parser.add_argument("--msp", required=True)
    parser.add_argument("--fb", required=True)
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()

    chr_id = args.chromosome

    with open(args.msp) as f:
        header = f.readline().lstrip("#").strip()
    codes = {}
    for token in header.replace("Subpopulation order/codes:", "").split("\t"):
        if "=" in token:
            key, value = token.split("=")
            codes[int(value)] = key.strip()

    with open(args.msp) as f:
        f.readline()
        cols = f.readline().rstrip("\n").split("\t")
        smp_cols = cols[6:]
        samples = [c.rsplit(".", 1)[0] for c in smp_cols[0::2]]
        tracts = {s: [] for s in samples}
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 6 + 2 * len(samples):
                continue
            spos, epos, nsnp = int(p[1]), int(p[2]), int(p[5])
            for i, s in enumerate(samples):
                h1 = int(p[6 + 2 * i])
                h2 = int(p[7 + 2 * i])
                t = tracts[s]
                if t and t[-1]["h1"] == h1 and t[-1]["h2"] == h2 and t[-1]["end"] + 1 >= spos:
                    t[-1]["end"] = epos
                    t[-1]["n_snps"] += nsnp
                    t[-1]["n_msp"] += 1
                else:
                    t.append({"h1": h1, "h2": h2, "start": spos, "end": epos,
                              "n_snps": nsnp, "n_msp": 1})

    for s in samples:
        for t in tracts[s]:
            for h in ("h1", "h2"):
                t[h + "_sum"] = 0.0
                t[h + "_n"] = 0
                t[h + "_min"] = None
                t[h + "_ge80"] = 0
                t[h + "_vals"] = []

    with open(args.fb) as f:
        f.readline()
        fbcols = f.readline().rstrip("\n").split("\t")
        base = {}
        ci = 4
        for s in samples:
            base[s] = ci
            ci += 14
        starts = {s: [t["start"] for t in tracts[s]] for s in samples}
        for line in f:
            p = line.rstrip("\n").split("\t")
            pos = int(p[1])
            for s in samples:
                st = starts[s]
                j = bisect.bisect_right(st, pos) - 1
                if j < 0:
                    continue
                t = tracts[s][j]
                if pos > t["end"]:
                    continue
                b = base[s]
                for h, off in (("h1", 0), ("h2", 7)):
                    code = t[h]
                    v = float(p[b + off + code])
                    t[h + "_sum"] += v
                    t[h + "_n"] += 1
                    t[h + "_ge80"] += (1 if v >= 0.80 else 0)
                    if t[h + "_min"] is None or v < t[h + "_min"]:
                        t[h + "_min"] = v
                    t[h + "_vals"].append(v)

    dip_out = gzip.open("%s/chr%s.diploid_tracts.tsv.gz" % (args.out_dir, chr_id), "wt")
    hap_out = gzip.open("%s/chr%s.haplotype_segments.tsv.gz" % (args.out_dir, chr_id), "wt")
    dip_hdr = ["chromosome", "sample", "breed", "diploid_tract_id", "start", "end", "length_bp",
               "hap1_code", "hap1_label", "hap2_code", "hap2_label", "ordered_diploid_pair",
               "unordered_pair", "tract_n_snps", "n_msp_intervals", "n_fb_markers",
               "hap1_mean_assignment_probability", "hap1_median", "hap1_min", "hap1_prop_ge080",
               "hap2_mean_assignment_probability", "hap2_median", "hap2_min", "hap2_prop_ge080",
               "posterior_pass", "snp30_pass", "high_confidence_local_ancestry_segment", "failure_reason"]
    hap_hdr = ["chromosome", "sample", "breed", "haplotype", "diploid_tract_id", "start", "end", "length_bp",
               "ancestry_numeric_code", "ancestry_label", "paired_haplotype_ancestry", "ordered_diploid_pair",
               "n_msp_intervals", "n_fb_markers", "n_rfmix_snps", "assigned_posterior_mean",
               "assigned_posterior_median", "assigned_posterior_min", "prop_markers_posterior_ge_080",
               "posterior_pass", "snp30_pass", "high_confidence_local_ancestry_segment", "failure_reason"]
    dip_out.write("\t".join(dip_hdr) + "\n")
    hap_out.write("\t".join(hap_hdr) + "\n")

    n_total = n_hc = 0
    for s in samples:
        br = breed_of(s)
        for k, t in enumerate(tracts[s], 1):
            n_total += 1
            tid = "chr%s_%s_%06d" % (chr_id, s, k)

            def stat(h):
                n = t[h + "_n"]
                if n == 0:
                    return ("NA", "NA", "NA", "NA")
                mean = t[h + "_sum"] / n
                med = statistics.median(t[h + "_vals"])
                mn = ("%.6f" % t[h + "_min"]) if t[h + "_min"] is not None else "NA"
                return ("%.6f" % mean, "%.6f" % med, mn, "%.6f" % (t[h + "_ge80"] / n))

            m1, md1, mn1, pr1 = stat("h1")
            m2, md2, mn2, pr2 = stat("h2")
            fb_n = max(t["h1_n"], t["h2_n"])
            post_pass = (t["h1_n"] > 0 and t["h2_n"] > 0
                         and t["h1_sum"] / t["h1_n"] >= 0.80
                         and t["h2_sum"] / t["h2_n"] >= 0.80)
            snp_pass = t["n_snps"] >= 30
            hc = post_pass and snp_pass
            if hc:
                n_hc += 1
            reasons = []
            if not post_pass:
                if t["h1_n"] == 0 or t["h2_n"] == 0:
                    reasons.append("no_fb_markers")
                else:
                    if t["h1_sum"] / t["h1_n"] < 0.80:
                        reasons.append("hap1_mean_lt_0.80")
                    if t["h2_sum"] / t["h2_n"] < 0.80:
                        reasons.append("hap2_mean_lt_0.80")
            if not snp_pass:
                reasons.append("n_rfmix_snps_lt_30")
            fr = ";".join(reasons) if reasons else ""
            l1 = codes[t["h1"]]
            l2 = codes[t["h2"]]
            ordered = "%s/%s" % (l1, l2)
            unordered = "/".join(sorted([l1, l2]))
            dip_out.write("\t".join(map(str, [chr_id, s, br, tid, t["start"], t["end"],
                    t["end"] - t["start"] + 1, t["h1"], l1, t["h2"], l2, ordered, unordered,
                    t["n_snps"], t["n_msp"], fb_n, m1, md1, mn1, pr1, m2, md2, mn2, pr2,
                    post_pass, snp_pass, hc, fr])) + "\n")
            for h, hkey, lab, plab, mean, med, mn, pr in (("hap1", "h1", l1, l2, m1, md1, mn1, pr1),
                                                           ("hap2", "h2", l2, l1, m2, md2, mn2, pr2)):
                hap_out.write("\t".join(map(str, [chr_id, s, br, h, tid, t["start"], t["end"],
                    t["end"] - t["start"] + 1, t[hkey], lab, plab, ordered, t["n_msp"], fb_n,
                    t["n_snps"], mean, med, mn, pr,
                    post_pass, snp_pass, hc, fr])) + "\n")
    dip_out.close()
    hap_out.close()

    with open("%s/chr%s.tract_parsing.qc.tsv" % (args.out_dir, chr_id), "w") as q:
        q.write("chromosome\tn_samples\tn_total_diploid_tracts\tn_high_confidence_tracts\tstatus\n")
        q.write("%s\t%d\t%d\t%d\tPASS\n" % (chr_id, len(samples), n_total, n_hc))
    print("chr%s: %d tracts, %d high-confidence" % (chr_id, n_total, n_hc))


if __name__ == "__main__":
    main()
