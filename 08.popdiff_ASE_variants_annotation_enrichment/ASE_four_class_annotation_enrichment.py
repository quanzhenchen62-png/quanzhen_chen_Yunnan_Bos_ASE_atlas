#!/usr/bin/env python3
"""
Enrichment analysis for four ASE classes across:

1. chromatin states
2. Animal QTLdb annotations
3. genomic-region / functional annotations from raw ANN consequence names

The script assumes the input label table already contains the four final
ASE-class labels:

    indicine_taurine_label
    sai_eai_label
    yak_shared_label
    gayal_shared_label

The `threshold` column can be used for repeated ASE sets such as
brd-ASE-1, brd-ASE-2 and brd-ASE-3; each threshold is analyzed
independently.
"""

from __future__ import annotations

import argparse
import math
import shlex
import subprocess
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact


CATEGORY_SPECS = {
    "indicine_taurine": {"label_col": "indicine_taurine_label", "flag_col": "Indicine.vs.taurine"},
    "sai_eai": {"label_col": "sai_eai_label", "flag_col": "SAI.vs.EAI"},
    "yak_shared": {"label_col": "yak_shared_label", "flag_col": "yak-related introgression"},
    "gayal_shared": {"label_col": "gayal_shared_label", "flag_col": "gayal-related introgression"},
}

DEFAULT_CHROMHMM_TRACKS = {
    "Liver": ("liver", ""),
    "Lung": ("lung", ""),
    "Muscle": ("muscle", ""),
    "Spleen": ("spleen", ""),
    "Hump": ("muscle", "Hump_as_muscle"),
}


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser()
    ap.add_argument("--labels", required=True, help="Final ASE label table.")
    ap.add_argument("--chromhmm-dir", required=True, help="Directory containing ChromHMM BED files.")
    ap.add_argument("--qtl-bed", required=True, help="Animal QTLdb BED-like file with trait/class/category columns.")
    ap.add_argument("--outdir", required=True, help="Output directory.")
    return ap.parse_args()


def mkdir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def run(cmd: str) -> None:
    proc = subprocess.run(cmd, shell=True, text=True, capture_output=True)
    if proc.returncode != 0:
        raise RuntimeError(f"command failed: {cmd}\n{proc.stderr}")


def quote(path: Path | str) -> str:
    return shlex.quote(str(path))


def normalize_chr(values: pd.Series) -> pd.Series:
    return values.astype(str).str.replace("^chr", "", regex=True)


def bh(pvals: pd.Series) -> np.ndarray:
    p = pd.to_numeric(pvals, errors="coerce").to_numpy(dtype=float)
    out = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    if ok.sum() == 0:
        return out
    idx = np.where(ok)[0]
    order = idx[np.argsort(p[idx])]
    ranked = p[order] * ok.sum() / np.arange(1, ok.sum() + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out[order] = np.minimum(ranked, 1.0)
    return out


def add_bh_by_group(df: pd.DataFrame, p_col: str, group_cols: tuple[str, ...], out_col: str) -> pd.DataFrame:
    df = df.copy()
    df[out_col] = np.nan
    if df.empty:
        return df
    for _, idx in df.groupby(list(group_cols)).groups.items():
        df.loc[idx, out_col] = bh(df.loc[idx, p_col])
    return df


def odds_ci(a: int, b: int, c: int, d: int) -> tuple[float, float, float]:
    aa, bb, cc, dd = (a, b, c, d)
    if min(a, b, c, d) == 0:
        aa, bb, cc, dd = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    odds_ratio = (aa * dd) / (bb * cc) if (bb * cc) else np.nan
    se = math.sqrt((1 / aa) + (1 / bb) + (1 / cc) + (1 / dd))
    lo = math.exp(math.log(odds_ratio) - 1.96 * se)
    hi = math.exp(math.log(odds_ratio) + 1.96 * se)
    return odds_ratio, lo, hi


def fisher_row(a: int, b: int, c: int, d: int) -> dict[str, float]:
    if (a + b) == 0 or (c + d) == 0:
        return {
            "target_prop": np.nan,
            "background_prop": np.nan,
            "odds_ratio": np.nan,
            "ci95_low": np.nan,
            "ci95_high": np.nan,
            "fold_enrichment": np.nan,
            "p_value": np.nan,
        }
    p_value = fisher_exact([[a, b], [c, d]], alternative="two-sided")[1]
    odds_ratio, ci_low, ci_high = odds_ci(a, b, c, d)
    target_prop = a / (a + b)
    background_prop = c / (c + d)
    fold = (target_prop / background_prop) if background_prop > 0 else np.nan
    return {
        "target_prop": target_prop,
        "background_prop": background_prop,
        "odds_ratio": odds_ratio,
        "ci95_low": ci_low,
        "ci95_high": ci_high,
        "fold_enrichment": fold,
        "p_value": p_value,
    }


def make_interval_bed_from_sites(df: pd.DataFrame, path: Path) -> None:
    bed = df[["chrom", "pos", "variant_id"]].drop_duplicates().copy()
    bed["chrom"] = normalize_chr(bed["chrom"])
    bed["start"] = pd.to_numeric(bed["pos"], errors="coerce").astype(int) - 1
    bed["end"] = pd.to_numeric(bed["pos"], errors="coerce").astype(int)
    bed[["chrom", "start", "end", "variant_id"]].to_csv(path, sep="\t", index=False, header=False)


def read_labels(path: str) -> pd.DataFrame:
    labels = pd.read_csv(path, sep="\t", compression="infer")
    labels["chrom"] = normalize_chr(labels["chrom"])
    labels["threshold"] = pd.to_numeric(labels["threshold"], errors="coerce").astype("Int64")
    return labels


def enabled_subset(df: pd.DataFrame, flag_col: str) -> pd.DataFrame:
    return df.loc[df[flag_col].astype(str).str.lower() == "yes"].copy()


def chromhmm_enrichment(labels: pd.DataFrame, chromhmm_dir: Path, outdir: Path) -> None:
    tmpdir = outdir / "tmp_chromhmm"
    mkdir(tmpdir)

    annotated_tables = []
    qc_rows = []
    stat_rows = []

    for tissue, (track, proxy_reason) in DEFAULT_CHROMHMM_TRACKS.items():
        bed = chromhmm_dir / f"{track}_14_holstein_model_14states_segments.bed"
        if not bed.exists():
            qc_rows.append({"tissue": tissue, "status": "missing_chromhmm_bed", "path": str(bed)})
            continue

        tissue_df = labels.loc[labels["tissue"].astype(str) == tissue].copy()
        if tissue_df.empty:
            continue

        site_bed = tmpdir / f"{tissue}.sites.bed"
        overlap_tsv = tmpdir / f"{tissue}.chromhmm.overlap.tsv"
        make_interval_bed_from_sites(tissue_df, site_bed)
        run(f"bedtools intersect -a {quote(site_bed)} -b {quote(bed)} -wa -wb > {quote(overlap_tsv)}")

        if overlap_tsv.exists() and overlap_tsv.stat().st_size > 0:
            overlap = pd.read_csv(overlap_tsv, sep="\t", header=None, dtype=str)
        else:
            overlap = pd.DataFrame()

        state_map = {}
        if not overlap.empty:
            for variant_id, group in overlap.groupby(3):
                states = sorted(set(group.iloc[:, 7].astype(str)))
                state_map[str(variant_id)] = states[0] if len(states) == 1 else "MULTI_STATE"

        tissue_df["chromhmm_state"] = tissue_df["variant_id"].map(state_map).fillna("NO_STATE")
        tissue_df["annotation_tissue"] = track.capitalize()
        tissue_df["proxy_reason"] = proxy_reason
        annotated_tables.append(
            tissue_df[["tissue", "threshold", "variant_id", "chrom", "pos", "chromhmm_state", "annotation_tissue", "proxy_reason"]]
            .drop_duplicates()
        )

        valid_states = sorted([x for x in tissue_df["chromhmm_state"].unique() if x not in {"NO_STATE", "MULTI_STATE"}])
        for threshold, threshold_df in tissue_df.groupby("threshold", sort=True):
            for category, spec in CATEGORY_SPECS.items():
                layer = enabled_subset(threshold_df, spec["flag_col"])
                if layer.empty:
                    continue
                agg = layer.groupby("variant_id", as_index=False).agg({
                    "chromhmm_state": "first",
                    spec["label_col"]: "max",
                })
                label_values = pd.to_numeric(agg[spec["label_col"]], errors="coerce")
                usable = agg.loc[label_values.notna() & agg["chromhmm_state"].isin(valid_states)].copy()
                target = usable.loc[pd.to_numeric(usable[spec["label_col"]], errors="coerce") == 1].copy()
                background = usable.loc[pd.to_numeric(usable[spec["label_col"]], errors="coerce") == 0].copy()

                for state in valid_states:
                    a = int((target["chromhmm_state"] == state).sum())
                    b = len(target) - a
                    c = int((background["chromhmm_state"] == state).sum())
                    d = len(background) - c
                    stats = fisher_row(a, b, c, d)
                    stat_rows.append({
                        "tissue": tissue,
                        "threshold": int(threshold),
                        "category": category,
                        "state": state,
                        "target_hits": a,
                        "target_nonhits": b,
                        "background_hits": c,
                        "background_nonhits": d,
                        "target_n": len(target),
                        "background_n": len(background),
                        **stats,
                    })

    if annotated_tables:
        pd.concat(annotated_tables, ignore_index=True).drop_duplicates().to_csv(
            outdir / "chromatin_state_variant_annotations.tsv.gz",
            sep="\t",
            index=False,
            compression="gzip",
        )
    chrom = pd.DataFrame(stat_rows)
    if not chrom.empty:
        chrom = add_bh_by_group(chrom, p_col="p_value", group_cols=("threshold",), out_col="BH_FDR")
    chrom.to_csv(outdir / "chromatin_state_enrichment.tsv", sep="\t", index=False)
    pd.DataFrame(qc_rows).to_csv(outdir / "chromatin_state_qc.tsv", sep="\t", index=False)


def qtl_enrichment(labels: pd.DataFrame, qtl_bed: Path, outdir: Path) -> None:
    tmpdir = outdir / "tmp_qtl"
    mkdir(tmpdir)

    qtl = pd.read_csv(
        qtl_bed,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "trait", "class", "category"],
        dtype=str,
    )
    qtl["chrom"] = normalize_chr(qtl["chrom"])
    qtl["start"] = pd.to_numeric(qtl["start"], errors="coerce").astype("Int64")
    qtl["end"] = pd.to_numeric(qtl["end"], errors="coerce").astype("Int64")
    qtl = qtl.dropna(subset=["start", "end", "trait", "class", "category"]).copy()
    qtl.to_csv(outdir / "qtl.cleaned_categories.tsv.gz", sep="\t", index=False, compression="gzip")

    site_bed = tmpdir / "all_label_sites.bed"
    qtl_bed_clean = tmpdir / "qtl.cleaned.bed"
    overlap_tsv = tmpdir / "qtl.overlap.tsv"
    make_interval_bed_from_sites(labels, site_bed)
    qtl[["chrom", "start", "end", "trait", "class", "category"]].to_csv(qtl_bed_clean, sep="\t", index=False, header=False)
    run(f"bedtools intersect -a {quote(site_bed)} -b {quote(qtl_bed_clean)} -wa -wb > {quote(overlap_tsv)}")

    overlap = pd.read_csv(overlap_tsv, sep="\t", header=None, dtype=str) if overlap_tsv.stat().st_size else pd.DataFrame()
    variant_to_qtl = defaultdict(set)
    qtl_variant_sets = defaultdict(set)
    if not overlap.empty:
        for row in overlap.itertuples(index=False):
            qkey = (row[7], row[8], row[9])
            variant_id = row[3]
            variant_to_qtl[variant_id].add(qkey)
            qtl_variant_sets[qkey].add(variant_id)

    qtl_units = [
        tuple(x) for x in qtl[["trait", "class", "category"]]
        .drop_duplicates()
        .sort_values(["category", "class", "trait"])
        .to_numpy()
    ]

    coverage = qtl.groupby(["category", "class", "trait"], as_index=False).agg(
        n_qtl=("trait", "size"),
        total_interval_bp=("end", lambda x: int((qtl.loc[x.index, "end"].astype(int) - qtl.loc[x.index, "start"].astype(int)).sum())),
    )
    coverage.to_csv(outdir / "qtl.category_coverage.tsv", sep="\t", index=False)

    stat_rows = []
    for (breed, threshold), layer0 in labels.groupby(["breed", "threshold"], sort=True):
        for category, spec in CATEGORY_SPECS.items():
            layer = enabled_subset(layer0, spec["flag_col"])
            if layer.empty:
                continue
            label_values = pd.to_numeric(layer[spec["label_col"]], errors="coerce")
            usable = layer.loc[label_values.notna()].drop_duplicates("variant_id").copy()
            target = usable.loc[pd.to_numeric(usable[spec["label_col"]], errors="coerce") == 1].copy()
            background = usable.loc[pd.to_numeric(usable[spec["label_col"]], errors="coerce") == 0].copy()
            target_variants = set(target["variant_id"])
            background_variants = set(background["variant_id"])

            for qtl_trait, qtl_class, qtl_category in qtl_units:
                qkey = (qtl_trait, qtl_class, qtl_category)
                qvars = qtl_variant_sets.get(qkey, set())
                a = len(target_variants & qvars)
                b = len(target) - a
                c = len(background_variants & qvars)
                d = len(background) - c
                stats = fisher_row(a, b, c, d)
                stat_rows.append({
                    "breed": breed,
                    "threshold": int(threshold),
                    "category": category,
                    "qtl_trait": qtl_trait,
                    "qtl_class": qtl_class,
                    "qtl_category": qtl_category,
                    "target_hits": a,
                    "target_nonhits": b,
                    "background_hits": c,
                    "background_nonhits": d,
                    "target_n": len(target),
                    "background_n": len(background),
                    **stats,
                })

    stats = pd.DataFrame(stat_rows)
    if not stats.empty:
        stats = add_bh_by_group(stats, p_col="p_value", group_cols=("threshold",), out_col="BH_FDR")
    stats.to_csv(outdir / "qtl_enrichment.tsv", sep="\t", index=False)

    variant_matrix = pd.DataFrame(
        [
            {
                "variant_id": variant_id,
                "qtl_trait": trait,
                "qtl_class": qclass,
                "qtl_category": qcat,
                "overlap_at_least_one": 1,
            }
            for variant_id, qset in variant_to_qtl.items()
            for trait, qclass, qcat in qset
        ]
    )
    variant_matrix.to_csv(outdir / "qtl_variant_overlap.long.tsv.gz", sep="\t", index=False, compression="gzip")


def parse_ann_effects(labels: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    effect_map = {}
    all_effects = set()

    for sig_file in labels["sig_file"].dropna().drop_duplicates():
        table = pd.read_csv(
            sig_file,
            sep="\t",
            dtype=str,
            usecols=lambda col: col in {"contig", "position", "refAllele", "altAllele", "INFO"},
        )
        if "INFO" not in table.columns:
            continue
        for row in table.itertuples(index=False):
            info = getattr(row, "INFO", None)
            if not isinstance(info, str) or "ANN=" not in info:
                continue
            variant_id = f"{str(row.contig).replace('chr', '')}:{int(row.position)}:{row.refAllele}:{row.altAllele}"
            ann = info.split("ANN=", 1)[1].split(";", 1)[0]
            effects = effect_map.setdefault(variant_id, set())
            for transcript in ann.split(","):
                fields = transcript.split("|")
                if len(fields) > 1:
                    for effect in fields[1].split("&"):
                        if effect:
                            effects.add(effect)
                            all_effects.add(effect)

    effect_names = sorted(all_effects)
    rows = []
    for variant_id in sorted(labels["variant_id"].drop_duplicates()):
        effects = effect_map.get(variant_id, set())
        row = {"variant_id": variant_id, "annotation_source": "existing_ANN" if effects else "missing_ANN"}
        for effect in effect_names:
            row[effect] = int(effect in effects)
        rows.append(row)
    return pd.DataFrame(rows), pd.DataFrame([{
        "n_unique_variants": labels["variant_id"].nunique(),
        "n_with_ann": len(effect_map),
        "n_raw_ann_categories": len(effect_names),
    }])


def genomic_region_enrichment(labels: pd.DataFrame, outdir: Path) -> None:
    matrix, qc = parse_ann_effects(labels)
    matrix.to_csv(outdir / "genomic_region_variant_matrix.tsv.gz", sep="\t", index=False, compression="gzip")
    qc.to_csv(outdir / "genomic_region_qc.tsv", sep="\t", index=False)

    effect_names = [x for x in matrix.columns if x not in {"variant_id", "annotation_source"}]
    stat_rows = []

    for threshold, layer0 in labels.groupby("threshold", sort=True):
        for category, spec in CATEGORY_SPECS.items():
            layer = enabled_subset(layer0, spec["flag_col"])
            if layer.empty:
                continue
            label_values = pd.to_numeric(layer[spec["label_col"]], errors="coerce")
            target_ids = set(layer.loc[label_values == 1, "variant_id"])
            background_ids = set(layer.loc[label_values == 0, "variant_id"]) - target_ids
            target = matrix.loc[(matrix["variant_id"].isin(target_ids)) & (matrix["annotation_source"] != "missing_ANN")].copy()
            background = matrix.loc[(matrix["variant_id"].isin(background_ids)) & (matrix["annotation_source"] != "missing_ANN")].copy()

            for effect in effect_names:
                a = int(target[effect].sum())
                b = len(target) - a
                c = int(background[effect].sum())
                d = len(background) - c
                stats = fisher_row(a, b, c, d)
                stat_rows.append({
                    "threshold": int(threshold),
                    "category": category,
                    "function_category": effect,
                    "target_hits": a,
                    "target_nonhits": b,
                    "background_hits": c,
                    "background_nonhits": d,
                    "target_n": len(target),
                    "background_n": len(background),
                    **stats,
                })

    stats = pd.DataFrame(stat_rows)
    if not stats.empty:
        stats = add_bh_by_group(stats, p_col="p_value", group_cols=("threshold",), out_col="BH_FDR")
    stats.to_csv(outdir / "genomic_region_enrichment.tsv", sep="\t", index=False)


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    mkdir(outdir)

    labels = read_labels(args.labels)
    chromhmm_enrichment(labels, Path(args.chromhmm_dir), outdir)
    qtl_enrichment(labels, Path(args.qtl_bed), outdir)
    genomic_region_enrichment(labels, outdir)


if __name__ == "__main__":
    main()
