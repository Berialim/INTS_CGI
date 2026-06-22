#!/usr/bin/env python3
"""MA-like plot for splicing error change after treatment."""

import argparse
import os
import re
import sys
import warnings
from bisect import bisect_right
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")

ATTR_RE = re.compile(r'([A-Za-z0-9_]+)\s+"([^"]+)"')


def setup_font():
    available = {f.name for f in fm.fontManager.ttflist}
    for font in ["Arial", "Helvetica", "Liberation Sans", "FreeSans", "DejaVu Sans"]:
        if font in available:
            plt.rcParams["font.family"] = font
            return
    plt.rcParams["font.family"] = "sans-serif"


setup_font()
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42


def parse_args():
    parser = argparse.ArgumentParser(description="Plot MA-like splicing error change after treatment.")
    parser.add_argument("-j", "--junction-dir", default="error_splicing")
    parser.add_argument("-d", "--database", required=True)
    parser.add_argument("-g", "--gtf", required=True)
    parser.add_argument("--ctrl", default="CTRL")
    parser.add_argument("--treat", default="IAA", help="Treatment condition to compare against CTRL")
    parser.add_argument("--metric", default="novel_junction_ratio", choices=["novel_junction_ratio", "novel_read_ratio"])
    parser.add_argument("--min-total-junctions", type=int, default=1)
    parser.add_argument("--min-total-reads", type=int, default=1)
    parser.add_argument("--refseq-map", default=None)
    parser.add_argument("--label-top-n", type=int, default=12)
    parser.add_argument("--output", default="error_splicing/splicing_error_MAplot.pdf")
    parser.add_argument("--table-output", default="error_splicing/splicing_error_MAplot.tsv")
    return parser.parse_args()


def parse_gtf_attributes(attr_text):
    return {key: value for key, value in ATTR_RE.findall(attr_text)}


def load_gtf_resources(gtf_path):
    gene_bounds = {}
    total_records = 0
    with open(gtf_path, "r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            chrom, _, feature, start, end, _, _, _, attrs = fields
            if feature != "transcript":
                continue
            total_records += 1
            if chrom.startswith("yr"):
                chrom = chrom[2:]
            strand = fields[6]
            attr_map = parse_gtf_attributes(attrs)
            gene_name = attr_map.get("gene_name") or attr_map.get("gene_id")
            if not gene_name:
                continue
            entry = gene_bounds.setdefault(chrom, {}).setdefault(gene_name, {"strand": strand, "intervals": []})
            entry["intervals"].append((int(start), int(end)))

    chrom_index = {}
    for chrom, genes in gene_bounds.items():
        entries = []
        for gene_name, meta in genes.items():
            intervals = meta["intervals"]
            entries.append((min(x[0] for x in intervals), max(x[1] for x in intervals), gene_name, meta["strand"]))
        entries.sort(key=lambda x: (x[0], x[1], x[2]))
        chrom_index[chrom] = {"entries": entries, "starts": [x[0] for x in entries]}
    print(f"[INFO] loaded {total_records} transcripts", file=sys.stderr)
    return chrom_index


def genes_covering_site(chrom_index, chrom, pos, required_strand=None):
    if chrom not in chrom_index:
        return set()
    data = chrom_index[chrom]
    idx = bisect_right(data["starts"], pos)
    genes = set()
    for i in range(idx):
        start, end, gene_name, gene_strand = data["entries"][i]
        if end < pos:
            continue
        if start <= pos <= end and (required_strand is None or gene_strand == required_strand):
            genes.add(gene_name)
    return genes


def load_sj_strand_map(sample):
    sj_path = Path(f"{sample}SJ.out.tab")
    if not sj_path.exists():
        return None
    strand_map = {}
    with open(sj_path, "r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.strip():
                continue
            fs = line.rstrip("\n").split("\t")
            if len(fs) >= 4:
                strand_map[(fs[0], int(fs[1]), int(fs[2]))] = fs[3]
    return strand_map


def sj_strand_to_gene_strand(strand_code):
    if str(strand_code) == "1":
        return "+"
    if str(strand_code) == "2":
        return "-"
    return None


def load_junction_table(path):
    df = pd.read_csv(path, sep="\t")
    df["chrom"] = df["chrom"].astype(str).str.replace(r"^yr", "", regex=True)
    df["donor"] = pd.to_numeric(df["donor"], errors="coerce")
    df["acceptor"] = pd.to_numeric(df["acceptor"], errors="coerce")
    df["read_count"] = pd.to_numeric(df["read_count"], errors="coerce").fillna(0)
    return df.dropna(subset=["donor", "acceptor"])


def build_gene_level_table(junction_dir, db, chrom_index):
    sample_to_condition = dict(zip(db["name"].astype(str), db["condition"]))
    all_rows = []
    for path in sorted(Path(junction_dir).glob("*.junctions.tsv")):
        sample = path.name.replace(".junctions.tsv", "")
        if sample not in sample_to_condition:
            continue
        condition = sample_to_condition[sample]
        strand_map = load_sj_strand_map(sample)
        if strand_map is None:
            continue
        df = load_junction_table(path)
        gene_counts = {}
        for row in df.itertuples(index=False):
            strand_code = strand_map.get((row.chrom, int(row.donor) + 1, int(row.acceptor)))
            required_strand = sj_strand_to_gene_strand(strand_code)
            if required_strand is None:
                continue
            donor_genes = genes_covering_site(chrom_index, row.chrom, int(row.donor), required_strand)
            acceptor_genes = genes_covering_site(chrom_index, row.chrom, int(row.acceptor), required_strand)
            same_genes = donor_genes & acceptor_genes
            if not same_genes:
                continue
            is_novel = str(row.type).lower() != "annotated"
            for gene in same_genes:
                stat = gene_counts.setdefault(
                    gene,
                    {"sample": sample, "condition": condition, "gene_name": gene,
                     "total_junctions": 0, "novel_junctions": 0, "total_reads": 0.0, "novel_reads": 0.0}
                )
                stat["total_junctions"] += 1
                stat["total_reads"] += float(row.read_count)
                if is_novel:
                    stat["novel_junctions"] += 1
                    stat["novel_reads"] += float(row.read_count)
        sample_df = pd.DataFrame(gene_counts.values())
        if not sample_df.empty:
            sample_df["novel_junction_ratio"] = np.where(
                sample_df["total_junctions"] > 0,
                sample_df["novel_junctions"] / sample_df["total_junctions"] * 100.0,
                np.nan,
            )
            sample_df["novel_read_ratio"] = np.where(
                sample_df["total_reads"] > 0,
                sample_df["novel_reads"] / sample_df["total_reads"] * 100.0,
                np.nan,
            )
            all_rows.append(sample_df)
    if not all_rows:
        raise ValueError("No gene-level results generated")
    return pd.concat(all_rows, ignore_index=True)


def adjust_bh(pvalues):
    pvalues = np.asarray(pvalues, dtype=float)
    out = np.full_like(pvalues, np.nan)
    mask = np.isfinite(pvalues)
    if not np.any(mask):
        return out
    vals = pvalues[mask]
    order = np.argsort(vals)
    ranked = vals[order]
    n = len(ranked)
    adj = np.minimum.accumulate((ranked * n / np.arange(1, n + 1))[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    restored = np.empty_like(adj)
    restored[order] = adj
    out[mask] = restored
    return out


def build_volcano_table(data, ctrl, treat, metric, min_total_junctions, min_total_reads):
    data = data.loc[(data["total_junctions"] >= min_total_junctions) & (data["total_reads"] >= min_total_reads)].copy()
    rows = []
    for gene, gene_df in data.groupby("gene_name", sort=False):
        ctrl_vals = gene_df.loc[gene_df["condition"] == ctrl, metric].dropna().values
        treat_vals = gene_df.loc[gene_df["condition"] == treat, metric].dropna().values
        if len(ctrl_vals) == 0 or len(treat_vals) == 0:
            continue
        mean_ctrl = float(np.mean(ctrl_vals))
        mean_treat = float(np.mean(treat_vals))
        diff = mean_treat - mean_ctrl
        pvalue = np.nan
        if len(ctrl_vals) >= 2 and len(treat_vals) >= 2:
            try:
                _, pvalue = stats.mannwhitneyu(ctrl_vals, treat_vals, alternative="two-sided")
            except Exception:
                pvalue = np.nan
        rows.append({
            "gene_name": gene,
            "ctrl_mean": mean_ctrl,
            "treat_mean": mean_treat,
            "difference": diff,
            "pvalue": pvalue,
        })
    result = pd.DataFrame(rows)
    result["padj"] = adjust_bh(result["pvalue"].values)
    safe_p = pd.to_numeric(result["pvalue"], errors="coerce").clip(lower=1e-300)
    result["neglog10_pvalue"] = -np.log10(safe_p)
    result["mean_signal"] = (result["ctrl_mean"] + result["treat_mean"]) / 2.0
    result["direction"] = np.where(result["difference"] > 0, "Up", np.where(result["difference"] < 0, "Down", "Flat"))
    return result


def load_refseq_name_map(path):
    df = pd.read_csv(path, sep="\t")
    if not {"name", "name2"}.issubset(df.columns):
        raise ValueError("RefSeq_all.txt must contain name and name2 columns")
    return (
        df.loc[:, ["name", "name2"]]
        .dropna()
        .drop_duplicates(subset=["name"])
        .rename(columns={"name": "transcript_id", "name2": "gene_symbol"})
    )


def add_gene_symbols(df, refseq_path):
    name_map = load_refseq_name_map(refseq_path)
    merged = df.merge(name_map, left_on="gene_name", right_on="transcript_id", how="left")
    merged["plot_label"] = merged["gene_symbol"].fillna(merged["gene_name"])
    return merged.drop(columns=["transcript_id"], errors="ignore")


def plot_ma_like(df, treat, metric, output_path, label_top_n):
    fig, ax = plt.subplots(figsize=(6.4, 6.0), facecolor="white")
    colors = df["direction"].map({"Up": "#D55E00", "Down": "#0072B2", "Flat": "#BDBDBD"}).fillna("#BDBDBD")
    ax.scatter(df["mean_signal"], df["difference"], c=colors, s=14, alpha=0.78, linewidths=0)
    ax.axhline(0, color="#333333", linestyle=":", linewidth=1.0)
    ax.set_xlabel(f"Mean splicing error ({metric})")
    ax.set_ylabel(f"{treat} - CTRL ({metric})")
    ax.set_title(f"Splicing error MA-like plot: {treat} vs CTRL", fontweight="700")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(alpha=0.2, linestyle="--", linewidth=0.5)

    top_df = (
        df.assign(abs_difference=df["difference"].abs())
        .sort_values(["abs_difference", "neglog10_pvalue"], ascending=[False, False])
        .head(label_top_n)
    )
    for row in top_df.itertuples(index=False):
        ax.text(
            row.mean_signal,
            row.difference,
            str(row.plot_label),
            fontsize=8,
            ha="left",
            va="bottom",
            color="#222222",
        )

    fig.tight_layout()
    fig.savefig(output_path, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)


def main():
    args = parse_args()
    db = pd.read_csv(args.database, dtype={"name": str})
    chrom_index = load_gtf_resources(args.gtf)
    gene_sample_table = build_gene_level_table(args.junction_dir, db, chrom_index)
    volcano_df = build_volcano_table(
        gene_sample_table,
        ctrl=args.ctrl,
        treat=args.treat,
        metric=args.metric,
        min_total_junctions=args.min_total_junctions,
        min_total_reads=args.min_total_reads,
    )
    volcano_df = add_gene_symbols(volcano_df, args.refseq_map)
    os.makedirs(os.path.dirname(args.table_output) or ".", exist_ok=True)
    volcano_df.to_csv(args.table_output, sep="\t", index=False)
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    plot_ma_like(volcano_df, args.treat, args.metric, args.output, args.label_top_n)
    print(f"[INFO] table saved: {args.table_output}", file=sys.stderr)
    print(f"[INFO] figure saved: {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
