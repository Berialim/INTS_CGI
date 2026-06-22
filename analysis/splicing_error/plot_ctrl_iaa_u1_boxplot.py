#!/usr/bin/env python3
"""Plot CTRL/IAA/U1 gene-level splicing error ratios for three gene sets.

绘制三种基因集合的 gene-level 错误 splicing 比例箱线图：
1. active genes
2. with CGI
3. no CGI

每个面板都使用 seaborn boxplot，并用 statannotations.Annotator 标注 full pvalue，
最终保存为 PDF。
"""

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
import seaborn as sns
from scipy import stats
from statannotations.Annotator import Annotator

warnings.filterwarnings("ignore")


ATTR_RE = re.compile(r'([A-Za-z0-9_]+)\s+"([^"]+)"')
PALETTE = {
    "CTRL": "#4D4D4D",
    "IAA": "#D55E00",
    "U1": "#0072B2",
}
BED_COMPARE_PALETTE = {
    "with CGI": "#D55E00",
    "no CGI": "#0072B2",
}


def setup_font():
    available = {f.name for f in fm.fontManager.ttflist}
    preferred = ["Arial", "Helvetica", "Liberation Sans", "FreeSans", "DejaVu Sans"]
    for font in preferred:
        if font in available:
            plt.rcParams["font.family"] = font
            print(f"[INFO] 使用字体: {font}", file=sys.stderr)
            return
    plt.rcParams["font.family"] = "sans-serif"
    print("[INFO] 使用字体: sans-serif (系统默认)", file=sys.stderr)


setup_font()
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42


def parse_args():
    parser = argparse.ArgumentParser(
        description="绘制 active/with CGI/no CGI 三类基因在 CTRL/IAA/U1 下的 splicing error boxplot。"
    )
    parser.add_argument("-j", "--junction-dir", default="error_splicing", help="存放 *.junctions.tsv 的目录")
    parser.add_argument("-d", "--database", required=True, help="database.csv，至少包含 name,condition 两列")
    parser.add_argument("-g", "--gtf", required=True, help="基因注释 GTF")
    parser.add_argument(
        "--active-bed", required=True,
        help="active genes 对应 BED，第4列 transcript name",
    )
    parser.add_argument(
        "--with-cgi-bed", required=True,
        help="with CGI 对应 BED，第4列 transcript name",
    )
    parser.add_argument(
        "--no-cgi-bed", required=True,
        help="no CGI 对应 BED，第4列 transcript name",
    )
    parser.add_argument("--condition-order", nargs="+", default=["CTRL", "IAA", "U1"], help="condition 顺序")
    parser.add_argument("--test", default="mannwhitney", choices=["mannwhitney", "ttest"], help="组间检验方法")
    parser.add_argument("--min-total-junctions", type=int, default=1, help="最小总 junction 数阈值")
    parser.add_argument("--min-total-reads", type=int, default=1, help="最小总 reads 阈值")
    parser.add_argument("-o", "--output", default="error_splicing/ctrl_iaa_u1_splicing_error_three_sets_boxplot.pdf", help="输出 PDF 路径")
    parser.add_argument("--table-output", default="error_splicing/ctrl_iaa_u1_splicing_error_three_sets_summary.tsv", help="输出 gene-level 表")
    parser.add_argument("--dpi", type=int, default=150, help="PNG 输出 DPI")
    return parser.parse_args()


def parse_gtf_attributes(attr_text):
    return {key: value for key, value in ATTR_RE.findall(attr_text)}


def load_transcript_set(bed_path):
    transcripts = set()
    with open(bed_path, "r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 4 and fields[3]:
                transcripts.add(fields[3])
    print(f"[INFO] 从 BED 载入 {len(transcripts)} 个 transcripts: {bed_path}", file=sys.stderr)
    return transcripts


def load_gtf_resources(gtf_path):
    gene_bounds = {}
    transcript_to_gene = {}
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
            transcript_id = attr_map.get("transcript_id")
            gene_name = attr_map.get("gene_name") or attr_map.get("gene_id")
            if not transcript_id or not gene_name:
                continue
            entry = gene_bounds.setdefault(chrom, {}).setdefault(
                gene_name, {"strand": strand, "intervals": []}
            )
            entry["intervals"].append((int(start), int(end)))
            transcript_to_gene[transcript_id] = gene_name

    chrom_index = {}
    for chrom, genes in gene_bounds.items():
        entries = []
        for gene_name, meta in genes.items():
            intervals = meta["intervals"]
            entries.append(
                (min(x[0] for x in intervals), max(x[1] for x in intervals), gene_name, meta["strand"])
            )
        entries.sort(key=lambda x: (x[0], x[1], x[2]))
        chrom_index[chrom] = {"entries": entries, "starts": [x[0] for x in entries]}

    print(
        f"[INFO] 从 GTF 读取 {total_records} 个 transcript，得到 {sum(len(v['entries']) for v in chrom_index.values())} 个 gene 区间",
        file=sys.stderr,
    )
    return chrom_index, transcript_to_gene


def genes_covering_site(chrom_index, chrom, pos, required_strand=None):
    if chrom not in chrom_index:
        return set()
    data = chrom_index[chrom]
    idx = bisect_right(data["starts"], pos)
    genes = set()
    for i in range(idx):
        gene_start, gene_end, gene_name, gene_strand = data["entries"][i]
        if gene_end < pos:
            continue
        if gene_start <= pos <= gene_end:
            if required_strand is None or gene_strand == required_strand:
                genes.add(gene_name)
    return genes


def load_database(database_path):
    db = pd.read_csv(database_path, dtype={"name": str})
    required = {"name", "condition"}
    missing = required - set(db.columns)
    if missing:
        raise ValueError(f"database 缺少列: {sorted(missing)}")
    return db


def load_junction_table(path):
    df = pd.read_csv(path, sep="\t")
    required = {"chrom", "donor", "acceptor", "read_count", "type"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{path} 缺少列: {sorted(missing)}")
    df["chrom"] = df["chrom"].astype(str).str.replace(r"^yr", "", regex=True)
    df["donor"] = pd.to_numeric(df["donor"], errors="coerce")
    df["acceptor"] = pd.to_numeric(df["acceptor"], errors="coerce")
    df["read_count"] = pd.to_numeric(df["read_count"], errors="coerce").fillna(0)
    return df.dropna(subset=["donor", "acceptor"])


def load_sj_strand_map(sample, base_dir):
    sj_path = Path(base_dir) / f"{sample}SJ.out.tab"
    if not sj_path.exists():
        print(f"[WARN] sample={sample} 缺少 {sj_path.name}，无法做 strand 过滤，将跳过该 sample", file=sys.stderr)
        return None
    strand_map = {}
    with open(sj_path, "r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                continue
            strand_map[(fields[0], int(fields[1]), int(fields[2]))] = fields[3]
    return strand_map


def sj_strand_to_gene_strand(strand_code):
    if str(strand_code) == "1":
        return "+"
    if str(strand_code) == "2":
        return "-"
    return None


def build_gene_level_table(junction_dir, db, chrom_index):
    sample_to_condition = dict(zip(db["name"].astype(str), db["condition"]))
    files = sorted(
        os.path.join(junction_dir, name)
        for name in os.listdir(junction_dir)
        if name.endswith(".junctions.tsv")
    )
    if not files:
        raise FileNotFoundError(f"{junction_dir} 中没有找到 *.junctions.tsv")

    all_rows = []
    for path in files:
        sample = os.path.basename(path).replace(".junctions.tsv", "")
        if sample not in sample_to_condition:
            continue
        condition = sample_to_condition[sample]
        if condition not in {"CTRL", "IAA", "U1"}:
            continue
        strand_map = load_sj_strand_map(sample, Path(junction_dir).parent)
        if strand_map is None:
            continue
        df = load_junction_table(path)
        gene_counts = {}
        unmapped = 0
        wrong_strand = 0
        not_same_gene = 0

        for row in df.itertuples(index=False):
            strand_code = strand_map.get((row.chrom, int(row.donor) + 1, int(row.acceptor)))
            required_gene_strand = sj_strand_to_gene_strand(strand_code)
            if required_gene_strand is None:
                wrong_strand += 1
                continue
            donor_genes = genes_covering_site(chrom_index, row.chrom, int(row.donor), required_gene_strand)
            acceptor_genes = genes_covering_site(chrom_index, row.chrom, int(row.acceptor), required_gene_strand)
            same_genes = donor_genes & acceptor_genes
            if not donor_genes or not acceptor_genes:
                unmapped += 1
                continue
            if not same_genes:
                not_same_gene += 1
                continue
            is_novel = str(row.type).lower() != "annotated"
            for gene in same_genes:
                stat = gene_counts.setdefault(
                    gene,
                    {
                        "sample": sample,
                        "condition": condition,
                        "gene_name": gene,
                        "total_junctions": 0,
                        "novel_junctions": 0,
                        "total_reads": 0.0,
                        "novel_reads": 0.0,
                    },
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

        print(
            f"[INFO] sample={sample} condition={condition} mapped_genes={len(gene_counts)} "
            f"unmapped_junctions={unmapped} not_same_gene={not_same_gene} wrong_or_unknown_strand={wrong_strand}",
            file=sys.stderr,
        )

    if not all_rows:
        raise ValueError("没有生成任何 gene-level 统计结果。")
    return pd.concat(all_rows, ignore_index=True)


def filter_gene_table(data, min_total_junctions, min_total_reads):
    filtered = data.loc[
        (data["total_junctions"] >= min_total_junctions) &
        (data["total_reads"] >= min_total_reads)
    ].copy()
    print(f"[INFO] gene-level 记录: 原始 {len(data)} 行，过滤后 {len(filtered)} 行", file=sys.stderr)
    return filtered


def average_replicates_by_condition(data):
    value_cols = [
        "total_junctions", "novel_junctions", "total_reads", "novel_reads",
        "novel_junction_ratio", "novel_read_ratio",
    ]
    return (
        data.groupby(["condition", "gene_name"], as_index=False)[value_cols]
        .mean()
        .sort_values(["condition", "gene_name"])
    )


def map_transcripts_to_genes(transcript_set, transcript_to_gene, label):
    genes = {transcript_to_gene[t] for t in transcript_set if t in transcript_to_gene}
    missing = len(transcript_set) - sum(1 for t in transcript_set if t in transcript_to_gene)
    print(f"[INFO] {label}: {len(transcript_set)} transcripts -> {len(genes)} genes，未映射 transcript={missing}", file=sys.stderr)
    return genes


def subset_gene_set(data, genes, label):
    subset = data.loc[data["gene_name"].isin(genes)].copy()
    subset["gene_set"] = label
    print(f"[INFO] {label}: 用于绘图的 gene-level 记录 {len(subset)} 行", file=sys.stderr)
    return subset


def format_pvalue(p_value):
    text = f"{p_value:.1e}"
    return text.replace(".0e", "e").replace("e-0", "e-").replace("e+0", "e+")


def run_test(vals1, vals2, method):
    vals1 = np.asarray(vals1, dtype=float)
    vals2 = np.asarray(vals2, dtype=float)
    vals1 = vals1[np.isfinite(vals1)]
    vals2 = vals2[np.isfinite(vals2)]
    if len(vals1) < 2 or len(vals2) < 2:
        return None, None
    try:
        if method == "mannwhitney":
            stat, p = stats.mannwhitneyu(vals1, vals2, alternative="two-sided")
        else:
            stat, p = stats.ttest_ind(vals1, vals2, equal_var=False)
        return stat, p
    except Exception:
        return None, None


def metric_display_name(metric):
    if metric == "novel_junction_ratio":
        return "novel_junction_ratio (%)"
    if metric == "novel_read_ratio":
        return "junction_reads_ratio (%)"
    return metric


def describe_groups(data, metric, condition_order, title):
    print(f"\n[STAT] {title} | {metric}", file=sys.stderr)
    for condition in condition_order:
        vals = data.loc[data["condition"] == condition, metric].dropna().values
        if len(vals) == 0:
            print(f"  {condition:<8} N=0", file=sys.stderr)
            continue
        print(f"  {condition:<8} N={len(vals):<6} median={np.median(vals):.4f} mean={np.mean(vals):.4f}", file=sys.stderr)


def describe_ctrl_compare(data, metric, labels):
    print(f"\n[STAT] CTRL compare | {metric}", file=sys.stderr)
    for label in labels:
        vals = data.loc[data["gene_set"] == label, metric].dropna().values
        if len(vals) == 0:
            print(f"  {label:<18} N=0", file=sys.stderr)
            continue
        print(f"  {label:<18} N={len(vals):<6} median={np.median(vals):.4f} mean={np.mean(vals):.4f}", file=sys.stderr)


def plot_panel(ax, data, metric, condition_order, title, test_method, show_ylabel=True):
    plot_data = data.loc[data["condition"].isin(condition_order), ["condition", metric]].dropna().copy()
    sns.boxplot(
        data=plot_data,
        x="condition",
        y=metric,
        order=condition_order,
        palette=PALETTE,
        width=0.58,
        showfliers=False,
        linewidth=1.2,
        ax=ax,
    )
    for artist in ax.artists:
        artist.set_alpha(0.82)
    ax.set_title(title, fontsize=12, fontweight="700", pad=8)
    ax.set_ylabel(metric_display_name(metric) if show_ylabel else "", fontsize=10)
    ax.set_xlabel("")
    ax.set_xticklabels(condition_order, fontsize=10, fontweight="600")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(1.2)
    ax.spines["bottom"].set_linewidth(1.2)
    ax.grid(axis="y", linestyle="--", linewidth=0.5, alpha=0.5, zorder=0)
    ax.tick_params(axis="both", labelsize=9)

    pairs = [("CTRL", "IAA"), ("CTRL", "U1"), ("IAA", "U1")]
    valid_pairs = []
    annotations = []
    for left, right in pairs:
        vals_left = plot_data.loc[plot_data["condition"] == left, metric].dropna().values
        vals_right = plot_data.loc[plot_data["condition"] == right, metric].dropna().values
        _, p_value = run_test(vals_left, vals_right, test_method)
        if p_value is None or np.isnan(p_value):
            continue
        valid_pairs.append((left, right))
        annotations.append(f"p = {format_pvalue(p_value)}")
        print(f"    {title}: {left} vs {right}: p={format_pvalue(p_value)}", file=sys.stderr)

    if valid_pairs:
        annotator = Annotator(
            ax,
            valid_pairs,
            data=plot_data,
            x="condition",
            y=metric,
            order=condition_order,
        )
        annotator.configure(
            test=None,
            loc="outside",
            text_format="full",
            line_height=0.02,
            line_width=1.2,
            text_offset=2,
            fontsize=8,
            comparisons_correction=None,
            verbose=0,
        )
        annotator.set_custom_annotations(annotations)
        annotator.annotate()


def plot_ctrl_compare_panel(ax, data, metric, labels, test_method):
    plot_data = data.loc[data["gene_set"].isin(labels), ["gene_set", metric]].dropna().copy()
    sns.boxplot(
        data=plot_data,
        x="gene_set",
        y=metric,
        order=labels,
        palette=BED_COMPARE_PALETTE,
        width=0.58,
        showfliers=False,
        linewidth=1.2,
        ax=ax,
    )
    
    
    for artist in ax.artists:
        artist.set_alpha(0.82)
    ax.set_title("CTRL: special beds", fontsize=12, fontweight="700", pad=8)
    ax.set_ylabel("")
    ax.set_xlabel("")
    ax.set_xticklabels(labels, fontsize=10, fontweight="600", rotation=12, ha="right")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(1.2)
    ax.spines["bottom"].set_linewidth(1.2)
    ax.grid(axis="y", linestyle="--", linewidth=0.5, alpha=0.5, zorder=0)
    ax.tick_params(axis="both", labelsize=9)

    vals_left = plot_data.loc[plot_data["gene_set"] == labels[0], metric].dropna().values
    vals_right = plot_data.loc[plot_data["gene_set"] == labels[1], metric].dropna().values
    _, p_value = run_test(vals_left, vals_right, test_method)
    if p_value is None or np.isnan(p_value):
        return
    print(f"    CTRL compare: {labels[0]} vs {labels[1]}: p={format_pvalue(p_value)}", file=sys.stderr)
    annotator = Annotator(
        ax,
        [(labels[0], labels[1])],
        data=plot_data,
        x="gene_set",
        y=metric,
        order=labels,
    )
    annotator.configure(
        test=None,
        loc="outside",
        text_format="full",
        line_height=0.02,
        line_width=1.2,
        text_offset=2,
        fontsize=8,
        comparisons_correction=None,
        verbose=0,
    )
    annotator.set_custom_annotations([f"p = {format_pvalue(p_value)}"])
    annotator.annotate()


def make_figure(data_by_set, condition_order, output_path, test_method, dpi):
    metrics = ["novel_junction_ratio", "novel_read_ratio"]
    fig, axes = plt.subplots(2, 4, figsize=(21, 10.5), facecolor="#FAFAFA")
    fig.patch.set_facecolor("#FAFAFA")
    ctrl_compare_labels = ["with CGI", "no CGI"]
    ctrl_compare_df = pd.concat(
        [subset.loc[subset["condition"] == "CTRL"].copy() for title, subset in data_by_set if title in ctrl_compare_labels],
        ignore_index=True,
    )
    for row_idx, metric in enumerate(metrics):
        for col_idx, (title, subset) in enumerate(data_by_set):
            ax = axes[row_idx, col_idx]
            ax.set_facecolor("#FAFAFA")
            describe_groups(subset, metric, condition_order, title)
            panel_title = title if row_idx == 0 else ""
            plot_panel(
                ax,
                subset,
                metric,
                condition_order,
                panel_title,
                test_method,
                show_ylabel=(col_idx == 0),
            )
            if row_idx == 0:
                ax.set_title(title, fontsize=12, fontweight="700", pad=8)
        ctrl_ax = axes[row_idx, 3]
        ctrl_ax.set_facecolor("#FAFAFA")
        describe_ctrl_compare(ctrl_compare_df, metric, ctrl_compare_labels)
        plot_ctrl_compare_panel(ctrl_ax, ctrl_compare_df, metric, ctrl_compare_labels, test_method)
        axes[row_idx, 0].annotate(
            metric_display_name(metric),
            xy=(-0.38, 0.5),
            xycoords="axes fraction",
            rotation=90,
            va="center",
            ha="center",
            fontsize=11,
            fontweight="700",
        )

    fig.suptitle("Splicing Error by Condition Across Gene Sets", fontsize=14, fontweight="800", y=0.995)
    plt.tight_layout(rect=[0.04, 0, 1, 0.975])
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f"[INFO] 图片已保存: {output_path}", file=sys.stderr)


def main():
    args = parse_args()
    db = load_database(args.database)
    chrom_index, transcript_to_gene = load_gtf_resources(args.gtf)
    gene_sample_table = build_gene_level_table(args.junction_dir, db, chrom_index)
    gene_sample_table = filter_gene_table(gene_sample_table, args.min_total_junctions, args.min_total_reads)
    if gene_sample_table.empty:
        raise ValueError("过滤后没有可用于分析的 gene-level 数据。")

    gene_condition_table = average_replicates_by_condition(gene_sample_table)
    condition_order = [cond for cond in args.condition_order if cond in set(gene_condition_table["condition"].unique())]
    required = {"CTRL", "IAA", "U1"}
    if required - set(condition_order):
        raise ValueError(f"缺少 condition: {sorted(required - set(condition_order))}")

    active_genes = map_transcripts_to_genes(load_transcript_set(args.active_bed), transcript_to_gene, "Active genes")
    with_cgi_genes = map_transcripts_to_genes(load_transcript_set(args.with_cgi_bed), transcript_to_gene, "with CGI")
    no_cgi_genes = map_transcripts_to_genes(load_transcript_set(args.no_cgi_bed), transcript_to_gene, "no CGI")

    active_df = subset_gene_set(gene_condition_table, active_genes, "Active genes")
    with_cgi_df = subset_gene_set(gene_condition_table, with_cgi_genes, "with CGI")
    no_cgi_df = subset_gene_set(gene_condition_table, no_cgi_genes, "no CGI")

    combined = pd.concat([active_df, with_cgi_df, no_cgi_df], ignore_index=True)
    os.makedirs(os.path.dirname(args.table_output) or ".", exist_ok=True)
    combined.to_csv(args.table_output, sep="\t", index=False)
    print(f"[INFO] gene-level 表已保存: {args.table_output}", file=sys.stderr)

    data_by_set = [
        ("Active genes", active_df),
        ("with CGI", with_cgi_df),
        ("no CGI", no_cgi_df),
    ]
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    make_figure(data_by_set, condition_order, args.output, args.test, args.dpi)


if __name__ == "__main__":
    main()
