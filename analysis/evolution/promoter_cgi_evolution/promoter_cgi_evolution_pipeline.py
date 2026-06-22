#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Promoter CGI evolution analysis on a selected cross-species panel.

This pipeline:
1. loads downloaded genomes and annotations from `data/<species>/`
2. extracts protein-coding gene TSSs
3. evaluates whether the promoter around each TSS contains a TSS-overlapping CGI
4. evaluates matched random genomic positions as controls
5. summarizes species-level CGI fractions and merges protein-system scores
6. generates phylogenetic heatmaps and coevolution-style regplots
"""

from __future__ import annotations

import argparse
import gzip
import math
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Optional

import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["font.family"] = "Arial"

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LinearSegmentedColormap, Normalize
from scipy.stats import pearsonr, spearmanr

sys.path.append(str(Path(__file__).resolve().parent.parent))
from taxonomy_groups_refined import (  # noqa: E402
    GROUP_COLORS,
    GROUP_LABELS,
    GROUP_ORDER,
    clone_tree_template,
    normalize_species_meta,
)
TREE_TEMPLATE = clone_tree_template()

CGI_CMAP = LinearSegmentedColormap.from_list(
    "promoter_cgi_amber",
    ["#fff7ed", "#fed7aa", "#fb923c", "#ea580c", "#9a3412"],
)
COMPLEX_CMAP = LinearSegmentedColormap.from_list(
    "promoter_cgi_blue",
    ["#f8fbff", "#dbeafe", "#93c5fd", "#3b82f6", "#1d4ed8"],
)

HEATMAP_COLUMNS = [
    ("promoter_cgi_fraction", "Promoter"),
    ("random_ctrl_cgi_fraction", "Background"),
    ("cgi_delta", "Prom-Bg"),
    ("DNA_methylation", "DNA meth."),
    ("H3.3", "H3.3"),
    ("RNA_polymerase_control", "RNAP II"),
    ("Ribosome_control", "Ribosome"),
]

HEATMAP_HEADER_COLORS = {
    "promoter_cgi_fraction": "#c2410c",
    "random_ctrl_cgi_fraction": "#ea580c",
    "cgi_delta": "#9a3412",
    "DNA_methylation": "#1d4ed8",
    "H3.3": "#2563eb",
    "RNA_polymerase_control": "#1e40af",
    "Ribosome_control": "#1e3a8a",
}

HEATMAP_COLUMN_SPECS = {
    "promoter_cgi_fraction": {"cmap": CGI_CMAP, "norm": Normalize(vmin=0.0, vmax=1.0), "group": "CGI metrics"},
    "random_ctrl_cgi_fraction": {"cmap": CGI_CMAP, "norm": Normalize(vmin=0.0, vmax=1.0), "group": "CGI metrics"},
    "cgi_delta": {"cmap": CGI_CMAP, "norm": Normalize(vmin=0.0, vmax=1.0), "group": "CGI metrics"},
    "DNA_methylation": {"cmap": COMPLEX_CMAP, "norm": Normalize(vmin=0.0, vmax=1.0), "group": "Complex scores"},
    "H3.3": {"cmap": COMPLEX_CMAP, "norm": Normalize(vmin=0.0, vmax=1.0), "group": "Complex scores"},
    "RNA_polymerase_control": {"cmap": COMPLEX_CMAP, "norm": Normalize(vmin=0.0, vmax=1.0), "group": "Complex scores"},
    "Ribosome_control": {"cmap": COMPLEX_CMAP, "norm": Normalize(vmin=0.0, vmax=1.0), "group": "Complex scores"},
}


@dataclass
class TreeNode:
    label: str
    children: list["TreeNode"] = field(default_factory=list)
    species_name: Optional[str] = None
    taxonomy_class: Optional[str] = None
    depth: int = 0
    y: float = 0.0
    ymin: float = 0.0
    ymax: float = 0.0

    @property
    def is_leaf(self) -> bool:
        return self.species_name is not None


@dataclass
class RegplotSpec:
    pair_df: pd.DataFrame
    x_col: str
    y_col: str
    x_label: str
    y_label: str
    title: str
    n_label: str
    n_value: int
    pearson_r: float
    pearson_p: float
    spearman_rho: float
    spearman_p: float
    group_col: Optional[str] = None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run promoter CGI evolution analysis.")
    parser.add_argument(
        "--manifest",
        default="species_manifest_refined_50.tsv",
        help="Species manifest TSV.",
    )
    parser.add_argument("--data-dir", default="data", help="Downloaded genome/annotation directory.")
    parser.add_argument(
        "--output-dir",
        default="results_refined_59_species",
        help="Output directory.",
    )
    parser.add_argument(
        "--plot-only",
        action="store_true",
        help="Reuse an existing species summary TSV and regenerate plot-related outputs only.",
    )
    parser.add_argument(
        "--species-summary-input",
        default="species_promoter_cgi_summary.tsv",
        help="Existing species summary TSV filename under output-dir when --plot-only is enabled.",
    )
    parser.add_argument(
        "--species",
        nargs="*",
        default=[],
        help="Optional subset of species names to analyze.",
    )
    parser.add_argument(
        "--dna-summary",
        default="../results/refined_taxonomy/completeness_summary_coevolution.tsv",
        help="DNA methylation protein summary TSV from the protein pipeline.",
    )
    parser.add_argument("--promoter-upstream", type=int, default=1000, help="Upstream promoter span.")
    parser.add_argument("--promoter-downstream", type=int, default=500, help="Downstream promoter span.")
    parser.add_argument("--cgi-window", type=int, default=200, help="Fixed window size for CGI detection.")
    parser.add_argument("--min-gc", type=float, default=0.50, help="Minimum GC fraction for CGI.")
    parser.add_argument("--min-cpg-oe", type=float, default=0.60, help="Minimum CpG Obs/Exp for CGI.")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for control selection.")
    parser.add_argument(
        "--max-genes",
        type=int,
        default=0,
        help="Optional cap on the number of genes per species for quick testing.",
    )
    return parser.parse_args()


def parse_gtf_attributes(attr_text: str) -> dict[str, str]:
    attrs = {}
    for part in attr_text.strip().split(";"):
        part = part.strip()
        if not part or " " not in part:
            continue
        key, value = part.split(" ", 1)
        attrs[key] = value.strip().strip('"')
    return attrs


def parse_gff3_attributes(attr_text: str) -> dict[str, str]:
    attrs = {}
    for part in attr_text.strip().split(";"):
        if "=" not in part:
            continue
        key, value = part.split("=", 1)
        attrs[key] = value
    return attrs


def is_protein_coding(attrs: dict[str, str]) -> bool:
    for key in ["gene_biotype", "gene_type", "biotype"]:
        if key in attrs:
            return attrs[key].lower().replace("-", "_") == "protein_coding"
    return True


def open_text_maybe_gzip(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="ignore")
    return path.open("r", encoding="utf-8", errors="ignore")


def load_gene_annotations(annotation_path: Path) -> pd.DataFrame:
    rows = []
    is_gtf = annotation_path.name.endswith(".gtf") or annotation_path.name.endswith(".gtf.gz")

    with open_text_maybe_gzip(annotation_path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue

            attrs = parse_gtf_attributes(fields[8]) if is_gtf else parse_gff3_attributes(fields[8])
            if not is_protein_coding(attrs):
                continue

            chrom = fields[0]
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]

            gene_id = (
                attrs.get("gene_id")
                or attrs.get("ID")
                or attrs.get("gene")
                or attrs.get("Name")
                or f"{chrom}:{start}-{end}"
            )
            gene_name = attrs.get("gene_name") or attrs.get("Name") or gene_id

            rows.append(
                {
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "gene_id": gene_id,
                    "gene_name": gene_name,
                }
            )

    if not rows:
        raise ValueError(f"No protein-coding gene records were parsed from {annotation_path}")

    gene_df = pd.DataFrame(rows).drop_duplicates(subset=["gene_id"]).copy()
    return gene_df


def reverse_complement(seq: str) -> str:
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1]


def fetch_promoter_sequence(
    fasta: pysam.FastaFile,
    chrom: str,
    tss: int,
    strand: str,
    chrom_len: int,
    upstream: int,
    downstream: int,
) -> tuple[str, int, int, int]:
    if strand == "+":
        start = max(1, tss - upstream)
        end = min(chrom_len, tss + downstream)
        seq = fasta.fetch(chrom, start - 1, end).upper()
        tss_idx = tss - start
    else:
        start = max(1, tss - downstream)
        end = min(chrom_len, tss + upstream)
        seq = reverse_complement(fasta.fetch(chrom, start - 1, end).upper())
        tss_idx = end - tss
    return seq, tss_idx, start, end


def gc_fraction(seq: str) -> float:
    gc = seq.count("G") + seq.count("C")
    return gc / len(seq) if seq else 0.0


def cpg_observed_expected(seq: str) -> float:
    c = seq.count("C")
    g = seq.count("G")
    if c == 0 or g == 0:
        return 0.0
    cg = sum(1 for i in range(len(seq) - 1) if seq[i : i + 2] == "CG")
    return (cg * len(seq)) / (c * g)


def scan_tss_overlapping_cgi(
    seq: str,
    tss_idx: int,
    cgi_window: int,
    min_gc: float,
    min_cpg_oe: float,
) -> dict[str, float]:
    if len(seq) < cgi_window:
        return {
            "cgi_present": 0.0,
            "best_gc": 0.0,
            "best_cpg_oe": 0.0,
            "best_score": 0.0,
        }

    start_min = max(0, tss_idx - cgi_window + 1)
    start_max = min(tss_idx, len(seq) - cgi_window)

    best_gc = 0.0
    best_oe = 0.0
    best_score = 0.0
    cgi_present = 0.0

    for start in range(start_min, start_max + 1):
        window = seq[start : start + cgi_window]
        if "N" in window:
            continue
        gc = gc_fraction(window)
        oe = cpg_observed_expected(window)
        score = min(gc / min_gc if min_gc > 0 else 0.0, oe / min_cpg_oe if min_cpg_oe > 0 else 0.0)
        if score > best_score:
            best_gc = gc
            best_oe = oe
            best_score = score
        if gc >= min_gc and oe >= min_cpg_oe:
            cgi_present = 1.0

    return {
        "cgi_present": cgi_present,
        "best_gc": best_gc,
        "best_cpg_oe": best_oe,
        "best_score": best_score,
    }


def sample_random_control(
    fasta: pysam.FastaFile,
    chrom: str,
    strand: str,
    chrom_len: int,
    upstream: int,
    downstream: int,
    cgi_window: int,
    min_gc: float,
    min_cpg_oe: float,
    rng: np.random.Generator,
    max_attempts: int = 50,
) -> dict[str, float]:
    margin = max(upstream, downstream) + 2
    if chrom_len <= margin * 2:
        return {
            "cgi_present": 0.0,
            "best_gc": 0.0,
            "best_cpg_oe": 0.0,
            "best_score": 0.0,
            "center": math.nan,
        }

    for _ in range(max_attempts):
        center = int(rng.integers(margin, chrom_len - margin))
        seq, tss_idx, _, _ = fetch_promoter_sequence(
            fasta, chrom, center, strand, chrom_len, upstream, downstream
        )
        if len(seq) < cgi_window or seq.count("N") / max(len(seq), 1) > 0.20:
            continue
        result = scan_tss_overlapping_cgi(seq, tss_idx, cgi_window, min_gc, min_cpg_oe)
        result["center"] = center
        return result

    return {
        "cgi_present": 0.0,
        "best_gc": 0.0,
        "best_cpg_oe": 0.0,
        "best_score": 0.0,
        "center": math.nan,
    }


def load_protein_system_scores(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame(columns=["DNA_methylation", "H3.3", "H2A.Z", "RNA_polymerase_control", "Ribosome_control"])
    df = pd.read_csv(path, sep="\t", index_col=0)
    keep_cols = [col for col in ["DNA_methylation", "H3.3", "H2A.Z", "RNA_polymerase_control", "Ribosome_control"] if col in df.columns]
    return df[keep_cols].copy()


def analyze_species(
    species_row: pd.Series,
    data_dir: Path,
    gene_detail_dir: Path,
    upstream: int,
    downstream: int,
    cgi_window: int,
    min_gc: float,
    min_cpg_oe: float,
    seed: int,
    max_genes: int,
) -> dict[str, float]:
    species = str(species_row["species"])
    species_dir = data_dir / species
    fasta_path = species_dir / f"{species}.fa"
    raw_dir = species_dir / "raw"

    if not fasta_path.exists():
        raise FileNotFoundError(f"Missing FASTA for {species}: {fasta_path}")

    annotation_candidates = [
        p
        for p in (sorted(raw_dir.glob("*.gtf.gz")) + sorted(raw_dir.glob("*.gff3.gz")))
        if not p.name.startswith("._")
    ]
    if not annotation_candidates:
        annotation_candidates = [
            p
            for p in (sorted(raw_dir.glob("*.gtf")) + sorted(raw_dir.glob("*.gff3")))
            if not p.name.startswith("._")
        ]
    if not annotation_candidates:
        raise FileNotFoundError(f"Missing annotation file under {raw_dir}")
    annotation_path = annotation_candidates[0]

    genes = load_gene_annotations(annotation_path)
    if max_genes > 0:
        genes = genes.head(max_genes).copy()

    fasta = pysam.FastaFile(str(fasta_path))
    chrom_sizes = dict(zip(fasta.references, fasta.lengths))
    genes = genes[genes["chrom"].isin(chrom_sizes)].copy()
    if genes.empty:
        raise ValueError(f"No annotation chromosomes matched FASTA references for {species}")

    rng = np.random.default_rng(seed + sum(ord(ch) for ch in species))
    detail_rows = []
    for row in genes.itertuples(index=False):
        chrom_len = chrom_sizes[row.chrom]
        tss = row.start if row.strand == "+" else row.end
        seq, tss_idx, promoter_start, promoter_end = fetch_promoter_sequence(
            fasta, row.chrom, tss, row.strand, chrom_len, upstream, downstream
        )
        promoter_eval = scan_tss_overlapping_cgi(seq, tss_idx, cgi_window, min_gc, min_cpg_oe)
        ctrl_eval = sample_random_control(
            fasta,
            row.chrom,
            row.strand,
            chrom_len,
            upstream,
            downstream,
            cgi_window,
            min_gc,
            min_cpg_oe,
            rng,
        )
        detail_rows.append(
            {
                "species": species,
                "gene_id": row.gene_id,
                "gene_name": row.gene_name,
                "chrom": row.chrom,
                "strand": row.strand,
                "tss": tss,
                "promoter_start": promoter_start,
                "promoter_end": promoter_end,
                "promoter_cgi": promoter_eval["cgi_present"],
                "promoter_best_gc": promoter_eval["best_gc"],
                "promoter_best_cpg_oe": promoter_eval["best_cpg_oe"],
                "promoter_best_score": promoter_eval["best_score"],
                "ctrl_center": ctrl_eval["center"],
                "ctrl_cgi": ctrl_eval["cgi_present"],
                "ctrl_best_gc": ctrl_eval["best_gc"],
                "ctrl_best_cpg_oe": ctrl_eval["best_cpg_oe"],
                "ctrl_best_score": ctrl_eval["best_score"],
            }
        )

    detail_df = pd.DataFrame(detail_rows)
    gene_detail_dir.mkdir(parents=True, exist_ok=True)
    detail_df.to_csv(gene_detail_dir / f"{species}_promoter_cgi_details.tsv.gz", sep="\t", index=False)

    promoter_fraction = float(detail_df["promoter_cgi"].mean())
    ctrl_fraction = float(detail_df["ctrl_cgi"].mean())
    cgi_delta = promoter_fraction - ctrl_fraction
    cgi_log2_enrichment = float(np.log2((promoter_fraction + 1e-6) / (ctrl_fraction + 1e-6)))

    return {
        "species": species,
        "display_name": species_row["display_name"],
        "taxonomy_class": species_row["taxonomy_class"],
        "n_genes": len(detail_df),
        "promoter_cgi_fraction": promoter_fraction,
        "random_ctrl_cgi_fraction": ctrl_fraction,
        "cgi_delta": cgi_delta,
        "cgi_log2_enrichment": cgi_log2_enrichment,
        "mean_promoter_best_gc": float(detail_df["promoter_best_gc"].mean()),
        "mean_promoter_best_cpg_oe": float(detail_df["promoter_best_cpg_oe"].mean()),
    }


def build_group_summary(summary_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for group in GROUP_ORDER:
        subset = summary_df[summary_df["taxonomy_class"] == group].copy()
        if subset.empty:
            continue
        numeric = subset.select_dtypes(include=[np.number]).mean(axis=0)
        row = numeric.to_dict()
        row["species"] = group
        row["display_name"] = GROUP_LABELS.get(group, group.replace("_", " "))
        row["taxonomy_class"] = group
        row["n_species"] = len(subset)
        rows.append(row)
    return pd.DataFrame(rows)


def sort_species_summary(summary_df: pd.DataFrame) -> pd.DataFrame:
    df = normalize_species_meta(summary_df).copy()
    df["group_rank"] = df["taxonomy_class"].map({group: idx for idx, group in enumerate(GROUP_ORDER)})
    df["display_name"] = df["display_name"].fillna(df["species"])
    sort_cols = ["group_rank", "DNA_methylation", "promoter_cgi_fraction", "display_name"]
    ascending = [True, True, False, True]
    for col in ["DNA_methylation", "promoter_cgi_fraction"]:
        if col not in df.columns:
            df[col] = np.nan
    df = df.sort_values(sort_cols, ascending=ascending).reset_index(drop=True)
    return df


def get_cgi_metric_specs() -> list[tuple[str, str]]:
    return [
        ("promoter_cgi_fraction", "Promoter CGI fraction"),
        ("random_ctrl_cgi_fraction", "Matched random CGI fraction"),
        ("cgi_delta", "Promoter CGI fraction minus random control"),
        ("cgi_log2_enrichment", "Promoter CGI log2 enrichment"),
    ]


def get_protein_target_specs() -> list[tuple[str, str]]:
    return [
        ("DNA_methylation", "DNA methylation protein score"),
        ("H3.3", "H3.3 protein score"),
        ("H2A.Z", "H2A.Z protein score"),
        ("RNA_polymerase_control", "RNAPII protein score"),
        ("Ribosome_control", "Ribosome protein score"),
    ]


def build_regplot_specs(
    summary_df: pd.DataFrame,
    *,
    title_prefix: str = "",
    min_points: int = 3,
    group_col: Optional[str] = None,
    pair_label: str = "n",
) -> list[RegplotSpec]:
    specs: list[RegplotSpec] = []
    for metric, label in get_cgi_metric_specs():
        for target_col, target_label in get_protein_target_specs():
            if target_col not in summary_df.columns:
                continue
            cols = [metric, target_col] + ([group_col] if group_col else [])
            pair = summary_df[cols].dropna().copy()
            if len(pair) < min_points:
                continue
            pearson_r, pearson_p = pearsonr(pair[metric], pair[target_col])
            spearman_rho, spearman_p = spearmanr(pair[metric], pair[target_col])
            specs.append(
                RegplotSpec(
                    pair_df=pair,
                    x_col=metric,
                    y_col=target_col,
                    x_label=label,
                    y_label=target_label,
                    title=f"{title_prefix}{label} vs {target_label}",
                    n_label=pair_label,
                    n_value=len(pair),
                    pearson_r=pearson_r,
                    pearson_p=pearson_p,
                    spearman_rho=spearman_rho,
                    spearman_p=spearman_p,
                    group_col=group_col,
                )
            )

    if {"H3.3", "H2A.Z"}.issubset(summary_df.columns):
        cols = ["H3.3", "H2A.Z"] + ([group_col] if group_col else [])
        pair = summary_df[cols].dropna().copy()
        if len(pair) >= min_points:
            pearson_r, pearson_p = pearsonr(pair["H3.3"], pair["H2A.Z"])
            spearman_rho, spearman_p = spearmanr(pair["H3.3"], pair["H2A.Z"])
            specs.append(
                RegplotSpec(
                    pair_df=pair,
                    x_col="H3.3",
                    y_col="H2A.Z",
                    x_label="H3.3 protein score",
                    y_label="H2A.Z protein score",
                    title=f"{title_prefix}H3.3 vs H2A.Z",
                    n_label=pair_label,
                    n_value=len(pair),
                    pearson_r=pearson_r,
                    pearson_p=pearson_p,
                    spearman_rho=spearman_rho,
                    spearman_p=spearman_p,
                    group_col=group_col,
                )
            )
    return specs


def build_correlation_rows(summary_df: pd.DataFrame) -> list[dict[str, float | str | int]]:
    rows = []
    for spec in build_regplot_specs(summary_df, title_prefix="", min_points=3):
        if spec.x_col == "H3.3" and spec.y_col == "H2A.Z":
            continue
        rows.append(
            {
                "cgi_metric": spec.x_col,
                "protein_system": spec.y_col,
                "pearson_r": spec.pearson_r,
                "pearson_p": spec.pearson_p,
                "spearman_rho": spec.spearman_rho,
                "spearman_p": spec.spearman_p,
                "n_species": spec.n_value,
            }
        )
    return rows


def make_correlation_table(summary_df: pd.DataFrame) -> pd.DataFrame:
    rows = build_correlation_rows(summary_df)
    corr_df = pd.DataFrame(rows)
    if corr_df.empty:
        return pd.DataFrame(
            columns=["cgi_metric", "protein_system", "pearson_r", "pearson_p", "spearman_rho", "spearman_p", "n_species"]
        )
    return corr_df.sort_values(["cgi_metric", "pearson_r"], ascending=[True, False])


def plot_regplots(specs: list[RegplotSpec], output_path: Path) -> None:
    if not specs:
        return
    with PdfPages(output_path) as pdf:
        for spec in specs:
            fig, ax = plt.subplots(figsize=(5, 5))
            sns.regplot(
                data=spec.pair_df,
                x=spec.x_col,
                y=spec.y_col,
                scatter_kws={"s": 40, "alpha": 0.8, "rasterized": True},
                line_kws={"color": "#c23b22", "linewidth": 2},
                ax=ax,
            )
            ax.set_title(spec.title)
            ax.set_xlabel(spec.x_label)
            ax.set_ylabel(spec.y_label)
            ax.text(
                0.03,
                0.97,
                (
                    f"{spec.n_label} = {spec.n_value}\n"
                    f"Pearson r = {spec.pearson_r:.2f}\n"
                    f"Pearson p = {spec.pearson_p:.2e}\n"
                    f"Spearman rho = {spec.spearman_rho:.2f}\n"
                    f"Spearman p = {spec.spearman_p:.2e}"
                ),
                transform=ax.transAxes,
                va="top",
                ha="left",
                fontsize=9,
                bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "alpha": 0.85},
            )
            fig.tight_layout()
            pdf.savefig(fig, dpi=300, bbox_inches="tight")
            plt.close(fig)

def plot_group_regplots(specs: list[RegplotSpec], output_path: Path, label_lookup: dict[str, str]) -> None:
    if not specs:
        return
    with PdfPages(output_path) as pdf:
        for spec in specs:
            fig, ax = plt.subplots(figsize=(5.4, 5.4))
            sns.regplot(
                data=spec.pair_df,
                x=spec.x_col,
                y=spec.y_col,
                scatter=False,
                line_kws={"color": "#c23b22", "linewidth": 2},
                ax=ax,
            )
            ax.scatter(
                spec.pair_df[spec.x_col],
                spec.pair_df[spec.y_col],
                s=58,
                c=[GROUP_COLORS.get(group, "#475569") for group in spec.pair_df[spec.group_col]],
                alpha=0.95,
                edgecolors="white",
                linewidths=0.7,
                zorder=3,
            )
            for _, row in spec.pair_df.iterrows():
                group_name = row[spec.group_col]
                ax.text(
                    row[spec.x_col] + 0.01,
                    row[spec.y_col],
                    label_lookup.get(group_name, group_name).replace("_", " "),
                    fontsize=7.8,
                    color=GROUP_COLORS.get(group_name, "#334155"),
                    va="center",
                    ha="left",
                )
            ax.set_title(spec.title)
            ax.set_xlabel(spec.x_label)
            ax.set_ylabel(spec.y_label)
            ax.text(
                0.03,
                0.97,
                (
                    f"{spec.n_label} = {spec.n_value}\n"
                    f"Pearson r = {spec.pearson_r:.2f}\n"
                    f"Pearson p = {spec.pearson_p:.2e}\n"
                    f"Spearman rho = {spec.spearman_rho:.2f}\n"
                    f"Spearman p = {spec.spearman_p:.2e}"
                ),
                transform=ax.transAxes,
                va="top",
                ha="left",
                fontsize=9,
                bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "alpha": 0.85},
            )
            fig.tight_layout()
            pdf.savefig(fig, dpi=300, bbox_inches="tight")
            plt.close(fig)


def make_regplots(summary_df: pd.DataFrame, output_path: Path) -> None:
    specs = build_regplot_specs(summary_df, min_points=3, pair_label="n")
    plot_regplots(specs, output_path)


def make_group_regplots(group_summary_df: pd.DataFrame, output_path: Path) -> None:
    specs = build_regplot_specs(
        group_summary_df,
        title_prefix="Broad groups: ",
        min_points=2,
        group_col="species",
        pair_label="n groups",
    )
    label_lookup = dict(zip(group_summary_df["species"], group_summary_df["display_name"]))
    plot_group_regplots(specs, output_path, label_lookup)


def render_outputs(species_summary: pd.DataFrame, output_dir: Path) -> None:
    species_summary = sort_species_summary(species_summary)
    species_summary.to_csv(output_dir / "species_promoter_cgi_summary.tsv", sep="\t", index=False)
    group_summary = build_group_summary(species_summary)
    group_summary.to_csv(output_dir / "broad_group_promoter_cgi_summary.tsv", sep="\t", index=False)

    corr_df = make_correlation_table(species_summary)
    corr_df.to_csv(output_dir / "promoter_cgi_vs_dna_methylation_correlations.tsv", sep="\t", index=False)
    make_regplots(species_summary, output_dir / "promoter_cgi_vs_dna_methylation_regplots.pdf")
    make_group_regplots(group_summary, output_dir / "promoter_cgi_vs_dna_methylation_broad_group_regplots.pdf")

    species_meta = species_summary[["species", "display_name", "taxonomy_class"]].copy()
    plot_phylogenetic_heatmap(
        summary_df=species_summary,
        row_meta=species_meta,
        output_path=output_dir / "species_promoter_cgi_phylogenetic_heatmap.pdf",
        title="Species-Level Phylogenetic Heatmap of Promoter CGI and Complex Evolution",
        subtitle="Protein-coding gene promoters are summarized for all species, with promoter CGI metrics and matched protein-system scores aligned on the same taxonomic scaffold.",
    )

    group_meta = group_summary[["species", "display_name", "taxonomy_class"]].copy()
    plot_phylogenetic_heatmap(
        summary_df=group_summary,
        row_meta=group_meta,
        output_path=output_dir / "broad_group_promoter_cgi_phylogenetic_heatmap.pdf",
        title="Broad-Group Phylogenetic Heatmap of Promoter CGI and Complex Evolution",
        subtitle="Broad-group means are shown for Promoter, Background, Prom-Bg, DNA meth., H3.3, RNAP II and Ribosome, with white-orange for CGI-related columns and white-blue for complex scores.",
    )

    print(f"Broad-group summary saved to: {output_dir / 'broad_group_promoter_cgi_summary.tsv'}")
    print(f"Correlation table saved to: {output_dir / 'promoter_cgi_vs_dna_methylation_correlations.tsv'}")
    print(f"Regplots saved to: {output_dir / 'promoter_cgi_vs_dna_methylation_regplots.pdf'}")
    print(f"Broad-group regplots saved to: {output_dir / 'promoter_cgi_vs_dna_methylation_broad_group_regplots.pdf'}")
    print(f"Species heatmap saved to: {output_dir / 'species_promoter_cgi_phylogenetic_heatmap.pdf'}")
    print(f"Broad-group heatmap saved to: {output_dir / 'broad_group_promoter_cgi_phylogenetic_heatmap.pdf'}")


def build_tree(template: dict, row_meta: pd.DataFrame) -> Optional[TreeNode]:
    group_key = template.get("group")
    if group_key:
        group_df = row_meta[row_meta["taxonomy_class"] == group_key].copy()
        if group_df.empty:
            return None
        if len(group_df) == 1:
            row = group_df.iloc[0]
            return TreeNode(
                label=str(row["display_name"]),
                species_name=str(row["species"]),
                taxonomy_class=group_key,
            )
        leaves = [
            TreeNode(
                label=str(row["display_name"]),
                species_name=str(row["species"]),
                taxonomy_class=group_key,
            )
            for _, row in group_df.iterrows()
        ]
        return TreeNode(label=str(template["label"]), children=leaves, taxonomy_class=group_key)

    children = []
    for child_template in template.get("children", []):
        child = build_tree(child_template, row_meta)
        if child is not None:
            children.append(child)
    if not children:
        return None
    return TreeNode(label=str(template["label"]), children=children)


def assign_layout(node: TreeNode, species_to_y: dict[str, float], depth: int = 0) -> None:
    node.depth = depth
    if node.is_leaf:
        node.y = species_to_y[node.species_name]
        node.ymin = node.y
        node.ymax = node.y
        return
    for child in node.children:
        assign_layout(child, species_to_y, depth + 1)
    node.y = float(np.mean([child.y for child in node.children]))
    node.ymin = min(child.ymin for child in node.children)
    node.ymax = max(child.ymax for child in node.children)


def iter_nodes(node: TreeNode) -> Iterable[TreeNode]:
    yield node
    for child in node.children:
        yield from iter_nodes(child)


def iter_leaves(node: TreeNode) -> Iterable[TreeNode]:
    if node.is_leaf:
        yield node
        return
    for child in node.children:
        yield from iter_leaves(child)


def leaf_count(node: TreeNode) -> int:
    return sum(1 for _ in iter_leaves(node))


def max_internal_depth(node: TreeNode) -> int:
    depths = [n.depth for n in iter_nodes(node) if not n.is_leaf]
    return max(depths) if depths else 0


def draw_tree(ax: plt.Axes, node: TreeNode, leaf_x: float) -> None:
    if node.is_leaf:
        return
    child_ys = [child.y for child in node.children]
    ax.plot([node.depth, node.depth], [min(child_ys), max(child_ys)], color="#334155", lw=1.15, zorder=3)
    for child in node.children:
        child_x = leaf_x if child.is_leaf else child.depth
        ax.plot([node.depth, child_x], [child.y, child.y], color="#334155", lw=1.15, zorder=3)
        if not child.is_leaf:
            draw_tree(ax, child, leaf_x)


def annotate_internal_nodes(ax: plt.Axes, node: TreeNode) -> None:
    if node.is_leaf:
        return
    if node.depth > 0:
        label = node.label
        if node.taxonomy_class:
            label = f"{label} (n={leaf_count(node)})"
        ax.text(
            node.depth + 0.06,
            node.y,
            label,
            fontsize=7.0,
            color="#475569",
            va="center",
            ha="left",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.65, "pad": 0.15},
        )
    for child in node.children:
        annotate_internal_nodes(ax, child)


def add_group_backgrounds(ax: plt.Axes, row_meta: pd.DataFrame) -> None:
    for group in GROUP_ORDER:
        group_df = row_meta[row_meta["taxonomy_class"] == group]
        if group_df.empty:
            continue
        y0 = group_df["y"].min() - 0.5
        y1 = group_df["y"].max() + 0.5
        ax.axhspan(y0, y1, color=GROUP_COLORS[group], alpha=0.08, lw=0, zorder=0)


def plot_phylogenetic_heatmap(
    summary_df: pd.DataFrame,
    row_meta: pd.DataFrame,
    output_path: Path,
    title: str,
    subtitle: str,
) -> None:
    row_meta = row_meta.copy()
    row_meta["y"] = np.arange(len(row_meta), dtype=float)
    species_to_y = dict(zip(row_meta["species"], row_meta["y"]))
    root = build_tree(TREE_TEMPLATE, row_meta)
    if root is None:
        raise ValueError("Tree construction failed for the selected rows.")
    assign_layout(root, species_to_y)

    plot_cols = [col for col, _ in HEATMAP_COLUMNS]
    plot_labels = [label for _, label in HEATMAP_COLUMNS]

    fig_height = max(6.2, len(row_meta) * 0.45 + 1.8)
    fig = plt.figure(figsize=(13.5, fig_height))
    gs = fig.add_gridspec(nrows=1, ncols=2, width_ratios=[4.5, 1.8], wspace=0.03)
    ax_tree = fig.add_subplot(gs[0, 0])
    ax_heat = fig.add_subplot(gs[0, 1], sharey=ax_tree)

    add_group_backgrounds(ax_tree, row_meta)
    add_group_backgrounds(ax_heat, row_meta)

    leaf_x = max_internal_depth(root) + 0.9
    draw_tree(ax_tree, root, leaf_x)
    annotate_internal_nodes(ax_tree, root)
    for leaf in iter_leaves(root):
        ax_tree.text(
            leaf_x + 0.08,
            leaf.y,
            leaf.label,
            fontsize=8.5 if len(row_meta) > 10 else 10.0,
            fontweight="bold" if len(row_meta) <= 10 else "normal",
            color="#0f172a",
            va="center",
            ha="left",
        )

    plot_df = summary_df.set_index("species").loc[row_meta["species"], plot_cols].copy()

    values = plot_df.to_numpy(dtype=float)
    rgba = np.zeros((values.shape[0], values.shape[1], 4), dtype=float)
    for col_idx, column in enumerate(plot_cols):
        spec = HEATMAP_COLUMN_SPECS[column]
        rgba[:, col_idx, :] = spec["cmap"](spec["norm"](values[:, col_idx]))
    ax_heat.imshow(
        rgba,
        aspect="auto",
        interpolation="nearest",
        origin="upper",
    )
    for xpos in np.arange(-0.5, len(plot_cols), 1.0):
        ax_heat.axvline(xpos, color="white", lw=1.0)
    ax_heat.axvline(2.5, color="#cbd5e1", lw=1.5)
    ax_heat.axvline(3.5, color="#cbd5e1", lw=1.5)

    ax_heat.set_xticks(np.arange(len(plot_cols)))
    ax_heat.set_xticklabels(plot_labels, rotation=35, ha="left", fontsize=9.5)
    ax_heat.xaxis.tick_top()
    ax_heat.tick_params(axis="x", pad=8, length=0)
    for tick, column in zip(ax_heat.get_xticklabels(), plot_cols):
        tick.set_color(HEATMAP_HEADER_COLORS[column])
        tick.set_fontweight("bold")

    ax_heat.set_yticks([])
    ax_tree.set_xticks([])
    ax_tree.set_yticks([])
    ax_tree.set_xlim(-0.3, leaf_x + 2.8)
    ax_tree.set_ylim(len(row_meta) - 0.5, -0.5)
    ax_heat.set_ylim(len(row_meta) - 0.5, -0.5)
    for ax in (ax_tree, ax_heat):
        for spine in ax.spines.values():
            spine.set_visible(False)

    fig.suptitle(title, fontsize=15, fontweight="bold", y=0.995)
    fig.text(0.5, 0.982, subtitle, ha="center", va="top", fontsize=9.5, color="#475569")
    cbar_x = 0.78
    cbar_width = 0.16
    cbar_height = 0.018
    cbar_specs = [
        ("Promoter / Background / Prom-Bg", CGI_CMAP, Normalize(vmin=0.0, vmax=1.0), 0.095),
        ("DNA meth. / H3.3 / RNAP II / Ribosome", COMPLEX_CMAP, Normalize(vmin=0.0, vmax=1.0), 0.050),
    ]
    for label, cmap, norm, y_pos in cbar_specs:
        cax = fig.add_axes([cbar_x, y_pos, cbar_width, cbar_height])
        mappable = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
        cbar = fig.colorbar(mappable, cax=cax, orientation="horizontal")
        cbar.ax.tick_params(labelsize=7.0, length=2)
        cbar.set_label(label, fontsize=7.8, labelpad=2)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    manifest_path = Path(args.manifest).resolve()
    base_dir = manifest_path.parent
    data_dir = (base_dir / args.data_dir).resolve()
    output_dir = (base_dir / args.output_dir).resolve()
    gene_detail_dir = output_dir / "gene_level_tables"
    output_dir.mkdir(parents=True, exist_ok=True)

    if args.plot_only:
        species_summary_path = (output_dir / args.species_summary_input).resolve()
        species_summary = pd.read_csv(species_summary_path, sep="\t")
        species_summary = normalize_species_meta(species_summary)
        render_outputs(species_summary, output_dir)
        print(f"Species summary reused from: {species_summary_path}")
        return

    manifest = pd.read_csv(manifest_path, sep="\t")
    manifest = normalize_species_meta(manifest)
    if args.species:
        manifest = manifest[manifest["species"].isin(set(args.species))].copy()
    protein_scores = load_protein_system_scores((base_dir / args.dna_summary).resolve())

    summaries = []
    for _, row in manifest.iterrows():
        print(f"Analyzing {row['species']}...")
        result = analyze_species(
            species_row=row,
            data_dir=data_dir,
            gene_detail_dir=gene_detail_dir,
            upstream=args.promoter_upstream,
            downstream=args.promoter_downstream,
            cgi_window=args.cgi_window,
            min_gc=args.min_gc,
            min_cpg_oe=args.min_cpg_oe,
            seed=args.seed,
            max_genes=args.max_genes,
        )
        if row["species"] in protein_scores.index:
            for col in protein_scores.columns:
                result[col] = float(protein_scores.loc[row["species"], col])
        else:
            for col in protein_scores.columns:
                result[col] = np.nan
        summaries.append(result)

    species_summary = pd.DataFrame(summaries)
    render_outputs(species_summary, output_dir)

    print(f"Species summary saved to: {output_dir / 'species_promoter_cgi_summary.tsv'}")


if __name__ == "__main__":
    main()
