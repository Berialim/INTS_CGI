#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Plot TSS CGI coverage proportions for four representative species using
previously summarized promoter CGI results.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["font.family"] = "Arial"

import matplotlib.pyplot as plt
import pandas as pd


SPECIES_ORDER = [
    ("saccharomyces_cerevisiae", "sc", "#7c3aed"),
    ("drosophila_melanogaster", "dm", "#f59e0b"),
    ("mus_musculus", "mm", "#2563eb"),
    ("homo_sapiens", "human", "#dc2626"),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot TSS CGI coverage proportions for sc, dm, mm and human."
    )
    parser.add_argument(
        "--summary",
        default="promoter_cgi_evolution/results/species_promoter_cgi_summary.tsv",
        help="Previously summarized species-level promoter CGI table.",
    )
    parser.add_argument(
        "--output",
        default="promoter_cgi_evolution/tss_cgi_coverage_sc_dm_mm_human_barplot.pdf",
        help="Output barplot PDF.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    summary = pd.read_csv(Path(args.summary), sep="\t")

    rows = []
    for species, short_label, color in SPECIES_ORDER:
        hit = summary.loc[summary["species"] == species]
        if hit.empty:
            raise ValueError(f"Species not found in summary: {species}")
        rows.append(
            {
                "species": species,
                "label": short_label,
                "display_name": str(hit["display_name"].iloc[0]),
                "promoter_cgi_fraction": float(hit["promoter_cgi_fraction"].iloc[0]),
                "color": color,
            }
        )

    plot_df = pd.DataFrame(rows)

    fig, ax = plt.subplots(figsize=(5.0, 4.6))
    bars = ax.bar(
        plot_df["label"],
        plot_df["promoter_cgi_fraction"],
        color=plot_df["color"],
        edgecolor="white",
        linewidth=0.8,
        width=0.72,
    )

    for bar, value in zip(bars, plot_df["promoter_cgi_fraction"]):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            value + 0.015,
            f"{value:.3f}",
            ha="center",
            va="bottom",
            fontsize=9,
            color="#334155",
        )

    ax.set_ylim(0, max(0.75, plot_df["promoter_cgi_fraction"].max() + 0.08))
    ax.set_ylabel("TSS CGI coverage proportion among protein-coding genes", fontsize=10)
    ax.set_xlabel("Species", fontsize=10)
    ax.set_title("TSS CGI Coverage in Protein-Coding Genes", fontsize=12, fontweight="bold")
    ax.grid(axis="y", color="#e2e8f0", linewidth=0.8)
    ax.set_axisbelow(True)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    fig.tight_layout()
    fig.savefig(args.output, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved barplot to: {args.output}")
    for _, row in plot_df.iterrows():
        print(f"{row['label']}: {row['promoter_cgi_fraction']:.6f}")


if __name__ == "__main__":
    main()
