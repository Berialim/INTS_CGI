#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Reuse the 26-species promoter CGI summary and add U1 snRNP / Integrator
protein-system scores for matched phylogenetic heatmaps and regplots.
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
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import pearsonr, spearmanr

import promoter_cgi_evolution_pipeline as cgi_pipe


U1_COMPONENTS = ["U1-70K", "U1A", "U1C"]
INTEGRATOR_COMPONENTS = [f"INTS{i}" for i in range(1, 16)]

LEGACY_GROUP_ORDER = [
    "Fungi",
    "Invertebrates",
    "Basal_chordates",
    "Fish_Amphibians",
    "Sauropsids",
    "Mammals_non_primate",
    "Primates_Hominids",
]

LEGACY_GROUP_LABELS = {
    "Fungi": "Fungi",
    "Invertebrates": "Invertebrates",
    "Basal_chordates": "Basal chordates",
    "Fish_Amphibians": "Fish / Amphibians",
    "Sauropsids": "Sauropsids",
    "Mammals_non_primate": "Other mammals",
    "Primates_Hominids": "Primates / Hominids",
}

LEGACY_GROUP_COLORS = {
    "Fungi": "#b07aa1",
    "Invertebrates": "#f28e2b",
    "Basal_chordates": "#76b7b2",
    "Fish_Amphibians": "#4e79a7",
    "Sauropsids": "#59a14f",
    "Mammals_non_primate": "#9c755f",
    "Primates_Hominids": "#e15759",
}

LEGACY_TREE_TEMPLATE = {
    "label": "root",
    "children": [
        {"group": "Fungi", "label": "Fungi"},
        {
            "label": "Metazoa",
            "children": [
                {"group": "Invertebrates", "label": "Invertebrates"},
                {
                    "label": "Chordates",
                    "children": [
                        {"group": "Basal_chordates", "label": "Basal chordates"},
                        {
                            "label": "Vertebrates",
                            "children": [
                                {"group": "Fish_Amphibians", "label": "Fish / Amphibians"},
                                {
                                    "label": "Amniotes",
                                    "children": [
                                        {"group": "Sauropsids", "label": "Sauropsids"},
                                        {
                                            "label": "Mammalia",
                                            "children": [
                                                {"group": "Mammals_non_primate", "label": "Other mammals"},
                                                {"group": "Primates_Hominids", "label": "Primates / Hominids"},
                                            ],
                                        },
                                    ],
                                },
                            ],
                        },
                    ],
                },
            ],
        },
    ],
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot 26-species promoter CGI evolution against U1 snRNP and Integrator protein scores."
    )
    parser.add_argument(
        "--species-summary",
        default="promoter_cgi_evolution/results/species_promoter_cgi_summary.tsv",
        help="26-species promoter CGI summary TSV.",
    )
    parser.add_argument(
        "--presence-matrix",
        default="presence_matrix_coevolution.tsv",
        help="Protein component support matrix TSV.",
    )
    parser.add_argument(
        "--output-dir",
        default="promoter_cgi_evolution/results_u1_integrator",
        help="Output directory for U1/Integrator CGI comparison.",
    )
    return parser.parse_args()


def add_u1_integrator_scores(species_summary: pd.DataFrame, presence_matrix_path: Path) -> pd.DataFrame:
    component_df = pd.read_csv(presence_matrix_path, sep="\t")
    species_cols = [species for species in species_summary["species"] if species in component_df.columns]
    component_lookup = component_df.set_index("component")[species_cols]

    def mean_signal(components: list[str]) -> pd.Series:
        available = [component for component in components if component in component_lookup.index]
        if not available:
            return pd.Series(0.0, index=species_cols, dtype=float)
        return component_lookup.loc[available].mean(axis=0)

    u1_signal = mean_signal(U1_COMPONENTS)
    integrator_signal = mean_signal(INTEGRATOR_COMPONENTS)

    out = species_summary.copy()
    out["U1_snRNP"] = out["species"].map(u1_signal.to_dict()).fillna(0.0)
    out["Integrator"] = out["species"].map(integrator_signal.to_dict()).fillna(0.0)
    return out


def convert_to_legacy_taxonomy(species_summary: pd.DataFrame) -> pd.DataFrame:
    out = species_summary.copy()
    legacy_map = {
        "Fish": "Fish_Amphibians",
        "Amphibians": "Fish_Amphibians",
        "Monotremes_Marsupials": "Mammals_non_primate",
    }
    out["taxonomy_class"] = out["taxonomy_class"].replace(legacy_map)
    return out


def sort_species_summary_legacy(species_summary: pd.DataFrame) -> pd.DataFrame:
    df = species_summary.copy()
    df["display_name"] = df["display_name"].fillna(df["species"])
    df["group_rank"] = df["taxonomy_class"].map({group: idx for idx, group in enumerate(LEGACY_GROUP_ORDER)})
    sort_cols = ["group_rank", "DNA_methylation", "promoter_cgi_fraction", "display_name"]
    ascending = [True, True, False, True]
    for col in ["DNA_methylation", "promoter_cgi_fraction"]:
        if col not in df.columns:
            df[col] = np.nan
    return df.sort_values(sort_cols, ascending=ascending).reset_index(drop=True)


def build_named_pair_correlation_rows(
    species_summary: pd.DataFrame,
    group_summary: pd.DataFrame,
) -> pd.DataFrame:
    pair_defs = [
        ("promoter_cgi_fraction", "Integrator", "Promoter CGI", "Integrator"),
        ("U1_snRNP", "Integrator", "U1 snRNP", "Integrator"),
    ]
    rows = []
    for level, df_in, n_label in [
        ("all_species", species_summary, "n_species"),
        ("broad_group", group_summary, "n_groups"),
    ]:
        for x_col, y_col, x_label, y_label in pair_defs:
            pair = df_in[[x_col, y_col]].dropna().copy()
            if len(pair) < 2:
                continue
            try:
                pearson_r, pearson_p = pearsonr(pair[x_col], pair[y_col])
            except ValueError:
                pearson_r, pearson_p = np.nan, np.nan
            try:
                spearman_rho, spearman_p = spearmanr(pair[x_col], pair[y_col])
            except ValueError:
                spearman_rho, spearman_p = np.nan, np.nan
            rows.append(
                {
                    "level": level,
                    "x_metric": x_col,
                    "y_metric": y_col,
                    "x_label": x_label,
                    "y_label": y_label,
                    "pearson_r": pearson_r,
                    "pearson_p": pearson_p,
                    "spearman_rho": spearman_rho,
                    "spearman_p": spearman_p,
                    n_label: len(pair),
                }
            )
    return pd.DataFrame(rows)


def plot_named_pair_correlations(
    species_summary: pd.DataFrame,
    group_summary: pd.DataFrame,
    output_path: Path,
) -> None:
    label_lookup = dict(zip(group_summary["species"], group_summary["display_name"]))
    pair_defs = [
        ("promoter_cgi_fraction", "Integrator", "Promoter CGI fraction", "Integrator protein score"),
        ("U1_snRNP", "Integrator", "U1 snRNP protein score", "Integrator protein score"),
    ]
    with PdfPages(output_path) as pdf:
        for x_col, y_col, x_label, y_label in pair_defs:
            fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.8))

            species_pair = species_summary[[x_col, y_col]].dropna().copy()
            sns.regplot(
                data=species_pair,
                x=x_col,
                y=y_col,
                scatter_kws={"s": 40, "alpha": 0.8, "rasterized": True},
                line_kws={"color": "#c23b22", "linewidth": 2},
                ax=axes[0],
            )
            try:
                pearson_r, pearson_p = pearsonr(species_pair[x_col], species_pair[y_col])
            except ValueError:
                pearson_r, pearson_p = np.nan, np.nan
            try:
                spearman_rho, spearman_p = spearmanr(species_pair[x_col], species_pair[y_col])
            except ValueError:
                spearman_rho, spearman_p = np.nan, np.nan
            axes[0].set_title(f"All species: {x_label} vs {y_label}")
            axes[0].set_xlabel(x_label)
            axes[0].set_ylabel(y_label)
            axes[0].text(
                0.03,
                0.97,
                (
                    f"n = {len(species_pair)}\n"
                    f"Pearson r = {pearson_r:.2f}\n"
                    f"Pearson p = {pearson_p:.2e}\n"
                    f"Spearman rho = {spearman_rho:.2f}\n"
                    f"Spearman p = {spearman_p:.2e}"
                ),
                transform=axes[0].transAxes,
                va="top",
                ha="left",
                fontsize=9,
                bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "alpha": 0.85},
            )

            group_pair = group_summary[[x_col, y_col, "species"]].dropna().copy()
            sns.regplot(
                data=group_pair,
                x=x_col,
                y=y_col,
                scatter=False,
                line_kws={"color": "#c23b22", "linewidth": 2},
                ax=axes[1],
            )
            axes[1].scatter(
                group_pair[x_col],
                group_pair[y_col],
                s=58,
                c=[cgi_pipe.GROUP_COLORS.get(group, "#475569") for group in group_pair["species"]],
                alpha=0.95,
                edgecolors="white",
                linewidths=0.7,
                zorder=3,
            )
            for _, row in group_pair.iterrows():
                group_name = row["species"]
                axes[1].text(
                    row[x_col] + 0.01,
                    row[y_col],
                    label_lookup.get(group_name, group_name).replace("_", " "),
                    fontsize=7.8,
                    color=cgi_pipe.GROUP_COLORS.get(group_name, "#334155"),
                    va="center",
                    ha="left",
                )
            try:
                pearson_r, pearson_p = pearsonr(group_pair[x_col], group_pair[y_col])
            except ValueError:
                pearson_r, pearson_p = np.nan, np.nan
            try:
                spearman_rho, spearman_p = spearmanr(group_pair[x_col], group_pair[y_col])
            except ValueError:
                spearman_rho, spearman_p = np.nan, np.nan
            axes[1].set_title(f"Broad groups: {x_label} vs {y_label}")
            axes[1].set_xlabel(x_label)
            axes[1].set_ylabel(y_label)
            axes[1].text(
                0.03,
                0.97,
                (
                    f"n groups = {len(group_pair)}\n"
                    f"Pearson r = {pearson_r:.2f}\n"
                    f"Pearson p = {pearson_p:.2e}\n"
                    f"Spearman rho = {spearman_rho:.2f}\n"
                    f"Spearman p = {spearman_p:.2e}"
                ),
                transform=axes[1].transAxes,
                va="top",
                ha="left",
                fontsize=9,
                bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "alpha": 0.85},
            )

            fig.tight_layout()
            pdf.savefig(fig, dpi=300, bbox_inches="tight")
            plt.close(fig)


def render_outputs(species_summary: pd.DataFrame, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    species_summary = convert_to_legacy_taxonomy(species_summary)

    original_group_order = cgi_pipe.GROUP_ORDER
    original_group_labels = cgi_pipe.GROUP_LABELS
    original_group_colors = cgi_pipe.GROUP_COLORS
    original_tree_template = cgi_pipe.TREE_TEMPLATE
    cgi_pipe.GROUP_ORDER = LEGACY_GROUP_ORDER
    cgi_pipe.GROUP_LABELS = LEGACY_GROUP_LABELS
    cgi_pipe.GROUP_COLORS = LEGACY_GROUP_COLORS
    cgi_pipe.TREE_TEMPLATE = LEGACY_TREE_TEMPLATE

    species_summary = sort_species_summary_legacy(species_summary)
    species_summary.to_csv(output_dir / "species_promoter_cgi_u1_integrator_summary.tsv", sep="\t", index=False)

    group_summary = cgi_pipe.build_group_summary(species_summary)
    group_summary.to_csv(output_dir / "broad_group_promoter_cgi_u1_integrator_summary.tsv", sep="\t", index=False)

    correlation_targets = ["U1_snRNP", "Integrator"]
    rows = []
    for metric, label in cgi_pipe.get_cgi_metric_specs():
        for target in correlation_targets:
            pair = species_summary[[metric, target]].dropna().copy()
            if len(pair) < 3:
                continue
            pearson_r = pair[metric].corr(pair[target], method="pearson")
            spearman_rho = pair[metric].corr(pair[target], method="spearman")
            rows.append(
                {
                    "cgi_metric": metric,
                    "protein_system": target,
                    "pearson_r": pearson_r,
                    "spearman_rho": spearman_rho,
                    "n_species": len(pair),
                }
            )
    corr_df = pd.DataFrame(rows)
    corr_df.to_csv(output_dir / "promoter_cgi_vs_u1_integrator_correlations.tsv", sep="\t", index=False)
    named_pair_df = build_named_pair_correlation_rows(species_summary, group_summary)
    named_pair_df.to_csv(output_dir / "integrator_cgi_u1_pairwise_correlations.tsv", sep="\t", index=False)
    plot_named_pair_correlations(
        species_summary,
        group_summary,
        output_dir / "integrator_cgi_u1_pairwise_correlations.pdf",
    )

    original_target_specs = cgi_pipe.get_protein_target_specs
    original_heatmap_columns = cgi_pipe.HEATMAP_COLUMNS
    original_header_colors = cgi_pipe.HEATMAP_HEADER_COLORS.copy()
    original_column_specs = cgi_pipe.HEATMAP_COLUMN_SPECS.copy()
    try:
        cgi_pipe.get_protein_target_specs = lambda: [
            ("U1_snRNP", "U1 snRNP protein score"),
            ("Integrator", "Integrator protein score"),
        ]
        cgi_pipe.HEATMAP_COLUMNS = [
            ("promoter_cgi_fraction", "Promoter"),
            ("random_ctrl_cgi_fraction", "Background"),
            ("cgi_delta", "Prom-Bg"),
            ("Integrator", "Integrator"),
            ("U1_snRNP", "U1 snRNP"),
            ("RNA_polymerase_control", "RNAP II"),
            ("Ribosome_control", "Ribosome"),
        ]
        cgi_pipe.HEATMAP_HEADER_COLORS.update(
            {
                "U1_snRNP": "#1d4ed8",
                "Integrator": "#0f172a",
                "RNA_polymerase_control": "#1e40af",
                "Ribosome_control": "#1e3a8a",
            }
        )
        cgi_pipe.HEATMAP_COLUMN_SPECS.update(
            {
                "U1_snRNP": {
                    "cmap": cgi_pipe.COMPLEX_CMAP,
                    "norm": cgi_pipe.Normalize(vmin=0.0, vmax=1.0),
                    "group": "Complex scores",
                },
                "Integrator": {
                    "cmap": cgi_pipe.COMPLEX_CMAP,
                    "norm": cgi_pipe.Normalize(vmin=0.0, vmax=1.0),
                    "group": "Complex scores",
                },
                "RNA_polymerase_control": {
                    "cmap": cgi_pipe.COMPLEX_CMAP,
                    "norm": cgi_pipe.Normalize(vmin=0.0, vmax=1.0),
                    "group": "Complex scores",
                },
                "Ribosome_control": {
                    "cmap": cgi_pipe.COMPLEX_CMAP,
                    "norm": cgi_pipe.Normalize(vmin=0.0, vmax=1.0),
                    "group": "Complex scores",
                },
            }
        )

        cgi_pipe.make_regplots(species_summary, output_dir / "promoter_cgi_vs_u1_integrator_regplots.pdf")
        cgi_pipe.make_group_regplots(group_summary, output_dir / "promoter_cgi_vs_u1_integrator_broad_group_regplots.pdf")

        species_meta = species_summary[["species", "display_name", "taxonomy_class"]].copy()
        cgi_pipe.plot_phylogenetic_heatmap(
            summary_df=species_summary,
            row_meta=species_meta,
            output_path=output_dir / "species_promoter_cgi_u1_integrator_phylogenetic_heatmap.pdf",
            title="Species-Level Phylogenetic Heatmap of Promoter CGI, U1 snRNP, Integrator and Control Complex Evolution",
            subtitle="Protein-coding gene promoters are summarized for all 26 species, with promoter CGI metrics aligned to U1 snRNP, Integrator, RNAP II and Ribosome protein-system scores.",
        )

        group_meta = group_summary[["species", "display_name", "taxonomy_class"]].copy()
        cgi_pipe.plot_phylogenetic_heatmap(
            summary_df=group_summary,
            row_meta=group_meta,
            output_path=output_dir / "broad_group_promoter_cgi_u1_integrator_phylogenetic_heatmap.pdf",
            title="Broad-Group Phylogenetic Heatmap of Promoter CGI, U1 snRNP, Integrator and Control Complex Evolution",
            subtitle="Broad-group means are shown for Promoter, Background, Prom-Bg, U1 snRNP, Integrator, RNAP II and Ribosome.",
        )
    finally:
        cgi_pipe.get_protein_target_specs = original_target_specs
        cgi_pipe.HEATMAP_COLUMNS = original_heatmap_columns
        cgi_pipe.HEATMAP_HEADER_COLORS.clear()
        cgi_pipe.HEATMAP_HEADER_COLORS.update(original_header_colors)
        cgi_pipe.HEATMAP_COLUMN_SPECS.clear()
        cgi_pipe.HEATMAP_COLUMN_SPECS.update(original_column_specs)
        cgi_pipe.GROUP_ORDER = original_group_order
        cgi_pipe.GROUP_LABELS = original_group_labels
        cgi_pipe.GROUP_COLORS = original_group_colors
        cgi_pipe.TREE_TEMPLATE = original_tree_template


def main() -> None:
    args = parse_args()
    species_summary = pd.read_csv(Path(args.species_summary), sep="\t")
    species_summary = add_u1_integrator_scores(species_summary, Path(args.presence_matrix))
    render_outputs(species_summary, Path(args.output_dir))
    print(f"Saved outputs to: {args.output_dir}")


if __name__ == "__main__":
    main()
