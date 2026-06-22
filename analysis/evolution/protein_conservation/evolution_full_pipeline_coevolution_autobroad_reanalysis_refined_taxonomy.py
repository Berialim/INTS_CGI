#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
evolution_full_pipeline_coevolution_autobroad_reanalysis_refined_taxonomy.py

输入：
    Compara.115.protein_default.homologies.tsv.gz
    ../gene_map.csv

输出：
    presence_matrix_coevolution.tsv
    completeness_summary_coevolution.tsv
    homolog_percentage_summary.tsv
    species_order_dnamethylation_coevolution.tsv
    class_complex_summary_coevolution.tsv
    broad_group_complex_summary.tsv
    broad_group_homolog_percentage_summary.tsv
    excluded_outlier_species.tsv
    excluded_unclassified_species.tsv
    coevolution_permutation_tests.tsv
    all_species_correlation_stats.tsv
    analysis_metadata_summary.md
    complex_correlation_pearson.tsv
    complex_correlation_spearman.tsv
    evolution_heatmap_coevolution.pdf
    broad_group_evolution_heatmap.pdf
    complex_pairwise_scatter_coevolution.pdf
    complex_pairwise_scatter_all_species.pdf
    complex_correlation_heatmap_coevolution.pdf
    complex_correlation_heatmap_all_species.pdf
    outlier_species_qc.pdf
    coevolution_rank_scatter.pdf

说明：
    1. 沿用当前脚本的 component 级别分析方式
    2. 蛋白支持度使用 identity + goc_score + wga_coverage 组合打分
    3. 使用分析文件中的全部可用物种，并在后续去除进化分类上的离群物种
    4. 对照组改为 Ribosome 与 RNA polymerase
    5. H3.3 分析中去掉 H3F3A/B，只保留伴侣与装载相关组分
    6. broad_group 使用物种名自动归类，不再依赖手写且不完整的 species map
    7. 新增多种图用于展示共进化与异常点
    8. 物种排序仍按 DNA methylation complex 的加权进化程度
    9. 新版会额外输出源文件的 species / homology_species 覆盖情况，便于确认缺失物种来自官方 dump 还是来自后续过滤
"""
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['pdf.fonttype'] = 42  # 使 PDF 中字体可编辑
matplotlib.rcParams['ps.fonttype'] = 42
plt.rcParams["font.family"] = "Arial"
from pathlib import Path
from itertools import combinations

import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import pearsonr, spearmanr
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LinearSegmentedColormap
from taxonomy_groups_refined import GROUP_COLORS, GROUP_ORDER, infer_broad_group


# ======================
# 0. 参数
# ======================
compara_file = "Compara.115.protein_default.human_related.normalized.tsv.gz"
map_file = "gene_map.csv"
REFERENCE_SPECIES = "homo_sapiens"
fallback_compara_file = "Compara.115.protein_default.homologies.tsv"
OUTPUT_DIR = Path("results/refined_taxonomy")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

PRIMARY_COMPLEXES = [
    "H3.3",
    "H2A.Z",
    "DNA_methylation",
    "Integrator",
    "H3K4me3",
    "snRNP",
    "SMN",
]
OUTLIER_ROBUST_Z_THRESHOLD = 3.5
PERMUTATION_N = 2000
LOW_CONFIDENCE_PENALTY = 0.5
WHITE_BLUE_CMAP = LinearSegmentedColormap.from_list(
    "white_blue_coevolution",
    ["#ffffff", "#e6f0ff", "#bfdbfe", "#93c5fd", "#60a5fa", "#2563eb", "#1e3a8a"],
)

# 为空时使用分析文件中的全部可用物种
FOCUS_SPECIES = []

phylo_hint = [
    "saccharomyces_cerevisiae",
    "schizosaccharomyces_pombe",
    "neurospora_crassa",
    "aspergillus_nidulans",
    "caenorhabditis_elegans",
    "caenorhabditis_briggsae",
    "pristionchus_pacificus",
    "drosophila_melanogaster",
    "anopheles_gambiae",
    "apis_mellifera",
    "tribolium_castaneum",
    "strongylocentrotus_purpuratus",
    "branchiostoma_floridae",
    "petromyzon_marinus",
    "danio_rerio",
    "takifugu_rubripes",
    "oryzias_latipes",
    "xenopus_tropicalis",
    "anolis_carolinensis",
    "chelonia_mydas",
    "gallus_gallus",
    "taeniopygia_guttata",
    "ornithorhynchus_anatinus",
    "monodelphis_domestica",
    "oryctolagus_cuniculus",
    "mus_musculus",
    "rattus_norvegicus",
    "canis_lupus_familiaris",
    "felis_catus",
    "sus_scrofa",
    "bos_taurus",
    "equus_caballus",
    "loxodonta_africana",
    "macaca_mulatta",
    "chlorocebus_sabaeus",
    "gorilla_gorilla",
    "pan_troglodytes",
    "homo_sapiens",
]
phylo_rank = {sp: i for i, sp in enumerate(phylo_hint)}

broad_group_order = GROUP_ORDER
broad_group_colors = GROUP_COLORS

species_display = {
    "saccharomyces_cerevisiae": "S. cerevisiae",
    "schizosaccharomyces_pombe": "S. pombe",
    "neurospora_crassa": "N. crassa",
    "aspergillus_nidulans": "A. nidulans",
    "caenorhabditis_elegans": "C. elegans",
    "caenorhabditis_briggsae": "C. briggsae",
    "pristionchus_pacificus": "P. pacificus",
    "drosophila_melanogaster": "D. melanogaster",
    "anopheles_gambiae": "A. gambiae",
    "apis_mellifera": "A. mellifera",
    "tribolium_castaneum": "T. castaneum",
    "strongylocentrotus_purpuratus": "S. purpuratus",
    "branchiostoma_floridae": "B. floridae",
    "petromyzon_marinus": "P. marinus",
    "danio_rerio": "D. rerio",
    "takifugu_rubripes": "T. rubripes",
    "oryzias_latipes": "O. latipes",
    "xenopus_tropicalis": "X. tropicalis",
    "anolis_carolinensis": "A. carolinensis",
    "chelonia_mydas": "C. mydas",
    "gallus_gallus": "G. gallus",
    "taeniopygia_guttata": "T. guttata",
    "ornithorhynchus_anatinus": "O. anatinus",
    "monodelphis_domestica": "M. domestica",
    "oryctolagus_cuniculus": "O. cuniculus",
    "mus_musculus": "M. musculus",
    "rattus_norvegicus": "R. norvegicus",
    "canis_lupus_familiaris": "C. familiaris",
    "felis_catus": "F. catus",
    "sus_scrofa": "S. scrofa",
    "bos_taurus": "B. taurus",
    "equus_caballus": "E. caballus",
    "loxodonta_africana": "L. africana",
    "macaca_mulatta": "M. mulatta",
    "chlorocebus_sabaeus": "C. sabaeus",
    "gorilla_gorilla": "G. gorilla",
    "pan_troglodytes": "P. troglodytes",
    "homo_sapiens": "H. sapiens",
}


# ======================
# 1. 定义 complexes / components
# ======================
complex_definitions = {
    "H3.3": [
        ("ASF1A/B", ["ASF1A", "ASF1B"]),
        ("HIRA", ["HIRA"]),
        ("UBN1", ["UBN1"]),
        ("CABIN1", ["CABIN1"]),
        ("DAXX", ["DAXX"]),
        ("ATRX", ["ATRX"]),
        ("DEK", ["DEK"]),
    ],
    "DNA_methylation": [
        ("DNMT1", ["DNMT1"]),
        ("UHRF1", ["UHRF1"]),
        ("UHRF2", ["UHRF2"]),
        ("DNMT3A", ["DNMT3A"]),
        ("DNMT3B", ["DNMT3B"]),
        ("DNMT3L", ["DNMT3L"]),
        ("TET1", ["TET1"]),
        ("TET2", ["TET2"]),
        ("TET3", ["TET3"]),
        ("TDG", ["TDG"]),
        ("HELLS", ["HELLS"]),
        ("CDCA7", ["CDCA7"]),
        ("ZBTB24", ["ZBTB24"]),
        ("MECP2", ["MECP2"]),
        ("MBD1", ["MBD1"]),
        ("MBD2", ["MBD2"]),
        ("MBD4", ["MBD4"]),
        ("MBD5/6", ["MBD5", "MBD6"]),
        ("ZBTB38", ["ZBTB38"]),
    ],
    "H3K4me3": [
        ("KMT2A/B", ["KMT2A", "KMT2B"]),
        ("KMT2C/D", ["KMT2C", "KMT2D"]),
        ("SETD1A/B", ["SETD1A", "SETD1B"]),
        ("WDR5", ["WDR5"]),
        ("RBBP5", ["RBBP5"]),
        ("ASH2L", ["ASH2L"]),
        ("DPY30", ["DPY30"]),
        ("CXXC1", ["CXXC1"]),
        ("WDR82", ["WDR82"]),
        ("HCFC1/2", ["HCFC1", "HCFC2"]),
        ("MEN1", ["MEN1"]),
        ("KDM5A", ["KDM5A"]),
        ("KDM5B", ["KDM5B"]),
        ("KDM5C", ["KDM5C"]),
        ("KDM5D", ["KDM5D"]),
        ("BPTF", ["BPTF"]),
        ("CHD1", ["CHD1"]),
        ("TAF3", ["TAF3"]),
        ("PHF8", ["PHF8"]),
        ("ING1/2", ["ING1", "ING2"]),
        ("ING4/5", ["ING4", "ING5"]),
        ("SPIN1/4", ["SPIN1", "SPIN4"]),
    ],
    "H2A.Z": [
        ("H2AFZ/V", ["H2AFZ", "H2AFV"]),
        ("SRCAP", ["SRCAP"]),
        ("VPS72", ["VPS72"]),
        ("ZNHIT1", ["ZNHIT1"]),
        ("DMAP1", ["DMAP1"]),
        ("RUVBL1/2", ["RUVBL1", "RUVBL2"]),
        ("EP400", ["EP400"]),
        ("YEATS4", ["YEATS4"]),
        ("KAT5", ["KAT5"]),
        ("TRRAP", ["TRRAP"]),
        ("EPC1/2", ["EPC1", "EPC2"]),
        ("ANP32E", ["ANP32E"]),
        ("ACTL6A/B", ["ACTL6A", "ACTL6B"]),
    ],
    "Integrator": [(f"INTS{i}", [f"INTS{i}"]) for i in range(1, 16)],
    "snRNP": [
        ("U1-70K", ["SNRNP70"]),
        ("U1A", ["SNRPA"]),
        ("U1C", ["SNRPC"]),
        ("SF3A1", ["SF3A1"]),
        ("SF3A2", ["SF3A2"]),
        ("SF3A3", ["SF3A3"]),
        ("SF3B1", ["SF3B1"]),
        ("SF3B2", ["SF3B2"]),
        ("SF3B3", ["SF3B3"]),
        ("PRPF8", ["PRPF8"]),
        ("EFTUD2", ["EFTUD2"]),
        ("SNRNP200", ["SNRNP200"]),
        ("PRPF6", ["PRPF6"]),
        ("PRPF31", ["PRPF31"]),
        ("SNRNP40", ["SNRNP40"]),
        ("TXNL4A", ["TXNL4A"]),
        ("SmB/B'", ["SNRPB"]),
        ("SmD2", ["SNRPD2"]),
        ("SmE", ["SNRPE"]),
    ],
    "SMN": [
        ("SMN1/2", ["SMN1", "SMN2"]),
        ("GEMIN2", ["GEMIN2"]),
        ("GEMIN3/DDX20", ["DDX20"]),
        ("GEMIN4", ["GEMIN4"]),
        ("GEMIN5", ["GEMIN5"]),
        ("GEMIN6", ["GEMIN6"]),
        ("GEMIN7", ["GEMIN7"]),
        ("GEMIN8", ["GEMIN8"]),
        ("STRAP", ["STRAP"]),
    ],
    "RNA_polymerase_control": [
        ("POLR2A", ["POLR2A"]),
        ("POLR2B", ["POLR2B"]),
        ("POLR2C", ["POLR2C"]),
        ("POLR2D", ["POLR2D"]),
        ("POLR2E", ["POLR2E"]),
        ("POLR2F", ["POLR2F"]),
        ("POLR2G", ["POLR2G"]),
        ("POLR2H", ["POLR2H"]),
        ("POLR2I", ["POLR2I"]),
        ("POLR2J", ["POLR2J"]),
        ("POLR2K", ["POLR2K"]),
        ("POLR2L", ["POLR2L"]),
    ],
    "Ribosome_control": [
        ("RPL3", ["RPL3"]),
        ("RPL5", ["RPL5"]),
        ("RPL7", ["RPL7"]),
        ("RPL11", ["RPL11"]),
        ("RPL23A", ["RPL23A"]),
        ("RPLP0", ["RPLP0"]),
        ("RPS3", ["RPS3"]),
        ("RPS6", ["RPS6"]),
        ("RPS8", ["RPS8"]),
        ("RPS14", ["RPS14"]),
        ("RPS18", ["RPS18"]),
    ],

}

target_genes = sorted(
    {
        gene
        for components in complex_definitions.values()
        for _, genes in components
        for gene in genes
    }
)

CORRELATION_COMPLEXES = list(complex_definitions.keys())
KEY_SPECIES_TO_TRACK = [
    "saccharomyces_cerevisiae",
    "schizosaccharomyces_pombe",
    "caenorhabditis_elegans",
    "drosophila_melanogaster",
    "ciona_intestinalis",
    "petromyzon_marinus",
    "danio_rerio",
    "xenopus_tropicalis",
    "oryzias_latipes",
    "gallus_gallus",
    "anolis_carolinensis",
    "mus_musculus",
    "mus_caroli",
    "canis_lupus_familiaris",
    "macaca_mulatta",
    "pan_troglodytes",
    "homo_sapiens",
]


# ======================
# 2. 工具函数
# ======================
def read_table_flexible(path):
    path = Path(path)
    guessed_sep = "\t" if path.suffix.lower() == ".tsv" else ","

    try:
        return pd.read_csv(path, sep=guessed_sep, low_memory=False)
    except Exception:
        return pd.read_csv(path, sep=None, engine="python")


def format_species_label(species_name):
    if species_name in species_display:
        return species_display[species_name]
    tokens = str(species_name).split("_")
    if len(tokens) >= 2:
        genus = tokens[0]
        species = tokens[1]
        return f"{genus[0].upper()}. {species}"
    return str(species_name)


def out_path(filename):
    return str(OUTPUT_DIR / filename)


def combine_evidence_score(df_in):
    metric_weights = {
        "homology_identity": 0.6,
        "goc_score": 0.2,
        "wga_coverage": 0.2,
    }

    weighted_sum = pd.Series(0.0, index=df_in.index, dtype=float)
    weight_sum = pd.Series(0.0, index=df_in.index, dtype=float)

    for col, weight in metric_weights.items():
        values = pd.to_numeric(df_in[col], errors="coerce").clip(0, 100) / 100.0
        valid = values.notna()
        weighted_sum.loc[valid] += values.loc[valid] * weight
        weight_sum.loc[valid] += weight

    return weighted_sum.div(weight_sum.where(weight_sum > 0)).fillna(0.0).clip(0.0, 1.0)


def apply_confidence_penalty(df_in, base_score, low_confidence_penalty):
    confidence = pd.to_numeric(df_in["is_high_confidence"], errors="coerce").fillna(0)
    factor = np.where(confidence == 1, 1.0, low_confidence_penalty)
    return (base_score * factor).clip(0.0, 1.0)


def restrict_to_focus_species(weighted_matrix, focus_species):
    matrix = weighted_matrix.copy()

    if "homo_sapiens" in focus_species and "homo_sapiens" not in matrix.columns:
        matrix["homo_sapiens"] = 1.0

    selected = [sp for sp in focus_species if sp in matrix.columns]
    return matrix[selected]


def select_all_available_species(gene_matrix):
    """
    使用分析文件中的全部可用物种。
    同时补入 human 参考物种列。
    """
    matrix = gene_matrix.copy()
    if "homo_sapiens" not in matrix.columns:
        matrix["homo_sapiens"] = 1.0

    support_count = (matrix > 0).sum(axis=0)
    selected = support_count[support_count > 0].index.tolist()
    return matrix[selected]


def build_component_matrix(gene_matrix, complex_defs, species_order):
    rows = []
    for system_name, components in complex_defs.items():
        for component_name, genes in components:
            available_genes = [gene for gene in genes if gene in gene_matrix.index]
            if available_genes:
                values = gene_matrix.loc[available_genes, species_order].max(axis=0)
            else:
                values = pd.Series(0.0, index=species_order, dtype=float)

            row = values.to_dict()
            row["component_id"] = f"{system_name}:{component_name}"
            row["component"] = component_name
            row["system"] = system_name
            rows.append(row)

    return pd.DataFrame(rows).set_index("component_id")


def build_component_group_matrix(component_df, species_meta_df, group_col):
    group_order = species_meta_df[group_col].drop_duplicates().tolist()
    grouped_values = []

    for _, row in component_df.iterrows():
        group_means = (
            species_meta_df.assign(value=[row[sp] for sp in species_meta_df["species"]])
            .groupby(group_col)["value"]
            .mean()
            .reindex(group_order)
        )
        group_row = group_means.to_dict()
        group_row["component_id"] = row.name
        group_row["component"] = row["component"]
        group_row["system"] = row["system"]
        grouped_values.append(group_row)

    return pd.DataFrame(grouped_values).set_index("component_id")


def calc_system_scores(component_df, complex_defs, species_order):
    summary = {}
    for system_name in complex_defs:
        sub = component_df[component_df["system"] == system_name]
        summary[system_name] = sub[species_order].mean(axis=0)
    return pd.DataFrame(summary, index=species_order).fillna(0.0)


def calc_system_homolog_percentages(component_df, complex_defs, species_order, threshold=0.0):
    summary = {}
    for system_name in complex_defs:
        sub = component_df[component_df["system"] == system_name]
        summary[system_name] = (sub[species_order] > threshold).mean(axis=0) * 100.0
    return pd.DataFrame(summary, index=species_order).fillna(0.0)


def make_taxonomy_summary(summary_df):
    df_out = summary_df.copy()
    df_out["taxonomy_class"] = [infer_broad_group(sp) for sp in df_out.index]
    df_out["taxonomy_rank"] = [
        broad_group_order.index(cls) if cls in broad_group_order else len(broad_group_order) + 100
        for cls in df_out["taxonomy_class"]
    ]
    class_summary = (
        df_out.groupby(["taxonomy_rank", "taxonomy_class"], sort=True)
        .mean(numeric_only=True)
        .reset_index()
        .sort_values(["taxonomy_rank", "taxonomy_class"])
        .drop(columns=["taxonomy_rank"])
        .set_index("taxonomy_class")
    )
    return class_summary


def make_group_summary(summary_df, group_map, group_order, group_col_name):
    df_out = summary_df.copy()
    df_out[group_col_name] = [group_map.get(sp, "Other") for sp in df_out.index]
    df_out["group_rank"] = [
        group_order.index(group) if group in group_order else len(group_order) + 100
        for group in df_out[group_col_name]
    ]
    group_summary = (
        df_out.groupby(["group_rank", group_col_name], sort=True)
        .mean(numeric_only=True)
        .reset_index()
        .sort_values(["group_rank", group_col_name])
        .drop(columns=["group_rank"])
        .set_index(group_col_name)
    )
    return group_summary


def filter_species_with_defined_group(species_list, group_map):
    retained = [sp for sp in species_list if group_map.get(sp) in broad_group_order]
    excluded = [sp for sp in species_list if group_map.get(sp) not in broad_group_order]
    return retained, excluded


def make_pairwise_correlation_table(summary_df, complexes, method):
    return summary_df[complexes].corr(method=method)


def make_pairwise_correlation_stats(summary_df, complexes):
    rows = []
    for x_col, y_col in combinations(complexes, 2):
        plot_df = summary_df[[x_col, y_col]].dropna().copy()
        pearson_r, pearson_p = pearsonr(plot_df[x_col], plot_df[y_col])
        spearman_rho, spearman_p = spearmanr(plot_df[x_col], plot_df[y_col])
        rows.append(
            {
                "x_complex": x_col,
                "y_complex": y_col,
                "n_species": len(plot_df),
                "pearson_r": pearson_r,
                "pearson_p": pearson_p,
                "spearman_rho": spearman_rho,
                "spearman_p": spearman_p,
            }
        )
    return pd.DataFrame(rows)


def permutation_test_correlations(summary_df, complexes, n_perm, random_seed=42):
    rng = np.random.default_rng(random_seed)
    rows = []
    for x_col, y_col in combinations(complexes, 2):
        observed_pearson = summary_df[x_col].corr(summary_df[y_col], method="pearson")
        observed_spearman = summary_df[x_col].corr(summary_df[y_col], method="spearman")

        perm_pearson = []
        perm_spearman = []
        y_values = summary_df[y_col].to_numpy(copy=True)

        for _ in range(n_perm):
            permuted = rng.permutation(y_values)
            perm_series = pd.Series(permuted, index=summary_df.index)
            perm_pearson.append(summary_df[x_col].corr(perm_series, method="pearson"))
            perm_spearman.append(summary_df[x_col].corr(perm_series, method="spearman"))

        perm_pearson = np.asarray(perm_pearson, dtype=float)
        perm_spearman = np.asarray(perm_spearman, dtype=float)

        p_pearson = (np.sum(np.abs(perm_pearson) >= abs(observed_pearson)) + 1) / (n_perm + 1)
        p_spearman = (np.sum(np.abs(perm_spearman) >= abs(observed_spearman)) + 1) / (n_perm + 1)

        rows.append(
            {
                "x_complex": x_col,
                "y_complex": y_col,
                "pearson_r": observed_pearson,
                "pearson_perm_p": p_pearson,
                "spearman_rho": observed_spearman,
                "spearman_perm_p": p_spearman,
                "n_species": len(summary_df),
            }
        )

    return pd.DataFrame(rows)


def detect_outlier_species(summary_df, primary_complexes, robust_z_threshold):
    """
    使用 robust z-score(MAD) 在主 complex 维度上识别离群物种。
    """
    work = summary_df[primary_complexes].copy()
    medians = work.median(axis=0)
    mad = (work - medians).abs().median(axis=0)
    mad = mad.replace(0, pd.NA)

    robust_z = ((work - medians).abs()).div(1.4826 * mad)
    robust_z = robust_z.fillna(0.0)

    outlier_mask = robust_z.max(axis=1) > robust_z_threshold
    outlier_df = robust_z.loc[outlier_mask].copy()
    if outlier_df.empty:
        return pd.DataFrame(columns=["species"] + primary_complexes + ["max_robust_z"])

    outlier_df["species"] = outlier_df.index
    outlier_df["max_robust_z"] = robust_z.loc[outlier_mask].max(axis=1)
    ordered_cols = ["species"] + primary_complexes + ["max_robust_z"]
    outlier_df = outlier_df[ordered_cols].sort_values(
        by=["max_robust_z", "species"], ascending=[False, True]
    )
    return outlier_df


def make_pairwise_regression_plot(summary_df, complexes, output_path, title_prefix=""):
    pairs = list(combinations(complexes, 2))
    with PdfPages(output_path) as pdf:
        for x_col, y_col in pairs:
            fig, ax = plt.subplots(figsize=(5, 5))
            plot_df = summary_df[[x_col, y_col]].dropna().copy()
            sns.regplot(
                data=plot_df,
                x=x_col,
                y=y_col,
                scatter_kws={"s": 36, "alpha": 0.75, "rasterized": True},
                line_kws={"color": "#c23b22", "linewidth": 2},
                ax=ax,
            )
            pearson_r, pearson_p = pearsonr(plot_df[x_col], plot_df[y_col])
            spearman_r, spearman_p = spearmanr(plot_df[x_col], plot_df[y_col])
            ax.set_title(f"{title_prefix}{x_col} vs {y_col}")
            ax.set_xlabel(x_col)
            ax.set_ylabel(y_col)
            ax.text(
                0.03,
                0.97,
                (
                    f"n = {len(plot_df)}\n"
                    f"Pearson r = {pearson_r:.2f}\n"
                    f"Pearson p = {pearson_p:.2e}\n"
                    f"Spearman rho = {spearman_r:.2f}\n"
                    f"Spearman p = {spearman_p:.2e}"
                ),
                transform=ax.transAxes,
                va="top",
                ha="left",
                fontsize=9,
                bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "alpha": 0.8},
            )
            fig.tight_layout()
            pdf.savefig(fig, dpi=300, bbox_inches="tight")
            plt.close(fig)


def make_pairwise_group_regression_plot(summary_df, complexes, output_path, title_prefix=""):
    pairs = list(combinations(complexes, 2))
    label_lookup = {group: group.replace("_", " ") for group in summary_df.index}
    with PdfPages(output_path) as pdf:
        for x_col, y_col in pairs:
            plot_df = summary_df[[x_col, y_col]].dropna().copy()
            if len(plot_df) < 2:
                continue
            fig, ax = plt.subplots(figsize=(5.4, 5.4))
            sns.regplot(
                data=plot_df,
                x=x_col,
                y=y_col,
                scatter=False,
                line_kws={"color": "#c23b22", "linewidth": 2},
                ax=ax,
            )
            ax.scatter(
                plot_df[x_col],
                plot_df[y_col],
                s=58,
                c=[broad_group_colors.get(group, "#475569") for group in plot_df.index],
                alpha=0.95,
                edgecolors="white",
                linewidths=0.7,
                zorder=3,
            )
            for group_name, row in plot_df.iterrows():
                ax.text(
                    row[x_col] + 0.01,
                    row[y_col],
                    label_lookup.get(group_name, group_name),
                    fontsize=7.8,
                    color=broad_group_colors.get(group_name, "#334155"),
                    va="center",
                    ha="left",
                )
            pearson_r, pearson_p = pearsonr(plot_df[x_col], plot_df[y_col])
            spearman_r, spearman_p = spearmanr(plot_df[x_col], plot_df[y_col])
            ax.set_title(f"{title_prefix}{x_col} vs {y_col}")
            ax.set_xlabel(x_col)
            ax.set_ylabel(y_col)
            ax.text(
                0.03,
                0.97,
                (
                    f"n groups = {len(plot_df)}\n"
                    f"Pearson r = {pearson_r:.2f}\n"
                    f"Pearson p = {pearson_p:.2e}\n"
                    f"Spearman rho = {spearman_r:.2f}\n"
                    f"Spearman p = {spearman_p:.2e}"
                ),
                transform=ax.transAxes,
                va="top",
                ha="left",
                fontsize=9,
                bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "alpha": 0.8},
            )
            fig.tight_layout()
            pdf.savefig(fig, dpi=300, bbox_inches="tight")
            plt.close(fig)


def make_outlier_qc_plot(summary_df, outlier_df, primary_complexes, output_path):
    medians = summary_df[primary_complexes].median(axis=0)
    mad = (summary_df[primary_complexes] - medians).abs().median(axis=0).replace(0, pd.NA)
    robust_z = ((summary_df[primary_complexes] - medians).abs()).div(1.4826 * mad).fillna(0.0)
    robust_z["max_robust_z"] = robust_z.max(axis=1)
    robust_z = robust_z.sort_values("max_robust_z", ascending=False)

    plt.figure(figsize=(max(8, len(robust_z) * 0.18), 4.6))
    colors = ["#c23b22" if sp in set(outlier_df["species"]) else "#4c78a8" for sp in robust_z.index]
    plt.bar(range(len(robust_z)), robust_z["max_robust_z"].values, color=colors)
    plt.axhline(OUTLIER_ROBUST_Z_THRESHOLD, color="black", linestyle="--", linewidth=1)
    plt.ylabel("Max robust z-score")
    plt.xlabel("Species")
    plt.title(f"Outlier QC Across Coevolution Complexes (n={len(robust_z)})")
    plt.xticks([])
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()


def make_rank_scatter_plot(summary_df, complexes, output_path):
    pairs = list(combinations(complexes, 2))
    with PdfPages(output_path) as pdf:
        for x_col, y_col in pairs:
            fig, ax = plt.subplots(figsize=(5, 5))
            rank_df = summary_df[[x_col, y_col]].rank(method="average")
            sns.regplot(
                data=rank_df,
                x=x_col,
                y=y_col,
                scatter_kws={"s": 32, "alpha": 0.75, "rasterized": True},
                line_kws={"color": "#1f6f8b", "linewidth": 2},
                ax=ax,
            )
            spearman_r, spearman_p = spearmanr(summary_df[x_col], summary_df[y_col])
            ax.set_title(f"Rank trend: {x_col} vs {y_col}")
            ax.set_xlabel(f"{x_col} rank")
            ax.set_ylabel(f"{y_col} rank")
            ax.text(
                0.03,
                0.97,
                (
                    f"n = {len(rank_df)}\n"
                    f"Spearman rho = {spearman_r:.2f}\n"
                    f"Spearman p = {spearman_p:.2e}"
                ),
                transform=ax.transAxes,
                va="top",
                ha="left",
                fontsize=9,
                bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "alpha": 0.8},
            )
            fig.tight_layout()
            pdf.savefig(fig, dpi=300, bbox_inches="tight")
            plt.close(fig)


def write_analysis_metadata_summary(
    output_path,
    complex_defs,
    species_meta_df,
    broad_group_summary_df,
    broad_group_homolog_percentage_summary_df,
    class_summary_df,
    outlier_df,
    excluded_unclassified_species,
):
    def df_block(df, index=True):
        if df.empty:
            return "(empty)"
        return df.to_string(index=index)

    lines = []
    lines.append("# Analysis Metadata Summary")
    lines.append("")
    lines.append("## Important Notes")
    lines.append("")
    lines.append(
        "- This analysis is based on a human-reference Ensembl Compara homology table; "
        "deeply diverged complexes in fungi and other distant taxa can be systematically underestimated."
    )
    lines.append(
        f"- Low-confidence orthologs are retained with a penalty factor of {LOW_CONFIDENCE_PENALTY:.2f} "
        "instead of being discarded completely."
    )
    lines.append("")
    lines.append("## Complex Definitions")
    lines.append("")
    for complex_name, components in complex_defs.items():
        component_strings = [f"- {component}: {', '.join(genes)}" for component, genes in components]
        lines.append(f"### {complex_name}")
        lines.extend(component_strings)
        lines.append("")

    lines.append("## Broad Group Species Counts")
    lines.append("")
    broad_counts = species_meta_df.groupby("broad_group").size().sort_index()
    for group_name, count in broad_counts.items():
        lines.append(f"- {group_name}: {count}")
    lines.append("")

    lines.append("## Taxonomy Class Species Counts")
    lines.append("")
    class_counts = species_meta_df.groupby("taxonomy_class").size().sort_index()
    for class_name, count in class_counts.items():
        lines.append(f"- {class_name}: {count}")
    lines.append("")

    lines.append("## Broad Group Summary Table")
    lines.append("")
    lines.append("```text")
    lines.append(df_block(broad_group_summary_df))
    lines.append("```")
    lines.append("")

    lines.append("## Broad Group Homolog Percentage Summary")
    lines.append("")
    lines.append("```text")
    lines.append(df_block(broad_group_homolog_percentage_summary_df))
    lines.append("```")
    lines.append("")

    lines.append("## Taxonomy Class Summary Table")
    lines.append("")
    lines.append("```text")
    lines.append(df_block(class_summary_df))
    lines.append("```")
    lines.append("")

    lines.append("## Species Metadata")
    lines.append("")
    lines.append("```text")
    lines.append(df_block(species_meta_df, index=False))
    lines.append("```")
    lines.append("")

    lines.append("## Excluded Outlier Species")
    lines.append("")
    if outlier_df.empty:
        lines.append("- None")
    else:
        lines.append("```text")
        lines.append(df_block(outlier_df, index=False))
        lines.append("```")
    lines.append("")

    lines.append("## Excluded Unclassified Species")
    lines.append("")
    if excluded_unclassified_species:
        for sp in excluded_unclassified_species:
            lines.append(f"- {sp}")
    else:
        lines.append("- None")
    lines.append("")

    Path(output_path).write_text("\n".join(lines), encoding="utf-8")


def write_source_species_reports(df_in, output_species_path, output_key_species_path):
    species_summary = (
        df_in.groupby("species")
        .size()
        .rename("n_rows")
        .reset_index()
        .sort_values(["n_rows", "species"], ascending=[False, True])
    )
    homology_summary = (
        df_in.groupby("homology_species")
        .size()
        .rename("n_rows")
        .reset_index()
        .rename(columns={"homology_species": "species"})
        .sort_values(["n_rows", "species"], ascending=[False, True])
    )
    species_summary["column_role"] = "species"
    homology_summary["column_role"] = "homology_species"
    combined = pd.concat([species_summary, homology_summary], ignore_index=True)
    combined.to_csv(output_species_path, sep="\t", index=False)

    observed_homology_species = set(df_in["homology_species"].dropna())
    observed_species = set(df_in["species"].dropna())
    rows = []
    for species in KEY_SPECIES_TO_TRACK:
        rows.append(
            {
                "species": species,
                "present_in_species_column": species in observed_species,
                "present_in_homology_species_column": species in observed_homology_species,
            }
        )
    pd.DataFrame(rows).to_csv(output_key_species_path, sep="\t", index=False)


# ======================
# 3. 读取 Compara
# ======================
print("Loading compara...")
compara_path = Path(compara_file)
if not compara_path.exists():
    print(f"Primary compara source not found yet: {compara_file}")
    print(f"Falling back to local copy: {fallback_compara_file}")
    compara_path = Path(fallback_compara_file)
try:
    df = pd.read_csv(compara_path, sep="\t", low_memory=False)
except Exception as exc:
    fallback_path = Path(fallback_compara_file)
    if compara_path != fallback_path and fallback_path.exists():
        print(f"Primary compara source could not be read ({exc}).")
        print(f"Falling back to local copy: {fallback_compara_file}")
        compara_path = fallback_path
        df = pd.read_csv(compara_path, sep="\t", low_memory=False)
    else:
        raise
write_source_species_reports(
    df_in=df,
    output_species_path=out_path("source_species_coverage.tsv"),
    output_key_species_path=out_path("key_species_presence_in_source.tsv"),
)
unique_species = sorted(df["species"].dropna().unique().tolist())
print(f"Source species column unique values: {unique_species}")
if REFERENCE_SPECIES in set(df["species"].dropna()):
    df = df[df["species"] == REFERENCE_SPECIES].copy()
elif len(unique_species) == 1:
    print(
        "Reference species filter skipped because the source dump already contains a single species column value."
    )
else:
    raise ValueError(
        f"Reference species {REFERENCE_SPECIES} not found in source species column: {unique_species[:10]}"
    )
df = df[df["homology_type"].str.contains("ortholog", na=False)].copy()
df = df[df["gene_stable_id"].notna()].copy()


# ======================
# 4. 构建 gene support matrix
# ======================
print("Building matrix...")
base_score = combine_evidence_score(df)
df["protein_weight"] = apply_confidence_penalty(
    df_in=df,
    base_score=base_score,
    low_confidence_penalty=LOW_CONFIDENCE_PENALTY,
)

matrix = df.pivot_table(
    index="gene_stable_id",
    columns="homology_species",
    values="protein_weight",
    aggfunc="max",
    fill_value=0,
)


# ======================
# 5. gene id -> gene name
# ======================
print("Mapping gene names...")
map_df = read_table_flexible(map_file)
required_cols = {"Gene stable ID", "Gene name"}
missing_cols = required_cols - set(map_df.columns)
if missing_cols:
    raise ValueError(
        f"gene_map 文件缺少必要列: {sorted(missing_cols)}; 当前列为: {list(map_df.columns)}"
    )

map_df = map_df[["Gene stable ID", "Gene name"]].rename(
    columns={"Gene stable ID": "gene_stable_id", "Gene name": "gene_name"}
)

matrix = matrix.reset_index().merge(map_df, on="gene_stable_id", how="left")
matrix = matrix.dropna(subset=["gene_name"]).copy()
species_cols = [c for c in matrix.columns if c not in {"gene_stable_id", "gene_name"}]
matrix = matrix.groupby("gene_name", as_index=True)[species_cols].max()


# ======================
# 6. 目标基因与物种
# ======================
matrix_targets = matrix.reindex(target_genes).fillna(0.0).copy()

if FOCUS_SPECIES:
    matrix_targets = restrict_to_focus_species(matrix_targets, FOCUS_SPECIES)
else:
    matrix_targets = select_all_available_species(matrix_targets)

all_species = matrix_targets.columns.tolist()
broad_group_map = {sp: infer_broad_group(sp) for sp in all_species}

if not all_species:
    raise ValueError("未找到任何目标物种列，无法继续分析。")

initial_component_matrix = build_component_matrix(matrix_targets, complex_definitions, all_species)


# ======================
# 7. complex score、离群物种过滤与排序
# ======================
print("Calculating weighted system scores...")
initial_summary = calc_system_scores(initial_component_matrix, complex_definitions, all_species)

print("Skipping outlier species removal...")
outlier_df = pd.DataFrame(columns=["species"] + PRIMARY_COMPLEXES + ["max_robust_z"])
outlier_df.to_csv(out_path("excluded_outlier_species.tsv"), sep="\t", index=False)

retained_species = all_species.copy()

retained_species, excluded_unclassified_species = filter_species_with_defined_group(
    retained_species,
    broad_group_map,
)
pd.DataFrame({"species": excluded_unclassified_species}).to_csv(
    out_path("excluded_unclassified_species.tsv"), sep="\t", index=False
)
if not retained_species:
    raise ValueError("去除未分类物种后没有剩余物种，无法继续分析。")

matrix_targets = matrix_targets[retained_species]
all_species = retained_species
component_matrix = build_component_matrix(matrix_targets, complex_definitions, all_species)
summary = calc_system_scores(component_matrix, complex_definitions, all_species)
homolog_percentage_summary = calc_system_homolog_percentages(
    component_matrix, complex_definitions, all_species
)

sort_df = summary.copy()
sort_df["phylo_rank"] = [phylo_rank.get(sp, len(phylo_rank) + 10000) for sp in sort_df.index]
sort_df["species_name"] = sort_df.index

species_order = (
    sort_df.sort_values(
        by=["DNA_methylation", "phylo_rank", "species_name"],
        ascending=[True, True, True],
    )
    .index.tolist()
)

summary = summary.loc[species_order]
summary.to_csv(out_path("completeness_summary_coevolution.tsv"), sep="\t")
homolog_percentage_summary = homolog_percentage_summary.loc[species_order]
homolog_percentage_summary.to_csv(out_path("homolog_percentage_summary.tsv"), sep="\t")

species_meta = pd.DataFrame(
    {
        "species": species_order,
        "display_name": [format_species_label(sp) for sp in species_order],
        "taxonomy_class": [infer_broad_group(sp) for sp in species_order],
        "broad_group": [broad_group_map.get(sp, infer_broad_group(sp)) for sp in species_order],
    }
)
species_meta.to_csv(out_path("species_order_dnamethylation_coevolution.tsv"), sep="\t", index=False)
broad_group_summary = make_group_summary(
    summary_df=summary,
    group_map=broad_group_map,
    group_order=broad_group_order,
    group_col_name="broad_group",
)
broad_group_summary.to_csv(out_path("broad_group_complex_summary.tsv"), sep="\t")
broad_group_homolog_percentage_summary = make_group_summary(
    summary_df=homolog_percentage_summary,
    group_map=broad_group_map,
    group_order=broad_group_order,
    group_col_name="broad_group",
)
broad_group_homolog_percentage_summary.to_csv(
    out_path("broad_group_homolog_percentage_summary.tsv"), sep="\t"
)


# ======================
# 8. component 画图矩阵
# ======================
matrix_plot = component_matrix.copy()
system_order = list(complex_definitions.keys())
matrix_plot["system"] = pd.Categorical(matrix_plot["system"], categories=system_order, ordered=True)
matrix_plot["component_order"] = range(len(matrix_plot))
matrix_plot = matrix_plot.sort_values(by=["system", "component_order"])
matrix_plot_export = matrix_plot.copy()
matrix_plot_export.to_csv(out_path("presence_matrix_coevolution.tsv"), sep="\t")

group_component_matrix = build_component_group_matrix(
    component_df=matrix_plot,
    species_meta_df=species_meta,
    group_col="broad_group",
)
group_component_matrix_only = group_component_matrix[broad_group_summary.index.tolist()].copy()
group_component_matrix_only.index = group_component_matrix["component"]


# ======================
# 9. 大分类 heatmap
# ======================
print("Plotting broad-group heatmap...")
heatmap_width = max(10, len(broad_group_summary.index) * 1.0)
heatmap_height = max(10, len(group_component_matrix_only.index) * 0.36)

plt.figure(figsize=(heatmap_width, heatmap_height))
ax = sns.heatmap(
    group_component_matrix_only,
    cmap=WHITE_BLUE_CMAP,
    vmin=0,
    vmax=1,
    cbar_kws={"label": "Evolutionary support score"},
    linewidths=0.3,
    linecolor="white",
)

boundary_positions = []
current = 0
for system_name in system_order[:-1]:
    current += (matrix_plot["system"] == system_name).sum()
    boundary_positions.append(current)
for pos in boundary_positions:
    ax.hlines(pos, *ax.get_xlim(), colors="black", linewidth=1.0)

ax.set_xticklabels(group_component_matrix_only.columns.tolist(), rotation=35, ha="right")
plt.title(f"Broad-Group Component Support Heatmap (n={len(species_meta)})")
plt.xlabel("Broad taxonomic group")
plt.ylabel("Components")
plt.tight_layout()
plt.savefig(out_path("evolution_heatmap_coevolution.pdf"), dpi=300, bbox_inches="tight")
plt.close()


# ======================
# 10. 汇总表，不再绘制小分类热图
# ======================
print("Preparing taxonomy-level summary...")
class_summary = make_taxonomy_summary(summary)
class_summary.to_csv(out_path("class_complex_summary_coevolution.tsv"), sep="\t")


# ======================
# 11. 额外输出 broad-group 图
# ======================
print("Saving broad-group summary figures...")

plt.figure(figsize=(max(10, len(broad_group_summary.index) * 1.05), 4.8))
sns.heatmap(
    broad_group_summary.T,
    cmap=WHITE_BLUE_CMAP,
    vmin=0,
    vmax=1,
    annot=True,
    fmt=".2f",
    linewidths=0.4,
    linecolor="white",
    cbar_kws={"label": "Mean evolutionary score"},
)
plt.title(f"Complex Evolution Across Broad Taxonomic Groups (n={len(species_meta)})")
plt.xlabel("Broad taxonomic group")
plt.ylabel("Complex")
plt.xticks(rotation=35, ha="right")
plt.tight_layout()
plt.savefig(out_path("broad_group_evolution_heatmap.pdf"), dpi=300, bbox_inches="tight")
plt.close()


# ======================
# 12. 新增：共进化相关性分析
# ======================
print("Calculating coevolution correlations...")
all_species_corr_stats = make_pairwise_correlation_stats(initial_summary, CORRELATION_COMPLEXES)
all_species_corr_stats.to_csv(out_path("all_species_correlation_stats.tsv"), sep="\t", index=False)

all_species_pearson_corr = make_pairwise_correlation_table(
    initial_summary, CORRELATION_COMPLEXES, method="pearson"
)
plt.figure(figsize=(max(5.2, len(CORRELATION_COMPLEXES) * 1.2), max(4.6, len(CORRELATION_COMPLEXES) * 1.0)))
sns.heatmap(
    all_species_pearson_corr,
    cmap="vlag",
    vmin=-1,
    vmax=1,
    annot=True,
    fmt=".2f",
    linewidths=0.4,
    linecolor="white",
    cbar_kws={"label": "Pearson correlation"},
)
plt.title(f"All-Species Complex Correlation Including Controls (n={len(initial_summary)})")
plt.tight_layout()
plt.savefig(out_path("complex_correlation_heatmap_all_species.pdf"), dpi=300, bbox_inches="tight")
plt.close()

make_pairwise_regression_plot(
    summary_df=initial_summary,
    complexes=CORRELATION_COMPLEXES,
    output_path=out_path("complex_pairwise_scatter_all_species.pdf"),
    title_prefix="All species: ",
)
initial_broad_group_summary = make_group_summary(
    summary_df=initial_summary,
    group_map={sp: infer_broad_group(sp) for sp in initial_summary.index},
    group_order=broad_group_order,
    group_col_name="broad_group",
)
make_pairwise_group_regression_plot(
    summary_df=initial_broad_group_summary,
    complexes=CORRELATION_COMPLEXES,
    output_path=out_path("complex_pairwise_scatter_all_species_broad_groups.pdf"),
    title_prefix="All species broad groups: ",
)

pearson_corr = make_pairwise_correlation_table(summary, CORRELATION_COMPLEXES, method="pearson")
spearman_corr = make_pairwise_correlation_table(summary, CORRELATION_COMPLEXES, method="spearman")
perm_test_df = permutation_test_correlations(
    summary_df=summary,
    complexes=CORRELATION_COMPLEXES,
    n_perm=PERMUTATION_N,
)

pearson_corr.to_csv(out_path("complex_correlation_pearson.tsv"), sep="\t")
spearman_corr.to_csv(out_path("complex_correlation_spearman.tsv"), sep="\t")
perm_test_df.to_csv(out_path("coevolution_permutation_tests.tsv"), sep="\t", index=False)

plt.figure(figsize=(max(5.2, len(CORRELATION_COMPLEXES) * 1.2), max(4.6, len(CORRELATION_COMPLEXES) * 1.0)))
sns.heatmap(
    pearson_corr,
    cmap="vlag",
    vmin=-1,
    vmax=1,
    annot=True,
    fmt=".2f",
    linewidths=0.4,
    linecolor="white",
    cbar_kws={"label": "Pearson correlation"},
)
plt.title(f"Complex Coevolution Correlation Including Controls (n={len(summary)})")
plt.tight_layout()
plt.savefig(out_path("complex_correlation_heatmap_coevolution.pdf"), dpi=300, bbox_inches="tight")
plt.close()

make_pairwise_regression_plot(
    summary_df=summary,
    complexes=CORRELATION_COMPLEXES,
    output_path=out_path("complex_pairwise_scatter_coevolution.pdf"),
    title_prefix="Filtered species: ",
)
make_pairwise_group_regression_plot(
    summary_df=broad_group_summary,
    complexes=CORRELATION_COMPLEXES,
    output_path=out_path("complex_pairwise_scatter_coevolution_broad_groups.pdf"),
    title_prefix="Filtered broad groups: ",
)
make_rank_scatter_plot(
    summary_df=summary,
    complexes=CORRELATION_COMPLEXES,
    output_path=out_path("coevolution_rank_scatter.pdf"),
)
make_outlier_qc_plot(
    summary_df=initial_summary,
    outlier_df=outlier_df,
    primary_complexes=PRIMARY_COMPLEXES,
    output_path=out_path("outlier_species_qc.pdf"),
)
write_analysis_metadata_summary(
    output_path=out_path("analysis_metadata_summary.md"),
    complex_defs=complex_definitions,
    species_meta_df=species_meta,
    broad_group_summary_df=broad_group_summary,
    broad_group_homolog_percentage_summary_df=broad_group_homolog_percentage_summary,
    class_summary_df=class_summary,
    outlier_df=outlier_df,
    excluded_unclassified_species=excluded_unclassified_species,
)

print("Done!")
print(f"Species included: {len(species_order)}")
print(f"Species order saved to: {out_path('species_order_dnamethylation_coevolution.tsv')}")
print(f"Class summary saved to: {out_path('class_complex_summary_coevolution.tsv')}")
print(f"Broad-group summary saved to: {out_path('broad_group_complex_summary.tsv')}")
print(f"Homolog percentage summary saved to: {out_path('homolog_percentage_summary.tsv')}")
print(
    f"Broad-group homolog percentage summary saved to: {out_path('broad_group_homolog_percentage_summary.tsv')}"
)
print(f"Excluded outlier species saved to: {out_path('excluded_outlier_species.tsv')}")
print(f"Excluded unclassified species saved to: {out_path('excluded_unclassified_species.tsv')}")
print(
    f"Correlation tables saved to: {out_path('complex_correlation_pearson.tsv')}, "
    f"{out_path('complex_correlation_spearman.tsv')}"
)
print(
    f"Broad-group scatter plots saved to: {out_path('complex_pairwise_scatter_all_species_broad_groups.pdf')}, "
    f"{out_path('complex_pairwise_scatter_coevolution_broad_groups.pdf')}"
)
print(f"Permutation tests saved to: {out_path('coevolution_permutation_tests.tsv')}")
print(f"All-species correlation stats saved to: {out_path('all_species_correlation_stats.tsv')}")
print(f"Metadata summary saved to: {out_path('analysis_metadata_summary.md')}")
print(f"Source species coverage saved to: {out_path('source_species_coverage.tsv')}")
print(f"Key species source presence saved to: {out_path('key_species_presence_in_source.tsv')}")
print(f"Compara source used: {compara_path}")
