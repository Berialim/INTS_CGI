#!/usr/bin/env python3

import sys
import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from statannotations.Annotator import Annotator

cgi_bed = sys.argv[1]
gene_bed = sys.argv[2]
prefix = sys.argv[3]

# =========================
# Step 1: 提取 TSS
# 同时保存原始 gene bed 信息
# =========================
tss_file = prefix + ".tss.bed"
gene_records = []

with open(gene_bed) as f, open(tss_file, "w") as out:
    for line in f:
        if line.startswith("#") or line.strip() == "":
            continue

        fields = line.strip().split()
        chrom, start, end, name, score, strand = fields[:6]
        start = int(start)
        end = int(end)

        if strand == "+":
            tss_start = start
            tss_end = start + 1
        else:
            tss_start = end - 1
            tss_end = end

        out.write(f"{chrom}\t{tss_start}\t{tss_end}\t{name}\t0\t{strand}\n")

        gene_records.append({
            "chr": chrom,
            "start": start,
            "end": end,
            "gene": name,
            "score": score,
            "strand": strand,
            "tss_start": tss_start,
            "tss_end": tss_end,
            "length": end - start
        })

gene_df = pd.DataFrame(gene_records)

# =========================
# Step 2: intersect
# =========================
intersect_file = prefix + ".tss_cgi_intersect.bed"
cmd = f"bedtools intersect -a {tss_file} -b {cgi_bed} -wa -wb > {intersect_file}"
subprocess.run(cmd, shell=True, check=True)

# =========================
# Step 3: 解析 intersect
# =========================
cols = [
    "chr", "tss_start", "tss_end", "gene", "score", "strand",
    "cgi_chr", "cgi_start", "cgi_end"
]

try:
    df = pd.read_csv(intersect_file, sep="\t", header=None)
    df = df.iloc[:, :9]
    df.columns = cols
except pd.errors.EmptyDataError:
    df = pd.DataFrame(columns=cols)

# 输出文件
out1 = open(prefix + ".tss_overlapped_cgi.bed", "w")
out2 = open(prefix + ".tss_to_cgi_edge.bed", "w")
out3 = open(prefix + ".upstream_to_cgi_edge_opposite.bed", "w")
out4 = open(prefix + ".genes_without_tss_cgi.bed", "w")

# =========================
# Step 4: 构建三个区域
# =========================
for _, row in df.iterrows():
    chrom = row["chr"]
    tss = int(row["tss_start"])
    strand = row["strand"]
    gene = row["gene"]

    cgi_start = int(row["cgi_start"])
    cgi_end = int(row["cgi_end"])

    # 1. 覆盖TSS的CGI
    out1.write(f"{chrom}\t{cgi_start}\t{cgi_end}\t{gene}\t0\t{strand}\n")

    # 2. TSS → CGI边界（基因方向）
    if strand == "+":
        start = tss
        end = cgi_end
    else:
        start = cgi_start
        end = tss

    if start < end:
        out2.write(f"{chrom}\t{start}\t{end}\t{gene}\t0\t{strand}\n")

    # 3. 上游 → CGI边界（反方向）
    if strand == "+":
        start = cgi_start
        end = tss
    else:
        start = tss
        end = cgi_end

    if start < end:
        opposite_strand = "-" if strand == "+" else "+"
        out3.write(f"{chrom}\t{start}\t{end}\t{gene}\t0\t{opposite_strand}\n")

# =========================
# Step 5: 输出 TSS 上没有 CGI 的基因
# =========================
overlapped_keys = set(
    zip(df["chr"], df["tss_start"], df["tss_end"], df["gene"], df["strand"])
)

for rec in gene_records:
    key = (rec["chr"], rec["tss_start"], rec["tss_end"], rec["gene"], rec["strand"])
    if key not in overlapped_keys:
        out4.write(
            f'{rec["chr"]}\t{rec["start"]}\t{rec["end"]}\t{rec["gene"]}\t{rec["score"]}\t{rec["strand"]}\n'
        )

out1.close()
out2.close()
out3.close()
out4.close()

# =========================
# Step 6: 统计长度并输出 PDF
# =========================
pdf_file = prefix + ".gene_length_compare.pdf"

gene_df["group"] = "No CGI"
gene_df.loc[
    gene_df.apply(
        lambda r: (r["chr"], r["tss_start"], r["tss_end"], r["gene"], r["strand"]) in overlapped_keys,
        axis=1
    ),
    "group"
] = "With CGI"

plot_df = gene_df[["gene", "length", "group"]].copy()

# 字体设置
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["font.family"] = "Arial"

sns.set(style="whitegrid", font="Arial")
fig, ax = plt.subplots(figsize=(10, 10))

order = ["With CGI", "No CGI"]

sns.boxplot(
    data=plot_df,
    x="group",
    y="length",
    order=order,
    ax=ax,
    width=0.8,
    showfliers=False
)



count_with = (plot_df["group"] == "With CGI").sum()
count_without = (plot_df["group"] == "No CGI").sum()

ax.set_xlabel("")
ax.set_ylabel("Gene length", fontsize=42)
ax.set_xticklabels([
    f"With CGI\nn={count_with}",
    f"No CGI\nn={count_without}"
], fontsize=42)
ax.tick_params(axis="y", labelsize=42)

# Annotator 标注 full pvalue
if count_with > 0 and count_without > 0:
    pairs = [("With CGI", "No CGI")]
    annotator = Annotator(
        ax,
        pairs,
        data=plot_df,
        x="group",
        y="length",
        order=order
    )
    annotator.configure(
        test="Mann-Whitney",
        text_format="full",
        loc="outside",
        fontsize=42,
        verbose=0
    )
    annotator.apply_and_annotate()

plt.tight_layout()
plt.savefig(pdf_file, format="pdf", bbox_inches="tight")
plt.close()

print("Done!")
