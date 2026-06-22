#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  4 19:47:19 2025

@author: berial
"""
# generate length histgram with bed file, and quantile with length, and annotate histgram with 4 quantile alpha
#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os 
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42  # 使 PDF 中字体可编辑
plt.rcParams["font.family"] = "Arial"
os.chdir("/Users/berial/reference/hg38/CpG/")

# =========================
# 参数设置
# =========================
bed_file = "CpGtoPROMPTs.bed"   # 修改为你的 bed 文件
output_file = "length_CpGtoPROMPTs_histogram.pdf"


def plot(bed_file, output_file, col):
    cols = ["chrom", "start", "end"]
    df = pd.read_csv(bed_file, sep="\t", header=None, usecols=[0,1,2], names=cols)
    
    # 计算长度
    df["length"] = df["end"] - df["start"]
    
    # =========================
    # 计算分位数
    # =========================
    quantiles = df["length"].quantile([0.25, 0.5, 0.75]).to_dict()
    q1, q2, q3 = quantiles[0.25], quantiles[0.5], quantiles[0.75]
    
    print("Q1 (25%):", q1)
    print("Q2 (50%):", q2)
    print("Q3 (75%):", q3)
    
    # =========================
    # 绘制横向直方图
    # =========================
    fig, ax = plt.subplots(figsize=(6,8))  # 交换宽高
    
    counts, bins, patches = ax.hist(
        df["length"], bins=800, color=col, edgecolor='none', alpha=0.6, orientation="horizontal"
    )
    
    # 给每个分位区间着色
    for patch, bottom_edge in zip(patches, bins[:-1]):
        if bottom_edge < q1:
            patch.set_alpha(0.25)
        elif bottom_edge < q2:
            patch.set_alpha(0.5)
        elif bottom_edge < q3:
            patch.set_alpha(0.75)
        else:
            patch.set_alpha(1)
    
    # 标注分位数横线
    for q, label in zip([q1,q2,q3], ["Q1","Q2","Q3"]):
        ax.axhline(q, color="black", linestyle="--")
        ax.text(ax.get_xlim()[1]*0.9, q, label, va="center", ha="right", fontsize=10)
    
    ax.set_title("Length distribution with quartiles", fontsize=14)
    ax.set_ylabel("Length")
    ax.set_xlabel("Count")
    plt.ylim(0,1300)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.show()
    
plot(bed_file, output_file, "#E68414")
plot("CpGtoGene.bed", "length_CpGtoGene_histogram.pdf", "#1E95D4")
