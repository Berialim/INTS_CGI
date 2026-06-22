#!/usr/bin/env python3

import argparse
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
plt.rcParams["font.family"] = "Arial"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot a CpG-region length histogram with quartile guides."
    )
    parser.add_argument("--bed", required=True, help="Input BED file")
    parser.add_argument("--output", required=True, help="Output PDF or PNG file")
    parser.add_argument("--color", default="#E68414", help="Histogram color")
    parser.add_argument("--max-length", type=int, default=1300, help="Y-axis upper limit")
    return parser.parse_args()


def main():
    args = parse_args()
    cols = ["chrom", "start", "end"]
    df = pd.read_csv(args.bed, sep="\t", header=None, usecols=[0, 1, 2], names=cols)
    df["length"] = df["end"] - df["start"]

    q1, q2, q3 = df["length"].quantile([0.25, 0.5, 0.75]).tolist()

    fig, ax = plt.subplots(figsize=(6, 8))
    counts, bins, patches = ax.hist(
        df["length"],
        bins=800,
        color=args.color,
        edgecolor="none",
        alpha=0.6,
        orientation="horizontal",
    )

    for patch, bottom_edge in zip(patches, bins[:-1]):
        if bottom_edge < q1:
            patch.set_alpha(0.25)
        elif bottom_edge < q2:
            patch.set_alpha(0.5)
        elif bottom_edge < q3:
            patch.set_alpha(0.75)
        else:
            patch.set_alpha(1.0)

    for q, label in zip([q1, q2, q3], ["Q1", "Q2", "Q3"]):
        ax.axhline(q, color="black", linestyle="--")
        ax.text(ax.get_xlim()[1] * 0.9, q, label, va="center", ha="right", fontsize=10)

    ax.set_title("Length distribution with quartiles", fontsize=14)
    ax.set_ylabel("Length")
    ax.set_xlabel("Count")
    ax.set_ylim(0, args.max_length)
    plt.tight_layout()
    plt.savefig(args.output, dpi=300)


if __name__ == "__main__":
    main()
