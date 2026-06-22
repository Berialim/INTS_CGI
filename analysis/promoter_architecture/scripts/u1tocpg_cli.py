#!/usr/bin/env python3

import argparse
import pandas as pd
import pybedtools


def parse_args():
    parser = argparse.ArgumentParser(
        description="Summarize strand-aware U1-to-CpG distances for overlapping intervals."
    )
    parser.add_argument("--u1-bed", required=True, help="BED6 file of predicted U1 sites")
    parser.add_argument("--cpg-bed", required=True, help="BED6 file of CpG-associated regions")
    parser.add_argument(
        "--output",
        default="cpg_u1_distance.tsv",
        help="Output TSV path",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    u1_binding_sites = pybedtools.BedTool(args.u1_bed)
    cpg_islands = pybedtools.BedTool(args.cpg_bed)
    intersections = u1_binding_sites.intersect(cpg_islands, wa=True, wb=True, s=True)

    data = []
    for feature in intersections:
        u1_start = int(feature[1])
        u1_end = int(feature[2])
        u1_strand = feature[5]
        cpg_start = int(feature[7])
        cpg_end = int(feature[8])
        cpg_name = feature[9]

        u1_center = (u1_start + u1_end) / 2
        if u1_strand == "+":
            distance = cpg_end - u1_center
        else:
            distance = u1_center - cpg_start
        data.append((cpg_name, distance))

    df = pd.DataFrame(data, columns=["name", "distance_tocpg"])
    if df.empty:
        df.to_csv(args.output, sep="\t", index=False)
        return

    df = df.groupby("name")["distance_tocpg"].max().reset_index()
    cpg_df = pd.read_csv(
        args.cpg_bed,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "name", "score", "strand"],
    )
    cpg_df["cpg_length"] = cpg_df["end"] - cpg_df["start"]
    df = df.merge(cpg_df[["name", "cpg_length"]], on="name", how="left")
    df.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
