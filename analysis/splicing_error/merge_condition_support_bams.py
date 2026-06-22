#!/usr/bin/env python3
"""Merge annotated/novel support BAMs by condition from database.csv."""

import argparse
from collections import defaultdict
from pathlib import Path

import pandas as pd
import pysam


def parse_args():
    parser = argparse.ArgumentParser(
        description="Merge annotated/novel support BAMs by condition."
    )
    parser.add_argument("-d", "--database", required=True, help="database.csv with name,condition")
    parser.add_argument(
        "-i", "--input-dir", default="error_splicing/igv_junction_bams",
        help="Directory containing per-sample support BAMs",
    )
    parser.add_argument(
        "-o", "--output-dir", default="error_splicing/igv_condition_bams",
        help="Directory for merged condition BAMs",
    )
    return parser.parse_args()


def merge_bams(input_paths, output_path):
    if not input_paths:
        return False
    unsorted_path = output_path.with_suffix(".unsorted.bam")
    with pysam.AlignmentFile(str(input_paths[0]), "rb") as template:
        with pysam.AlignmentFile(str(unsorted_path), "wb", template=template) as out_bam:
            for bam_path in input_paths:
                with pysam.AlignmentFile(str(bam_path), "rb") as bam_in:
                    for read in bam_in.fetch(until_eof=True):
                        out_bam.write(read)
    pysam.sort("-o", str(output_path), str(unsorted_path))
    unsorted_path.unlink(missing_ok=True)
    pysam.index(str(output_path))
    return True


def main():
    args = parse_args()
    db = pd.read_csv(args.database, dtype={"name": str})
    if not {"name", "condition"}.issubset(db.columns):
        raise ValueError("database.csv must contain name and condition columns")

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    grouped = defaultdict(lambda: {"annotated": [], "novel": []})
    for row in db.itertuples(index=False):
        sample = str(row.name)
        condition = str(row.condition)
        ann_path = input_dir / f"{sample}.annotated.support.bam"
        nov_path = input_dir / f"{sample}.novel.support.bam"
        if ann_path.exists():
            grouped[condition]["annotated"].append(ann_path)
        if nov_path.exists():
            grouped[condition]["novel"].append(nov_path)

    for condition, bam_map in grouped.items():
        for kind, paths in bam_map.items():
            output_path = output_dir / f"{condition}.{kind}.merged.bam"
            if merge_bams(paths, output_path):
                print(f"[OK] {condition} {kind}: merged {len(paths)} BAMs -> {output_path.name}")
            else:
                print(f"[SKIP] {condition} {kind}: no input BAMs found")


if __name__ == "__main__":
    main()
