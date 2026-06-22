#!/usr/bin/env python3
"""Export junction TSV files to IGV-friendly junction tracks.

For each *.junctions.tsv file, generate two BEDPE files:
  - annotated junctions (blue track)
  - novel junctions (red track)

Also generate BED12 files with color + score derived from read_count, which are
more useful for showing junction strength in IGV.
"""

import argparse
import csv
import os
from pathlib import Path


BLUE = "0,0,255"
RED = "255,0,0"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Export junction TSV files as colored BEDPE tracks for IGV."
    )
    parser.add_argument(
        "-i", "--input-dir", default="error_splicing",
        help="Directory containing *.junctions.tsv files",
    )
    parser.add_argument(
        "-o", "--output-dir", default="error_splicing/igv_bedpe",
        help="Directory for junction track outputs",
    )
    return parser.parse_args()


def normalize_type(raw_type):
    value = str(raw_type).strip().lower()
    return "annotated" if value == "annotated" else "novel"


def junction_to_bedpe_row(sample, row, index):
    chrom = str(row["chrom"])
    donor = int(float(row["donor"]))
    acceptor = int(float(row["acceptor"]))
    read_count = int(float(row["read_count"]))
    junction_type = normalize_type(row["type"])

    # BEDPE start coordinates are 0-based; use 1-bp anchors around splice sites.
    start1 = donor - 1
    end1 = donor
    start2 = acceptor - 1
    end2 = acceptor
    name = f"{sample}_{junction_type}_{index}"

    return [
        chrom, start1, end1,
        chrom, start2, end2,
        name, read_count, ".", ".", junction_type, read_count,
    ]


def junction_to_strength_bedpe_row(sample, row, index):
    chrom = str(row["chrom"])
    donor = int(float(row["donor"]))
    acceptor = int(float(row["acceptor"]))
    read_count = int(float(row["read_count"]))
    junction_type = normalize_type(row["type"])

    start1 = donor - 1
    end1 = donor
    start2 = acceptor - 1
    end2 = acceptor
    name = f"{sample}|{junction_type}|reads={read_count}|id={index}"

    # Keep read_count in score-related columns so IGV can use it for interaction strength.
    return [
        chrom, start1, end1,
        chrom, start2, end2,
        name, read_count, ".", ".", junction_type, read_count,
    ]


def junction_to_bed12_row(sample, row, index):
    chrom = str(row["chrom"])
    donor = int(float(row["donor"]))
    acceptor = int(float(row["acceptor"]))
    read_count = int(float(row["read_count"]))
    junction_type = normalize_type(row["type"])

    left = min(donor, acceptor)
    right = max(donor, acceptor)
    chrom_start = left - 1
    chrom_end = right
    score = min(read_count, 1000)
    color = BLUE if junction_type == "annotated" else RED
    name = f"{sample}|{junction_type}|reads={read_count}|id={index}"

    # Two 1-bp blocks mark donor and acceptor.
    block_sizes = "1,1"
    block_starts = f"0,{chrom_end - chrom_start - 1}"

    return [
        chrom,
        chrom_start,
        chrom_end,
        name,
        score,
        ".",
        chrom_start,
        chrom_end,
        color,
        2,
        block_sizes,
        block_starts,
    ]


def write_track(path, track_name, color, rows):
    with open(path, "w", newline="", encoding="utf-8") as handle:
        handle.write(f'track name="{track_name}" description="{track_name}" color={color}\n')
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        for row in rows:
            writer.writerow(row)


def export_file(tsv_path, output_dir):
    sample = tsv_path.name.replace(".junctions.tsv", "")
    annotated_rows = []
    novel_rows = []
    annotated_strength_bedpe_rows = []
    novel_strength_bedpe_rows = []
    annotated_bed12_rows = []
    novel_bed12_rows = []

    with open(tsv_path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for index, row in enumerate(reader, start=1):
            bedpe_row = junction_to_bedpe_row(sample, row, index)
            strength_bedpe_row = junction_to_strength_bedpe_row(sample, row, index)
            bed12_row = junction_to_bed12_row(sample, row, index)
            if bedpe_row[10] == "annotated":
                annotated_rows.append(bedpe_row)
                annotated_strength_bedpe_rows.append(strength_bedpe_row)
                annotated_bed12_rows.append(bed12_row)
            else:
                novel_rows.append(bedpe_row)
                novel_strength_bedpe_rows.append(strength_bedpe_row)
                novel_bed12_rows.append(bed12_row)

    annotated_path = output_dir / f"{sample}.annotated.bedpe"
    novel_path = output_dir / f"{sample}.novel.bedpe"
    annotated_strength_bedpe_path = output_dir / f"{sample}.annotated.strength.bedpe"
    novel_strength_bedpe_path = output_dir / f"{sample}.novel.strength.bedpe"
    annotated_bed12_path = output_dir / f"{sample}.annotated.strength.bed"
    novel_bed12_path = output_dir / f"{sample}.novel.strength.bed"

    write_track(annotated_path, f"{sample} annotated junctions", BLUE, annotated_rows)
    write_track(novel_path, f"{sample} novel junctions", RED, novel_rows)
    write_track(annotated_strength_bedpe_path, f"{sample} annotated junction strength (bedpe)", BLUE, annotated_strength_bedpe_rows)
    write_track(novel_strength_bedpe_path, f"{sample} novel junction strength (bedpe)", RED, novel_strength_bedpe_rows)
    write_track(annotated_bed12_path, f"{sample} annotated junction strength", BLUE, annotated_bed12_rows)
    write_track(novel_bed12_path, f"{sample} novel junction strength", RED, novel_bed12_rows)

    return {
        "sample": sample,
        "annotated_path": annotated_path,
        "annotated_count": len(annotated_rows),
        "novel_path": novel_path,
        "novel_count": len(novel_rows),
        "annotated_strength_bedpe_path": annotated_strength_bedpe_path,
        "novel_strength_bedpe_path": novel_strength_bedpe_path,
        "annotated_bed12_path": annotated_bed12_path,
        "novel_bed12_path": novel_bed12_path,
    }


def main():
    args = parse_args()
    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    files = sorted(input_dir.glob("*.junctions.tsv"))
    if not files:
        raise FileNotFoundError(f"No *.junctions.tsv files found in {input_dir}")

    summaries = []
    for tsv_path in files:
        summary = export_file(tsv_path, output_dir)
        summaries.append(summary)
        print(
            f"[OK] {summary['sample']}: "
            f"annotated={summary['annotated_count']} -> "
            f"{summary['annotated_path'].name} / {summary['annotated_strength_bedpe_path'].name} / {summary['annotated_bed12_path'].name}, "
            f"novel={summary['novel_count']} -> "
            f"{summary['novel_path'].name} / {summary['novel_strength_bedpe_path'].name} / {summary['novel_bed12_path'].name}"
        )


if __name__ == "__main__":
    main()
