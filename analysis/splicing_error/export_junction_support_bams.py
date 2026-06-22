#!/usr/bin/env python3
"""Export BAM subsets supporting annotated or novel junctions.

For each sample with both BAM and *.junctions.tsv present:
  - collect splice junction coordinates from the TSV
  - scan the BAM CIGAR for spliced alignments (N operations)
  - keep reads supporting annotated junctions in one BAM
  - keep reads supporting novel junctions in another BAM
  - create BAM indexes

This is intended for IGV visualization of real splice-supporting read strength.
"""

import argparse
import re
from bisect import bisect_right
from pathlib import Path

import pysam

ATTR_RE = re.compile(r'([A-Za-z0-9_]+)\s+"([^"]+)"')


def parse_args():
    parser = argparse.ArgumentParser(
        description="Export annotated/novel junction-supporting BAM subsets for IGV."
    )
    parser.add_argument(
        "-i", "--input-dir", default="error_splicing",
        help="Directory containing *.junctions.tsv files",
    )
    parser.add_argument(
        "-b", "--bam-dir", default=".",
        help="Directory containing source BAM files",
    )
    parser.add_argument(
        "-o", "--output-dir", default="error_splicing/igv_junction_bams",
        help="Directory for BAM subset outputs",
    )
    parser.add_argument(
        "-g", "--gtf", required=True,
        help="GTF used to require both splice sites fall within the same gene",
    )
    return parser.parse_args()


def normalize_type(raw_type):
    value = str(raw_type).strip().lower()
    return "annotated" if value == "annotated" else "novel"


def parse_gtf_attributes(attr_text):
    return {key: value for key, value in ATTR_RE.findall(attr_text)}


def load_gene_index(gtf_path):
    gene_bounds = {}
    with open(gtf_path, "r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            chrom, _, feature, start, end, _, _, _, attrs = fields
            if feature != "transcript":
                continue
            if chrom.startswith("yr"):
                chrom = chrom[2:]
            attr_map = parse_gtf_attributes(attrs)
            gene_name = attr_map.get("gene_name") or attr_map.get("gene_id")
            if not gene_name:
                continue
            gene_bounds.setdefault(chrom, {}).setdefault(gene_name, []).append((int(start), int(end)))

    chrom_index = {}
    for chrom, genes in gene_bounds.items():
        entries = []
        for gene_name, intervals in genes.items():
            entries.append((min(x[0] for x in intervals), max(x[1] for x in intervals), gene_name))
        entries.sort(key=lambda x: (x[0], x[1], x[2]))
        chrom_index[chrom] = {"entries": entries, "starts": [x[0] for x in entries]}
    return chrom_index


def genes_covering_pos(chrom_index, chrom, pos):
    if chrom not in chrom_index:
        return set()
    data = chrom_index[chrom]
    idx = bisect_right(data["starts"], pos)
    genes = set()
    for i in range(idx):
        start, end, gene_name = data["entries"][i]
        if end < pos:
            continue
        if start <= pos <= end:
            genes.add(gene_name)
    return genes


def load_junction_sets(tsv_path, chrom_index):
    annotated = set()
    novel = set()
    kept = 0
    skipped = 0
    with open(tsv_path, "r", encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        index = {name: i for i, name in enumerate(header)}
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            chrom = fields[index["chrom"]]
            donor = int(fields[index["donor"]])
            acceptor = int(fields[index["acceptor"]])
            donor_genes = genes_covering_pos(chrom_index, chrom, donor)
            acceptor_genes = genes_covering_pos(chrom_index, chrom, acceptor)
            if not (donor_genes & acceptor_genes):
                skipped += 1
                continue
            key = (chrom, donor, acceptor)
            if normalize_type(fields[index["type"]]) == "annotated":
                annotated.add(key)
            else:
                novel.add(key)
            kept += 1
    return annotated, novel, kept, skipped


def extract_splice_junctions(read):
    if read.is_unmapped or read.cigartuples is None:
        return []
    ref_pos = read.reference_start
    junctions = []
    for op, length in read.cigartuples:
        # M, =, X
        if op in (0, 7, 8):
            ref_pos += length
        # D, N
        elif op in (2, 3):
            if op == 3:
                donor = ref_pos
                acceptor = ref_pos + length
                junctions.append((read.reference_name, donor, acceptor))
            ref_pos += length
        # I, S, H, P do not consume reference
        else:
            continue
    return junctions


def export_sample(sample, bam_path, junction_path, output_dir, chrom_index):
    annotated_set, novel_set, kept_junctions, skipped_junctions = load_junction_sets(junction_path, chrom_index)
    annotated_out = output_dir / f"{sample}.annotated.support.bam"
    novel_out = output_dir / f"{sample}.novel.support.bam"

    written_annotated = set()
    written_novel = set()
    ann_count = 0
    nov_count = 0

    with pysam.AlignmentFile(bam_path, "rb") as bam_in:
        with pysam.AlignmentFile(annotated_out, "wb", template=bam_in) as bam_ann, \
             pysam.AlignmentFile(novel_out, "wb", template=bam_in) as bam_nov:
            for read in bam_in.fetch(until_eof=True):
                junctions = extract_splice_junctions(read)
                if not junctions:
                    continue

                has_annotated = any(junc in annotated_set for junc in junctions)
                has_novel = any(junc in novel_set for junc in junctions)

                key = (read.query_name, read.flag, read.reference_id, read.reference_start)
                if has_annotated and key not in written_annotated:
                    bam_ann.write(read)
                    written_annotated.add(key)
                    ann_count += 1
                if has_novel and key not in written_novel:
                    bam_nov.write(read)
                    written_novel.add(key)
                    nov_count += 1

    pysam.index(str(annotated_out))
    pysam.index(str(novel_out))
    return ann_count, nov_count, annotated_out, novel_out, kept_junctions, skipped_junctions


def main():
    args = parse_args()
    input_dir = Path(args.input_dir)
    bam_dir = Path(args.bam_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    chrom_index = load_gene_index(args.gtf)

    junction_files = sorted(input_dir.glob("*.junctions.tsv"))
    if not junction_files:
        raise FileNotFoundError(f"No *.junctions.tsv files found in {input_dir}")

    for junction_path in junction_files:
        sample = junction_path.name.replace(".junctions.tsv", "")
        bam_path = bam_dir / f"{sample}.bam"
        if not bam_path.exists():
            print(f"[SKIP] {sample}: missing BAM {bam_path.name}")
            continue
        ann_count, nov_count, ann_out, nov_out, kept_junctions, skipped_junctions = export_sample(
            sample, bam_path, junction_path, output_dir, chrom_index
        )
        print(
            f"[OK] {sample}: kept_junctions={kept_junctions}, skipped_junctions={skipped_junctions}, "
            f"annotated_reads={ann_count} -> {ann_out.name}, "
            f"novel_reads={nov_count} -> {nov_out.name}"
        )


if __name__ == "__main__":
    main()
