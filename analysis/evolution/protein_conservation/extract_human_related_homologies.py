#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Extract all human-related orthology rows from the full Ensembl Compara table.

This script keeps rows where:
- species == homo_sapiens, or
- homology_species == homo_sapiens

Rows from the second case are normalized so that the output always uses
`species == homo_sapiens` on the left side.
"""

from __future__ import annotations

import argparse
import csv
import gzip
from collections import Counter
from pathlib import Path


REFERENCE_SPECIES = "homo_sapiens"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Extract human-related homology rows from a full Compara dump.")
    parser.add_argument(
        "--input",
        default="Compara.115.protein_default.homologies.tsv",
        help="Input full Compara TSV.",
    )
    parser.add_argument(
        "--output",
        default="Compara.115.protein_default.human_related.normalized.tsv.gz",
        help="Output normalized TSV.gz.",
    )
    parser.add_argument(
        "--report",
        default="human_related_extraction_report.tsv",
        help="Output extraction summary TSV.",
    )
    parser.add_argument(
        "--progress-every",
        type=int,
        default=2_000_000,
        help="Progress log interval in rows.",
    )
    return parser.parse_args()


def normalize_row(row: dict[str, str]) -> dict[str, str]:
    if row["species"] == REFERENCE_SPECIES:
        return row

    # Human is on the right; swap columns so output is always human-left.
    normalized = row.copy()
    normalized["gene_stable_id"], normalized["homology_gene_stable_id"] = (
        row["homology_gene_stable_id"],
        row["gene_stable_id"],
    )
    normalized["protein_stable_id"], normalized["homology_protein_stable_id"] = (
        row["homology_protein_stable_id"],
        row["protein_stable_id"],
    )
    normalized["species"], normalized["homology_species"] = (
        REFERENCE_SPECIES,
        row["species"],
    )
    normalized["identity"], normalized["homology_identity"] = (
        row["homology_identity"],
        row["identity"],
    )
    return normalized


def main() -> None:
    args = parse_args()
    input_path = Path(args.input)
    output_path = Path(args.output)
    report_path = Path(args.report)

    kept_from_left = 0
    kept_from_right = 0
    total_rows = 0
    kept_rows = 0
    partner_counts = Counter()

    with input_path.open("r", encoding="utf-8", newline="") as src, gzip.open(
        output_path, "wt", encoding="utf-8", newline=""
    ) as dst:
        reader = csv.DictReader(src, delimiter="\t")
        fieldnames = reader.fieldnames
        if not fieldnames:
            raise ValueError("Input file is missing a header.")
        writer = csv.DictWriter(dst, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()

        for row in reader:
            total_rows += 1

            if "ortholog" not in (row.get("homology_type") or ""):
                if args.progress_every and total_rows % args.progress_every == 0:
                    print(f"Scanned {total_rows:,} rows; kept {kept_rows:,}")
                continue

            left_is_human = row.get("species") == REFERENCE_SPECIES
            right_is_human = row.get("homology_species") == REFERENCE_SPECIES

            if not (left_is_human or right_is_human):
                if args.progress_every and total_rows % args.progress_every == 0:
                    print(f"Scanned {total_rows:,} rows; kept {kept_rows:,}")
                continue

            normalized = normalize_row(row) if right_is_human and not left_is_human else row
            writer.writerow(normalized)
            kept_rows += 1

            if left_is_human:
                kept_from_left += 1
                partner_counts[row["homology_species"]] += 1
            else:
                kept_from_right += 1
                partner_counts[row["species"]] += 1

            if args.progress_every and total_rows % args.progress_every == 0:
                print(f"Scanned {total_rows:,} rows; kept {kept_rows:,}")

    report_rows = [
        {"metric": "total_rows_scanned", "value": total_rows},
        {"metric": "kept_rows_total", "value": kept_rows},
        {"metric": "kept_rows_human_left", "value": kept_from_left},
        {"metric": "kept_rows_human_right", "value": kept_from_right},
        {"metric": "n_partner_species", "value": len(partner_counts)},
    ]
    top_partners = [
        {"metric": f"partner_species::{species}", "value": count}
        for species, count in partner_counts.most_common()
    ]

    with report_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["metric", "value"], delimiter="\t")
        writer.writeheader()
        writer.writerows(report_rows + top_partners)

    print(f"Finished. Scanned {total_rows:,} rows and kept {kept_rows:,}.")
    print(f"Output written to: {output_path}")
    print(f"Report written to: {report_path}")


if __name__ == "__main__":
    main()
