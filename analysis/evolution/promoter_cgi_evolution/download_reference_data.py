#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Resolve and download genome FASTA + gene annotation files for promoter CGI analysis.

The script uses a curated species manifest and official Ensembl / Ensembl Genomes
FTP directory listings to discover the current filenames for each species.
"""

from __future__ import annotations

import argparse
import gzip
import re
import shutil
import sys
from contextlib import nullcontext
from pathlib import Path
from urllib.parse import urljoin
from urllib.request import Request, urlopen

import pandas as pd
from pyfaidx import Fasta


SOURCE_BASES = {
    "ensembl": "https://ftp.ensembl.org/pub/release-{release}/",
    "metazoa": "http://ftp.ensemblgenomes.org/pub/metazoa/release-{release}/",
    "fungi": "http://ftp.ensemblgenomes.org/pub/fungi/release-{release}/",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Download genomes and annotations for CGI evolution analysis.")
    parser.add_argument("--manifest", default="species_manifest.tsv", help="Species manifest TSV.")
    parser.add_argument("--data-dir", default="data", help="Directory for downloaded data.")
    parser.add_argument("--timeout", type=int, default=60, help="HTTP timeout in seconds.")
    parser.add_argument("--limit", type=int, default=0, help="Optional limit on the number of species.")
    parser.add_argument(
        "--species",
        nargs="*",
        default=[],
        help="Optional subset of species names to resolve/download.",
    )
    parser.add_argument("--dry-run", action="store_true", help="Only resolve remote URLs without downloading.")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing local files.")
    return parser.parse_args()


def fetch_links(url: str, timeout: int) -> list[str]:
    req = Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urlopen(req, timeout=timeout) as handle:
        html = handle.read().decode("utf-8", errors="ignore")
    links = re.findall(r'href="([^"]+)"', html)
    return [link for link in links if not link.startswith("?")]


def choose_first(candidates: list[str], suffixes: list[str]) -> str:
    for suffix in suffixes:
        matches = [name for name in candidates if name.endswith(suffix)]
        if matches:
            return sorted(matches)[0]
    raise FileNotFoundError(f"No remote file matched suffixes: {suffixes}")


def resolve_species_assets(row: pd.Series, timeout: int) -> dict[str, str]:
    source_db = str(row["source_db"])
    release = str(row["release"])
    species = str(row["species"])
    base = SOURCE_BASES[source_db].format(release=release)

    fasta_dir = urljoin(base, f"fasta/{species}/dna/")
    fasta_links = fetch_links(fasta_dir, timeout)
    fasta_name = choose_first(
        fasta_links,
        [
            "dna.primary_assembly.fa.gz",
            "dna.toplevel.fa.gz",
            "dna_sm.primary_assembly.fa.gz",
            "dna_sm.toplevel.fa.gz",
        ],
    )

    ann_name = ""
    ann_format = ""
    ann_dir = ""
    for subdir, suffixes, fmt in [
        (f"gtf/{species}/", [".gtf.gz"], "gtf"),
        (f"gff3/{species}/", [".gff3.gz"], "gff3"),
    ]:
        try:
            links = fetch_links(urljoin(base, subdir), timeout)
        except Exception:
            continue
        filtered = [
            name
            for name in links
            if any(name.endswith(suffix) for suffix in suffixes)
            and "abinitio" not in name
            and ".CHECKSUMS" not in name
        ]
        if filtered:
            ann_name = sorted(filtered)[0]
            ann_format = fmt
            ann_dir = urljoin(base, subdir)
            break

    if not ann_name:
        raise FileNotFoundError(f"No GTF/GFF3 annotation file found for {species}")

    return {
        "species": species,
        "source_db": source_db,
        "release": release,
        "fasta_url": urljoin(fasta_dir, fasta_name),
        "annotation_url": urljoin(ann_dir, ann_name),
        "annotation_format": ann_format,
    }


def download_file(url: str, destination: Path, timeout: int, overwrite: bool) -> None:
    if destination.exists() and not overwrite:
        return
    destination.parent.mkdir(parents=True, exist_ok=True)
    req = Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urlopen(req, timeout=timeout) as src, destination.open("wb") as dst:
        shutil.copyfileobj(src, dst)


def gunzip_file(source_gz: Path, destination: Path, overwrite: bool) -> None:
    if destination.exists() and not overwrite:
        return
    with gzip.open(source_gz, "rb") as src, destination.open("wb") as dst:
        shutil.copyfileobj(src, dst)


def rewrite_fasta_with_fixed_wrapping(source: Path, line_width: int = 60) -> None:
    temp_path = source.with_suffix(source.suffix + ".rewrapped")
    with source.open("r", encoding="utf-8", errors="ignore") as src, temp_path.open(
        "w", encoding="utf-8"
    ) as dst:
        seq_chunks: list[str] = []
        for line in src:
            if line.startswith(">"):
                if seq_chunks:
                    seq = "".join(seq_chunks)
                    for i in range(0, len(seq), line_width):
                        dst.write(seq[i : i + line_width] + "\n")
                    seq_chunks = []
                dst.write(line if line.endswith("\n") else line + "\n")
                continue
            seq_chunks.append(line.strip())
        if seq_chunks:
            seq = "".join(seq_chunks)
            for i in range(0, len(seq), line_width):
                dst.write(seq[i : i + line_width] + "\n")
    temp_path.replace(source)


def ensure_fasta_index(fasta_path: Path) -> None:
    fai_path = fasta_path.with_suffix(fasta_path.suffix + ".fai")
    if fai_path.exists():
        fai_path.unlink()
    try:
        Fasta(str(fasta_path), as_raw=True, sequence_always_upper=True)
    except Exception as exc:
        if "Line length of fasta file is not consistent" not in str(exc):
            raise
        print(f"Rewrapping FASTA before indexing: {fasta_path.name}")
        rewrite_fasta_with_fixed_wrapping(fasta_path)
        if fai_path.exists():
            fai_path.unlink()
        Fasta(str(fasta_path), as_raw=True, sequence_always_upper=True)


def main() -> None:
    args = parse_args()
    manifest_path = Path(args.manifest)
    base_dir = manifest_path.parent.resolve()
    data_dir = (base_dir / args.data_dir).resolve()
    manifest = pd.read_csv(manifest_path, sep="\t")

    if args.species:
        manifest = manifest[manifest["species"].isin(set(args.species))].copy()

    if args.limit > 0:
        manifest = manifest.head(args.limit).copy()

    resolved_rows = []
    for _, row in manifest.iterrows():
        resolved = resolve_species_assets(row, timeout=args.timeout)
        resolved_rows.append(resolved)
        print(f"Resolved {resolved['species']}")

    resolved_df = pd.DataFrame(resolved_rows)
    resolved_df.to_csv(base_dir / "resolved_downloads.tsv", sep="\t", index=False)

    if args.dry_run:
        print("Dry run only. Resolved URLs saved to: resolved_downloads.tsv")
        return

    for row in resolved_rows:
        species = row["species"]
        species_dir = data_dir / species
        raw_dir = species_dir / "raw"
        raw_dir.mkdir(parents=True, exist_ok=True)

        fasta_gz = raw_dir / Path(row["fasta_url"]).name
        annotation_gz = raw_dir / Path(row["annotation_url"]).name
        fasta_path = species_dir / f"{species}.fa"

        download_file(row["fasta_url"], fasta_gz, timeout=args.timeout, overwrite=args.overwrite)
        download_file(row["annotation_url"], annotation_gz, timeout=args.timeout, overwrite=args.overwrite)
        gunzip_file(fasta_gz, fasta_path, overwrite=args.overwrite)
        ensure_fasta_index(fasta_path)
        print(f"Downloaded and indexed {species}")

    print(f"All downloads completed under: {data_dir}")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
