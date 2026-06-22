#!/usr/bin/env python3

import pandas as pd
import argparse

# =========================================================
# attribute parser
# =========================================================
def parse_attribute(attr, key):

    for item in attr.split(';'):

        item = item.strip()

        if item.startswith(key):

            return (
                item.replace(key, '')
                .replace('"', '')
                .strip()
            )

    return "Unknown"

# =========================================================
# extract first exon
# =========================================================
def extract_first_exons(gtf_file, out_prefix):

    print(f"\n读取 GTF 文件: {gtf_file}")

    # =====================================================
    # read GTF
    # =====================================================
    cols = [
        'seqname',
        'source',
        'feature',
        'start',
        'end',
        'score',
        'strand',
        'frame',
        'attribute'
    ]

    df = pd.read_csv(
        gtf_file,
        sep='\t',
        comment='#',
        header=None,
        names=cols
    )

    print(f"总记录数: {len(df):,}")

    # =====================================================
    # exon only
    # =====================================================
    exons = df[
        df['feature'] == 'exon'
    ].copy()

    print(f"Exon 数量: {len(exons):,}")

    # =====================================================
    # extract gene_id / transcript_id
    # =====================================================
    exons['gene_id'] = exons['attribute'].apply(
        lambda x: parse_attribute(x, 'gene_id')
    )

    exons['transcript_id'] = exons['attribute'].apply(
        lambda x: parse_attribute(x, 'transcript_id')
    )

    # =====================================================
    # sort
    # =====================================================
    exons = exons.sort_values(
        by=['seqname', 'start']
    )

    # =====================================================
    # gene-level first exon
    # =====================================================

    # + strand -> smallest start
    gene_plus = exons[
        exons['strand'] == '+'
    ].groupby('gene_id').head(1)

    # - strand -> largest end
    gene_minus = exons[
        exons['strand'] == '-'
    ].groupby('gene_id').tail(1)

    gene_first = pd.concat([
        gene_plus,
        gene_minus
    ])

    # =====================================================
    # transcript-level first exon
    # =====================================================

    tx_plus = exons[
        exons['strand'] == '+'
    ].groupby('transcript_id').head(1)

    tx_minus = exons[
        exons['strand'] == '-'
    ].groupby('transcript_id').tail(1)

    tx_first = pd.concat([
        tx_plus,
        tx_minus
    ])

    # =====================================================
    # BED conversion
    # BED start = 0-based
    # =====================================================

    # gene BED
    gene_bed = gene_first[
        [
            'seqname',
            'start',
            'end',
            'gene_id',
            'score',
            'strand'
        ]
    ].copy()

    gene_bed['start'] = gene_bed['start'] - 1

    # transcript BED
    tx_bed = tx_first[
        [
            'seqname',
            'start',
            'end',
            'transcript_id',
            'score',
            'strand'
        ]
    ].copy()

    tx_bed['start'] = tx_bed['start'] - 1

    # =====================================================
    # output
    # =====================================================
    gene_out = f"{out_prefix}.gene_id.bed"
    tx_out   = f"{out_prefix}.transcript_id.bed"

    gene_bed.to_csv(
        gene_out,
        sep='\t',
        header=False,
        index=False
    )

    tx_bed.to_csv(
        tx_out,
        sep='\t',
        header=False,
        index=False
    )

    # =====================================================
    # summary
    # =====================================================
    print("\n===== 输出完成 =====")

    print(
        f"Gene-level first exon BED : "
        f"{gene_out}"
    )

    print(
        f"Transcript-level first exon BED : "
        f"{tx_out}"
    )

    print(
        f"\nGene 数量: "
        f"{len(gene_bed):,}"
    )

    print(
        f"Transcript 数量: "
        f"{len(tx_bed):,}"
    )

# =========================================================
# main
# =========================================================
if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="""
Extract first exons from GTF.

Outputs:
1. gene_id BED
2. transcript_id BED
"""
    )

    parser.add_argument(
        '-g',
        '--gtf',
        required=True,
        help='Input GTF file'
    )

    parser.add_argument(
        '-o',
        '--out',
        default='first_exons',
        help='Output prefix'
    )

    args = parser.parse_args()

    extract_first_exons(
        args.gtf,
        args.out
    )