# Indel eQTL Analysis

This directory contains the GTEx v11 indel/eQTL scripts used to analyze insertion and deletion effects in U1-CGI-related genomic regions.

Main preserved scripts:

- `insert_del_eachtissue.py`
- `insert_2.py`
- `CGIinsert_slop_1.py`
- `summary_inster.py`

These scripts were consolidated from a working analysis directory and kept largely intact so the original logic is preserved.

Inputs typically include:

- region BED6 files defining U1-CGI or related intervals
- GTEx `*.signif_pairs.parquet` files
- GTEx `*.eGenes.txt.gz` annotation files

Because some preserved scripts still contain original hard-coded paths, treat them as analysis templates and update path arguments or variables before running.
