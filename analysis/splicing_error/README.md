# Splicing Error Analysis

This directory contains scripts for gene-level splicing error quantification and follow-up visualization.

Main workflow:

1. Run `splicingerro.py` on indexed BAM files plus a reference GTF.
2. Use `plot_ctrl_iaa_u1_boxplot.py` to compare splicing error ratios across conditions and gene classes.
3. Use `plot_splicing_error_volcano.py` for gene-level differential splicing-error summaries.
4. Use `export_junctions_bedpe.py`, `export_junction_support_bams.py`, and `merge_condition_support_bams.py` for IGV-ready inspection outputs.

Expected inputs:

- indexed BAM files
- GTF annotation
- `database.csv` with at least `name` and `condition`
- optional BED files defining active genes and CGI-related classes
