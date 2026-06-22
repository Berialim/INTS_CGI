# CGI-INTS11 Analysis Code

This repository collects the core analysis code used for a multi-layer study of transcription elongation, promoter CpG island architecture, U1-associated regulation, splicing fidelity, and evolutionary conservation. The code was originally developed across several working directories and has been reorganized here into a cleaner GitHub-ready project structure.

The main purpose of this repository is to preserve the runnable analytical logic behind the methods, especially for:

- TT-seq and PRO-seq data processing
- promoter CpG island and U1-window analyses
- cross-species protein conservation analysis
- promoter CGI evolution analysis
- splicing error quantification from TT-seq alignments
- indel/eQTL analysis in U1-CGI-related genomic regions

This is a code-centric repository. Large raw datasets, most intermediate files, and large final outputs are intentionally not included.

## Scientific scope

The analyses in this repository revolve around a common biological question: how promoter architecture, especially CpG island extent and its spatial relationship to predicted U1-recognition sites, relates to transcriptional output, splicing fidelity, regulatory perturbation, and evolutionary conservation.

The repository integrates several analysis layers:

1. Primary sequencing processing
2. Promoter architecture quantification
3. U1-CGI distance and CpG-related feature analysis
4. Splicing error measurement under perturbation conditions
5. GTEx indel/eQTL effect analysis in CpG-rich regulatory contexts
6. Cross-species protein conservation and promoter CGI evolution

## Methods-to-code map

The repository was organized to match the major method sections.

### TT-seq and PRO-seq data processing

Mapped primarily to:

- `analysis/sequencing_processing/TT-seq.sh`
- `analysis/sequencing_processing/PRO-seq.sh`
- `analysis/sequencing_processing/featureCounts.sh`

These scripts cover:

- adapter trimming
- basic FASTQ cleanup
- rRNA depletion handling
- alignment with `STAR`
- filtering for unique/high-confidence alignments
- BAM sorting and indexing
- spike-in-based scaling for TT-seq
- BigWig generation
- featureCounts-based quantification

### ChIP-seq-like processing

Mapped primarily to:

- `analysis/sequencing_processing/chipseq_pe.sh`

This script preserves the original paired-end ChIP-like processing logic using:

- `trim_galore`
- `bowtie2`
- `samtools`
- `picard`
- `bamCoverage`

### NFR / TSS-associated CpG island / promoter architecture analysis

Mapped primarily to:

- `analysis/promoter_architecture/scripts/1stexon.py`
- `analysis/promoter_architecture/scripts/CGIcovTSS_sep.py`
- `analysis/promoter_architecture/scripts/CpG_togeneP_lengthdistribution.py`
- `analysis/promoter_architecture/scripts/u1tocpg.py`
- `analysis/promoter_architecture/scripts/u1tocpg_cli.py`
- `analysis/promoter_architecture/scripts/CpG_togeneP_lengthdistribution_cli.py`

These scripts support:

- first-exon extraction from GTF
- TSS-overlapping CGI identification
- generation of TSS-to-CGI and antisense-oriented CGI intervals
- CGI length distribution summaries
- U1-to-CpG distance calculations

### U1-CGI window, expression stratification, and downstream gene set interpretation

Mapped primarily to:

- `analysis/promoter_architecture/scripts/gsea_human.py`

This script supports gene set enrichment analysis on ranked genes after transcription-related stratification or perturbation analyses.

### Evolutionary conservation analysis

Mapped primarily to:

- `analysis/evolution/protein_conservation/extract_human_related_homologies.py`
- `analysis/evolution/protein_conservation/evolution_full_pipeline_coevolution_autobroad_reanalysis_refined_taxonomy.py`
- `analysis/evolution/taxonomy_groups_refined.py`

These scripts implement:

- extraction of human-centered orthology records from Ensembl Compara
- normalization of human-left / human-right rows
- weighted conservation scoring
- complex-level conservation summaries
- phylogenetic grouping and visualization
- species-level and broad-group correlation analyses

### Promoter CGI evolution analysis

Mapped primarily to:

- `analysis/evolution/promoter_cgi_evolution/download_reference_data.py`
- `analysis/evolution/promoter_cgi_evolution/promoter_cgi_evolution_pipeline.py`
- `analysis/evolution/promoter_cgi_evolution/plot_promoter_cgi_vs_u1_integrator.py`
- `analysis/evolution/promoter_cgi_evolution/plot_tss_cgi_coverage_four_species.py`

These scripts support:

- downloading genome and annotation references for selected species
- promoter definition around TSSs
- scanning for TSS-overlapping CGI-like windows
- matched background comparisons
- CGI conservation summaries
- integration with protein conservation profiles

### Splicing error quantification and gene-level comparison

Mapped primarily to:

- `analysis/splicing_error/splicingerro.py`
- `analysis/splicing_error/plot_ctrl_iaa_u1_boxplot.py`
- `analysis/splicing_error/plot_splicing_error_volcano.py`
- `analysis/splicing_error/export_junctions_bedpe.py`
- `analysis/splicing_error/export_junction_support_bams.py`
- `analysis/splicing_error/merge_condition_support_bams.py`

These scripts cover:

- junction extraction from BAM alignments
- annotation-aware classification of annotated versus novel junctions
- gene-level novel junction ratio and novel read ratio calculation
- condition-wise statistical comparison
- IGV-ready junction track export
- support BAM export for visual validation

### Indel/eQTL analysis in U1-CGI-associated regions

Mapped primarily to:

- `analysis/indel_eqtl/insert_del_eachtissue.py`
- `analysis/indel_eqtl/insert_2.py`
- `analysis/indel_eqtl/CGIinsert_slop_1.py`
- `analysis/indel_eqtl/summary_inster.py`

These scripts preserve the logic for:

- parsing GTEx v11 significant variant-gene pairs
- restricting to protein-coding genes when appropriate
- extracting insertion and deletion events based on allele-length differences
- intersecting variants with U1-CGI or related BED intervals
- summarizing effect sizes globally and by tissue

## Repository layout

```text
.
├── README.md
├── requirements.txt
├── .gitignore
└── analysis
    ├── sequencing_processing
    │   ├── README.md
    │   ├── TT-seq.sh
    │   ├── PRO-seq.sh
    │   ├── chipseq_pe.sh
    │   └── featureCounts.sh
    ├── promoter_architecture
    │   ├── README.md
    │   └── scripts
    │       ├── 1stexon.py
    │       ├── CGIcovTSS_sep.py
    │       ├── CpG_togeneP_lengthdistribution.py
    │       ├── CpG_togeneP_lengthdistribution_cli.py
    │       ├── gsea_human.py
    │       ├── u1tocpg.py
    │       └── u1tocpg_cli.py
    ├── evolution
    │   ├── README.md
    │   ├── taxonomy_groups_refined.py
    │   ├── protein_conservation
    │   │   ├── evolution_full_pipeline_coevolution_autobroad_reanalysis_refined_taxonomy.py
    │   │   ├── extract_human_related_homologies.py
    │   │   └── taxonomy_groups_refined.py
    │   └── promoter_cgi_evolution
    │       ├── download_reference_data.py
    │       ├── promoter_cgi_evolution_pipeline.py
    │       ├── plot_promoter_cgi_vs_u1_integrator.py
    │       └── plot_tss_cgi_coverage_four_species.py
    ├── splicing_error
    │   ├── README.md
    │   ├── splicingerro.py
    │   ├── plot_ctrl_iaa_u1_boxplot.py
    │   ├── plot_splicing_error_volcano.py
    │   ├── export_junctions_bedpe.py
    │   ├── export_junction_support_bams.py
    │   └── merge_condition_support_bams.py
    └── indel_eqtl
        ├── README.md
        ├── insert_del_eachtissue.py
        ├── insert_2.py
        ├── CGIinsert_slop_1.py
        └── summary_inster.py
```

## Directory guide

### `analysis/sequencing_processing`

This directory contains shell-based templates for primary sequencing processing. These scripts largely preserve the original lab execution style and are best viewed as workflow templates rather than turn-key portable pipelines.

Typical tasks supported here:

- single-end TT-seq preprocessing and alignment
- single-end PRO-seq preprocessing and alignment
- paired-end ChIP-like data alignment and BigWig generation
- featureCounts-based expression quantification

Typical assumptions:

- input FASTQ files are placed in an `fq/` directory
- reference indices already exist
- command-line tools are installed and available
- sample naming follows the original conventions expected by the scripts

### `analysis/promoter_architecture`

This directory contains scripts for constructing and analyzing promoter-related interval definitions.

Representative tasks:

- extract first exons from a GTF
- identify genes with TSS-overlapping CpG islands
- derive promoter-side and antisense-side CGI intervals
- summarize length distributions
- measure U1-to-CGI distances
- perform GSEA-style downstream interpretation

Two helper CLI copies were added during consolidation:

- `u1tocpg_cli.py`
- `CpG_togeneP_lengthdistribution_cli.py`

These provide simpler argument-based entry points for common tasks.

### `analysis/evolution`

This directory contains the evolutionary modules and a shared taxonomy helper.

It is split into:

- `protein_conservation/`
- `promoter_cgi_evolution/`

The top-level `taxonomy_groups_refined.py` defines the shared phylogenetic grouping scheme used across these analyses.

### `analysis/splicing_error`

This directory contains the gene-level splicing fidelity workflow.

The intended progression is:

1. quantify per-sample junction statistics from BAM files
2. aggregate results at the gene level
3. compare gene sets or conditions
4. export tracks or BAM subsets for IGV inspection

### `analysis/indel_eqtl`

This directory contains the GTEx indel/eQTL analysis code. The scripts were preserved from the working analysis directory with minimal restructuring so that the original analysis logic remains visible.

## Typical input files

The repository does not ship data, but most analyses expect some combination of the following:

- FASTQ files
- BAM files and BAM indexes
- STAR or Bowtie/Bowtie2 reference indices
- GTF or GFF3 annotations
- BED6 region files
- BigWig files
- Ensembl Compara TSV tables
- GTEx v11 `*.signif_pairs.parquet` files
- GTEx `*.eGenes.txt.gz` files
- sample metadata tables such as `database.csv`

## Typical outputs

Depending on the module, outputs may include:

- BAM and indexed BAM files
- BigWig tracks
- count matrices
- BED / BEDPE / BED12 interval files
- per-gene summary tables
- species-level or tissue-level summary tables
- boxplots, KDE plots, volcano-like plots, heatmaps, phylogenetic plots, and scatter/regression summaries

## Software dependencies

The project uses both Python and common command-line genomics tools.

### Python packages

The main Python dependencies are listed in `requirements.txt` and commonly include:

- `pandas`
- `numpy`
- `matplotlib`
- `seaborn`
- `scipy`
- `pysam`
- `pyarrow`
- `gseapy`
- `statannotations`
- `pybedtools`
- `sympy`

One preserved exploratory script also imports `torch`, so it remains listed in the requirements file.

### Command-line tools

Several shell pipelines or Python wrappers assume the availability of:

- `trim_galore`
- `STAR`
- `bowtie`
- `bowtie2`
- `samtools`
- `picard`
- `deepTools`
- `multiqc`
- `featureCounts`
- `bedtools`

## Installation

### Python environment

Create a Python environment and install the listed packages:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Depending on your platform, `pybedtools`, `pysam`, and `pyarrow` may require additional system libraries or conda-based installation.

### External genomics tools

Install the command-line tools separately using your preferred environment manager, conda environment, module system, or cluster software stack.

## Recommended usage pattern

Because the original code was developed inside an active research workflow, the most practical way to reuse this repository is:

1. Prepare your own working directories and references.
2. Start from the module that matches your analysis goal.
3. Read that module's `README.md`.
4. Adjust input paths or command-line arguments as needed.
5. Treat the shell scripts as templates for batch execution on HPC or a lab server.
6. Treat the Python scripts as analysis entry points that may still need local path cleanup.

## Example workflows

### A. Splicing error workflow

Example logic:

1. Run `splicingerro.py` on one or more indexed BAM files with a reference GTF.
2. Generate per-sample junction summary tables.
3. Compare `CTRL`, `IAA`, or `U1` conditions using `plot_ctrl_iaa_u1_boxplot.py`.
4. Export BEDPE or support BAMs for IGV inspection.

Example command:

```bash
python analysis/splicing_error/splicingerro.py \
  -b sample1.bam sample2.bam \
  -g annotation.gtf \
  -o results/splicing_error
```

### B. Promoter architecture workflow

Example logic:

1. Extract first exons from a GTF.
2. Intersect TSS positions with CpG islands.
3. Derive TSS-to-CGI and antisense-oriented intervals.
4. Quantify U1-to-CGI relationships.

Example command:

```bash
python analysis/promoter_architecture/scripts/1stexon.py \
  -g annotation.gtf \
  -o first_exons
```

### C. Protein conservation workflow

Example logic:

1. Extract human-related orthology rows from a full Compara dump.
2. Run the refined taxonomy coevolution pipeline.
3. Generate system-level conservation tables and summary plots.

Example command:

```bash
python analysis/evolution/protein_conservation/extract_human_related_homologies.py \
  --input Compara.115.protein_default.homologies.tsv \
  --output Compara.115.protein_default.human_related.normalized.tsv.gz
```

### D. Promoter CGI evolution workflow

Example logic:

1. Download species references.
2. Run the promoter CGI evolution pipeline.
3. Integrate the results with protein conservation summaries.

Example command:

```bash
python analysis/evolution/promoter_cgi_evolution/download_reference_data.py --dry-run
python analysis/evolution/promoter_cgi_evolution/promoter_cgi_evolution_pipeline.py
```

### E. GTEx indel/eQTL workflow

Example logic:

1. Prepare a BED file describing U1-CGI or related regions.
2. Collect GTEx `signif_pairs.parquet` and matching `eGenes` files.
3. Run the indel/eQTL summary scripts.
4. Compare insertion and deletion effect-size distributions across tissues.

## Reproducibility notes

This repository is a cleaned aggregation of previously distributed project code. As a result:

- some scripts are already portable and argument-driven
- some scripts still preserve the original lab assumptions
- some modules are cleaner than others because the original code quality varied by analysis phase

This is expected. The consolidation goal was to preserve analytical logic first, then improve portability where feasible without changing the scientific behavior of the code.

## What was intentionally not included

The following were intentionally left out:

- large raw sequencing datasets
- large GTEx or Compara data files
- most intermediate analysis products
- PDF figure outputs
- temporary files
- platform-specific resource-fork files

## Current limitations

Users of this repository should be aware of a few practical limitations:

- Several shell pipelines still assume a specific folder structure and local environment.
- Some preserved Python scripts still contain hard-coded paths and may need editing before direct reuse.
- Not every module has been converted into a fully standardized command-line package.
- The repository does not yet include unit tests or example datasets.

## Recommended future improvements

If you want to make this repository more fully public-facing and easier for other groups to reuse, the most valuable next steps would be:

1. Replace all remaining hard-coded paths with command-line arguments or config files.
2. Add a small example dataset for each major module.
3. Add exact input/output schema documentation per script.
4. Standardize plotting scripts to use one argument style.
5. Convert the shell workflows into Snakemake or Nextflow.
6. Add environment files such as `environment.yml`.
7. Add validation notebooks or example result snapshots.

## Citation and attribution

If you use or adapt this code, please cite the associated study or repository release as appropriate. Since this repository was assembled from active research code, attribution should follow the publication or project context in which the analysis was developed.

## Contact

For reuse, adaptation, or clarification, the most useful way to approach this repository is to start from the analysis module closest to your own question and review the script-level assumptions before execution.
