# Cas13 MASH Project Documentation

This repository contains scripts and data for analyzing bulk RNA-seq, single-cell RNA-seq, and GWAS data related to liver diseases and Cas13 experiments.

## Project Structure

The project is organized into the following main directories:

- `bulk_RNA/`: Contains bulk RNA-seq data analysis scripts and results
- `GWAS/`: Contains GWAS data analysis results
- `liver_atlas/`: Contains single-cell RNA-seq data analysis scripts and results for liver atlas

## bulk_RNA Directory

The bulk RNA-seq analysis directory contains scripts for downloading and processing RNA-seq data from GEO for NAFLD/NASH studies.

### Scripts in `bulk_RNA/scripts/`

- `bulk_RNA.R`: Downloads and processes GEO datasets for NAFLD/NASH studies, creates directories for data storage, saves metadata to CSV, and retrieves SRA information for downloading FASTQ files.
- `download_GSE130970.py`, `download_GSE126848.py`, `download_GSE135251.py`: Python scripts to download specific GEO datasets (Govaere et al. 2020, Suppli et al. 2019, Liu et al. 2020).
- `fastq_download.py`, `fastq_download.R`, `download_fastq.R`: Scripts to download FASTQ files from SRA.
- `fastq_download.sh`: Shell script to download FASTQ files using SRA-toolkit.
- `run_fastqc.sh`: Runs FastQC on downloaded FASTQ files for quality control.
- `bulk.ipynb`: Jupyter notebook for bulk RNA-seq data analysis.

### Data and Results

- `fastq/`: Contains the raw FASTQ files for each sample
- `alignment/`: Contains gene counts for each sample
- `metadata/`: Contains metadata for each sample
- `NAFLD_NASH_RNAseq/`: Contains processed RNA-seq data for NAFLD/NASH studies
- `figures/`: Contains visualization results

### Pipeline Scripts

- `run_rnaseq_pipeline.sh`: Main script to run the nf-core/rnaseq pipeline (v3.14.0), which performs quality control, adapter trimming, alignment, and quantification.
- `nextflow.config`: Configuration file for the Nextflow pipeline.
- `samplesheet.csv`: Sample sheet for the nf-core/rnaseq pipeline.

## GWAS Directory

Contains GWAS analysis results:

- `variant_to_closest.csv`: Mapping of variants to closest genes
- `closest_tss_list.txt`: List of closest transcription start sites
- `closest_gene_list.txt`: List of closest genes to GWAS variants

## liver_atlas Directory

Contains single-cell RNA-seq data analysis for liver cell atlas.

### Bash Scripts in `liver_atlas/bash/`

- `cellranger_human.sh`, `cellranger_human2.sh`: Run Cell Ranger for human single-cell RNA-seq data
- `cellranger_mouse_v3.sh`, `cellranger_mouse_v3.1.sh`, `cellranger_mouse_v3.1_2.sh`: Run Cell Ranger for mouse single-cell RNA-seq data with different versions
- `cellranger_mouse_v3 (using library.csv).sh`, `cellranger_mouse_v3.1 (using library.csv).sh`: Run Cell Ranger using library CSV files
- `rename_fastq.sh`, `organize_fastq.sh`: Organize and rename FASTQ files
- `convert_sra_to_fastq.sh`, `convert_sra_to_fastq_parallel.sh`: Convert SRA files to FASTQ format
- `fastq_download.sh`, `fastq_download_skipped.sh`: Download FASTQ files from SRA

### Directories

- `Seurat/`: Likely contains Seurat R scripts and objects for single-cell analysis
- `notebooks/`: Jupyter notebooks for data analysis
- `data/`: Raw and processed data files
- `figures/`: Generated figures including `Western_diet_volcano_plot.pdf` and related plots
- `metadata/`: Sample metadata files

## Data Processing Pipelines

### Bulk RNA-seq Pipeline

1. Quality control with FastQC
2. Adapter trimming with Trim Galore
3. Alignment with STAR/HISAT2
4. Quantification with Salmon/featureCounts
5. Quality control reports with MultiQC

### Single-cell RNA-seq Pipeline

1. FASTQ file organization and renaming
2. Cell Ranger processing for gene expression quantification
3. Seurat analysis for clustering, differential expression, and visualization 