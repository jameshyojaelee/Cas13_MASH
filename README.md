# Cas13 MASH Project Documentation

This repository contains scripts and data for analyzing bulk RNA-seq, single-cell RNA-seq, and GWAS data related to liver diseases.

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


*Bulk RNA-seq data*

This directory contains the bulk RNA-seq data for the Cas13 experiments. The data is organized by the following subdirectories:

- `fastq/`: contains the raw fastq files for each sample
- `alignment/`: contains the gene counts for each sample
- `metadata/`: contains the metadata for each sample
- `differential_expression/`: contains the differential expression analysis results for each sample
- `plots/`: contains the plots generated from the differential expression analysis
- `tables/`: contains the tables generated from the differential expression analysis
- `scripts/`: contains the scripts used to generate the differential expression analysis results
- `results/`: contains the results of the differential expression analysis
- `README.md`: this file

### Data description
  "GSE135251", # Liu et al. (2020) - 206 snap-frozen biopsy samples from 206 patients diagnosed with NAFLD in France, Germany, Italy, and the UK and enrolled in the European NAFLD Registry
  Patient samples were grouped: NAFL (n = 51) and NASH with fibrosis stages of F0/1 (n = 34), F2 (n = 53), F3 (n = 54) and F4 (n = 14)
  - NAFLD https://www.science.org/doi/10.1126/scitranslmed.aba4448?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200pubmed
  - NASH limits anti-tumor surveillance in HCC https://www.nature.com/articles/s41586-021-03362-0
  - Steatohepatits https://www.cell.com/cell-reports-medicine/fulltext/S2666-3791(24)00591-3?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS2666379124005913%3Fshowall%3Dtrue

  "GSE130970", # Suppli et al. (2019) - liver biopsies obtained from healthy normal weight (n=14) and obese (n=12) individuals, NAFL (n=15) and NASH (n=16) patients
  - Transcriptomic Profiling Across NAFLD Spectrum - https://journals.physiology.org/doi/full/10.1152/ajpgi.00358.2018?rfr_dat=cr_pub++0pubmed&url_ver=Z39.88-2003&rfr_id=ori%3Arid%3Acrossref.org 
  
  "GSE126848", # Hoang SA et al. (2019) - This dataset contains transcriptomic profiles of 78 distinct human liver biopsies. Of these, 6 are histologically normal, and 72 cover the full spectrum of nonalcoholic fatty liver disease
  - Transcriptomic Profiling of Human Liver Biopsies - https://www.nature.com/articles/s41598-019-48746-5 

### Data processing with nf-core/rnaseq pipeline

The bulk RNA-seq data is processed using the nf-core/rnaseq pipeline (v3.14.0). This pipeline performs:

1. Quality control with FastQC
2. Adapter trimming with Trim Galore
3. Alignment with STAR/HISAT2
4. Quantification with Salmon/featureCounts
5. Quality control reports with MultiQC

#### Running the pipeline

1. Ensure Nextflow is installed:
   ```
   curl -s https://get.nextflow.io | bash
   ```

2. Run nf-core/rnaseq pipeline:
   ```
   nextflow run nf-core/rnaseq -r 3.14.0 \
   -profile <profile> \
   --input samplesheet.csv \
   --outdir results/nf-core_rnaseq \
   --fasta <path_to_genome_fasta> \
   --gtf <path_to_genome_gtf> \
   --star_index <path_to_star_index> \
   -resume
   ```

See the `run_rnaseq_pipeline.sh` script for full implementation details.

### Data generation pipeline
Bulk RNA-seq data was generated using the following pipeline:

1. Quality control of raw fastq files using FastQC
2. Trimming of raw fastq files using Trimmomatic
3. Alignment of trimmed fastq files to the reference genome using STAR
4. Quantification of gene counts using featureCounts
5. Differential expression analysis using DESeq2
6. Visualization of differential expression analysis results using ggplot2

