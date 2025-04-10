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

