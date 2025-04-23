#!/bin/bash

# Full run script for the bulk RNA-seq pipeline on GSE130970
# This script runs the pipeline with all samples using STAR

# Load Nextflow module
module purge
module load Nextflow/24.04.2

# Set paths to required files
SAMPLESHEET="/gpfs/commons/home/jameslee/Cas13/bulk_RNA/samplesheets/GSE130970_samplesheet.csv"

# Reference genome files
REF_DIR="/gpfs/commons/home/jameslee/reference_genome/refdata-gex-GRCh38-2024-A"
GENOME="${REF_DIR}/fasta/genome.fa"
GTF="${REF_DIR}/genes/genes.gtf.gz"
STAR_INDEX="${REF_DIR}/star_2.7.10b"  # Using the STAR index compatible with STAR 2.7.10b

# Set output and work directories
OUTDIR="${HOME}/Cas13/bulk_RNA/RNAseq_results_GSE130970"

# Create output directory
mkdir -p $OUTDIR

# Run the pipeline with all samples
echo "Starting full run with all samples using STAR..."
nextflow run main.nf \
    --input ${SAMPLESHEET} \
    --genome ${GENOME} \
    --gtf ${GTF} \
    --star_index ${STAR_INDEX} \
    --outdir ${OUTDIR} \
    -profile slurm_no_singularity \
    -resume

echo "Full run completed. Check results in ${OUTDIR}"