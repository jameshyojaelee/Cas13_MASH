#!/bin/bash

# Script to regenerate STAR index with version 2.7.11b
# Author: James Lee
# Date: 2025-04-10

# Load modules
module purge
module load star/2.7.10b-GCC-11.3.0

# Input files
REF_DIR="/gpfs/commons/home/jameslee/reference_genome/refdata-gex-GRCh38-2024-A"
GENOME="${REF_DIR}/fasta/genome.fa"
GTF="${REF_DIR}/genes/genes.gtf.gz"

# Output directory
NEW_STAR_INDEX="${REF_DIR}/star_2.7.10b"

# Create output directory
mkdir -p ${NEW_STAR_INDEX}

# SLURM settings
#SBATCH --account=nslab
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=6:00:00
#SBATCH --job-name=STAR_index_gen
#SBATCH --output=%x-%j.out

# Generate STAR index
echo "Generating STAR index with version 2.7.10b..."
STAR --runMode genomeGenerate \
     --genomeDir ${NEW_STAR_INDEX} \
     --genomeFastaFiles ${GENOME} \
     --sjdbGTFfile ${GTF} \
     --runThreadN 16 \
     --sjdbOverhang 100

echo "STAR index generation completed. New index in ${NEW_STAR_INDEX}" 