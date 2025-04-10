#!/bin/bash

# Configure SLURM
#SBATCH --account=nslab
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=6:00:00
#SBATCH --job-name=star_index_fix
#SBATCH --output=star_index_fix-%j.out

# Set up environment
module purge
module load star/2.7.10b-GCC-11.3.0

# Set paths
REF_DIR="/gpfs/commons/home/jameslee/reference_genome/refdata-gex-GRCh38-2024-A"
GENOME="${REF_DIR}/fasta/genome.fa"
GTF_GZ="${REF_DIR}/genes/genes.gtf.gz"
GTF="${REF_DIR}/genes/genes.gtf"
NEW_INDEX="${REF_DIR}/star_2.7.10b"

# Create output directory
mkdir -p $NEW_INDEX

# Decompress the GTF file
echo "Decompressing GTF file at $(date)"
gunzip -c $GTF_GZ > $GTF

# Check if the GTF file contains exon lines
echo "Checking if GTF file contains exon annotations..."
grep -m 5 "exon" $GTF

# Generate STAR index
echo "Starting STAR index generation at $(date)"
STAR --runMode genomeGenerate \
     --genomeDir $NEW_INDEX \
     --genomeFastaFiles $GENOME \
     --sjdbGTFfile $GTF \
     --runThreadN 16 \
     --sjdbOverhang 100

echo "Completed STAR index generation at $(date)"

# Clean up
echo "Cleaning up..."
rm -f $GTF 