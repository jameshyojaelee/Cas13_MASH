#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH --output=fastqc_%A_%a.out
#SBATCH --error=fastqc_%A_%a.err
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jameslee@nygenome.org

module load fastqc/0.12.1

FASTQ_DIR="/gpfs/commons/home/jameslee/Cas13/bulk_RNA/fastq"
OUT_DIR="/gpfs/commons/home/jameslee/Cas13/bulk_RNA/fastqc_results"
mkdir -p "$OUT_DIR"

# Create an array of fastq files
files=("$FASTQ_DIR"/*.fastq)
# Select the file corresponding to the array index
file="${files[$SLURM_ARRAY_TASK_ID]}"
echo "Processing file: $file"

# Run FastQC on the selected file
fastqc -t 4 -o "$OUT_DIR" "$file"
