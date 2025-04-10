#!/bin/bash
#SBATCH --job-name=fastq_download
#SBATCH --output=fastq_download_%A_%a.out
#SBATCH --error=fastq_download_%A_%a.err
#SBATCH --time=90:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --array=0-2
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jameslee@nygenome.org

source /gpfs/commons/home/jameslee/miniconda3/etc/profile.d/conda.sh
conda activate RNA
module load SRA-Toolkit/3.2.0-gompi-2023b

# Define array of script names
scripts=(
  "/gpfs/commons/home/jameslee/Cas13/bulk_RNA/scripts/download_GSE135251.py"
  "/gpfs/commons/home/jameslee/Cas13/bulk_RNA/scripts/download_GSE130970.py"
  "/gpfs/commons/home/jameslee/Cas13/bulk_RNA/scripts/download_GSE126848.py"
)

# Run the selected Python script as an executable (ensure scripts are executable: chmod +x)
"${scripts[$SLURM_ARRAY_TASK_ID]}"