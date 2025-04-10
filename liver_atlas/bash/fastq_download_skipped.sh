#!/bin/bash
#SBATCH --job-name=download_skipped_files
#SBATCH --partition=pe2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=120G
#SBATCH --time=120:00:00
#SBATCH --output=download_skipped_%j.log
#SBATCH --error=download_skipped_%j.log

# Download skipped files with an increased max size
prefetch --option-file skipped_srr_list.txt --max-size 50G --output-directory ./sra_files
