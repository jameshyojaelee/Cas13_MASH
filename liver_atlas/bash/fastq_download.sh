#!/bin/bash
#SBATCH --job-name=fastq_download
#SBATCH --partition=pe2
#SBATCH --nodes=1                                                              
#SBATCH --cpus-per-task=4                                                                  
#SBATCH --mem=120G                                                                       
#SBATCH --time=120:00:00
#SBATCH --output=fastq_download_%j.log
#SBATCH --error=fastq_download_%j.log
# Author: James Lee
# Description: Script to download FASTQ files for Liver Atlas project
# Date: $(date +%Y-%m-%d)

prefetch --option-file SRR_Acc_List.txt --output-directory ./sra_files

cd ./sra_files

# Convert .sra files to FASTQ in the current directory
while read -r srr; do
    fastq-dump --split-files --gzip --outdir /gpfs/commons/home/jameslee/Cas13/Liver_Atlas/fastq_files $srr
done < ../SRR_Acc_List.txt

# Validate successful FASTQ downloads
while read -r srr; do
    if [ ! -f "/gpfs/commons/home/jameslee/Cas13/Liver_Atlas/fastq_files/${srr}_1.fastq.gz" ] || [ ! -f "/gpfs/commons/home/jameslee/Cas13/Liver_Atlas/fastq_files/${srr}_2.fastq.gz" ]; then
        echo "Download failed for $srr"
    fi
done

# Return to the original directory
cd ..
