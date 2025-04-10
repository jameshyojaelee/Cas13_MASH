#!/bin/sh
#SBATCH --job-name=convert_sra_to_fastq
#SBATCH --partition=pe2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --time=120:00:00
#SBATCH --output=convert_sra_to_fastq_%j.log
#SBATCH --error=convert_sra_to_fastq_%j.log

# File containing the list of SRR numbers
SRR_LIST=/gpfs/commons/groups/sanjana_lab/Cas13/Liver_Atlas/SRR_Acc_List.txt

# Directory containing SRA folders
SRA_DIR=/gpfs/commons/groups/sanjana_lab/Cas13/Liver_Atlas/sra_files

# Output directory for FASTQ files
FASTQ_DIR=/gpfs/commons/groups/sanjana_lab/Cas13/Liver_Atlas/fastq

# Ensure the output directory exists
mkdir -p "$FASTQ_DIR"

# Loop through each SRR number in the list
while read -r srr; do
    # Construct the path to the .sra file
    sra_file="$SRA_DIR/$srr/$srr.sra"

    # Check if the .sra file exists
    if [ -f "$sra_file" ]; then
        echo "Processing $srr"
        
        # Convert .sra to FASTQ
        fastq-dump --split-files --gzip --outdir "$FASTQ_DIR" "$sra_file"
    else
        echo "SRA file not found for $srr in $SRA_DIR/$srr"
    fi
done < "$SRR_LIST"

# Validate successful FASTQ conversion
while read -r srr; do
    if [ ! -f "$FASTQ_DIR/${srr}_1.fastq.gz" ] || [ ! -f "$FASTQ_DIR/${srr}_2.fastq.gz" ]; then
        echo "FASTQ conversion failed for $srr"
    else
        echo "FASTQ conversion succeeded for $srr"
    fi
done < "$SRR_LIST"
