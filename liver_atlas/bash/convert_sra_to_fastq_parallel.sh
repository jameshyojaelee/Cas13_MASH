#!/bin/sh
#SBATCH --job-name=parallel_fastq_conversion
#SBATCH --partition=pe2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=256G  # 16 CPUs * 16GB per CPU
#SBATCH --time=120:00:00
#SBATCH --output=parallel_fastq_%j.log
#SBATCH --error=parallel_fastq_%j.log

# Directory paths
SRA_DIR=/gpfs/commons/groups/sanjana_lab/Cas13/Liver_Atlas/sra_files
FASTQ_DIR=/gpfs/commons/groups/sanjana_lab/Cas13/Liver_Atlas/fastq
SRR_LIST=/gpfs/commons/groups/sanjana_lab/Cas13/Liver_Atlas/SRR_Acc_List.txt

# Ensure the output directory exists
mkdir -p "$FASTQ_DIR"

# Export necessary variables for GNU Parallel
export SRA_DIR FASTQ_DIR

# Use GNU Parallel to process 16 files concurrently
cat "$SRR_LIST" | parallel -j 16 \
    "fastq-dump --split-files --gzip --outdir $FASTQ_DIR $SRA_DIR/{}/{}.sra"

# Validation Step: Check if all FASTQ files are created
while read -r srr; do
    if [ ! -f "$FASTQ_DIR/${srr}_1.fastq.gz" ] || [ ! -f "$FASTQ_DIR/${srr}_2.fastq.gz" ]; then
        echo "FASTQ conversion failed for $srr"
    else
        echo "FASTQ conversion succeeded for $srr"
    fi
done < "$SRR_LIST"

