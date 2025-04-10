#!/bin/bash
#SBATCH --job-name=v3.1_human_cellranger
#SBATCH --nodes=1 
#SBATCH --cpus-per-task=4
#SBATCH --mem=180G
#SBATCH --time=120:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/Liver_Atlas/cellranger/human/logs/%x_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/Liver_Atlas/cellranger/human/logs/%x_%j.err

CELLRANGER_PATH="/gpfs/commons/home/jameslee/cellranger-9.0.0/cellranger"
TRANSCRIPTOME="/gpfs/commons/home/jameslee/reference_genome/refdata-gex-GRCh38-2024-A"
OUTPUT_BASE_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/Liver_Atlas/cellranger/human"

# Define the list of samples
SAMPLES=(
    "SRR17375057"
    "SRR17375058"
    "SRR17375051"
    "SRR17375052"
)

# Base directory for the fastq files
FASTQ_BASE_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/Liver_Atlas/fastq"

# Loop through each sample
for SAMPLE in "${SAMPLES[@]}"; do
    echo "Processing sample: ${SAMPLE}"

    # Define paths
    SAMPLE_FASTQ="${FASTQ_BASE_DIR}/${SAMPLE}"
    OUTPUT_DIR="${OUTPUT_BASE_DIR}/${SAMPLE}"

    # Run Cell Ranger for the sample
    "${CELLRANGER_PATH}" count \
        --id="${SAMPLE}" \
        --transcriptome="${TRANSCRIPTOME}" \
        --fastqs="$(dirname "${SAMPLE_FASTQ}")" \
        --sample="${SAMPLE}" \
        --create-bam=false \
        --output-dir="${OUTPUT_DIR}" \
        --nosecondary
done