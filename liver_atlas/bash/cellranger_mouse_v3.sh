#!/bin/bash
#SBATCH --job-name=v3_mouse_cellranger
#SBATCH --nodes=1 
#SBATCH --cpus-per-task=8
#SBATCH --mem=220G 
#SBATCH --time=120:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/Liver_Atlas/cellranger/mouse/logs/%x_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/Liver_Atlas/cellranger/mouse/logs/%x_%j.err

CELLRANGER_PATH="/gpfs/commons/home/jameslee/cellranger-9.0.0/cellranger"
TRANSCRIPTOME="/gpfs/commons/home/jameslee/reference_genome/refdata-gex-GRCm39-2024-A"
OUTPUT_BASE_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/Liver_Atlas/cellranger/mouse"

# Create logs directory if it doesn't exist
mkdir -p ${OUTPUT_BASE_DIR}/logs

# Define the list of samples
SAMPLES=(
    "SRR17374925"
    "SRR17374977"
    "SRR17374981"
    "SRR17374985"
    "SRR17374989"
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
    ${CELLRANGER_PATH} count \
        --id=${SAMPLE} \
        --transcriptome=${TRANSCRIPTOME} \
        --fastqs=$(dirname ${SAMPLE_FASTQ}) \
        --sample=${SAMPLE} \
        --chemistry=SC3Pv3 \
        --create-bam=false \
        --output-dir=${OUTPUT_DIR} \
        --nosecondary
done