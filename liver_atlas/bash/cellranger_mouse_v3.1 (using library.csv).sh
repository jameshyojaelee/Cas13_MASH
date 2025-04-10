#!/bin/bash
#SBATCH --job-name=mouse_Cellranger
#SBATCH --partition=pe2
#SBATCH --nodes=1 
#SBATCH --cpus-per-task=8 
#SBATCH --mem=220G 
#SBATCH --time=120:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/Liver_Atlas/cellranger/mouse/logs/%x_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/Liver_Atlas/cellranger/mouse/logs/%x_%j.err

CELLRANGER_PATH="/gpfs/commons/home/jameslee/cellranger-9.0.0/cellranger"
TRANSCRIPTOME="/gpfs/commons/home/jameslee/reference_genome/refdata-gex-GRCm39-2024-A"
OUTPUT_BASE_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/Liver_Atlas/cellranger/mouse"
LIBRARY_CSV="/gpfs/commons/home/jameslee/Cas13/liver_atlas/metadata/mouse_library_v3.1.csv"

# Create logs directory if it doesn't exist
# mkdir -p ${OUTPUT_BASE_DIR}/logs

# Get the sample list for Mus musculus
SAMPLE_CSV="/gpfs/commons/home/jameslee/Cas13/liver_atlas/metadata/mouse_SRR_list_v3.1.txt"
samples=($(awk '{print $1}' $SAMPLE_CSV))

for sample in "${samples[@]}"; do
  # Define the output directory for this sample (shortened for --id)
  OUTPUT_DIR="${OUTPUT_BASE_DIR}/cellranger_${sample}"
  
  # Use a shorter id (only the sample name)
  ID="${sample}"

  # Run Cell Ranger for the sample using the library CSV
  ${CELLRANGER_PATH} count \
      --id=${ID} \
      --description="Output for ${sample}" \
      --transcriptome=${TRANSCRIPTOME} \
      --libraries=${LIBRARY_CSV} \
      --create-bam=true \
      --output-dir=${OUTPUT_DIR}
done