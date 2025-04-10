#!/bin/bash

# Directory containing FASTQ files
FASTQ_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/Liver_Atlas/fastq"

# Loop over each FASTQ file
for file in "${FASTQ_DIR}"/*.fastq.gz; do
  # Extract the SRR number (text before the first underscore)
  srr_number=$(basename "$file" | cut -d_ -f1)
  
  # Create a directory named after the SRR number if it doesn't exist
  mkdir -p "${FASTQ_DIR}/${srr_number}"
  
  # Move the FASTQ file into the corresponding directory
  mv "$file" "${FASTQ_DIR}/${srr_number}/"
done

echo "FASTQ files organized into folders by SRR number."