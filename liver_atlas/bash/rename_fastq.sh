#!/bin/bash

# Directory containing the FASTQ files
FASTQ_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/Liver_Atlas/fastq"

# Loop through each SRR directory
for srr_dir in "${FASTQ_DIR}"/*/; do
  # Extract the SRR number (folder name)
  srr_number=$(basename "$srr_dir")
  
  # Rename `_1.fastq.gz` to `SampleName_S1_L001_R1_001.fastq.gz`
  if [ -f "${srr_dir}/${srr_number}_1.fastq.gz" ]; then
    mv "${srr_dir}/${srr_number}_1.fastq.gz" "${srr_dir}/${srr_number}_S1_L001_R1_001.fastq.gz"
  fi

  # Rename `_2.fastq.gz` to `SampleName_S1_L001_R2_001.fastq.gz`
  if [ -f "${srr_dir}/${srr_number}_2.fastq.gz" ]; then
    mv "${srr_dir}/${srr_number}_2.fastq.gz" "${srr_dir}/${srr_number}_S1_L001_R2_001.fastq.gz"
  fi
done

echo "FASTQ files renamed using the new naming convention."