#!/bin/bash

# Set the directory where your fastq files are located
FASTQ_DIR="/gpfs/commons/home/jameslee/Cas13/bulk_RNA/fastq"
OUTPUT_DIR="/gpfs/commons/home/jameslee/Cas13/bulk_RNA/results/nf-core_rnaseq"

# Generate samplesheet if not already created
if [ ! -f samplesheet.csv ]; then
  echo "sample,fastq_1,fastq_2,strandedness" > samplesheet.csv
  
  # Find all R1 fastq files and create samplesheet entries
  for R1 in ${FASTQ_DIR}/*_R1_*.fastq.gz ${FASTQ_DIR}/*_1.fastq.gz; do
    if [ -f "$R1" ]; then
      # Get the corresponding R2 file
      R2="${R1/_R1_/_R2_}"
      R2="${R2/_1.fastq/_2.fastq}"
      
      # Get the sample name
      SAMPLE=$(basename "$R1" | sed 's/_R1_.*$//' | sed 's/_1\.fastq.*$//')
      
      # Add to samplesheet if R2 exists
      if [ -f "$R2" ]; then
        echo "${SAMPLE},${R1},${R2},auto" >> samplesheet.csv
      else
        echo "${SAMPLE},${R1},,auto" >> samplesheet.csv
      fi
    fi
  done
fi

# Run the nf-core/rnaseq pipeline
nextflow run nf-core/rnaseq -r 3.14.0 \
  --input samplesheet.csv \
  --outdir ${OUTPUT_DIR} \
  --genome GRCh38 \
  --aligner star_salmon \
  --star_index "" \
  --salmon_index "" \
  --gencode \
  --email your.email@example.com \
  -profile standard \
  --max_memory '64.GB' \
  --max_cpus 8 \
  -resume

# Note: Remove the double quotes from --star_index and --salmon_index 
# if you have pre-built indices, or leave empty to build them during execution.
# Adjust -profile, --max_memory, and --max_cpus according to your computing environment.
