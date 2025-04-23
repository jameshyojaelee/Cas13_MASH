#!/bin/bash
#SBATCH --job-name=cellranger_batch
#SBATCH --output=cellranger_batch_%j.out
#SBATCH --error=cellranger_batch_%j.err
#SBATCH --time=48:00:00 # Adjusted time for potentially multiple samples
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=196G
#SBATCH --partition=standard

# Path to Cell Ranger installation
CELLRANGER_PATH="/gpfs/commons/home/jameslee/cellranger-9.0.0"

# Function to print usage information
function print_usage {
    echo "Usage:"
    echo "  sbatch [SBATCH OPTIONS] $0 <fastq_dir> <transcriptome_ref> [feature_ref]"
    echo ""
    echo "Arguments:"
    echo "  fastq_dir        - Directory containing all FASTQ files (e.g., /path/to/FASTQ_output)"
    echo "  transcriptome_ref - Path to transcriptome reference"
    echo "  feature_ref      - Path to feature reference (optional, required for CRISPR samples)"
    echo ""
    echo "SLURM Options:"
    echo "  Use standard sbatch options like --cpus-per-task=N and --mem=XG to specify resources."
    echo "  The script will find samples in fastq_dir and run cellranger for each."
    exit 1
}

# Argument parsing
if [ $# -lt 2 ] || [ $# -gt 3 ]; then
    print_usage
fi

FASTQ_DIR=$1
TRANSCRIPTOME_REF=$2
FEATURE_REF=${3:-} # Optional feature reference

# Ensure running within a SLURM job (optional, but good practice)
# if [ -z "$SLURM_JOB_ID" ]; then
#     echo "Warning: Not running within a SLURM job. Using default resources." >&2
#     # Uncomment the line below to enforce running via sbatch
#     # exit 1 
# fi

# Use SLURM allocated resources, provide defaults if not running via sbatch
LOCAL_CORES=${SLURM_CPUS_PER_TASK:-16}
# SLURM_MEM_PER_NODE is in MB, convert to GB for cellranger. Default to 196G if not set.
SLURM_MEM_MB=${SLURM_MEM_PER_NODE:-196000} 
LOCAL_MEM_GB=$((SLURM_MEM_MB / 1024))
# Ensure at least 1GB memory
if [ "$LOCAL_MEM_GB" -lt 1 ]; then
    LOCAL_MEM_GB=1
fi

# echo "Running Cell Ranger Batch Processing"
# echo "FASTQ Directory: $FASTQ_DIR"
# echo "Transcriptome Reference: $TRANSCRIPTOME_REF"
# if [ -n "$FEATURE_REF" ]; then
#     echo "Feature Reference: $FEATURE_REF"
# fi
# echo "Using Cores: $LOCAL_CORES"
# echo "Using Memory: ${LOCAL_MEM_GB}GB"

# Find unique sample prefixes in the FASTQ directory
# Extracts prefixes like GEX1, CRISPR1 from files like GEX1_S1_L001_R1_001.fastq.gz
# Ignores 'Undetermined'
sample_prefixes=$(ls "$FASTQ_DIR" | grep -Eo '^[A-Za-z0-9]+_S[0-9]+' | grep -v '^Undetermined_S[0-9]+' | sed -E 's/_S[0-9]+$//' | sort -u)

if [ -z "$sample_prefixes" ]; then
    echo "Error: No valid sample prefixes found in $FASTQ_DIR" >&2
    exit 1
fi

# echo "Found sample prefixes:"
# for sample_id in $sample_prefixes; do
#     echo "  - $sample_id"
# done

# Loop through each sample and run cellranger count
for sample_id in $sample_prefixes; do
    # echo "----------------------------------------"
    # echo "Processing sample: $sample_id"
    # echo "----------------------------------------"

    # Determine sample type and construct command
    cmd_args=(
        "--id=$sample_id" \
        "--fastqs=$FASTQ_DIR" \
        "--sample=$sample_id" \
        "--transcriptome=$TRANSCRIPTOME_REF" \
        "--localcores=$LOCAL_CORES" \
        "--localmem=$LOCAL_MEM_GB"
    )

    if [[ "$sample_id" == CRISPR* ]]; then
        # echo "Sample type detected: CRISPR"
        if [ -z "$FEATURE_REF" ]; then
            echo "Error: Feature reference is required for CRISPR sample '$sample_id' but not provided." >&2
            continue # Skip this sample
        elif [ ! -f "$FEATURE_REF" ]; then
             echo "Error: Feature reference file '$FEATURE_REF' not found for CRISPR sample '$sample_id'." >&2
             continue # Skip this sample
        fi
        cmd_args+=(--feature-ref="$FEATURE_REF")
    elif [[ "$sample_id" == GEX* ]]; then
        # echo "Sample type detected: GEX"
        # No extra args needed for GEX
    else
        echo "Warning: Unknown sample type for prefix '$sample_id'. Skipping." >&2
        continue
    fi

    # Construct the full command
    full_cmd="$CELLRANGER_PATH/cellranger count ${cmd_args[*]}"
    
    # echo "Running command:"
    # echo "$full_cmd"

    # Execute the command
    eval $full_cmd

    # Check exit status
    if [ $? -ne 0 ]; then
        echo "Error processing sample $sample_id. Check logs in the '$sample_id' output directory." >&2
        # Decide whether to continue or exit. For now, we continue.
        # exit 1 
    fi

done

# echo "----------------------------------------"
# echo "All samples processed."
