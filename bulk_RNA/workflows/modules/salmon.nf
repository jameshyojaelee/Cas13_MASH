process SALMON {
    tag "$sample"
    label 'high_memory'
    
    container 'quay.io/biocontainers/salmon:1.10.1--h7e5ed60_0'
    
    input:
    tuple val(sample), path(read1), path(read2), val(group), val(disease), val(fibrosis), val(nas), val(stage)
    path transcriptome
    path gtf
    
    output:
    path "${sample}", emit: salmon_out
    tuple val(sample), path("${sample}/quant.sf"), val(group), val(disease), val(fibrosis), val(nas), val(stage), emit: quant
    
    script:
    def index_cmd = ""
    def index_path = ""
    
    // Check if the provided transcriptome is an index or a FASTA
    if(transcriptome.isDirectory()) {
        // Use existing index
        index_path = transcriptome
    } else {
        // Create index
        index_cmd = """
        echo "Creating Salmon index at \$(date)"
        salmon index \
            -t ${transcriptome} \
            -i salmon_index \
            -p ${task.cpus} \
            --gencode
        echo "Salmon index created at \$(date)"
        """
        index_path = "salmon_index"
    }
    
    // Detect if reads are gzipped
    def is_gzipped = read1.getName().endsWith('.gz') ? true : false
    def read_opts = is_gzipped ? "--libType=A" : "--libType=A"
    
    """
    ${index_cmd}
    
    echo "Starting Salmon quantification for ${sample} at \$(date)"
    echo "Input files: ${read1} ${read2}"
    echo "Is gzipped: ${is_gzipped}"
    
    salmon quant \
        -i ${index_path} \
        -l A \
        -1 ${read1} \
        -2 ${read2} \
        -p ${task.cpus} \
        --validateMappings \
        --gcBias \
        --seqBias \
        -o ${sample}
    
    echo "Salmon quantification completed at \$(date)"
    """
} 