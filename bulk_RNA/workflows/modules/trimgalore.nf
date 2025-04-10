process TRIMGALORE {
    tag "$sample"
    label 'low_memory'
    
    container 'quay.io/biocontainers/trim-galore:0.6.7--hdfd78af_0'
    
    input:
    tuple val(sample), path(read1), path(read2), val(group), val(disease), val(fibrosis), val(nas), val(stage)
    
    output:
    tuple val(sample), path("*_val_1.fq.gz"), path("*_val_2.fq.gz"), val(group), val(disease), val(fibrosis), val(nas), val(stage), emit: trimmed_reads
    path "*_trimming_report.txt", emit: reports
    
    script:
    """
    trim_galore \
        --paired \
        --quality ${params.min_quality} \
        --length ${params.min_length} \
        --cores ${task.cpus} \
        --gzip \
        ${read1} ${read2}
    """
} 