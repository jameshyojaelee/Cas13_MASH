process FASTQC {
    tag "$sample"
    label 'low_memory'
    
    container 'quay.io/biocontainers/fastqc:0.11.9--0'
    
    input:
    tuple val(sample), path(read1), path(read2), val(group), val(disease), val(fibrosis), val(nas), val(stage)
    
    output:
    path "fastqc_${sample}_logs/*.zip", emit: zip
    path "fastqc_${sample}_logs/*.html"
    
    script:
    """
    mkdir fastqc_${sample}_logs
    fastqc -o fastqc_${sample}_logs -q ${read1} ${read2}
    """
} 