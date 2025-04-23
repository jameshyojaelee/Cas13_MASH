process MULTIQC {
    label 'low_memory'
    publishDir "${params.outdir}/multiqc", mode: 'copy'
    
    container 'quay.io/biocontainers/multiqc:1.12--pyhdfd78af_0'
    
    input:
    path '*'
    
    output:
    path 'multiqc_report.html'
    path 'multiqc_data'
    
    script:
    """
    multiqc . -f
    """
} 