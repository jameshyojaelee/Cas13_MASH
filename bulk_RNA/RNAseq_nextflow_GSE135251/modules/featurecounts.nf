process FEATURECOUNTS {
    tag "$sample"
    label 'medium_memory'
    
    container 'quay.io/biocontainers/subread:2.0.1--h7132678_2'
    
    input:
    tuple val(sample), path(bam)
    path gtf
    
    output:
    path "${sample}.counts.txt", emit: counts
    path "${sample}.counts.txt.summary", emit: summary
    
    script:
    """
    featureCounts \
        -p \
        -T ${task.cpus} \
        -a ${gtf} \
        -o ${sample}.counts.txt \
        ${bam}
    """
} 