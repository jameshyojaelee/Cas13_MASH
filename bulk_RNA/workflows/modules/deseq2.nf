process DESEQ2 {
    publishDir "${params.outdir}/deseq2", mode: 'copy'
    
    input:
    path(counts_files)
    val(sample_metadata)
    
    output:
    path("*.html")
    path("*.txt")
    path("*.tsv")
    path("*.csv")
    path("*.pdf")
    path("*.png")
    
    script:
    """
    run_deseq2.R --metadata='${sample_metadata}'
    """
} 