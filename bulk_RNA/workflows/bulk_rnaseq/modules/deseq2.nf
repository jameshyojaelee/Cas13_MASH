process DESEQ2 {
    tag "differential_expression"
    label 'medium_memory'
    
    container 'quay.io/biocontainers/bioconductor-deseq2:1.34.0--r41h399db7b_0'
    
    input:
    path counts_files
    tuple val(sample), val(group), val(disease), val(fibrosis), val(nas), val(stage)
    
    output:
    path "deseq2_results/*.csv"
    path "deseq2_results/*.rds"
    path "deseq2_results/plots/*.pdf"
    
    script:
    """
    #!/usr/bin/env Rscript
    
    library(DESeq2)
    library(ggplot2)
    
    # Create output directory
    dir.create("deseq2_results/plots", recursive=TRUE)
    
    # Read count data
    counts_data <- read.table("${counts_files}", header=TRUE, row.names=1)
    
    # Create metadata
    metadata <- data.frame(
        sample = "${sample}",
        group = "${group}",
        disease = "${disease}",
        fibrosis = "${fibrosis}",
        nas = "${nas}",
        stage = "${stage}"
    )
    
    # Create DESeq2 object
    dds <- DESeqDataSetFromMatrix(
        countData = counts_data,
        colData = metadata,
        design = ~ disease + stage
    )
    
    # Run DESeq2
    dds <- DESeq(dds)
    
    # Get results for different comparisons
    contrasts <- list(
        disease = c("disease", "NAFLD", "Control"),
        stage = c("stage", "moderate", "early")
    )
    
    for (contrast in names(contrasts)) {
        res <- results(dds, contrast=contrasts[[contrast]])
        res <- res[!is.na(res\$padj),]
        res <- res[res\$padj < ${params.deseq2_fdr},]
        res <- res[abs(res\$log2FoldChange) > ${params.min_fold_change},]
        
        # Save results
        write.csv(
            as.data.frame(res),
            file=paste0("deseq2_results/", contrast, "_results.csv")
        )
        
        # Create volcano plot
        pdf(paste0("deseq2_results/plots/", contrast, "_volcano.pdf"))
        plot(
            res\$log2FoldChange,
            -log10(res\$padj),
            main=paste0(contrast, " comparison"),
            xlab="log2 fold change",
            ylab="-log10 adjusted p-value"
        )
        dev.off()
    }
    
    # Save DESeq2 object
    saveRDS(dds, "deseq2_results/dds.rds")
    """
} 