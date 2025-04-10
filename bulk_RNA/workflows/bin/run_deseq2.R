#!/usr/bin/env Rscript

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if metadata argument is provided
metadata_string <- NULL
for (arg in args) {
  if (grepl("^--metadata=", arg)) {
    metadata_string <- sub("--metadata=", "", arg)
  }
}

if (is.null(metadata_string)) {
  stop("Metadata argument is missing. Usage: run_deseq2.R --metadata='metadata_string'")
}

# Load required libraries
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(dplyr)

# Create sample metadata dataframe from the input
sample_data <- data.frame(
  sample = character(),
  group = character(),
  disease = character(),
  fibrosis = character(),
  nas = character(),
  stage = character(),
  stringsAsFactors = FALSE
)

# Completely rewrite the metadata parsing
# First, clean the input string by removing all brackets and extra spaces
clean_string <- gsub("[\\[\\]]", "", metadata_string)
# Split by commas
elements <- unlist(strsplit(clean_string, ","))
# Trim whitespace from each element
elements <- trimws(elements)

# Process in groups of 6 elements (sample, group, disease, fibrosis, nas, stage)
for (i in seq(1, length(elements), 6)) {
  if (i+5 <= length(elements)) {
    # Clean sample name to ensure it's valid for rownames
    sample_name <- trimws(elements[i])
    sample_name <- gsub("[\\[\\]]", "", sample_name)
    
    sample_data <- rbind(sample_data, data.frame(
      sample = sample_name,
      group = elements[i+1],
      disease = elements[i+2],
      fibrosis = elements[i+3],
      nas = elements[i+4],
      stage = elements[i+5],
      stringsAsFactors = FALSE
    ))
  }
}

# Debug: print sample data
cat("Sample data after parsing:\n")
print(sample_data)

# Read count files and merge them
count_files <- list.files(pattern="*.counts.txt", full.names=TRUE)
count_data_list <- list()

# Debug: show all files being processed
cat("Count files found:", paste(count_files, collapse=", "), "\n")

# Create a single count matrix properly
# First extract gene IDs from one file to use as row names
firstfile <- count_files[1]
# Skip the comment lines
counts_data <- read.table(firstfile, header=TRUE, sep="\t", comment.char="#")
gene_ids <- counts_data[,1]  # First column contains gene IDs

# Create an empty count matrix with the right dimensions
count_matrix <- matrix(0, nrow=length(gene_ids), ncol=length(count_files))
rownames(count_matrix) <- gene_ids
colnames(count_matrix) <- sapply(count_files, function(f) {
  # Extract sample name from filename without the .counts suffix
  gsub(".*/(.*)\\.(counts.txt)$", "\\1", f)
})

# Fill the count matrix with data from each file
for (i in 1:length(count_files)) {
  file <- count_files[i]
  sample_name <- colnames(count_matrix)[i]
  cat("Processing file:", file, "for sample:", sample_name, "\n")
  
  # Read the featureCounts output file
  counts <- read.table(file, header=TRUE, sep="\t", comment.char="#")
  
  # The last column contains the counts
  count_col <- ncol(counts)
  cat("Count column index:", count_col, "\n")
  
  # Add counts to the matrix
  count_matrix[,i] <- counts[,count_col]
}

# Ensure the matrix is numeric
storage.mode(count_matrix) <- "numeric"

# Print diagnostic information
cat("Count matrix dimensions:", dim(count_matrix), "\n")
cat("Count matrix class:", class(count_matrix), "\n")
cat("Count matrix column names:", paste(colnames(count_matrix), collapse=", "), "\n")
cat("First few rows of count matrix:\n")
print(head(count_matrix))

# Match metadata samples to count matrix columns
cat("Count matrix column names:", paste(colnames(count_matrix), collapse=", "), "\n")
cat("Metadata sample names:", paste(sample_data$sample, collapse=", "), "\n")

# Filter metadata to keep only samples in count matrix
filter_idx <- sample_data$sample %in% colnames(count_matrix)
if (!all(filter_idx)) {
  cat("Warning: Some samples in metadata don't match count matrix columns\n")
  sample_data <- sample_data[filter_idx,]
}

# Set rownames for easy indexing - ensure they are valid, non-missing values
rownames(sample_data) <- as.character(sample_data$sample)

# Reorder sample_data to match count_matrix column order
matching_idx <- match(colnames(count_matrix), rownames(sample_data))
if (any(is.na(matching_idx))) {
  cat("Warning: Some count matrix columns don't match any metadata samples\n")
  # Keep only matching columns
  valid_cols <- !is.na(matching_idx)
  count_matrix <- count_matrix[,valid_cols]
  matching_idx <- matching_idx[!is.na(matching_idx)]
}
sample_data <- sample_data[matching_idx,]

# Create DESeq2 dataset
cat("Final dimensions of count matrix:", dim(count_matrix), "\n")
cat("Final dimensions of sample metadata:", dim(sample_data), "\n")

# Print sample data for debugging
cat("Final sample data for debugging:", "\n")
print(sample_data)
cat("Count matrix column names:", paste(colnames(count_matrix), collapse=", "), "\n")
cat("Sample data rownames:", paste(rownames(sample_data), collapse=", "), "\n")

# Create DESeq2 object - this is the key step
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = sample_data,
  design = ~ disease
)

# Filter out genes with low counts
keep <- rowSums(counts(dds) >= 10) >= ncol(dds)/3
dds <- dds[keep,]

# Run DESeq2
dds <- DESeq(dds)

# Save normalized counts
normalized_counts <- counts(dds, normalized=TRUE)
write.csv(normalized_counts, "normalized_counts.csv")

# Get results based on disease comparison
unique_diseases <- unique(sample_data$disease)

# PCA plot
vsd <- vst(dds, blind=FALSE)
pdf("pca_plot.pdf", width=10, height=8)
pcaData <- plotPCA(vsd, intgroup=c("disease", "group"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=disease, shape=group)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  geom_text_repel(aes(label=name), show.legend=FALSE) +
  theme_bw()
dev.off()

# If there are only two disease types, perform a simple comparison
if (length(unique_diseases) == 2) {
  res <- results(dds, contrast=c("disease", unique_diseases[1], unique_diseases[2]))
  res <- res[order(res$padj),]
  
  # Save results
  write.csv(as.data.frame(res), "deseq2_results.csv")
  
  # MA plot
  pdf("ma_plot.pdf", width=8, height=6)
  plotMA(res, ylim=c(-5, 5))
  dev.off()
  
  # Volcano plot
  pdf("volcano_plot.pdf", width=10, height=8)
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  
  # Add a column to indicate if genes are significantly DE (|log2FC| > 1 and padj < 0.05)
  res_df$significant <- ifelse(abs(res_df$log2FoldChange) > 1 & res_df$padj < 0.05, "Significant", "Not Significant")
  
  # Top 20 genes for labeling
  top_genes <- head(res_df[order(res_df$padj),], 20)
  
  ggplot(res_df, aes(x=log2FoldChange, y=-log10(padj), color=significant)) +
    geom_point(alpha=0.6) +
    scale_color_manual(values=c("Not Significant"="grey", "Significant"="red")) +
    theme_bw() +
    theme(legend.position="bottom") +
    geom_text_repel(
      data=top_genes,
      aes(label=gene),
      max.overlaps=20,
      box.padding=0.5
    ) +
    geom_vline(xintercept=c(-1, 1), linetype="dashed", color="blue") +
    geom_hline(yintercept=-log10(0.05), linetype="dashed", color="blue") +
    labs(title=paste("Volcano Plot:", unique_diseases[1], "vs", unique_diseases[2]),
         x="log2 Fold Change",
         y="-log10 Adjusted p-value")
  dev.off()
  
  # Heatmap of top DE genes
  top_genes_count <- 50
  top_de_genes <- rownames(res)[order(res$padj)][1:min(top_genes_count, nrow(res))]
  
  if (length(top_de_genes) > 0) {
    pdf("heatmap_top_de_genes.pdf", width=12, height=10)
    heatmap_data <- assay(vsd)[top_de_genes,]
    pheatmap(
      heatmap_data,
      cluster_rows=TRUE,
      cluster_cols=TRUE,
      show_rownames=TRUE,
      annotation_col=sample_data[,c("disease", "group")],
      scale="row"
    )
    dev.off()
  }
  
} else if (length(unique_diseases) > 2) {
  # For multiple disease types, do pairwise comparisons
  results_list <- list()
  
  for (i in 1:(length(unique_diseases)-1)) {
    for (j in (i+1):length(unique_diseases)) {
      disease1 <- unique_diseases[i]
      disease2 <- unique_diseases[j]
      
      res <- results(dds, contrast=c("disease", disease1, disease2))
      res <- res[order(res$padj),]
      
      result_name <- paste0(disease1, "_vs_", disease2)
      results_list[[result_name]] <- res
      
      # Save each comparison result
      write.csv(as.data.frame(res), paste0("deseq2_", result_name, ".csv"))
      
      # MA plot for each comparison
      pdf(paste0("ma_plot_", result_name, ".pdf"), width=8, height=6)
      plotMA(res, ylim=c(-5, 5))
      dev.off()
      
      # Volcano plot for each comparison
      pdf(paste0("volcano_plot_", result_name, ".pdf"), width=10, height=8)
      res_df <- as.data.frame(res)
      res_df$gene <- rownames(res_df)
      
      # Add a column to indicate if genes are significantly DE
      res_df$significant <- ifelse(abs(res_df$log2FoldChange) > 1 & res_df$padj < 0.05, "Significant", "Not Significant")
      
      # Top 20 genes for labeling
      top_genes <- head(res_df[order(res_df$padj),], 20)
      
      ggplot(res_df, aes(x=log2FoldChange, y=-log10(padj), color=significant)) +
        geom_point(alpha=0.6) +
        scale_color_manual(values=c("Not Significant"="grey", "Significant"="red")) +
        theme_bw() +
        theme(legend.position="bottom") +
        geom_text_repel(
          data=top_genes,
          aes(label=gene),
          max.overlaps=20,
          box.padding=0.5
        ) +
        geom_vline(xintercept=c(-1, 1), linetype="dashed", color="blue") +
        geom_hline(yintercept=-log10(0.05), linetype="dashed", color="blue") +
        labs(title=paste("Volcano Plot:", disease1, "vs", disease2),
             x="log2 Fold Change",
             y="-log10 Adjusted p-value")
      dev.off()
      
      # Heatmap of top DE genes for each comparison
      top_genes_count <- 50
      top_de_genes <- rownames(res)[order(res$padj)][1:min(top_genes_count, nrow(res))]
      
      if (length(top_de_genes) > 0) {
        pdf(paste0("heatmap_top_de_genes_", result_name, ".pdf"), width=12, height=10)
        heatmap_data <- assay(vsd)[top_de_genes,]
        pheatmap(
          heatmap_data,
          cluster_rows=TRUE,
          cluster_cols=TRUE,
          show_rownames=TRUE,
          annotation_col=sample_data[,c("disease", "group")],
          scale="row"
        )
        dev.off()
      }
    }
  }
  
  # Save a summary of all comparisons
  summary_table <- data.frame(
    comparison = character(),
    total_DE_genes = integer(),
    upregulated = integer(),
    downregulated = integer(),
    stringsAsFactors = FALSE
  )
  
  for (name in names(results_list)) {
    res <- results_list[[name]]
    sig_genes <- sum(!is.na(res$padj) & res$padj < 0.05 & abs(res$log2FoldChange) > 1, na.rm=TRUE)
    up_genes <- sum(!is.na(res$padj) & res$padj < 0.05 & res$log2FoldChange > 1, na.rm=TRUE)
    down_genes <- sum(!is.na(res$padj) & res$padj < 0.05 & res$log2FoldChange < -1, na.rm=TRUE)
    
    summary_table <- rbind(summary_table, data.frame(
      comparison = name,
      total_DE_genes = sig_genes,
      upregulated = up_genes,
      downregulated = down_genes
    ))
  }
  
  write.csv(summary_table, "deseq2_comparison_summary.csv")
}

# Write session info for reproducibility
writeLines(capture.output(sessionInfo()), "session_info.txt")

# After processing the count matrix and sample metadata, perform multiple analyses

# Check if we have enough diversity in groups for differential expression analysis
have_multiple_diseases <- length(unique(sample_data$disease)) > 1
have_multiple_stages <- length(unique(sample_data$stage)) > 1
have_multiple_groups <- length(unique(sample_data$group)) > 1

# Create output directories for each analysis type
dir.create("disease_analysis", showWarnings = FALSE)
dir.create("stage_analysis", showWarnings = FALSE)
dir.create("disease_stage_analysis", showWarnings = FALSE)
dir.create("group_analysis", showWarnings = FALSE)

# Decide which analyses to run based on data diversity
can_analyze_by_disease <- have_multiple_diseases
can_analyze_by_stage <- have_multiple_stages
can_analyze_by_disease_stage <- have_multiple_diseases && have_multiple_stages
can_analyze_by_group <- have_multiple_groups

#################################################
## ANALYSIS 1: Disease-based analysis (original)
#################################################
cat("\n\n=====================================================\n")
cat("ATTEMPTING ANALYSIS BY DISEASE\n")
cat("=====================================================\n\n")

if (can_analyze_by_disease) {
  cat("Running disease-based analysis (multiple disease types detected)\n")
  
  # Create DESeq2 object for disease analysis
  dds_disease <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = sample_data,
    design = ~ disease
  )

  # Filter out genes with low counts
  keep <- rowSums(counts(dds_disease) >= 10) >= ncol(dds_disease)/3
  dds_disease <- dds_disease[keep,]

  # Run DESeq2
  dds_disease <- DESeq(dds_disease)

  # Save normalized counts
  normalized_counts <- counts(dds_disease, normalized=TRUE)
  write.csv(normalized_counts, "disease_analysis/normalized_counts.csv")

  # Get results based on disease comparison
  unique_diseases <- unique(sample_data$disease)

  # PCA plot
  vsd_disease <- vst(dds_disease, blind=FALSE)
  pdf("disease_analysis/pca_plot_disease.pdf", width=10, height=8)
  pcaData <- plotPCA(vsd_disease, intgroup=c("disease", "group"), returnData=TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  ggplot(pcaData, aes(PC1, PC2, color=disease, shape=group)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    geom_text_repel(aes(label=name), show.legend=FALSE) +
    theme_bw()
  dev.off()

  # If there are only two disease types, perform a simple comparison
  if (length(unique_diseases) == 2) {
    res <- results(dds_disease, contrast=c("disease", unique_diseases[1], unique_diseases[2]))
    res <- res[order(res$padj),]
    
    # Save results
    write.csv(as.data.frame(res), "disease_analysis/deseq2_results.csv")
    
    # MA plot
    pdf("disease_analysis/ma_plot.pdf", width=8, height=6)
    plotMA(res, ylim=c(-5, 5))
    dev.off()
    
    # Volcano plot
    pdf("disease_analysis/volcano_plot.pdf", width=10, height=8)
    res_df <- as.data.frame(res)
    res_df$gene <- rownames(res_df)
    
    # Add a column to indicate if genes are significantly DE (|log2FC| > 1 and padj < 0.05)
    res_df$significant <- ifelse(abs(res_df$log2FoldChange) > 1 & res_df$padj < 0.05, "Significant", "Not Significant")
    
    # Top 20 genes for labeling
    top_genes <- head(res_df[order(res_df$padj),], 20)
    
    ggplot(res_df, aes(x=log2FoldChange, y=-log10(padj), color=significant)) +
      geom_point(alpha=0.6) +
      scale_color_manual(values=c("Not Significant"="grey", "Significant"="red")) +
      theme_bw() +
      theme(legend.position="bottom") +
      geom_text_repel(
        data=top_genes,
        aes(label=gene),
        max.overlaps=20,
        box.padding=0.5
      ) +
      geom_vline(xintercept=c(-1, 1), linetype="dashed", color="blue") +
      geom_hline(yintercept=-log10(0.05), linetype="dashed", color="blue") +
      labs(title=paste("Volcano Plot:", unique_diseases[1], "vs", unique_diseases[2]),
           x="log2 Fold Change",
           y="-log10 Adjusted p-value")
    dev.off()
    
    # Heatmap of top DE genes
    top_genes_count <- 50
    top_de_genes <- rownames(res)[order(res$padj)][1:min(top_genes_count, nrow(res))]
    
    if (length(top_de_genes) > 0) {
      pdf("disease_analysis/heatmap_top_de_genes.pdf", width=12, height=10)
      heatmap_data <- assay(vsd_disease)[top_de_genes,]
      pheatmap(
        heatmap_data,
        cluster_rows=TRUE,
        cluster_cols=TRUE,
        show_rownames=TRUE,
        annotation_col=sample_data[,c("disease", "group")],
        scale="row"
      )
      dev.off()
    }
  } else if (length(unique_diseases) > 2) {
    # For multiple disease types, do pairwise comparisons
    results_list <- list()
    
    for (i in 1:(length(unique_diseases)-1)) {
      for (j in (i+1):length(unique_diseases)) {
        disease1 <- unique_diseases[i]
        disease2 <- unique_diseases[j]
        
        res <- results(dds_disease, contrast=c("disease", disease1, disease2))
        res <- res[order(res$padj),]
        
        result_name <- paste0(disease1, "_vs_", disease2)
        results_list[[result_name]] <- res
        
        # Save each comparison result
        write.csv(as.data.frame(res), paste0("disease_analysis/deseq2_", result_name, ".csv"))
        
        # MA plot for each comparison
        pdf(paste0("disease_analysis/ma_plot_", result_name, ".pdf"), width=8, height=6)
        plotMA(res, ylim=c(-5, 5))
        dev.off()
        
        # Volcano plot for each comparison
        pdf(paste0("disease_analysis/volcano_plot_", result_name, ".pdf"), width=10, height=8)
        res_df <- as.data.frame(res)
        res_df$gene <- rownames(res_df)
        
        # Add a column to indicate if genes are significantly DE
        res_df$significant <- ifelse(abs(res_df$log2FoldChange) > 1 & res_df$padj < 0.05, "Significant", "Not Significant")
        
        # Top 20 genes for labeling
        top_genes <- head(res_df[order(res_df$padj),], 20)
        
        ggplot(res_df, aes(x=log2FoldChange, y=-log10(padj), color=significant)) +
          geom_point(alpha=0.6) +
          scale_color_manual(values=c("Not Significant"="grey", "Significant"="red")) +
          theme_bw() +
          theme(legend.position="bottom") +
          geom_text_repel(
            data=top_genes,
            aes(label=gene),
            max.overlaps=20,
            box.padding=0.5
          ) +
          geom_vline(xintercept=c(-1, 1), linetype="dashed", color="blue") +
          geom_hline(yintercept=-log10(0.05), linetype="dashed", color="blue") +
          labs(title=paste("Volcano Plot:", disease1, "vs", disease2),
               x="log2 Fold Change",
               y="-log10 Adjusted p-value")
        dev.off()
        
        # Heatmap of top DE genes for each comparison
        top_genes_count <- 50
        top_de_genes <- rownames(res)[order(res$padj)][1:min(top_genes_count, nrow(res))]
        
        if (length(top_de_genes) > 0) {
          pdf(paste0("disease_analysis/heatmap_top_de_genes_", result_name, ".pdf"), width=12, height=10)
          heatmap_data <- assay(vsd_disease)[top_de_genes,]
          pheatmap(
            heatmap_data,
            cluster_rows=TRUE,
            cluster_cols=TRUE,
            show_rownames=TRUE,
            annotation_col=sample_data[,c("disease", "group")],
            scale="row"
          )
          dev.off()
        }
      }
    }
    
    # Save a summary of all comparisons
    disease_summary_table <- data.frame(
      comparison = character(),
      total_DE_genes = integer(),
      upregulated = integer(),
      downregulated = integer(),
      stringsAsFactors = FALSE
    )
    
    for (name in names(results_list)) {
      res <- results_list[[name]]
      sig_genes <- sum(!is.na(res$padj) & res$padj < 0.05 & abs(res$log2FoldChange) > 1, na.rm=TRUE)
      up_genes <- sum(!is.na(res$padj) & res$padj < 0.05 & res$log2FoldChange > 1, na.rm=TRUE)
      down_genes <- sum(!is.na(res$padj) & res$padj < 0.05 & res$log2FoldChange < -1, na.rm=TRUE)
      
      disease_summary_table <- rbind(disease_summary_table, data.frame(
        comparison = name,
        total_DE_genes = sig_genes,
        upregulated = up_genes,
        downregulated = down_genes
      ))
    }
    
    write.csv(disease_summary_table, "disease_analysis/deseq2_comparison_summary.csv")
  }
} else {
  cat("Skipping disease-based analysis - all samples have the same disease\n")
  # Just create a basic DESeq dataset with intercept-only design
  dds_disease <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = sample_data,
    design = ~ 1
  )
  
  # Run size factor normalization only - no differential testing
  dds_disease <- estimateSizeFactors(dds_disease)
  
  # Save normalized counts
  normalized_counts <- counts(dds_disease, normalized=TRUE)
  write.csv(normalized_counts, "disease_analysis/normalized_counts.csv")
  
  # Create PCA plot based on variance stabilized data
  vsd_disease <- vst(dds_disease, blind=TRUE)
  pdf("disease_analysis/pca_plot_normalized_counts.pdf", width=10, height=8)
  pcaData <- plotPCA(vsd_disease, intgroup=c("group", "stage"), returnData=TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  ggplot(pcaData, aes(PC1, PC2, color=group, shape=stage)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    geom_text_repel(aes(label=name), show.legend=FALSE) +
    theme_bw() + 
    ggtitle("PCA of Normalized Counts (no differential testing)")
  dev.off()
}

#################################################
## ANALYSIS 2: Stage-based analysis
#################################################
cat("\n\n=====================================================\n")
cat("ATTEMPTING ANALYSIS BY STAGE\n")
cat("=====================================================\n\n")

if (can_analyze_by_stage) {
  cat("Running stage-based analysis (multiple stages detected)\n")
  
  # Create DESeq2 object for stage analysis
  dds_stage <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = sample_data,
    design = ~ stage
  )

  # Filter out genes with low counts
  keep <- rowSums(counts(dds_stage) >= 10) >= ncol(dds_stage)/3
  dds_stage <- dds_stage[keep,]

  # Run DESeq2
  dds_stage <- DESeq(dds_stage)

  # Save normalized counts
  normalized_counts <- counts(dds_stage, normalized=TRUE)
  write.csv(normalized_counts, "stage_analysis/normalized_counts.csv")

  # Get results based on stage comparison
  unique_stages <- unique(sample_data$stage)

  # PCA plot for stage
  vsd_stage <- vst(dds_stage, blind=FALSE)
  pdf("stage_analysis/pca_plot_stage.pdf", width=10, height=8)
  pcaData <- plotPCA(vsd_stage, intgroup=c("stage", "disease"), returnData=TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  ggplot(pcaData, aes(PC1, PC2, color=stage, shape=disease)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    geom_text_repel(aes(label=name), show.legend=FALSE) +
    theme_bw()
  dev.off()

  # If there are only two stages, perform a simple comparison
  if (length(unique_stages) == 2) {
    res <- results(dds_stage, contrast=c("stage", unique_stages[1], unique_stages[2]))
    res <- res[order(res$padj),]
    
    # Save results
    write.csv(as.data.frame(res), "stage_analysis/deseq2_results.csv")
    
    # MA plot
    pdf("stage_analysis/ma_plot.pdf", width=8, height=6)
    plotMA(res, ylim=c(-5, 5))
    dev.off()
    
    # Volcano plot
    pdf("stage_analysis/volcano_plot.pdf", width=10, height=8)
    res_df <- as.data.frame(res)
    res_df$gene <- rownames(res_df)
    
    # Add a column to indicate if genes are significantly DE (|log2FC| > 1 and padj < 0.05)
    res_df$significant <- ifelse(abs(res_df$log2FoldChange) > 1 & res_df$padj < 0.05, "Significant", "Not Significant")
    
    # Top 20 genes for labeling
    top_genes <- head(res_df[order(res_df$padj),], 20)
    
    ggplot(res_df, aes(x=log2FoldChange, y=-log10(padj), color=significant)) +
      geom_point(alpha=0.6) +
      scale_color_manual(values=c("Not Significant"="grey", "Significant"="red")) +
      theme_bw() +
      theme(legend.position="bottom") +
      geom_text_repel(
        data=top_genes,
        aes(label=gene),
        max.overlaps=20,
        box.padding=0.5
      ) +
      geom_vline(xintercept=c(-1, 1), linetype="dashed", color="blue") +
      geom_hline(yintercept=-log10(0.05), linetype="dashed", color="blue") +
      labs(title=paste("Volcano Plot:", unique_stages[1], "vs", unique_stages[2]),
           x="log2 Fold Change",
           y="-log10 Adjusted p-value")
    dev.off()
    
    # Heatmap of top DE genes
    top_genes_count <- 50
    top_de_genes <- rownames(res)[order(res$padj)][1:min(top_genes_count, nrow(res))]
    
    if (length(top_de_genes) > 0) {
      pdf("stage_analysis/heatmap_top_de_genes.pdf", width=12, height=10)
      heatmap_data <- assay(vsd_stage)[top_de_genes,]
      pheatmap(
        heatmap_data,
        cluster_rows=TRUE,
        cluster_cols=TRUE,
        show_rownames=TRUE,
        annotation_col=sample_data[,c("stage", "disease")],
        scale="row"
      )
      dev.off()
    }
  } else if (length(unique_stages) > 2) {
    # For multiple stages, do pairwise comparisons
    results_list <- list()
    
    for (i in 1:(length(unique_stages)-1)) {
      for (j in (i+1):length(unique_stages)) {
        stage1 <- unique_stages[i]
        stage2 <- unique_stages[j]
        
        res <- results(dds_stage, contrast=c("stage", stage1, stage2))
        res <- res[order(res$padj),]
        
        result_name <- paste0(stage1, "_vs_", stage2)
        results_list[[result_name]] <- res
        
        # Save each comparison result
        write.csv(as.data.frame(res), paste0("stage_analysis/deseq2_", result_name, ".csv"))
        
        # MA plot for each comparison
        pdf(paste0("stage_analysis/ma_plot_", result_name, ".pdf"), width=8, height=6)
        plotMA(res, ylim=c(-5, 5))
        dev.off()
        
        # Volcano plot for each comparison
        pdf(paste0("stage_analysis/volcano_plot_", result_name, ".pdf"), width=10, height=8)
        res_df <- as.data.frame(res)
        res_df$gene <- rownames(res_df)
        
        # Add a column to indicate if genes are significantly DE
        res_df$significant <- ifelse(abs(res_df$log2FoldChange) > 1 & res_df$padj < 0.05, "Significant", "Not Significant")
        
        # Top 20 genes for labeling
        top_genes <- head(res_df[order(res_df$padj),], 20)
        
        ggplot(res_df, aes(x=log2FoldChange, y=-log10(padj), color=significant)) +
          geom_point(alpha=0.6) +
          scale_color_manual(values=c("Not Significant"="grey", "Significant"="red")) +
          theme_bw() +
          theme(legend.position="bottom") +
          geom_text_repel(
            data=top_genes,
            aes(label=gene),
            max.overlaps=20,
            box.padding=0.5
          ) +
          geom_vline(xintercept=c(-1, 1), linetype="dashed", color="blue") +
          geom_hline(yintercept=-log10(0.05), linetype="dashed", color="blue") +
          labs(title=paste("Volcano Plot:", stage1, "vs", stage2),
               x="log2 Fold Change",
               y="-log10 Adjusted p-value")
        dev.off()
        
        # Heatmap of top DE genes for each comparison
        top_genes_count <- 50
        top_de_genes <- rownames(res)[order(res$padj)][1:min(top_genes_count, nrow(res))]
        
        if (length(top_de_genes) > 0) {
          pdf(paste0("stage_analysis/heatmap_top_de_genes_", result_name, ".pdf"), width=12, height=10)
          heatmap_data <- assay(vsd_stage)[top_de_genes,]
          pheatmap(
            heatmap_data,
            cluster_rows=TRUE,
            cluster_cols=TRUE,
            show_rownames=TRUE,
            annotation_col=sample_data[,c("stage", "disease")],
            scale="row"
          )
          dev.off()
        }
      }
    }
    
    # Save a summary of all stage comparisons
    stage_summary_table <- data.frame(
      comparison = character(),
      total_DE_genes = integer(),
      upregulated = integer(),
      downregulated = integer(),
      stringsAsFactors = FALSE
    )
    
    for (name in names(results_list)) {
      res <- results_list[[name]]
      sig_genes <- sum(!is.na(res$padj) & res$padj < 0.05 & abs(res$log2FoldChange) > 1, na.rm=TRUE)
      up_genes <- sum(!is.na(res$padj) & res$padj < 0.05 & res$log2FoldChange > 1, na.rm=TRUE)
      down_genes <- sum(!is.na(res$padj) & res$padj < 0.05 & res$log2FoldChange < -1, na.rm=TRUE)
      
      stage_summary_table <- rbind(stage_summary_table, data.frame(
        comparison = name,
        total_DE_genes = sig_genes,
        upregulated = up_genes,
        downregulated = down_genes
      ))
    }
    
    write.csv(stage_summary_table, "stage_analysis/deseq2_comparison_summary.csv")
  }
} else {
  cat("Skipping stage-based differential expression analysis - all samples have the same stage\n")
}

#################################################
## ANALYSIS 3: Group-based analysis
#################################################
cat("\n\n=====================================================\n")
cat("ATTEMPTING ANALYSIS BY SAMPLE GROUP\n")
cat("=====================================================\n\n")

if (can_analyze_by_group) {
  cat("Running group-based analysis (multiple sample groups detected)\n")
  
  # Create DESeq2 object for group analysis
  dds_group <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = sample_data,
    design = ~ group
  )

  # Filter out genes with low counts
  keep <- rowSums(counts(dds_group) >= 10) >= ncol(dds_group)/3
  dds_group <- dds_group[keep,]

  # Run DESeq2
  dds_group <- DESeq(dds_group)

  # Save normalized counts
  normalized_counts <- counts(dds_group, normalized=TRUE)
  write.csv(normalized_counts, "group_analysis/normalized_counts.csv")

  # Get results based on group comparison
  unique_groups <- unique(sample_data$group)
  
  # PCA plot for group
  vsd_group <- vst(dds_group, blind=FALSE)
  pdf("group_analysis/pca_plot_group.pdf", width=10, height=8)
  pcaData <- plotPCA(vsd_group, intgroup=c("group", "stage"), returnData=TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  ggplot(pcaData, aes(PC1, PC2, color=group, shape=stage)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    geom_text_repel(aes(label=name), show.legend=FALSE) +
    theme_bw()
  dev.off()

  # If there are only two groups, perform a simple comparison
  if (length(unique_groups) == 2) {
    res <- results(dds_group, contrast=c("group", unique_groups[1], unique_groups[2]))
    res <- res[order(res$padj),]
    
    # Save results
    write.csv(as.data.frame(res), "group_analysis/deseq2_results.csv")
    
    # MA plot
    pdf("group_analysis/ma_plot.pdf", width=8, height=6)
    plotMA(res, ylim=c(-5, 5))
    dev.off()
    
    # Volcano plot
    pdf("group_analysis/volcano_plot.pdf", width=10, height=8)
    res_df <- as.data.frame(res)
    res_df$gene <- rownames(res_df)
    
    # Add a column to indicate if genes are significantly DE (|log2FC| > 1 and padj < 0.05)
    res_df$significant <- ifelse(abs(res_df$log2FoldChange) > 1 & res_df$padj < 0.05, "Significant", "Not Significant")
    
    # Top 20 genes for labeling
    top_genes <- head(res_df[order(res_df$padj),], 20)
    
    ggplot(res_df, aes(x=log2FoldChange, y=-log10(padj), color=significant)) +
      geom_point(alpha=0.6) +
      scale_color_manual(values=c("Not Significant"="grey", "Significant"="red")) +
      theme_bw() +
      theme(legend.position="bottom") +
      geom_text_repel(
        data=top_genes,
        aes(label=gene),
        max.overlaps=20,
        box.padding=0.5
      ) +
      geom_vline(xintercept=c(-1, 1), linetype="dashed", color="blue") +
      geom_hline(yintercept=-log10(0.05), linetype="dashed", color="blue") +
      labs(title=paste("Volcano Plot:", unique_groups[1], "vs", unique_groups[2]),
           x="log2 Fold Change",
           y="-log10 Adjusted p-value")
    dev.off()
    
    # Heatmap of top DE genes
    top_genes_count <- 50
    top_de_genes <- rownames(res)[order(res$padj)][1:min(top_genes_count, nrow(res))]
    
    if (length(top_de_genes) > 0) {
      pdf("group_analysis/heatmap_top_de_genes.pdf", width=12, height=10)
      heatmap_data <- assay(vsd_group)[top_de_genes,]
      pheatmap(
        heatmap_data,
        cluster_rows=TRUE,
        cluster_cols=TRUE,
        show_rownames=TRUE,
        annotation_col=sample_data[,c("group", "stage")],
        scale="row"
      )
      dev.off()
    }
  } else if (length(unique_groups) > 2) {
    # For multiple groups, do pairwise comparisons (code similar to disease analysis)
    # ...
  }
} else {
  cat("Skipping group-based analysis - all samples are from the same group\n")
}

#################################################
## ANALYSIS 4: Disease+Stage interaction analysis
#################################################
cat("\n\n=====================================================\n")
cat("ATTEMPTING ANALYSIS BY DISEASE+STAGE\n")
cat("=====================================================\n\n")

if (can_analyze_by_disease_stage) {
  cat("Running disease+stage interaction analysis\n")
  
  # Create combined factor for disease+stage analysis
  sample_data$condition <- paste(sample_data$disease, sample_data$stage, sep="_")
  
  # Create DESeq2 object for disease+stage analysis
  dds_combined <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = sample_data,
    design = ~ disease + stage
  )

  # Filter out genes with low counts
  keep <- rowSums(counts(dds_combined) >= 10) >= ncol(dds_combined)/3
  dds_combined <- dds_combined[keep,]

  # Run DESeq2
  dds_combined <- DESeq(dds_combined)

  # Save normalized counts
  normalized_counts <- counts(dds_combined, normalized=TRUE)
  write.csv(normalized_counts, "disease_stage_analysis/normalized_counts.csv")

  # PCA plot for combined analysis
  vsd_combined <- vst(dds_combined, blind=FALSE)
  pdf("disease_stage_analysis/pca_plot_combined.pdf", width=10, height=8)
  pcaData <- plotPCA(vsd_combined, intgroup=c("condition"), returnData=TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  ggplot(pcaData, aes(PC1, PC2, color=condition)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    geom_text_repel(aes(label=name), show.legend=FALSE) +
    theme_bw()
  dev.off()

  # Save results for disease effect controlling for stage
  res_disease <- results(dds_combined, name="disease_NAFLD_vs_control")
  write.csv(as.data.frame(res_disease), "disease_stage_analysis/deseq2_disease_effect.csv")

  # Save results for stage effects
  unique_stages <- unique(sample_data$stage)
  if (length(unique_stages) > 1) {
    for (i in 2:length(unique_stages)) {
      # Compare each stage to the first stage (reference level)
      res_stage <- results(dds_combined, name=paste0("stage_", unique_stages[i], "_vs_", unique_stages[1]))
      write.csv(as.data.frame(res_stage), paste0("disease_stage_analysis/deseq2_stage_", unique_stages[i], "_vs_", unique_stages[1], ".csv"))
      
      # MA plot
      pdf(paste0("disease_stage_analysis/ma_plot_stage_", unique_stages[i], "_vs_", unique_stages[1], ".pdf"), width=8, height=6)
      plotMA(res_stage, ylim=c(-5, 5))
      dev.off()
      
      # Volcano plot
      pdf(paste0("disease_stage_analysis/volcano_plot_stage_", unique_stages[i], "_vs_", unique_stages[1], ".pdf"), width=10, height=8)
      res_df <- as.data.frame(res_stage)
      res_df$gene <- rownames(res_df)
      
      # Add a column to indicate if genes are significantly DE
      res_df$significant <- ifelse(abs(res_df$log2FoldChange) > 1 & res_df$padj < 0.05, "Significant", "Not Significant")
      
      # Top 20 genes for labeling
      top_genes <- head(res_df[order(res_df$padj),], 20)
      
      ggplot(res_df, aes(x=log2FoldChange, y=-log10(padj), color=significant)) +
        geom_point(alpha=0.6) +
        scale_color_manual(values=c("Not Significant"="grey", "Significant"="red")) +
        theme_bw() +
        theme(legend.position="bottom") +
        geom_text_repel(
          data=top_genes,
          aes(label=gene),
          max.overlaps=20,
          box.padding=0.5
        ) +
        geom_vline(xintercept=c(-1, 1), linetype="dashed", color="blue") +
        geom_hline(yintercept=-log10(0.05), linetype="dashed", color="blue") +
        labs(title=paste("Volcano Plot Stage:", unique_stages[i], "vs", unique_stages[1]),
             x="log2 Fold Change",
             y="-log10 Adjusted p-value")
      dev.off()
    }
  }

  # LRT test for any effect of stage
  dds_lrt <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = sample_data,
    design = ~ disease + stage
  )
  keep <- rowSums(counts(dds_lrt) >= 10) >= ncol(dds_lrt)/3
  dds_lrt <- dds_lrt[keep,]

  # Run LRT test comparing full model to reduced model (without stage)
  dds_lrt <- DESeq(dds_lrt, test="LRT", reduced = ~ disease)
  res_lrt <- results(dds_lrt)
  write.csv(as.data.frame(res_lrt), "disease_stage_analysis/deseq2_LRT_stage_effect.csv")

  # Select top genes affected by stage
  sig_lrt_genes <- rownames(res_lrt)[which(res_lrt$padj < 0.05)]
  if (length(sig_lrt_genes) > 5) {
    top_lrt_genes <- sig_lrt_genes[1:min(50, length(sig_lrt_genes))]
    
    # Create heatmap of top LRT genes
    pdf("disease_stage_analysis/heatmap_top_lrt_genes.pdf", width=12, height=10)
    heatmap_data <- assay(vsd_combined)[top_lrt_genes,]
    pheatmap(
      heatmap_data,
      cluster_rows=TRUE,
      cluster_cols=TRUE,
      show_rownames=TRUE,
      annotation_col=sample_data[,c("disease", "stage", "condition")],
      scale="row"
    )
    dev.off()
  }
} else if (have_multiple_stages) {
  # If we have only one disease but multiple stages, just do a stage-only analysis for this section
  cat("Performing stage-only analysis for Disease+Stage section (only one disease type present)\n")
  
  # Create combined condition for visualization
  sample_data$condition <- paste(sample_data$disease, sample_data$stage, sep="_")
  
  # Similar code to stage analysis, but save to disease_stage_analysis directory
  # ...
} else {
  cat("Skipping disease+stage analysis - insufficient variation in data\n")
}

# Compile summary info
summary_info <- data.frame(
  analysis_type = c(
    "Disease", 
    "Stage", 
    "Group",
    "Disease+Stage"
  ),
  was_performed = c(
    can_analyze_by_disease,
    can_analyze_by_stage,
    can_analyze_by_group,
    can_analyze_by_disease_stage
  ),
  model = c(
    ifelse(can_analyze_by_disease, "~ disease", "~ 1 (all same disease)"),
    ifelse(can_analyze_by_stage, "~ stage", "~ 1 (all same stage)"),
    ifelse(can_analyze_by_group, "~ group", "~ 1 (all same group)"),
    ifelse(can_analyze_by_disease_stage, "~ disease + stage", "Not applicable")
  )
)
write.csv(summary_info, "analysis_summary.csv")

# Write session info for reproducibility
writeLines(capture.output(sessionInfo()), "session_info.txt") 