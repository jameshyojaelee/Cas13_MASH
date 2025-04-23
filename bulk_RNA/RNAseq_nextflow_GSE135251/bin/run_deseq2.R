#!/usr/bin/env Rscript

# Basic error handling
tryCatch({
  # Parse arguments
  args <- commandArgs(trailingOnly = TRUE)
  metadata_string <- NULL
  for (arg in args) {
    if (grepl("^--metadata=", arg)) {
      metadata_string <- sub("--metadata=", "", arg)
    }
  }
  
  # Load libraries with minimal dependencies
  library(DESeq2, warn.conflicts = FALSE, quietly = TRUE)
  
  # Parse metadata
  cat("Parsing metadata...\n")
  clean_string <- gsub("[\\[\\]]", "", metadata_string)
  elements <- unlist(strsplit(clean_string, ","))
  elements <- trimws(elements)
  
  # Create sample data frame
  sample_data <- data.frame(
    sample = character(),
    group = character(),
    disease = character(),
    fibrosis = numeric(),
    nas = numeric(),
    stage = character(),
    stringsAsFactors = FALSE
  )
  
  # Process in groups of 6 elements
  for (i in seq(1, length(elements), 6)) {
    if (i+5 <= length(elements)) {
      sample_name <- trimws(elements[i])
      sample_name <- gsub("[\\[\\]]", "", sample_name)
      
      sample_data <- rbind(sample_data, data.frame(
        sample = sample_name,
        group = elements[i+1],
        disease = elements[i+2],
        fibrosis = as.numeric(elements[i+3]),
        nas = as.numeric(elements[i+4]),
        stage = elements[i+5],
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Read count files
  cat("Reading count files...\n")
  count_files <- list.files(pattern="*.counts.txt", full.names=TRUE)
  
  if (length(count_files) == 0) {
    stop("No count files found. Exiting.")
  }
  
  cat("Found", length(count_files), "count files\n")
  
  # Read first file to get gene IDs
  first_file <- count_files[1]
  counts_first <- read.table(first_file, header=TRUE, sep="\t", comment.char="#")
  genes <- counts_first$Geneid
  
  # Create count matrix
  count_matrix <- matrix(0, nrow=length(genes), ncol=length(count_files))
  rownames(count_matrix) <- genes
  colnames(count_matrix) <- sapply(count_files, function(f) {
    gsub(".*/(.*)\\.(counts.txt)$", "\\1", f)
  })
  
  # Read counts in batches
  batch_size <- 20
  for (batch_start in seq(1, length(count_files), batch_size)) {
    batch_end <- min(batch_start + batch_size - 1, length(count_files))
    cat(sprintf("Processing files %d to %d\n", batch_start, batch_end))
    
    for (i in batch_start:batch_end) {
      file <- count_files[i]
      counts <- read.table(file, header=TRUE, sep="\t", comment.char="#")
      count_matrix[,i] <- counts[,ncol(counts)]
    }
    
    # Force garbage collection between batches
    gc()
  }
  
  # Convert to numeric
  storage.mode(count_matrix) <- "numeric"
  
  # Match metadata and count matrix
  cat("Matching metadata with count matrix...\n")
  
  # Keep only samples present in both
  common_samples <- intersect(sample_data$sample, colnames(count_matrix))
  cat("Found", length(common_samples), "samples in both metadata and count files\n")
  
  if (length(common_samples) == 0) {
    cat("ERROR: No matching samples between metadata and count files\n")
    cat("Metadata samples:", paste(head(sample_data$sample), collapse=", "), "...\n")
    cat("Count file samples:", paste(head(colnames(count_matrix)), collapse=", "), "...\n")
    stop("No matching samples")
  }
  
  # Subset count matrix to common samples
  count_matrix <- count_matrix[, common_samples]
  
  # Subset and reorder sample_data
  rownames(sample_data) <- sample_data$sample
  sample_data <- sample_data[common_samples,]
  
  # Create DESeq2 dataset
  cat("Creating DESeq2 dataset...\n")
  dds <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = sample_data,
    design = ~ disease
  )
  
  # Filter out low count genes
  keep <- rowSums(counts(dds) >= 10) >= 3
  dds <- dds[keep,]
  
  # Run DESeq
  cat("Running DESeq2...\n")
  dds <- DESeq(dds)
  
  # Save normalized counts
  cat("Saving normalized counts...\n")
  normalized_counts <- counts(dds, normalized=TRUE)
  write.csv(normalized_counts, "normalized_counts.csv")
  
  # Get results if there are multiple disease types
  unique_diseases <- unique(sample_data$disease)
  
  if (length(unique_diseases) >= 2) {
    cat("Performing comparison between disease types...\n")
    res <- results(dds, contrast=c("disease", unique_diseases[1], unique_diseases[2]))
    res <- res[order(res$padj),]
    write.csv(as.data.frame(res), "deseq2_results.csv")
  } else {
    cat("Only one disease type found. Saving counts only.\n")
    write.csv(data.frame(disease=unique_diseases), "disease_types.csv")
  }
  
  cat("Analysis complete.\n")
  
}, error = function(e) {
  cat("ERROR:", conditionMessage(e), "\n")
  print(sessionInfo())
  quit(status = 1)
})

# Exit successfully
quit(status = 0) 