#!/usr/bin/env Rscript

#########################################################################
# Interactive DESeq2 Analysis Script for GSE135251
# This script has been optimized for interactive use in RStudio
#########################################################################

# ===== CONFIGURABLE PARAMETERS (modify these for your analysis) =====
# Input files and directories
INPUT_DIR <- "~/bulk_rnaseq_full_results/counts"  # Change to your input directory
OUTPUT_DIR <- "~/Cas13/bulk_RNA/deseq2_results_GSE135251"  # Change to desired output location
COUNT_FILE <- "counts_matrix.csv"          # Count matrix filename
METADATA_FILE <- "sample_metadata.csv"     # Metadata filename

# Analysis parameters
LOG2FC_CUTOFF <- 1                         # Log2 fold change cutoff for significance
PADJ_CUTOFF <- 0.05                        # Adjusted p-value cutoff
MIN_COUNT <- 10                            # Minimum count threshold for filtering
MIN_SAMPLES <- 3                           # Minimum samples passing threshold

# Visualization parameters
TOP_GENES_HEATMAP <- 50                    # Number of top genes to include in heatmaps
TOP_GENES_LABEL <- 20                      # Number of top genes to label in volcano plots
CREATE_PLOTS <- TRUE                       # Set to FALSE to skip plot generation

# Variables to analyze
ANALYSIS_VARIABLES <- c("disease", "stage", "nas") # Variables for differential expression analysis

# ===== LOAD LIBRARIES =====
# Load required libraries with error checking
required_packages <- c("DESeq2", "ggplot2", "pheatmap", "RColorBrewer", "ggrepel", "dplyr", "biomaRt")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(paste0("Package '", pkg, "' is not installed. Would you like to install it? (y/n): "))
    answer <- readline()
    if (tolower(answer) == "y") {
      if (pkg %in% c("DESeq2", "biomaRt")) {
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
          install.packages("BiocManager")
        }
        BiocManager::install(pkg)
      } else {
        install.packages(pkg)
      }
    } else {
      stop(paste0("Package '", pkg, "' is required but not installed."))
    }
  }
  library(pkg, character.only = TRUE)
}

# Function to check directory and create if needed
ensure_dir <- function(dir_path) {
  if (!dir.exists(dir_path)) {
    cat(paste0("Creating directory: ", dir_path, "\n"))
    dir.create(dir_path, recursive = TRUE)
  }
}

# Create output directory
ensure_dir(OUTPUT_DIR)

# Function to create full path
path <- function(...) {
  file.path(...)
}

# ===== INTERACTIVE SESSION SETUP =====
# Allow user to change directories interactively in RStudio
if (interactive()) {
  cat("Running in interactive mode.\n")
  
  # Ask user if they want to change default settings
  cat("Would you like to change the default settings? (y/n): ")
  answer <- readline()
  if (tolower(answer) == "y") {
    cat("Current input directory:", INPUT_DIR, "\n")
    cat("Would you like to change the input directory? (y/n): ")
    answer <- readline()
    if (tolower(answer) == "y") {
      cat("Enter new input directory: ")
      INPUT_DIR <- readline()
      ensure_dir(INPUT_DIR)
    }
    
    cat("Current output directory:", OUTPUT_DIR, "\n")
    cat("Would you like to change the output directory? (y/n): ")
    answer <- readline()
    if (tolower(answer) == "y") {
      cat("Enter new output directory: ")
      OUTPUT_DIR <- readline()
      ensure_dir(OUTPUT_DIR)
    }
    
    cat("Current count file:", COUNT_FILE, "\n")
    cat("Would you like to change the count file? (y/n): ")
    answer <- readline()
    if (tolower(answer) == "y") {
      cat("Enter new count file name: ")
      COUNT_FILE <- readline()
    }
    
    cat("Current metadata file:", METADATA_FILE, "\n")
    cat("Would you like to change the metadata file? (y/n): ")
    answer <- readline()
    if (tolower(answer) == "y") {
      cat("Enter new metadata file name: ")
      METADATA_FILE <- readline()
    }
  }
  
  # Set working directory for interactive session
  setwd(INPUT_DIR)
  cat("Working directory set to:", getwd(), "\n")
} else {
  # For non-interactive sessions, use the provided directory
  setwd(INPUT_DIR)
}

# ===== DATA LOADING =====
# Error handling for file loading
tryCatch({
  # Read input files
  message("Reading count matrix and metadata...")
  counts_file_path <- path(INPUT_DIR, COUNT_FILE)
  metadata_file_path <- path(INPUT_DIR, METADATA_FILE)
  
  if (!file.exists(counts_file_path)) {
    stop(paste0("Count file not found: ", counts_file_path))
  }
  if (!file.exists(metadata_file_path)) {
    stop(paste0("Metadata file not found: ", metadata_file_path))
  }
  
  counts_matrix <- read.csv(counts_file_path, row.names = 1)
  sample_metadata <- read.csv(metadata_file_path, stringsAsFactors = FALSE)
  
  # Explicitly convert all analysis variables to factors
  message("Converting analysis variables to factors...")
  for (variable in ANALYSIS_VARIABLES) {
    if (variable %in% colnames(sample_metadata)) {
      sample_metadata[[variable]] <- as.factor(sample_metadata[[variable]])
      message(paste("Converted", variable, "to factor with levels:", 
                    paste(levels(sample_metadata[[variable]]), collapse=", ")))
    }
  }
  
  # Initialize gene ID conversion settings
  convert_gene_ids <- TRUE
  
  # Create a mapping of Ensembl IDs to gene symbols
  ensembl_to_symbol <- NULL
  if (convert_gene_ids) {
    message("Converting Ensembl IDs to gene symbols...")
    
    # Extract the Ensembl IDs from rownames
    ensembl_ids <- rownames(counts_matrix)
    
    # Check if they look like Ensembl IDs
    if (all(grepl("^ENS", ensembl_ids[1:min(10, length(ensembl_ids))]))) {
      message("Detected Ensembl IDs. Getting gene symbols from biomaRt...")
      
      # Extract the Ensembl ID part before any version number (e.g., ENSG00000123456.1 -> ENSG00000123456)
      ensembl_ids_clean <- sub("\\.[0-9]+$", "", ensembl_ids)
      
      # Connect to Ensembl database
      message("Connecting to Ensembl database...")
      mart <- tryCatch({
        useMart("ensembl", dataset = "hsapiens_gene_ensembl")
      }, error = function(e) {
        message("Error connecting to Ensembl database. Will attempt to use archived version.")
        # Try an archived version if current fails
        useMart("ENSEMBL_MART_ENSEMBL", host = "https://mar2016.archive.ensembl.org", 
                dataset = "hsapiens_gene_ensembl")
      })
      
      # Get gene symbols
      message("Fetching gene symbols...")
      ensembl_to_symbol <- tryCatch({
        gene_info <- getBM(
          attributes = c("ensembl_gene_id", "external_gene_name", "description"),
          filters = "ensembl_gene_id",
          values = ensembl_ids_clean,
          mart = mart
        )
        
        # Create a map from Ensembl ID to gene symbol
        id_map <- setNames(gene_info$external_gene_name, gene_info$ensembl_gene_id)
        
        # For Ensembl IDs without a symbol, keep the original ID
        missing_symbols <- ensembl_ids_clean[!ensembl_ids_clean %in% gene_info$ensembl_gene_id]
        if (length(missing_symbols) > 0) {
          message(paste0(length(missing_symbols), " Ensembl IDs could not be mapped to gene symbols"))
          # Add the missing IDs back to the map using the original ID
          missing_map <- setNames(missing_symbols, missing_symbols)
          id_map <- c(id_map, missing_map)
        }
        
        id_map
      }, error = function(e) {
        message("Error fetching gene symbols: ", conditionMessage(e))
        message("Will continue with original IDs")
        NULL
      })
      
      if (!is.null(ensembl_to_symbol)) {
        message("Successfully mapped ", length(ensembl_to_symbol), " genes to symbols")
        
        # Create a data frame with both IDs for reference
        gene_id_mapping <- data.frame(
          ensembl_id = names(ensembl_to_symbol),
          gene_symbol = unname(ensembl_to_symbol),
          stringsAsFactors = FALSE
        )
        
        # Save the mapping for reference
        write.csv(gene_id_mapping, path(OUTPUT_DIR, "gene_id_mapping.csv"), row.names = FALSE)
        message("Saved gene ID mapping to ", path(OUTPUT_DIR, "gene_id_mapping.csv"))
      }
    } else {
      message("Row names don't appear to be Ensembl IDs. Skipping conversion.")
      convert_gene_ids <- FALSE
    }
  }
  
  # Basic quality checks
  message("Performing basic quality checks...")
  message(paste("Count matrix dimensions:", nrow(counts_matrix), "genes x", ncol(counts_matrix), "samples"))
  message(paste("Metadata dimensions:", nrow(sample_metadata), "samples x", ncol(sample_metadata), "features"))
  
  # Display first few rows of each for verification
  cat("\nFirst 5 rows of count matrix:\n")
  print(head(counts_matrix[, 1:min(5, ncol(counts_matrix))], 5))
  
  cat("\nFirst 5 rows of metadata:\n")
  print(head(sample_metadata, 5))
  
  # Check for sample mismatches - this is important and often a source of errors
  cat("\nChecking sample matching between count matrix and metadata...\n")
  count_samples <- colnames(counts_matrix)
  metadata_samples <- sample_metadata$sample
  
  missing_in_counts <- setdiff(metadata_samples, count_samples)
  missing_in_metadata <- setdiff(count_samples, metadata_samples)
  
  if (length(missing_in_counts) > 0) {
    warning(paste0(length(missing_in_counts), " samples in metadata are missing from count matrix:\n", 
                  paste(head(missing_in_counts, 5), collapse=", "), 
                  if(length(missing_in_counts) > 5) "..." else ""))
  }
  
  if (length(missing_in_metadata) > 0) {
    warning(paste0(length(missing_in_metadata), " samples in count matrix are missing from metadata:\n", 
                  paste(head(missing_in_metadata, 5), collapse=", "), 
                  if(length(missing_in_metadata) > 5) "..." else ""))
  }
  
  # Ensure sample order matches between metadata and count matrix
  message("Matching samples in metadata and count matrix...")
  common_samples <- intersect(count_samples, metadata_samples)
  if (length(common_samples) == 0) {
    stop("No matching samples between count matrix and metadata!")
  } else {
    message(paste(length(common_samples), "matching samples found"))
  }
  
  # Subset to common samples
  counts_matrix <- counts_matrix[, common_samples]
  rownames(sample_metadata) <- sample_metadata$sample
  sample_metadata <- sample_metadata[common_samples, ]
  
  # Verify matching after subsetting
  if(!all(rownames(sample_metadata) == colnames(counts_matrix))) {
    stop("Sample names in metadata still do not match column names in count matrix after subsetting")
  }
  
  # ===== DESEQ2 ANALYSIS =====
  # Convert count matrix to integer matrix (DESeq2 requirement)
  message("Preparing count matrix for DESeq2...")
  counts_matrix <- round(as.matrix(counts_matrix))
  mode(counts_matrix) <- "integer"
  
  # Function to run DESeq2 analysis for a given variable
  run_deseq2_analysis <- function(variable_name, counts_matrix, sample_metadata, ensembl_to_symbol, convert_gene_ids) {
    message(paste0("\n===== RUNNING DESEQ2 ANALYSIS FOR: ", variable_name, " ====="))
    
    # Create output directory for this variable
    variable_output_dir <- path(OUTPUT_DIR, variable_name)
    ensure_dir(variable_output_dir)
    
    # Create design formula based on variable
    design_formula <- as.formula(paste0("~ ", variable_name))
    message(paste("Using design formula:", deparse(design_formula)))
    
    # Create DESeq2 object
    dds <- DESeqDataSetFromMatrix(
      countData = counts_matrix,
      colData = sample_metadata,
      design = design_formula
    )
    
    # Pre-filter low count genes
    message("Filtering low count genes...")
    keep <- rowSums(counts(dds) >= MIN_COUNT) >= MIN_SAMPLES
    dds <- dds[keep, ]
    message(paste(sum(keep), "genes remaining after filtering"))
    
    # Run DESeq2
    message("Running DESeq2 analysis. This may take a while...")
    dds <- DESeq(dds)
    message("DESeq2 analysis complete!")
    
    # Save normalized counts to output directory
    message("Saving normalized counts...")
    normalized_counts <- counts(dds, normalized=TRUE)
    
    # Convert rownames to gene symbols if mapping is available
    if (!is.null(ensembl_to_symbol) && convert_gene_ids) {
      # Get original rownames (Ensembl IDs)
      original_ids <- rownames(normalized_counts)
      
      # Clean Ensembl IDs by removing version numbers
      original_ids_clean <- sub("\\.[0-9]+$", "", original_ids)
      
      # Map to gene symbols
      gene_symbols <- ensembl_to_symbol[original_ids_clean]
      
      # For any NA values, use the original Ensembl ID
      gene_symbols[is.na(gene_symbols)] <- original_ids[is.na(gene_symbols)]
      
      # Store both IDs in the output to maintain information
      normalized_counts_with_symbols <- normalized_counts
      normalized_counts_with_symbols <- cbind(
        data.frame(
          ensembl_id = original_ids,
          gene_symbol = gene_symbols,
          stringsAsFactors = FALSE
        ),
        as.data.frame(normalized_counts_with_symbols)
      )
      
      write.csv(normalized_counts_with_symbols, path(variable_output_dir, "deseq2_normalized_counts.csv"), row.names = FALSE)
    } else {
      # Save without conversion
      write.csv(normalized_counts, path(variable_output_dir, "deseq2_normalized_counts.csv"))
    }
    
    # Get unique values of the variable and check if we can make comparisons
    unique_values <- unique(sample_metadata[[variable_name]])
    message(paste("Found", length(unique_values), variable_name, "types:", paste(unique_values, collapse=", ")))
    
    # ===== VISUALIZATIONS =====
    # Only create plots if specified
    if (CREATE_PLOTS) {
      # Create PCA plot
      message("Generating PCA plot...")
      tryCatch({
        vsd <- vst(dds, blind=FALSE)
        
        # Create the PCA plot and save to file
        pcaData <- plotPCA(vsd, intgroup=c(variable_name, "group"), returnData=TRUE)
        percentVar <- round(100 * attr(pcaData, "percentVar"))
        pca_plot <- ggplot(pcaData, aes_string(x="PC1", y="PC2", color=variable_name, shape="group")) +
          geom_point(size=3) +
          xlab(paste0("PC1: ",percentVar[1],"% variance")) +
          ylab(paste0("PC2: ",percentVar[2],"% variance")) +
          geom_text_repel(aes(label=name), size=3, max.overlaps=50) +
          theme_bw() +
          ggtitle(paste("PCA Plot of RNA-seq Samples by", variable_name))
        
        # Save to file
        pdf_file <- path(variable_output_dir, "deseq2_pca_plot.pdf")
        message("Saving PCA plot to: ", pdf_file)
        pdf(pdf_file, width=10, height=8)
        print(pca_plot)
        dev.off()
        
        message("PCA plot completed successfully")
      }, error = function(e) {
        message("Error generating PCA plot: ", conditionMessage(e))
      })
      
      # Create sample-to-sample distance heatmap
      message("Generating sample distance heatmap...")
      tryCatch({
        sampleDists <- dist(t(assay(vsd)))
        sampleDistMatrix <- as.matrix(sampleDists)
        colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
        
        # Create annotation data
        annotation_data <- sample_metadata[, c(variable_name, "group"), drop=FALSE]
        colnames(annotation_data)[1] <- variable_name  # Ensure column name matches variable
        
        # Save to file
        pdf(path(variable_output_dir, "deseq2_sample_heatmap.pdf"), width=12, height=10)
        pheatmap(
          sampleDistMatrix,
          clustering_distance_rows=sampleDists,
          clustering_distance_cols=sampleDists,
          col=colors,
          annotation_row=annotation_data,
          main=paste("Sample-to-Sample Distance by", variable_name)
        )
        dev.off()
        message("Sample heatmap completed successfully")
      }, error = function(e) {
        message("Error generating sample heatmap: ", conditionMessage(e))
      })
      
      # ===== DIFFERENTIAL EXPRESSION ANALYSIS =====
      # Perform differential expression if multiple values exist
      if(length(unique_values) >= 2) {
        # Make pairwise comparisons
        message("Performing pairwise comparisons between", variable_name, "types...")
        
        # Create results directory within output directory
        results_dir <- path(variable_output_dir, "deseq2_results")
        ensure_dir(results_dir)
        
        # Special handling for "nas" variable
        if (variable_name == "nas") {
          # Make sure the levels are correctly ordered with "0" as the reference level
          message("Setting '0' as the reference level for nas variable")
          if (!"0" %in% levels(sample_metadata$nas)) {
            warning("Level '0' not found in nas variable. Please check your metadata.")
          } else {
            # Reorder levels to ensure "0" is first (reference level)
            current_levels <- levels(sample_metadata$nas)
            new_levels <- c("0", current_levels[current_levels != "0"])
            sample_metadata$nas <- factor(sample_metadata$nas, levels = new_levels)
            # Update the dds object with the new factor levels
            dds$nas <- factor(dds$nas, levels = new_levels)
            message("Reference level set to '0' for nas variable")
          }
        }
        
        # Comparisons data frame to store summary
        comparisons_summary <- data.frame(
          comparison = character(),
          total_DE_genes = integer(),
          upregulated = integer(),
          downregulated = integer(),
          stringsAsFactors = FALSE
        )
        
        # Let user select comparisons if in interactive mode
        if (interactive() && length(unique_values) > 2) {
          cat("\nMultiple", variable_name, "types found. You can choose specific comparisons.\n")
          cat(variable_name, "types:\n")
          for (i in 1:length(unique_values)) {
            cat(i, ": ", unique_values[i], "\n", sep="")
          }
          
          cat("\nWould you like to select specific", variable_name, "types to compare? (y/n): ")
          answer <- readline()
          
          if (tolower(answer) == "y") {
            compare_all <- FALSE
            cat("Enter first", variable_name, "type (number): ")
            value1_idx <- as.numeric(readline())
            cat("Enter second", variable_name, "type (number): ")
            value2_idx <- as.numeric(readline())
            
            if (is.na(value1_idx) || is.na(value2_idx) || 
                value1_idx < 1 || value1_idx > length(unique_values) ||
                value2_idx < 1 || value2_idx > length(unique_values)) {
              stop("Invalid", variable_name, "type selection")
            }
            
            # Force a specific comparison
            comparison_pairs <- list(c(value1_idx, value2_idx))
          } else {
            compare_all <- TRUE
          }
        } else {
          compare_all <- TRUE
        }
        
        # Generate comparison pairs
        if (compare_all) {
          comparison_pairs <- list()
          
          # Special handling for "nas" variable - only compare 0 to other values
          if (variable_name == "nas") {
            # Find the index of "0"
            zero_idx <- which(unique_values == "0")
            if (length(zero_idx) == 0) {
              message("Warning: Could not find '0' value for nas variable. Using first value as reference.")
              zero_idx <- 1
            }
            
            # Only create pairs comparing 0 to other values
            for (j in 1:length(unique_values)) {
              if (j != zero_idx) {
                comparison_pairs <- c(comparison_pairs, list(c(zero_idx, j)))
              }
            }
            message("For 'nas' variable: only comparing '0' (control) against other values")
          } else {
            # Normal pairwise comparisons for other variables
            for (i in 1:(length(unique_values)-1)) {
              for (j in (i+1):length(unique_values)) {
                comparison_pairs <- c(comparison_pairs, list(c(i, j)))
              }
            }
          }
        }
        
        # Process each comparison
        for (pair in comparison_pairs) {
          i <- pair[1]
          j <- pair[2]
          value1 <- unique_values[i]
          value2 <- unique_values[j]
          
          comparison_name <- paste0(value1, "_vs_", value2)
          message(paste("Analyzing:", comparison_name))
          
          # Create results directory for this comparison
          comp_dir <- path(results_dir, comparison_name)
          ensure_dir(comp_dir)
          
          # Get results
          # First check the actual levels in the DESeq2 object to ensure we're using valid levels
          dds_levels <- levels(dds[[variable_name]])
          message(paste("Factor levels in DESeq2 object for", variable_name, ":", paste(dds_levels, collapse=", ")))
          
          # Try to get results using error handling
          res <- tryCatch({
            # Try standard contrast
            results(dds, contrast=c(variable_name, as.character(value1), as.character(value2)))
          }, error = function(e) {
            message(paste("Error with standard contrast method:", conditionMessage(e)))
            message("Trying alternative approaches...")
            
            # Check if the values exist in the levels (exact match)
            if(as.character(value1) %in% dds_levels && as.character(value2) %in% dds_levels) {
              message(paste("Using exact levels:", value1, "vs", value2))
              return(results(dds, contrast=c(variable_name, as.character(value1), as.character(value2))))
            }
            
            # Try to get results using name parameter instead (default way DESeq2 names results)
            result_names <- resultsNames(dds)
            message(paste("Available result names:", paste(result_names, collapse=", ")))
            
            # Look for best matching result name
            expected_name1 <- paste0(variable_name, "_", value2, "_vs_", value1)
            expected_name2 <- paste0(variable_name, "_", value2, "_vs_", dds_levels[1])
            
            if(expected_name1 %in% result_names) {
              message(paste("Using result name:", expected_name1))
              return(results(dds, name=expected_name1))
            } else if(expected_name2 %in% result_names) {
              message(paste("Using result name:", expected_name2))
              return(results(dds, name=expected_name2))
            } else {
              # Last resort - use first result that contains the variable name
              var_results <- result_names[grepl(variable_name, result_names)]
              if(length(var_results) > 0) {
                message(paste("Using first available result for", variable_name, ":", var_results[1]))
                return(results(dds, name=var_results[1]))
              } else {
                stop(paste("No valid results found for", variable_name, "comparison between", value1, "and", value2))
              }
            }
          })
          
          res <- res[order(res$padj),]
          
          # Convert rownames to gene symbols if mapping is available
          if (!is.null(ensembl_to_symbol) && convert_gene_ids) {
            # Create a results data frame with both IDs
            res_df <- as.data.frame(res)
            
            # Get original rownames (Ensembl IDs)
            original_ids <- rownames(res_df)
            
            # Clean Ensembl IDs by removing version numbers
            original_ids_clean <- sub("\\.[0-9]+$", "", original_ids)
            
            # Map to gene symbols
            gene_symbols <- ensembl_to_symbol[original_ids_clean]
            
            # For any NA values, use the original Ensembl ID
            gene_symbols[is.na(gene_symbols)] <- original_ids[is.na(gene_symbols)]
            
            # Add the IDs as columns
            res_df$ensembl_id <- original_ids
            res_df$gene_symbol <- gene_symbols
            
            # Save results with both IDs
            write.csv(res_df, path(comp_dir, "differential_expression.csv"), row.names = FALSE)
          } else {
            # Save without conversion
            write.csv(as.data.frame(res), path(comp_dir, "differential_expression.csv"))
          }
          
          if (CREATE_PLOTS) {
            # MA plot
            pdf(path(comp_dir, "ma_plot.pdf"), width=8, height=6)
            plotMA(res, ylim=c(-5, 5), main=paste("MA Plot:", value1, "vs", value2))
            dev.off()
            
            # Volcano plot
            message(paste("Creating volcano plot for:", comparison_name))
            tryCatch({
              res_df <- as.data.frame(res)
              
              # Add gene name/ID column for labeling
              if (!is.null(ensembl_to_symbol) && convert_gene_ids) {
                # Clean Ensembl IDs by removing version numbers
                original_ids_clean <- sub("\\.[0-9]+$", "", rownames(res_df))
                
                # Map to gene symbols 
                gene_symbols <- ensembl_to_symbol[original_ids_clean]
                
                # For any NA values, use the original Ensembl ID
                gene_symbols[is.na(gene_symbols)] <- rownames(res_df)[is.na(gene_symbols)]
                
                res_df$gene <- gene_symbols
              } else {
                res_df$gene <- rownames(res_df)
              }
              
              # Add significance information
              res_df$significant <- ifelse(
                abs(res_df$log2FoldChange) > LOG2FC_CUTOFF & res_df$padj < PADJ_CUTOFF, 
                "Significant", 
                "Not Significant"
              )
              
              # Top genes for labeling
              top_genes <- head(res_df[order(res_df$padj),], TOP_GENES_LABEL)
              
              # Create volcano plot
              volcano_plot <- ggplot(res_df, aes(x=log2FoldChange, y=-log10(padj), color=significant)) +
                geom_point(alpha=0.6) +
                scale_color_manual(values=c("Not Significant"="grey", "Significant"="red")) +
                theme_bw() +
                theme(legend.position="bottom") +
                geom_text_repel(
                  data=top_genes,
                  aes(label=gene),
                  max.overlaps=50,
                  box.padding=0.5
                ) +
                geom_vline(xintercept=c(-LOG2FC_CUTOFF, LOG2FC_CUTOFF), linetype="dashed", color="blue") +
                geom_hline(yintercept=-log10(PADJ_CUTOFF), linetype="dashed", color="blue") +
                labs(title=paste("Volcano Plot:", value1, "vs", value2),
                     x="log2 Fold Change",
                     y="-log10 Adjusted p-value")
              
              # Save to file
              pdf_file <- path(comp_dir, "volcano_plot.pdf")
              message("Saving volcano plot to: ", pdf_file)
              pdf(pdf_file, width=10, height=8)
              print(volcano_plot)
              dev.off()
              message("Volcano plot completed successfully")
            }, error = function(e) {
              message("Error generating volcano plot: ", conditionMessage(e))
            })
            
            # Heatmap of top DE genes
            message(paste("Creating heatmap for:", comparison_name))
            tryCatch({
              sig_genes <- which(res$padj < PADJ_CUTOFF & abs(res$log2FoldChange) > LOG2FC_CUTOFF)
              top_genes_count <- min(TOP_GENES_HEATMAP, length(sig_genes))
              
              if(top_genes_count > 0) {
                top_de_genes <- rownames(res)[sig_genes][1:top_genes_count]
                
                # Get heatmap data
                heatmap_data <- assay(vsd)[top_de_genes,]
                
                # Convert rownames for the heatmap
                if (!is.null(ensembl_to_symbol) && convert_gene_ids) {
                  # Clean Ensembl IDs by removing version numbers
                  original_ids_clean <- sub("\\.[0-9]+$", "", rownames(heatmap_data))
                  
                  # Map to gene symbols
                  gene_symbols <- ensembl_to_symbol[original_ids_clean]
                  
                  # For any NA values, use the original Ensembl ID
                  gene_symbols[is.na(gene_symbols)] <- rownames(heatmap_data)[is.na(gene_symbols)]
                  
                  # Set new rownames
                  rownames(heatmap_data) <- gene_symbols
                }
                
                # Create heatmap
                heatmap_title <- paste("Top DE Genes:", value1, "vs", value2)
                
                # Save to file
                pdf_file <- path(comp_dir, "heatmap_top_de_genes.pdf")
                message("Saving heatmap to: ", pdf_file)
                pdf(pdf_file, width=12, height=10)
                pheatmap(
                  heatmap_data,
                  cluster_rows=TRUE,
                  cluster_cols=TRUE,
                  show_rownames=TRUE,
                  annotation_col=sample_metadata[,c(variable_name, "group")],
                  scale="row",
                  main=heatmap_title
                )
                dev.off()
                message("Heatmap completed successfully")
              } else {
                message("No significant genes found for this comparison - skipping heatmap")
              }
            }, error = function(e) {
              message("Error generating heatmap: ", conditionMessage(e))
            })
          }
          
          # Calculate summary statistics
          sig_genes <- sum(!is.na(res$padj) & res$padj < PADJ_CUTOFF & abs(res$log2FoldChange) > LOG2FC_CUTOFF, na.rm=TRUE)
          up_genes <- sum(!is.na(res$padj) & res$padj < PADJ_CUTOFF & res$log2FoldChange > LOG2FC_CUTOFF, na.rm=TRUE)
          down_genes <- sum(!is.na(res$padj) & res$padj < PADJ_CUTOFF & res$log2FoldChange < -LOG2FC_CUTOFF, na.rm=TRUE)
          
          # Add to summary
          comparisons_summary <- rbind(comparisons_summary, data.frame(
            comparison = comparison_name,
            total_DE_genes = sig_genes,
            upregulated = up_genes,
            downregulated = down_genes
          ))
        }
        
        # Save comparison summary
        write.csv(comparisons_summary, path(variable_output_dir, "deseq2_comparison_summary.csv"), row.names=FALSE)
        
        # Print a summary of the results
        message("\nDifferential Expression Summary for", variable_name, ":")
        print(comparisons_summary)
        
      } else {
        message("Only one", variable_name, "type found. Skipping differential expression analysis.")
      }
    } else {
      message("Skipping plot generation as specified by CREATE_PLOTS=FALSE")
    }
    
    # Return the DESeq2 object
    return(dds)
  }
  
  # Run DESeq2 analysis for each variable
  dds_list <- list()
  for (variable in ANALYSIS_VARIABLES) {
    # Check if the variable exists in the metadata
    if (variable %in% colnames(sample_metadata)) {
      dds_list[[variable]] <- run_deseq2_analysis(variable, counts_matrix, sample_metadata, ensembl_to_symbol, convert_gene_ids)
    } else {
      message(paste("Variable", variable, "not found in metadata. Skipping analysis."))
    }
  }
  
  # Save session info for reproducibility
  writeLines(capture.output(sessionInfo()), path(OUTPUT_DIR, "deseq2_session_info.txt"))
  
  message("\nAll DESeq2 analyses complete!")
  message(paste("Results saved to:", OUTPUT_DIR))
  
}, error = function(e) {
  message("\nERROR: ", conditionMessage(e))
  message("\nAn error occurred during analysis. Check the error message above.")
  
  # Print debug info in case of error
  if (exists("dds")) {
    message("\nDESeq2 object information:")
    print(dds)
  }
  
  message("\nSession info for debugging:")
  print(sessionInfo())
})

# If running interactively, print a final message
if (interactive()) {
  cat("\nScript execution completed.\n")
  cat("You can review the results in:", OUTPUT_DIR, "\n")
}

