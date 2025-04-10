# Install required packages if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("GEOquery", quietly = TRUE)) {
  BiocManager::install("GEOquery")
}

if (!requireNamespace("SRAdb", quietly = TRUE)) {
  BiocManager::install("SRAdb")
}

if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}

if (!requireNamespace("readr", quietly = TRUE)) {
  install.packages("readr")
}

# Load libraries
library(GEOquery)
library(SRAdb)
library(dplyr)
library(readr)

# Create directories for data storage
main_dir <- "NAFLD_NASH_RNAseq"
if (!dir.exists(main_dir)) {
  dir.create(main_dir)
}

metadata_dir <- file.path(main_dir, "metadata")
if (!dir.exists(metadata_dir)) {
  dir.create(metadata_dir)
}

raw_data_dir <- file.path(main_dir, "raw_data")
if (!dir.exists(raw_data_dir)) {
  dir.create(raw_data_dir)
}

# Define GEO accession numbers for the studies mentioned
geo_accessions <- c(
  # Human studies
  "GSE135251", # Liu et al. (2020)
  "GSE130970", # Govaere et al. (2020)
  "GSE126848", # Suppli et al. (2019)
  
  # Mouse studies
  "GSE63175",  # Teufel et al. (2016)
  "GSE119723", # Krishnan et al. (2018)
  "GSE179257", # Brunt et al. (2022)
  
  # Comparative studies
  "GSE137449", # Tran et al. (2020)
  "GSE83452"   # Cazanave et al. (2017)
)

# Function to download and process GEO dataset
download_geo_dataset <- function(geo_accession) {
  cat("Processing", geo_accession, "...\n")
  
  # Set output directory for this dataset
  dataset_dir <- file.path(raw_data_dir, geo_accession)
  if (!dir.exists(dataset_dir)) {
    dir.create(dataset_dir)
  }
  
  # Download GEO metadata
  gse <- getGEO(geo_accession, GSEMatrix = TRUE, getGPL = FALSE)
  
  # Extract metadata
  if (length(gse) > 1) {
    # If there are multiple platforms, use the first one
    metadata <- pData(phenoData(gse[[1]]))
  } else {
    metadata <- pData(phenoData(gse))
  }
  
  # Save metadata
  metadata_file <- file.path(metadata_dir, paste0(geo_accession, "_metadata.csv"))
  write.csv(metadata, metadata_file, row.names = FALSE)
  
  # Try to get SRA information if available
  sra_info <- NULL
  try({
    # Create or connect to the SRAdb SQLite database
    sra_dbname <- file.path(main_dir, "SRAmetadb.sqlite")
    if (!file.exists(sra_dbname)) {
      sra_dbname <- getSRAdbFile(destdir = main_dir)
    }
    con <- dbConnect(SQLite(), sra_dbname)
    
    # Find SRX IDs from the metadata (if available)
    srx_pattern <- "SRX[0-9]+"
    source_name_col <- which(grepl("source_name", colnames(metadata), ignore.case = TRUE))
    if (length(source_name_col) > 0) {
      potential_srx <- unlist(regmatches(metadata[, source_name_col], 
                                        gregexpr(srx_pattern, metadata[, source_name_col])))
    } else {
      # Look through all columns for SRX patterns
      all_text <- apply(metadata, 2, paste, collapse = " ")
      potential_srx <- unlist(regmatches(all_text, gregexpr(srx_pattern, all_text)))
    }
    
    if (length(potential_srx) > 0) {
      # Get SRA run information
      sra_info <- listSRAfile(search_terms = paste(potential_srx, collapse = " OR "), 
                              con = con, sraType = "sra")
    }
    
    # Close the SRA database connection
    dbDisconnect(con)
  }, silent = TRUE)
  
  # Extract expression data if available in GEO
  if (geo_accession %in% names(gse)) {
    expr_data <- exprs(gse[[geo_accession]])
    expr_file <- file.path(dataset_dir, paste0(geo_accession, "_expression_matrix.csv"))
    write.csv(expr_data, expr_file)
  }
  
  # Return information about the dataset
  return(list(
    accession = geo_accession,
    sample_count = nrow(metadata),
    sra_info = sra_info
  ))
}

# Process all GEO datasets
results <- lapply(geo_accessions, download_geo_dataset)

# Create a summary of downloaded datasets
summary_df <- data.frame(
  GEO_Accession = sapply(results, function(x) x$accession),
  Sample_Count = sapply(results, function(x) x$sample_count),
  SRA_Available = sapply(results, function(x) !is.null(x$sra_info)),
  SRA_Run_Count = sapply(results, function(x) ifelse(is.null(x$sra_info), 0, nrow(x$sra_info)))
)

summary_file <- file.path(main_dir, "dataset_summary.csv")
write.csv(summary_df, summary_file, row.names = FALSE)

# Define a function to download fastq files for a specific dataset
download_fastq_files <- function(geo_accession, max_samples = 2) {
  cat("Downloading FASTQ files for", geo_accession, "...\n")
  
  # Connect to SRA database
  sra_dbname <- file.path(main_dir, "SRAmetadb.sqlite")
  if (!file.exists(sra_dbname)) {
    sra_dbname <- getSRAdbFile(destdir = main_dir)
  }
  con <- dbConnect(SQLite(), sra_dbname)
  
  # Get SRA info for this GEO accession
  sra_info <- getSRAinfo(search_terms = geo_accession, con = con, sraType = "sra")
  
  # Limit to max_samples for demo purposes
  if (nrow(sra_info) > max_samples) {
    sra_info <- sra_info[1:max_samples, ]
  }
  
  # Download FASTQ files
  if (nrow(sra_info) > 0) {
    output_dir <- file.path(raw_data_dir, geo_accession)
    if (!dir.exists(output_dir)) {
      dir.create(output_dir)
    }
    
    # Get download URLs
    fastq_urls <- getFASTQinfo(sra_info$run, con = con)
    
    # Download files
    for (i in 1:nrow(fastq_urls)) {
      url <- fastq_urls$ftpURL[i]
      filename <- basename(url)
      output_file <- file.path(output_dir, filename)
      
      if (!file.exists(output_file)) {
        try({
          cat("Downloading", filename, "...\n")
          download.file(url, output_file)
        }, silent = TRUE)
      }
    }
  }
  
  # Close SRA database connection
  dbDisconnect(con)
}