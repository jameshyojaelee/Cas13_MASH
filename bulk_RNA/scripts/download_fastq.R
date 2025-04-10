# Load libraries
library(GEOquery)
library(SRAdb)
library(dplyr)
library(readr)
library(RSQLite) 

# Create directories for data storage
main_dir <- "NAFLD_NASH_RNAseq"
if (!dir.exists(main_dir)) {
  dir.create(main_dir)
  cat("Created directory:", main_dir, "\n")
}

metadata_dir <- file.path(main_dir, "metadata")
if (!dir.exists(metadata_dir)) {
  dir.create(metadata_dir)
  cat("Created directory:", metadata_dir, "\n")
}

raw_data_dir <- file.path(main_dir, "raw_data")
if (!dir.exists(raw_data_dir)) {
  dir.create(raw_data_dir)
  cat("Created directory:", raw_data_dir, "\n")
}

geo_accessions <- c(
  "GSE135251", # Liu et al. (2020)
  "GSE130970", # Govaere et al. (2020) - Transcriptomic Profiling Across NAFLD Spectrum
  "GSE126848", # Suppli et al. (2019)
)

# Function to check for SRAmetadb.sqlite and prompt for manual download if not found
check_sra_db <- function(destdir) {
  sra_dbname <- file.path(destdir, "SRAmetadb.sqlite")
  if (!file.exists(sra_dbname)) {
    stop("SRAmetadb.sqlite not found. Please download it manually from https://gbnci.cancer.gov/backup/SRAmetadb.sqlite.gz, extract it, and place it in the ", destdir, " directory.")
  }
  return(sra_dbname)
}

# Function to download FASTQ files for a specific dataset
download_fastq_files <- function(geo_accession, max_samples = 2) {
  cat("Processing", geo_accession, "...\n")
  
  # Connect to SRA database
  tryCatch({
    sra_dbname <- check_sra_db(destdir = main_dir)
    cat("Using SRA database:", sra_dbname, "\n")
    con <- dbConnect(SQLite(), sra_dbname)
    
    # First, get SRA study accession from GEO accession using GEOquery
    cat("Fetching GEO metadata for", geo_accession, "\n")
    geo_data <- tryCatch({
      getGEO(geo_accession, GSEMatrix = FALSE)
    }, error = function(e) {
      cat("Error fetching GEO data:", e$message, "\n")
      return(NULL)
    })
    
    if (is.null(geo_data)) {
      cat("Could not get GEO data for", geo_accession, "\n")
      dbDisconnect(con)
      return()
    }
    
    # Extract SRA identifiers
    sra_links <- lapply(GSMList(geo_data), function(x) {
      Meta(x)$relation
    })
    
    # Get SRA study identifiers
    sra_ids <- c()
    for (link_list in sra_links) {
      if (!is.null(link_list)) {
        # Extract SRP or similar identifiers from relation metadata
        ids <- grep("SRP|ERP|DRP", link_list, value = TRUE)
        if (length(ids) > 0) {
          # Extract actual SRP IDs
          extracted_ids <- gsub(".*=(SRP[0-9]+).*", "\\1", ids)
          sra_ids <- c(sra_ids, extracted_ids)
        }
      }
    }
    
    # Remove duplicates
    sra_ids <- unique(sra_ids)
    
    if (length(sra_ids) == 0) {
      cat("No SRA identifiers found for", geo_accession, "\n")
      dbDisconnect(con)
      return()
    }
    
    cat("Found SRA identifiers for", geo_accession, ":", paste(sra_ids, collapse=", "), "\n")
    
    # For each SRA study ID, get the runs
    for (sra_study in sra_ids) {
      cat("Processing SRA study:", sra_study, "\n")
      
      # Query directly using SQL
      query <- paste0("SELECT * FROM sra WHERE study_accession='", sra_study, "'")
      sra_info <- dbGetQuery(con, query)
      
      if (nrow(sra_info) == 0) {
        cat("No SRA entries found for study", sra_study, "\n")
        next
      }
      
      cat("Found", nrow(sra_info), "SRA entries for study", sra_study, "\n")
      
      # Limit to max_samples for demo purposes
      if (nrow(sra_info) > max_samples) {
        sra_info <- sra_info[1:max_samples, ]
        cat("Limited to", max_samples, "samples\n")
      }
      
      # Create output directory for this GEO accession
      output_dir <- file.path(raw_data_dir, geo_accession)
      if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
        cat("Created directory:", output_dir, "\n")
      }
      
      # Get FASTQ download URLs for each run
      for (i in 1:nrow(sra_info)) {
        run_id <- sra_info$run_accession[i]
        cat("Processing run:", run_id, "\n")
        
        # Query fastq table directly
        fastq_query <- paste0("SELECT * FROM fastq WHERE run_accession='", run_id, "'")
        fastq_urls <- dbGetQuery(con, fastq_query)
        
        if (nrow(fastq_urls) == 0) {
          cat("No FASTQ URLs found for run", run_id, "\n")
          next
        }
        
        # Download files
        for (j in 1:nrow(fastq_urls)) {
          url <- fastq_urls$url[j]
          if (is.na(url) || url == "") {
            cat("No download URL for run", run_id, "\n")
            next
          }
          
          filename <- basename(url)
          output_file <- file.path(output_dir, filename)
          
          if (!file.exists(output_file)) {
            cat("Downloading", filename, "from", url, "...\n")
            tryCatch({
              download.file(url, output_file, mode = "wb")
              cat("Download complete:", output_file, "\n")
            }, error = function(e) {
              cat("Error downloading", filename, ":", e$message, "\n")
            })
          } else {
            cat("File already exists:", output_file, "\n")
          }
        }
      }
    }
    
    # Close SRA database connection
    dbDisconnect(con)
    
  }, error = function(e) {
    cat("Error processing", geo_accession, ":", e$message, "\n")
  })
}

# Process each GEO accession
cat("Starting to process", length(geo_accessions), "GEO accessions\n")
for (geo_accession in geo_accessions) {
  download_fastq_files(geo_accession)
}
cat("Processing complete\n")