# Load libraries
library(GEOquery)
library(RSQLite)
library(dplyr)
library(readr)

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

# Define GEO accession numbers for all studies
geo_accessions <- c(
  # Previous studies
  "GSE135251", # Liu et al. (2020)
  "GSE130970", # Govaere et al. (2020) - Transcriptomic Profiling Across NAFLD Spectrum
  "GSE126848", # Suppli et al. (2019)
  "GSE63175",  # Teufel et al. (2016) - Comparative Hepatic Transcriptome Signatures
  "GSE119723", # Krishnan et al. (2018)
  "GSE179257", # Brunt et al. (2022)
  "GSE137449", # Tran et al. (2020) - Cross-Species Transcriptome Analysis
  "GSE83452",  # Cazanave et al. (2017)
  
  # Additional studies
  "GSE115469", # Single-Cell RNA Sequencing in Human NASH
  "GSE106737", # Additional Cross-Species Transcriptome Analysis
  "GSE48452",  # Additional Comparative Hepatic Transcriptome Signatures
  "GSE49541",  # Additional Transcriptomic Profiling Across NAFLD Spectrum
  "GSE151158", # Additional human NASH transcriptome analysis
  "GSE123876"  # Mouse model NAFLD/NASH progression
)

# Function to check for SRAmetadb.sqlite
check_sra_db <- function(destdir) {
  sra_dbname <- file.path(destdir, "SRAmetadb.sqlite")
  if (!file.exists(sra_dbname)) {
    stop("SRAmetadb.sqlite not found. Please download it manually from https://gbnci.cancer.gov/backup/SRAmetadb.sqlite.gz, extract it, and place it in the ", destdir, " directory.")
  }
  return(sra_dbname)
}

# Improved function to extract SRA IDs from GEO
extract_sra_from_geo <- function(geo_accession) {
  cat("Extracting SRA IDs for", geo_accession, "...\n")
  
  # Get GEO data
  geo_data <- tryCatch({
    getGEO(geo_accession, GSEMatrix = FALSE)
  }, error = function(e) {
    cat("Error fetching GEO data:", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(geo_data)) {
    cat("Could not get GEO data for", geo_accession, "\n")
    return(NULL)
  }
  
  # Try multiple approaches to find SRA IDs
  sra_ids <- c()
  
  # Approach 1: Extract from relations in samples
  cat("Trying to extract SRA IDs from sample relations...\n")
  sra_links <- lapply(GSMList(geo_data), function(x) {
    Meta(x)$relation
  })
  
  for (link_list in sra_links) {
    if (!is.null(link_list)) {
      # Look for SRP IDs
      srp_matches <- grep("SRP|ERP|DRP", link_list, value = TRUE)
      if (length(srp_matches) > 0) {
        srp_ids <- gsub(".*=(SRP[0-9]+).*", "\\1", srp_matches)
        sra_ids <- c(sra_ids, srp_ids)
      }
    }
  }
  
  # Approach 2: Look in the overall dataset metadata
  if (length(sra_ids) == 0) {
    cat("Trying to extract SRA IDs from dataset metadata...\n")
    dataset_meta <- Meta(geo_data)
    
    # Check relation field in dataset metadata
    if (!is.null(dataset_meta$relation)) {
      relation_matches <- grep("SRP|ERP|DRP", dataset_meta$relation, value = TRUE)
      if (length(relation_matches) > 0) {
        additional_srp_ids <- gsub(".*=(SRP[0-9]+).*", "\\1", relation_matches)
        sra_ids <- c(sra_ids, additional_srp_ids)
      }
    }
    
    # Check supplementary fields that might contain SRA references
    if (!is.null(dataset_meta$supplementary_file)) {
      supp_matches <- grep("SRP|ERP|DRP", dataset_meta$supplementary_file, value = TRUE)
      if (length(supp_matches) > 0) {
        additional_srp_ids <- gsub(".*/(SRP[0-9]+).*", "\\1", supp_matches)
        sra_ids <- c(sra_ids, additional_srp_ids)
      }
    }
  }
  
  # Approach 3: Try to find accessions in any text field
  if (length(sra_ids) == 0) {
    cat("Trying to extract SRA IDs from any text field...\n")
    # Convert all metadata to text and search
    meta_text <- capture.output(print(geo_data))
    srp_matches <- gregexpr("SRP[0-9]+", paste(meta_text, collapse=" "))
    srp_ids <- regmatches(paste(meta_text, collapse=" "), srp_matches)[[1]]
    sra_ids <- c(sra_ids, srp_ids)
  }
  
  # Approach 4: If GSE123876 specifically is passed, hardcode the known SRA study ID
  if (geo_accession == "GSE123876" && length(sra_ids) == 0) {
    cat("Using hardcoded SRA ID for GSE123876...\n")
    # This is the SRA study ID for GSE123876 (can be verified on NCBI GEO page)
    sra_ids <- c("SRP187707")
  }
  
  # Remove duplicates and invalid entries
  sra_ids <- unique(sra_ids)
  sra_ids <- sra_ids[grepl("^SRP|^ERP|^DRP", sra_ids)]
  
  if (length(sra_ids) == 0) {
    cat("No SRA study IDs found for", geo_accession, "\n")
    return(NULL)
  }
  
  cat("Found SRA study IDs for", geo_accession, ":", paste(sra_ids, collapse=", "), "\n")
  return(sra_ids)
}

# Increase R's default file download timeout
options(timeout = max(1000, getOption("timeout")))

# Function to directly download SRA files from ENA using a more reliable method
direct_download_fastq <- function(run_id, output_dir) {
  # Create URLs for different possible sources
  urls <- c(
    # ENA direct FTP URL pattern for FASTQ files
    paste0("ftp://ftp.sra.ebi.ac.uk/vol1/fastq/", 
           substr(run_id, 1, 6), "/",
           ifelse(nchar(run_id) >= 10, substr(run_id, 10, 10), ""), "/",
           run_id, "/", run_id, "_1.fastq.gz"),
    
    paste0("ftp://ftp.sra.ebi.ac.uk/vol1/fastq/", 
           substr(run_id, 1, 6), "/",
           ifelse(nchar(run_id) >= 10, substr(run_id, 10, 10), ""), "/",
           run_id, "/", run_id, "_2.fastq.gz"),
    
    paste0("ftp://ftp.sra.ebi.ac.uk/vol1/fastq/", 
           substr(run_id, 1, 6), "/",
           ifelse(nchar(run_id) >= 10, substr(run_id, 10, 10), ""), "/",
           run_id, "/", run_id, ".fastq.gz"),
    
    # Alternative ENA pattern
    paste0("ftp://ftp.sra.ebi.ac.uk/vol1/fastq/", 
           substr(run_id, 1, 6), "/00",
           ifelse(nchar(run_id) >= 10, substr(run_id, 10, 10), ""), "/",
           run_id, "/", run_id, "_1.fastq.gz"),
    
    paste0("ftp://ftp.sra.ebi.ac.uk/vol1/fastq/", 
           substr(run_id, 1, 6), "/00",
           ifelse(nchar(run_id) >= 10, substr(run_id, 10, 10), ""), "/",
           run_id, "/", run_id, "_2.fastq.gz"),
    
    paste0("ftp://ftp.sra.ebi.ac.uk/vol1/fastq/", 
           substr(run_id, 1, 6), "/00",
           ifelse(nchar(run_id) >= 10, substr(run_id, 10, 10), ""), "/",
           run_id, "/", run_id, ".fastq.gz"),
    
    # NCBI SRA format
    paste0("https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=", run_id)
  )
  
  success <- FALSE
  
  # Try each URL
  for (url in urls) {
    filename <- basename(url)
    output_file <- file.path(output_dir, filename)
    
    cat("Attempting to download", filename, "from", url, "...\n")
    
    download_result <- tryCatch({
      download.file(url, output_file,
                    mode = "wb",
                    method = "curl",
                    quiet = FALSE)
      TRUE
    }, error = function(e) {
      cat("ERROR: Failed to download from", url, "-", e$message, "\n")
      FALSE
    }, warning = function(w) {
      cat("WARNING:", w$message, "\n")
      FALSE
    })
    
    if (download_result) {
      cat("SUCCESS: Downloaded", filename, "to", output_file, "\n")
      # Check if the file has actual content
      file_info <- file.info(output_file)
      if (file_info$size > 1000) {  # Reasonable size for a fastq file
        success <- TRUE
        cat("File size seems reasonable:", file_info$size, "bytes\n")
      } else {
        cat("WARNING: File seems too small (", file_info$size, "bytes). It may be empty or corrupted.\n")
        file.remove(output_file)
      }
    }
  }
  
  # Last resort: try to download directly from NCBI SRA using prefetch/fastq-dump
  if (!success) {
    cat("Attempting to use prefetch/fastq-dump for", run_id, "...\n")
    
    # Check if SRA toolkit is installed
    has_prefetch <- system2("which", "prefetch", stdout = TRUE, stderr = TRUE, error = TRUE)
    has_fastqdump <- system2("which", "fastq-dump", stdout = TRUE, stderr = TRUE, error = TRUE)
    
    if (!is.null(has_prefetch) && !is.null(has_fastqdump)) {
      # Create a temporary directory for prefetch
      temp_dir <- tempdir()
      
      # Run prefetch
      prefetch_cmd <- paste("prefetch", run_id, "-O", temp_dir)
      cat("Running:", prefetch_cmd, "\n")
      system(prefetch_cmd)
      
      # Run fastq-dump
      fastqdump_cmd <- paste("fastq-dump --outdir", output_dir, "--gzip --split-3", 
                            file.path(temp_dir, run_id, paste0(run_id, ".sra")))
      cat("Running:", fastqdump_cmd, "\n")
      system(fastqdump_cmd)
      
      # Check if files were created
      fastq_files <- list.files(output_dir, pattern = paste0(run_id, ".*fastq.gz"))
      if (length(fastq_files) > 0) {
        cat("SUCCESS: Generated", length(fastq_files), "FASTQ files using SRA toolkit\n")
        success <- TRUE
      } else {
        cat("WARNING: Failed to generate FASTQ files using SRA toolkit\n")
      }
    } else {
      cat("WARNING: SRA toolkit (prefetch/fastq-dump) not found in PATH\n")
    }
  }
  
  return(success)
}

# Function to download FASTQ files
download_fastq_files <- function(geo_accession, max_samples = Inf) {
  cat("\n======== Processing", geo_accession, "========\n")
  
  # Connect to SRA database
  tryCatch({
    # Check SRA database file
    sra_dbname <- check_sra_db(destdir = main_dir)
    cat("Using SRA database:", sra_dbname, "\n")
    con <- dbConnect(SQLite(), sra_dbname)
    
    # Get SRA study IDs from GEO
    sra_studies <- extract_sra_from_geo(geo_accession)
    if (is.null(sra_studies) || length(sra_studies) == 0) {
      cat("No SRA studies found for", geo_accession, "\n")
      dbDisconnect(con)
      return()
    }
    
    # Create output directory
    output_dir <- file.path(raw_data_dir, geo_accession)
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
      cat("Created directory:", output_dir, "\n")
    }
    
    # Process each SRA study
    for (study_id in sra_studies) {
      cat("\nProcessing SRA study:", study_id, "\n")
      
      # Check table schema first
      tables <- dbListTables(con)
      cat("Available tables in SRA database:", paste(tables, collapse=", "), "\n")
      
      # Get run info for this study
      query <- paste0("SELECT * FROM sra WHERE study_accession='", study_id, "' OR study_accession='", 
                     gsub("^SRP", "SRA", study_id), "'")
      run_info <- dbGetQuery(con, query)
      
      if (nrow(run_info) == 0) {
        cat("No runs found for study", study_id, "\n")
        # Try to get column names
        col_info <- dbGetQuery(con, "PRAGMA table_info(sra)")
        cat("SRA table columns:", paste(col_info$name, collapse=", "), "\n")
        next
      }
      
      cat("Found", nrow(run_info), "runs for study", study_id, "\n")
      cat("Run accessions:", paste(head(run_info$run_accession, 5), collapse=", "), 
          ifelse(nrow(run_info) > 5, "... and more", ""), "\n")
      
      cat("run_info contents:\n")
      print(run_info)  # For debugging

      if (!("fastq" %in% tables) && !("sra_ft" %in% tables)) {
        cat("No fastq or sra_ft tables in SRA database. Skipping database URL checks.\n")
        cat("Falling back to direct_download_fastq only.\n")
        # No database URL fallback
      } else {
        # No limit on samples - process all of them
        cat("Processing all available samples\n")
        
        # Check which column to use for run_accession
        run_acc_col <- ifelse("run_accession" %in% names(run_info), "run_accession", "run")
        
        # Process each run - this is where we download the FASTQ files
        successful_downloads <- 0
        for (i in 1:nrow(run_info)) {
          run_id <- run_info[[run_acc_col]][i]
          cat("\n===> Processing run:", run_id, " (", i, "of", nrow(run_info), ")\n")
          
          # Skip if we already have files for this run
          existing_files <- list.files(output_dir, pattern = paste0(run_id, ".*fastq"))
          if (length(existing_files) > 0) {
            cat("Files for run", run_id, "already exist. Skipping.\n")
            successful_downloads <- successful_downloads + 1
            next
          }
          
          # Try direct download first - this is more reliable than using database URLs
          download_success <- direct_download_fastq(run_id, output_dir)
          
          if (!download_success) {
            cat("Direct download failed. Trying database URLs...\n")
            
            # Try database method as fallback
            # Try different tables for FASTQ URLs
            fastq_info <- NULL
            
            # Check fastq table
            if ("fastq" %in% tables) {
              fastq_query <- paste0("SELECT * FROM fastq WHERE run_accession='", run_id, "'")
              fastq_info <- dbGetQuery(con, fastq_query)
            }
            
            # If no results, try sra_ft table
            if (is.null(fastq_info) || nrow(fastq_info) == 0) {
              if ("sra_ft" %in% tables) {
                alt_query <- paste0("SELECT * FROM sra_ft WHERE run_accession='", run_id, "'")
                fastq_info <- dbGetQuery(con, alt_query)
              }
            }
            
            # If still no results, try to construct the URL using NCBI/EBI pattern
            if (is.null(fastq_info) || nrow(fastq_info) == 0) {
              cat("No FASTQ info found in database for run", run_id, "\n")
              cat("Constructing URL based on NCBI/EBI pattern...\n")
              
              # Create a dataframe with the standard NCBI/EBI FTP URLs
              fastq_info <- data.frame(
                run_accession = run_id,
                url = c(
                  # NCBI format
                  paste0("ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/",
                         substr(run_id, 1, 3), "/",
                         ifelse(nchar(run_id) > 6, substr(run_id, 4, 6), substr(run_id, 4, nchar(run_id))), "/",
                         run_id, "/", run_id, ".sra"),
                  # EBI format
                  paste0("ftp://ftp.sra.ebi.ac.uk/vol1/fastq/",
                         substr(run_id, 1, 6), "/",
                         ifelse(substr(run_id, 10, 10) == "", "00", substr(run_id, 10, 10)), "/",
                         run_id, "/", run_id, ".fastq.gz")
                )
              )
            }
            
            # Download each FASTQ file
            for (j in 1:nrow(fastq_info)) {
              # Try different column names for URL
              url <- NULL
              for (col_name in c("url", "ftp", "fastq_ftp", "ftp_path", "download_path")) {
                if (col_name %in% names(fastq_info) && !is.na(fastq_info[[col_name]][j]) && fastq_info[[col_name]][j] != "") {
                  url <- fastq_info[[col_name]][j]
                  break
                }
              }
              
              if (is.null(url)) {
                cat("No valid URL found for run", run_id, "\n")
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
            
            # If we reach this point and still have no success, log it clearly
            cat("WARNING: Failed to download any files for run", run_id, "\n")
          } else {
            successful_downloads <- successful_downloads + 1
          }
        }
        
        cat("\nSummary for study", study_id, ":", successful_downloads, "out of", 
            nrow(run_info), "runs successfully downloaded\n")
      }
    }
    
    # Close database connection
    dbDisconnect(con)
    
  }, error = function(e) {
    cat("ERROR processing", geo_accession, ":", e$message, "\n")
    print(traceback())  # Print stack trace for better debugging
  })
}

# Process each GEO accession
cat("Starting to process", length(geo_accessions), "GEO accessions\n")
for (geo_accession in geo_accessions) {
  download_fastq_files(geo_accession, max_samples = Inf)  # Use Infinity to download all samples
}
cat("Processing complete\n")
