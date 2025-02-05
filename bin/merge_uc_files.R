#!/usr/bin/env Rscript

## Merge UC files from different steps (dereplication, pre-clustering, clustering) into a single file

cat("Joining parquet files\n\n")

## Check time
start_time <- Sys.time()

## Function to load packages
load_pckg <- function(pkg = "data.table"){
    suppressPackageStartupMessages( library(package = pkg, character.only = TRUE) )
    cat(".. ", paste(pkg, packageVersion(pkg), "\n"))
}

cat("Loading packages:\n")
load_pckg("DBI")
load_pckg("duckdb")
load_pckg("qs")
# load_pckg("data.table")
# load_pckg("arrow")
# load_pckg("dplyr")

cat("\nParsing input options and arguments...\n")

suppressPackageStartupMessages(require(optparse))

## Parse arguments
option_list <- list(
  make_option("--ucderep",    action="store", default=NA, type='character', help="UC file from global dereplication"),
  make_option("--ucpreclust", action="store", default=NA, type='character', help="UC file from pre-clustering (optional)"),
  make_option("--ucclust",    action="store", default=NA, type='character', help="UC file from clustering"),
  make_option("--output",     action="store", default="UC_Pooled.parquet",  type='character', help="Output file name"),
  make_option(c("-t", "--threads"), action="store", default=4L, type='integer', help="Number of CPU threads for arrow, default 4")
)
opt <- parse_args(OptionParser(option_list=option_list))

## Function to convert text "NA"s to NA
to_na <- function(x){
  if(x %in% c("NA", "null", "Null")){ x <- NA }
  return(x)
}

## Replaces "null"s from Nextflow with NA
opt <- lapply(X = opt, FUN = to_na)


## Validation of the required argiments
if(is.na(opt$ucderep)){
  cat("Input file is not specified: UC file from global dereplication.\n", file=stderr())
  stop()
}
if(is.na(opt$ucclust)){
  cat("Input file is not specified: UC file from clustering.\n", file=stderr())
  stop()
}


## Assign variables
UCDEREP    <- opt$ucderep
UCPRECLUST <- opt$ucpreclust
UCCLUST    <- opt$ucclust
OUTPUT     <- opt$output
CPUTHREADS <- as.numeric( opt$threads )

## Log assigned variables
cat("\nParameters specified:\n")
cat(paste("  UC file from global dereplication: ",        UCDEREP,    "\n", sep=""))
cat(paste("  UC file from pre-clustering or denoising: ", UCPRECLUST, "\n", sep=""))
cat(paste("  UC file from clustering: ",                  UCCLUST,    "\n", sep=""))
cat(paste("  Output file name: ",                         OUTPUT,     "\n", sep=""))
cat(paste("  Number of CPU threads to use: ",             CPUTHREADS, "\n", sep=""))

cat("\n")

## Data for debugging
# UCDEREP    <- "UC_derep.parquet"
# UCPRECLUST <- "UC_preclust.parquet"
# UCCLUST    <- "UC_clust.parquet"
# OUTPUT     <- "UC_Pooled.parquet"
# CPUTHREADS <- 4


######################################
###################################### Load and process the data [duckdb]
######################################

## Initialize DuckDB connection
cat("..Initializing DuckDB connection\n")
con <- dbConnect(duckdb::duckdb())

## Register parquet files as tables
cat("..Registering parquet files as tables\n")
cat("...Dereplication UC\n")
dbExecute(con, sprintf("CREATE VIEW derep_seqs AS SELECT * FROM parquet_scan('%s')", UCDEREP))

if(!is.na(UCPRECLUST)) {
  cat("...Pre-clustering UC\n")
  dbExecute(con, sprintf("CREATE VIEW preclust_seqs AS SELECT * FROM parquet_scan('%s')", UCPRECLUST))
}

cat("...Clustering UC\n")
dbExecute(con, sprintf("CREATE VIEW clust_seqs AS SELECT * FROM parquet_scan('%s')", UCCLUST))

## Set number of threads
cat("..Setting number of threads for DuckDB\n")
dbExecute(con, sprintf("SET threads TO %d;", CPUTHREADS))


## Process and merge the data
if(is.na(UCPRECLUST)) {

  ## Two-file merge (no pre-clustering)
  cat("..Merging UC files [no pre-clustering or denoising]\n")
  dbExecute(con, sprintf("
    COPY (
      WITH derep AS (
        SELECT DISTINCT
          query as  SeqID,
          target as DerepID
        FROM derep_seqs
        QUALIFY ROW_NUMBER() OVER (PARTITION BY query ORDER BY target) = 1
      ),
      clust AS (
        SELECT DISTINCT
          query as  DerepID,
          target as OTU
        FROM clust_seqs
        QUALIFY ROW_NUMBER() OVER (PARTITION BY query ORDER BY target) = 1
      )
      SELECT 
        d.SeqID,
        d.DerepID,
        c.OTU
      FROM derep d
      LEFT JOIN clust c ON d.DerepID = c.DerepID
    ) TO '%s'
    (FORMAT 'parquet', COMPRESSION 'ZSTD', COMPRESSION_LEVEL 8)",
    OUTPUT))

} else {
  
  ## Three-file merge (with pre-clustering)
  cat("..Merging UC files [with pre-clustering or denoising]\n")

  dbExecute(con, sprintf("
    COPY (
      WITH derep AS (
        SELECT DISTINCT
          query as  SeqID,
          target as DerepID
        FROM derep_seqs
        QUALIFY ROW_NUMBER() OVER (PARTITION BY query ORDER BY target) = 1
      ),
      preclust AS (
        SELECT DISTINCT
          query as  DerepID,
          target as PreclusterID
        FROM preclust_seqs
        QUALIFY ROW_NUMBER() OVER (PARTITION BY query ORDER BY target) = 1
      ),
      clust AS (
        SELECT DISTINCT
          query as  PreclusterID,
          target as OTU
        FROM clust_seqs
        QUALIFY ROW_NUMBER() OVER (PARTITION BY query ORDER BY target) = 1
      )
      SELECT 
        d.SeqID,
        d.DerepID,
        p.PreclusterID,
        c.OTU
      FROM derep d
      LEFT JOIN preclust p ON d.DerepID = p.DerepID
      LEFT JOIN clust c ON p.PreclusterID = c.PreclusterID
    ) TO '%s'
    (FORMAT 'parquet', COMPRESSION 'ZSTD', COMPRESSION_LEVEL 8)",
    OUTPUT))

}

## Clean up
cat("..Disconnecting from DuckDB\n")
dbDisconnect(con, shutdown = TRUE)


cat("..Done!\n")



######################################
###################################### Load and process the data  [arrow + dplyr + data.table]
######################################

# ## Set number of threads for data.table
# cat("..Setting number of threads\n")
# setDTthreads(threads = CPUTHREADS)   # for data.table
# set_cpu_count(CPUTHREADS)            # for arrow
# 
# ## Globally dereplicated sequences (remove multi-target matches)
# cat("..Loading globally dereplicated sequences\n")
# UCA <- open_dataset(UCDEREP) %>% 
#   rename(SeqID = query, DerepID = target) %>% 
#   to_duckdb() %>%
#   distinct(SeqID, .keep_all = TRUE) %>%
#   collect() %>%
#   setDT()
# 
# ## Pre-clustered sequences (remove multi-target matches)
# if(!is.na(UCPRECLUST)){
#   cat("..Loading pre-clustered sequences\n")
#   UCP <- open_dataset(UCPRECLUST) %>% 
#     rename(DerepID = query, PreclusterID = target) %>% 
#     to_duckdb() %>%
#     distinct(DerepID, .keep_all = TRUE) %>%
#     collect() %>%
#     setDT()
# }
# 
# ## Clustered sequences (remove multi-target matches)
# cat("..Loading clustering UC file\n")
# if(!is.na(UCPRECLUST)){
#   UCO <- open_dataset(UCCLUST) %>% 
#     rename(PreclusterID = query, OTU = target) %>% 
#     to_duckdb() %>%
#     distinct(PreclusterID, .keep_all = TRUE) %>%
#     collect() %>%
#     setDT()
# } else {
#   UCO <- open_dataset(UCCLUST) %>% 
#     rename(DerepID = query, OTU = target) %>% 
#     to_duckdb() %>%
#     distinct(DerepID, .keep_all = TRUE) %>%
#     collect() %>%
#     setDT()
# }
# 
# ## Merge UC files
# if(is.na(UCPRECLUST)){
#   
#   ## No pre-clustering or denoising
#   cat("..Merging UC files [no pre-clustering or denoising ]\n")
#   RES <- merge(x = UCA, y = UCO, by = "DerepID", all.x = TRUE)
# 
# } else {
# 
#   ## Merge UC files with pre-clustering or denoising
#   cat("..Merging UC files [with pre-clustering or denoising ]\n")
#   
#   cat("... Adding pre-cluster or denoised IDs\n")
#   RES <- merge(x = UCA, y = UCP, by = "DerepID", all.x = TRUE)
# 
#   cat("... Adding clustering IDs\n")
#   RES <- merge(x = RES, y = UCO, by = "PreclusterID", all.x = TRUE)
# 
# }
# 

