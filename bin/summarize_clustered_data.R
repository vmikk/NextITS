#!/usr/bin/env Rscript

## Script to pool remove low-quality singletons and summarize sequence abundance at OTU level (per sample)

# Input:
#   1. Sequence tables in long format with de novo chimeras removed (`Seqs.parquet`)
#   2. UC file                         (`UC_Pooled.parquet`)
#   3. FASTA file with OTU sequences   (`Clustered.fa.gz`)
#   4. Max MEEP score

# Outputs:
#  - OTU table in long format    (`OTU_table_long.txt.gz` & `OTU_table_long.RData`)
#  - OTU table in wide format    (`OTU_table_wide.txt.gz` & `OTU_table_wide.RData`)
#  - FASTA file with sequences   (`OTUs.fa.gz`)

## Usage:
# ./summarize_clustered_data.R \
#    --seqtab  "Seqs.parquet" \
#    --uc      "UC_Pooled.parquet" \
#    --otus    "Clustered.fa.gz" \
#    --maxmeep 0.6 \
#    --recoversinglet TRUE \
#    --mergesamples TRUE \
#    --wide-rdata-format sparse \
#    --threads 4


## Quality threshold:
# MEEP score of 0.6 corresponds approximately to the average Phred score of 22.2

## Singleton recovery:
# If enabled, then singleton OTUs with MEEP score <= 0.6 & will be preserved
# Otherwise, singleton OTUs will be removed


############################################## Parse input parameters

## Check time
start_time <- Sys.time()


cat("Parsing input options and arguments:\n")

suppressPackageStartupMessages(require(optparse))

## Parse arguments
option_list <- list(
  make_option("--seqtab",  action = "store", default = NA,  type = "character", help = "Sequence tables in long format with de novo chimeras removed (Parquet format)"),
  make_option("--uc",      action = "store", default = NA,  type = "character", help = "UC file (Parquet format)"),
  make_option("--otus",    action = "store", default = NA,  type = "character", help = "FASTA file with OTU sequences"),
  make_option("--maxmeep", action = "store", default = 0.5, type = "double", help = "Max MEEP score"),
  make_option("--recoversinglet", action = "store", default = TRUE, type = "logical", help = "Recover singletons"),
  make_option(c("-m", "--mergesamples"), action = "store", default = FALSE, type = "logical", help = "Merge sample replicates (default, false)"),
  make_option("--wide-rdata-format", action = "store", default = "sparse", type = "character", help = "Format for the wide OTU table: sparse or dense"),
  make_option(c("-t", "--threads"), action = "store", default = 4L, type = "integer", help = "Number of CPU threads for data processing, default 4")
)
opt <- parse_args(OptionParser(option_list = option_list))

## Function to convert text "NA"s to NA
to_na <- function(x){
  if(x %in% c("NA", "null", "Null")){ x <- NA }
  return(x)
}

## Replaces "null"s from Nextflow with NA
opt <- lapply(X = opt, FUN = to_na)


## Validation of the required argiments
if(is.na(opt$seqtab)){
  cat("Input file is not specified: sequence tables in Parquet format.\n", file = stderr())
  stop()
}
if(is.na(opt$uc)){
  cat("Input file is not specified: UC file is required.\n", file = stderr())
  stop()
}
if(is.na(opt$otus)){
  cat("Input file is not specified: FASTA file with OTU sequences.\n", file = stderr())
  stop()
}
if(opt$recoversinglet == TRUE && is.na(opt$maxmeep)){
  cat("For singleton recovery, the max MEEP score must be specified.\n", file = stderr())
  stop()
}

## Assign variables
SEQTAB            <- opt$seqtab
UCF               <- opt$uc
MAXMEEP           <- as.numeric(opt$maxmeep)
RECOV_SINGLET     <- as.logical(opt$recoversinglet)
MERGE_SAMPLES     <- as.logical(opt$mergesamples)
OTUS              <- opt$otus
WIDE_RDATA_FORMAT <- tolower(as.character(opt$`wide-rdata-format`))
CPUTHREADS        <- as.numeric(opt$threads)

if(is.na(WIDE_RDATA_FORMAT) || !WIDE_RDATA_FORMAT %in% c("sparse", "dense")){
  stop("Invalid --wide-rdata-format value. Use 'sparse' or 'dense'.", call. = FALSE)
}

## Log assigned variables
cat(paste("Sequence tables (Parquet format): ", SEQTAB,            "\n", sep = ""))
cat(paste("UC file (Parquet format): ",         UCF,               "\n", sep = ""))
cat(paste("Max MEEP score: ",                   MAXMEEP,           "\n", sep = ""))
cat(paste("Low-quality singleton recovery: ",   RECOV_SINGLET,     "\n", sep = ""))
cat(paste("Merge sample replicates: ",          MERGE_SAMPLES,     "\n", sep = ""))
cat(paste("OTU sequences: ",                    OTUS,              "\n", sep = ""))
cat(paste("Wide RData format: ",                WIDE_RDATA_FORMAT, "\n", sep = ""))
cat(paste("Number of CPU threads to use: ",     CPUTHREADS,        "\n", sep = ""))

cat("\n")


############################################## Data for debugging

# SEQTAB            <- "Seqs.parquet"
# UCF               <- "UC_Pooled.parquet"
# MAXMEEP           <- 0.5
# RECOV_SINGLET     <- TRUE
# MERGE_SAMPLES     <- TRUE
# OTUS              <- "Clustered.fa.gz"
# WIDE_RDATA_FORMAT <- "sparse"
# CPUTHREADS        <- 4


############################################## Load packages and helpers

cat("Loading R packages:\n")

load_pckg <- function(pkg = "data.table"){
  suppressPackageStartupMessages(library(package = pkg, character.only = TRUE))
  cat(paste(pkg, packageVersion(pkg), "\n"))
}

load_pckg("data.table")
load_pckg("DBI")
load_pckg("duckdb")
load_pckg("Matrix")
load_pckg("Biostrings")

cat("\n")

## Set CPU thread number
cat("Setting number of CPU threads to: ", CPUTHREADS, "\n")
setDTthreads(threads = CPUTHREADS)

## Helper function to compress RData using pigz (in parallel)
saveRDS.gz <- function(object, file, threads = parallel::detectCores()) {
  con <- pipe(paste0("pigz -p", threads, " > ", file), "wb")
  suppressWarnings(saveRDS(object, file = con))
  close(con)
}

sql_num <- function(x){
  format(as.numeric(x), scientific = FALSE, trim = TRUE)
}

## Helper function to write wide table to text file in chunks
## Chunks are OTU-based, not sample-based
write_wide_text_chunks <- function(res, otu_order, sample_order, file, chunk_target_cells = 5e7){
    
  if(file.exists(file)){
    file.remove(file)
  }

  if(length(otu_order) == 0L){
    fwrite(
      x = data.table(OTU = character()),
      file = file,
      sep = "\t",
      compress = "gzip")
    return(invisible(NULL))
  }

  chunk_otu_n <- max(1L, as.integer(floor(chunk_target_cells / max(1L, length(sample_order)))))
  n_chunks <- as.integer(ceiling(length(otu_order) / chunk_otu_n))

  cat("...Writing wide table by OTU chunks\n")
  cat(".... Target cells per chunk: ", format(chunk_target_cells, scientific = FALSE), "\n", sep = "")
  cat(".... OTUs per chunk: ",        chunk_otu_n, "\n", sep = "")
  cat(".... Number of chunks: ",      n_chunks,    "\n", sep = "")

  for(idx in seq_len(n_chunks)){
    start_idx <- (idx - 1L) * chunk_otu_n + 1L
    end_idx   <- min(length(otu_order), idx * chunk_otu_n)
    otu_chunk <- otu_order[start_idx:end_idx]

    cat(".... Processing chunk ", idx, "/", n_chunks, " [", start_idx, "-", end_idx, "]\n", sep = "")

    res_chunk <- res[data.table(OTU = otu_chunk), on = "OTU", nomatch = 0L]

    res_chunk[, OTU := factor(OTU, levels = otu_chunk)]
    res_chunk[, SampleID := factor(SampleID, levels = sample_order)]

    wide_chunk <- dcast(
      data = res_chunk,
      formula = OTU ~ SampleID,
      fun.aggregate = sum,
      fill = 0,
      drop = FALSE,
      value.var = "Abundance")

    wide_chunk[, OTU := as.character(OTU)]

    fwrite(
      x = wide_chunk,
      file = file,
      sep = "\t",
      append = idx > 1L,
      col.names = idx == 1L,
      compress = "gzip")

    rm(res_chunk, wide_chunk)
    gc(verbose = FALSE)
  }

  invisible(NULL)
}

build_sparse_wide <- function(res, otu_order, sample_order){
  if(nrow(res) == 0L){
    return(Matrix::sparseMatrix(
      i = integer(),
      j = integer(),
      x = numeric(),
      dims = c(length(otu_order), length(sample_order)),
      dimnames = list(otu_order, sample_order)))
  }

  Matrix::sparseMatrix(
    i = match(res$OTU, otu_order),
    j = match(res$SampleID, sample_order),
    x = as.numeric(res$Abundance),
    dims = c(length(otu_order), length(sample_order)),
    dimnames = list(otu_order, sample_order))
}

build_dense_wide <- function(res, otu_order, sample_order, max_cells = 5e7){
  n_cells <- as.double(length(otu_order)) * as.double(length(sample_order))

  if(n_cells > max_cells){
    stop(
      "Dense wide RData requested, but the wide table has ",
      format(n_cells, scientific = FALSE),
      " cells. Use --wide-rdata-format sparse instead.",
      call. = FALSE)
  }

  dense <- matrix(
    0,
    nrow = length(otu_order),
    ncol = length(sample_order),
    dimnames = list(otu_order, sample_order))

  if(nrow(res) > 0L){
    dense[cbind(match(res$OTU, otu_order), match(res$SampleID, sample_order))] <- as.numeric(res$Abundance)
  }

  dense
}


######################################
###################################### Load and process the data with DuckDB
######################################

cat("\n..Initializing DuckDB\n")
con <- DBI::dbConnect(duckdb::duckdb())
on.exit({
  if(exists("con") && !is.null(con)){
    try(dbDisconnect(con, shutdown = TRUE), silent = TRUE)
  }
}, add = TRUE)

invisible(dbExecute(con, sprintf("SET threads TO %d", as.integer(CPUTHREADS))))

seqtab_sql <- as.character(DBI::dbQuoteString(con, normalizePath(SEQTAB, mustWork = TRUE, winslash = "/")))
uc_sql     <- as.character(DBI::dbQuoteString(con, normalizePath(UCF,    mustWork = TRUE, winslash = "/")))

cat("\n..Registering Parquet inputs\n")
invisible(dbExecute(con, sprintf("
  CREATE OR REPLACE TEMP VIEW seqtab AS
  SELECT
    SeqID,
    SampleID,
    Abundance,
    MEEP
  FROM parquet_scan(%s)
", seqtab_sql)))

invisible(dbExecute(con, sprintf("
  CREATE OR REPLACE TEMP VIEW uc AS
  SELECT
    SeqID,
    OTU
  FROM parquet_scan(%s)
", uc_sql)))

cat("..Inspecting sequence tables\n")
seqtab_stats <- dbGetQuery(con, "
  SELECT
    COUNT(*) AS n_records,
    COUNT(DISTINCT SeqID) AS n_seq,
    COUNT(DISTINCT SampleID) AS n_samples
  FROM seqtab
")

cat("... Total number of records: ",             seqtab_stats$n_records[[1]], "\n")
cat("... Total number unique sequences: ",       seqtab_stats$n_seq[[1]],     "\n")
cat("... Total number unique samples (files): ", seqtab_stats$n_samples[[1]], "\n")

cat("..Joining sequence table with OTU membership\n")
invisible(dbExecute(con, "
  CREATE OR REPLACE TEMP VIEW joined AS
  SELECT
    s.SeqID,
    s.SampleID,
    s.Abundance,
    s.MEEP,
    u.OTU
  FROM seqtab s
  LEFT JOIN uc u USING (SeqID)
"))

missing_otu <- dbGetQuery(con, "
  SELECT
    COUNT(*) AS n_missing,
    COALESCE(SUM(Abundance), 0) AS abundance_missing
  FROM joined
  WHERE OTU IS NULL
")

if(missing_otu$n_missing[[1]] > 0){
  cat("WARNING: not all sequences were assigned to OTUs\n")
  cat("..Removing missing/excluded sequences\n")
  cat(".. ", missing_otu$n_missing[[1]], " sequences with total abundance ",
      missing_otu$abundance_missing[[1]], " reads will be excluded\n", sep = "")
}

invisible(dbExecute(con, "
  CREATE OR REPLACE TEMP VIEW joined_filtered AS
  SELECT *
  FROM joined
  WHERE OTU IS NOT NULL
"))

cat("\n..Finding singleton OTUs\n")
invisible(dbExecute(con, "
  CREATE OR REPLACE TEMP VIEW singleton_otus AS
  SELECT
    OTU,
    SUM(Abundance) AS Abundance
  FROM joined_filtered
  GROUP BY OTU
  HAVING SUM(Abundance) < 2
"))

singleton_count <- dbGetQuery(con, "SELECT COUNT(*) AS n_singletons FROM singleton_otus")
cat("... Number of singleton OTUs: ", singleton_count$n_singletons[[1]], "\n")

if(RECOV_SINGLET == TRUE){
  invisible(dbExecute(con, sprintf("
    CREATE OR REPLACE TEMP VIEW excluded_singletons AS
    SELECT DISTINCT
      s.OTU
    FROM singleton_otus s
    JOIN joined_filtered j
      ON j.SeqID = s.OTU
    WHERE j.MEEP > %s
  ", sql_num(MAXMEEP))))

  singleton_filtered <- dbGetQuery(con, "SELECT COUNT(*) AS n_filtered FROM excluded_singletons")
  cat("... Number of singleton OTUs after filtering by MEEP score: ", singleton_filtered$n_filtered[[1]], "\n")
} else {
  invisible(dbExecute(con, "
    CREATE OR REPLACE TEMP VIEW excluded_singletons AS
    SELECT OTU
    FROM singleton_otus
  "))
}

excluded_singletons <- dbGetQuery(con, "SELECT COUNT(*) AS n_excluded FROM excluded_singletons")

if(excluded_singletons$n_excluded[[1]] > 0){
  cat("..Removing singleton OTUs\n")
}

invisible(dbExecute(con, "
  CREATE OR REPLACE TEMP VIEW joined_pruned AS
  SELECT
    j.*
  FROM joined_filtered j
  LEFT JOIN excluded_singletons e USING (OTU)
  WHERE e.OTU IS NULL
"))

records_after_singletons <- dbGetQuery(con, "SELECT COUNT(*) AS n_records FROM joined_pruned")
cat("... Number of records in sequence table after removing singleton OTUs: ", records_after_singletons$n_records[[1]], "\n")

cat("\n..Summarizing OTU abundance\n")
sample_expr_sql <- if(MERGE_SAMPLES == TRUE){
  cat("\n... Merging sample replicates (e.g., re-sequenced samples)\n")
  "CASE WHEN strpos(SampleID, '__') > 0 THEN split_part(SampleID, '__', 2) ELSE NULL END"
} else {
  "SampleID"
}

if(MERGE_SAMPLES == TRUE){
  cat(".... Summarizing abundance by sample and OTU\n")
} else {
  cat("... Summarizing abundance by sample and OTU\n")
}

RES <- dbGetQuery(con, sprintf("
  SELECT
    OTU,
    %s AS SampleID,
    SUM(Abundance) AS Abundance
  FROM joined_pruned
  GROUP BY 1, 2
  ORDER BY 1, 2
", sample_expr_sql))
setDT(RES)

if(nrow(RES) == 0L){
  RES <- data.table(OTU = character(), SampleID = character(), Abundance = numeric())
} else {
  RES[, Abundance := as.numeric(Abundance)]
}

cat("... Number of rows in summarized long table: ", nrow(RES), "\n")

cat("\n..Disconnecting DuckDB\n")
dbDisconnect(con, shutdown = TRUE)
con <- NULL


######################################
###################################### Prepare long and wide outputs
######################################

setorder(RES, OTU, SampleID)
setattr(x = RES, name = "Samples_merged", value = MERGE_SAMPLES)

cat("\nPreparing wide output metadata\n")
n_otu <- uniqueN(RES$OTU)
n_smp <- uniqueN(RES$SampleID)
n_cll <- as.double(n_otu) * as.double(n_smp)

cat("...In total, there are ", n_otu, " OTUs and ", n_smp, " samples\n")
cat("...The total number of cells in the wide table will be ", format(n_cll, scientific = FALSE), "\n")

if(nrow(RES) > 0L){
  otu_order <- RES[, .(TotalAbundance = sum(Abundance, na.rm = TRUE)), by = "OTU"][
    order(-TotalAbundance, OTU), OTU]
  sample_order <- sort(unique(RES$SampleID), na.last = TRUE)
  setindex(RES, OTU)
} else {
  otu_order    <- character()
  sample_order <- character()
}


cat("\nExporting results\n")

cat("..Exporting long table [R]\n")
saveRDS.gz(
  object = RES,
  file = "OTU_table_long.RData",
  threads = CPUTHREADS)

cat("..Exporting long table [tab-delimited]\n")
fwrite(x = RES, file = "OTU_table_long.txt.gz", sep = "\t", compress = "gzip")

cat("..Exporting wide table [tab-delimited]\n")
write_wide_text_chunks(
  res = RES,
  otu_order = otu_order,
  sample_order = sample_order,
  file = "OTU_table_wide.txt.gz")

cat("..Preparing wide table [R]\n")
if(WIDE_RDATA_FORMAT == "dense"){
  WIDE_OBJ <- build_dense_wide(
    res = RES,
    otu_order = otu_order,
    sample_order = sample_order)
} else {
  WIDE_OBJ <- build_sparse_wide(
    res = RES,
    otu_order = otu_order,
    sample_order = sample_order)
}
attr(WIDE_OBJ, "Samples_merged") <- MERGE_SAMPLES

cat("..Exporting wide table [R]\n")
saveRDS.gz(
  object = WIDE_OBJ,
  file = "OTU_table_wide.RData",
  threads = CPUTHREADS)


######################################
###################################### Export OTU sequences
######################################

cat("\nExporting OTU sequences to FASTA\n")
cat("..Preparing sequences\n")

cat("... Loading FASTA file\n")
SQS <- readDNAStringSet(filepath = OTUS, format = "fasta")
cat("... Extracting sequence IDs\n")
names(SQS) <- tstrsplit(x = names(SQS), split = ";", keep = 1)[[1]]

if(any(duplicated(names(SQS)))){
  cat("WARNING: duplicated OTU names detected!\n")
}

cat("... Subsetting OTUs\n")
otu_export_ids <- otu_order[otu_order %in% names(SQS)]
SQF <- SQS[otu_export_ids]

cat("....Total number of OTUs in input sequences: ", length(SQS),       "\n")
cat("....Number of OTUs to export: ",                length(SQF),       "\n")
cat("....Number of OTUs in the OTU table: ",         length(otu_order), "\n")

if(length(SQF) != length(otu_order)){
  cat("WARNING: not all OTUs from the OTU table were found in the input FASTA\n")
}

cat("... Writing FASTA file\n")
writeXStringSet(
  x = SQF,
  filepath = "OTUs.fa.gz",
  compress = TRUE,
  format = "fasta",
  width = 9999)


cat("\nAll done.\n")


##################### Session info

## Check time
end_time <- Sys.time()

tmm <- as.numeric(difftime(end_time, start_time, units = "min"))
cat("\nElapsed time: ", tmm, " minutes\n")

cat("\n")
cat("Session info:\n")
sessionInfo()
cat("\n")
