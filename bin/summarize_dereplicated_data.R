#!/usr/bin/env Rscript

## Script to pool dereplicated sequences (non-clustered or denoised),
## remove low-quality singletons, and summarize sequence abundance per sample

# Input:
#   1. Sequence tables in long format with de novo chimeras removed (`Seqs.parquet`)
#   2. UC file from dereplication (`UC_Pooled.parquet`)
#   3. Sequence table in Parquet format (`Dereplicated.parquet`)
#   4. Max MEEP score

# Outputs:
#  - OTU table in long format (`OTU_table_long.txt.gz` & `OTU_table_long.RData`)
#  - OTU table in wide format (`OTU_table_wide.txt.gz` & `OTU_table_wide.RData`)
#  - FASTA file with sequences (`OTUs.fa.gz`)

## Usage:
# ./summarize_dereplicated_data.R \
#    --seqtab  "Seqs.parquet" \
#    --uc      "UC_Pooled.parquet" \
#    --seqs    "Dereplicated.parquet" \
#    --maxmeep 0.6 \
#    --recoversinglet TRUE \
#    --mergesamples TRUE \
#    --wide-rdata-format sparse \
#    --threads 4


## Quality threshold:
# MEEP score of 0.6 corresponds approximately to the average Phred score of 22.2

## Singleton recovery:
# If enabled, then singleton dereplicated sequences with MEEP score <= 0.6
# will be preserved; otherwise, singleton dereplicated sequences will be removed

## TODO:
# - add option to keep singletons


############################################## Parse input parameters

## Check time
start_time <- Sys.time()

cat("Parsing input options and arguments:\n")

suppressPackageStartupMessages(require(optparse))

## Parse arguments
option_list <- list(
  make_option("--seqtab", action = "store", default = NA, type = "character", help = "Sequence tables in long format with de novo chimeras removed (Parquet format)"),
  make_option("--uc", action = "store", default = NA, type = "character", help = "UC file (Parquet format)"),
  make_option("--seqs", action = "store", default = NA, type = "character", help = "Dereplicated sequences in Parquet format"),
  make_option("--maxmeep", action = "store", default = 0.6, type = "double", help = "Max MEEP score"),
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


## Validation of the required arguments
if(is.na(opt$seqtab)){
  cat("Input file is not specified: sequence tables in Parquet format.\n", file = stderr())
  stop()
}
if(is.na(opt$uc)){
  cat("Input file is not specified: UC file is required.\n", file = stderr())
  stop()
}
if(is.na(opt$seqs)){
  cat("Input file is not specified: dereplicated sequence table in Parquet format.\n", file = stderr())
  stop()
}
if(opt$recoversinglet == TRUE && is.na(opt$maxmeep)){
  cat("For singleton recovery, the max MEEP score must be specified.\n", file = stderr())
  stop()
}

## Assign variables
SEQTAB            <- opt$seqtab
UCF               <- opt$uc
SEQS_PQ           <- opt$seqs
MAXMEEP           <- as.numeric(opt$maxmeep)
RECOV_SINGLET     <- as.logical(opt$recoversinglet)
MERGE_SAMPLES     <- as.logical(opt$mergesamples)
WIDE_RDATA_FORMAT <- tolower(as.character(opt$`wide-rdata-format`))
CPUTHREADS        <- as.numeric(opt$threads)

if(is.na(WIDE_RDATA_FORMAT) || !WIDE_RDATA_FORMAT %in% c("sparse", "dense")){
  stop("Invalid --wide-rdata-format value. Use 'sparse' or 'dense'.", call. = FALSE)
}

## Log assigned variables
cat(paste("Sequence tables (Parquet format): ", SEQTAB,            "\n", sep = ""))
cat(paste("UC file (Parquet format): ",         UCF,               "\n", sep = ""))
cat(paste("Dereplicated sequences (Parquet format): ", SEQS_PQ,    "\n", sep = ""))
cat(paste("Max MEEP score: ",                   MAXMEEP,           "\n", sep = ""))
cat(paste("Low-quality singleton recovery: ",   RECOV_SINGLET,     "\n", sep = ""))
cat(paste("Merge sample replicates: ",          MERGE_SAMPLES,     "\n", sep = ""))
cat(paste("Wide RData format: ",                WIDE_RDATA_FORMAT, "\n", sep = ""))
cat(paste("Number of CPU threads to use: ",     CPUTHREADS,        "\n", sep = ""))
cat("\n")


############################################## Data for debugging

# SEQTAB            <- "Seqs.parquet"
# UCF               <- "UC_Pooled.parquet"
# SEQS_PQ           <- "Dereplicated.parquet"
# MAXMEEP           <- 0.6
# RECOV_SINGLET     <- TRUE
# MERGE_SAMPLES     <- TRUE
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
cat("\n")

## Set CPU thread number
cat("Setting number of CPU threads to: ", CPUTHREADS, "\n")
setDTthreads(threads = CPUTHREADS)

## Helper function to compress RData using pigz (in parallel)
saveRDS.gz <- function(object, file, threads = parallel::detectCores()){
  con <- pipe(paste0("pigz -p", threads, " > ", file), "wb")
  suppressWarnings(saveRDS(object, file = con))
  close(con)
}

sql_num <- function(x){
  format(as.numeric(x), scientific = FALSE, trim = TRUE)
}

## Helper function to write wide table to text file in chunks
## Chunks are dereplicated-sequence based, not sample-based
write_wide_text_chunks <- function(res, derep_order, sample_order, file, chunk_target_cells = 5e7){

  if(file.exists(file)){
    file.remove(file)
  }

  if(length(derep_order) == 0L){
    fwrite(
      x = data.table(DerepID = character()),
      file = file,
      sep = "\t",
      compress = "gzip")
    return(invisible(NULL))
  }

  chunk_derep_n <- max(1L, as.integer(floor(chunk_target_cells / max(1L, length(sample_order)))))
  n_chunks <- as.integer(ceiling(length(derep_order) / chunk_derep_n))

  cat("...Writing wide table by dereplicated-sequence chunks\n")
  cat(".... Target cells per chunk: ", format(chunk_target_cells, scientific = FALSE), "\n", sep = "")
  cat(".... Dereplicated sequences per chunk: ", chunk_derep_n, "\n", sep = "")
  cat(".... Number of chunks: ", n_chunks, "\n", sep = "")

  for(idx in seq_len(n_chunks)){
    start_idx <- (idx - 1L) * chunk_derep_n + 1L
    end_idx   <- min(length(derep_order), idx * chunk_derep_n)
    derep_chunk <- derep_order[start_idx:end_idx]

    cat(".... Processing chunk ", idx, "/", n_chunks, " [", start_idx, "-", end_idx, "]\n", sep = "")

    res_chunk <- res[data.table(DerepID = derep_chunk), on = "DerepID", nomatch = 0L]

    res_chunk[, DerepID := factor(DerepID, levels = derep_chunk)]
    res_chunk[, SampleID := factor(SampleID, levels = sample_order)]

    wide_chunk <- dcast(
      data = res_chunk,
      formula = DerepID ~ SampleID,
      fun.aggregate = sum,
      fill = 0,
      drop = FALSE,
      value.var = "Abundance")

    wide_chunk[, DerepID := as.character(DerepID)]

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

build_sparse_wide <- function(res, derep_order, sample_order){
  if(nrow(res) == 0L){
    return(Matrix::sparseMatrix(
      i = integer(),
      j = integer(),
      x = numeric(),
      dims = c(length(derep_order), length(sample_order)),
      dimnames = list(derep_order, sample_order)))
  }

  Matrix::sparseMatrix(
    i = match(res$DerepID, derep_order),
    j = match(res$SampleID, sample_order),
    x = as.numeric(res$Abundance),
    dims = c(length(derep_order), length(sample_order)),
    dimnames = list(derep_order, sample_order))
}

build_dense_wide <- function(res, derep_order, sample_order, max_cells = 5e7){
  n_cells <- as.double(length(derep_order)) * as.double(length(sample_order))

  if(n_cells > max_cells){
    stop(
      "Dense wide RData requested, but the wide table has ",
      format(n_cells, scientific = FALSE),
      " cells. Use --wide-rdata-format sparse instead.",
      call. = FALSE)
  }

  dense <- matrix(
    0,
    nrow = length(derep_order),
    ncol = length(sample_order),
    dimnames = list(derep_order, sample_order))

  if(nrow(res) > 0L){
    dense[cbind(match(res$DerepID, derep_order), match(res$SampleID, sample_order))] <- as.numeric(res$Abundance)
  }

  dense
}

export_derep_fasta <- function(seqs_parquet, derep_order, file, threads = 1L, chunk_size = 10000L){
  con <- DBI::dbConnect(duckdb::duckdb())
  result <- NULL
  out_con <- NULL

  on.exit({
    if(!is.null(result)){
      try(DBI::dbClearResult(result), silent = TRUE)
    }
    if(!is.null(out_con)){
      try(close(out_con), silent = TRUE)
    }
    if(!is.null(con)){
      try(DBI::dbDisconnect(con, shutdown = TRUE), silent = TRUE)
    }
  }, add = TRUE)

  invisible(DBI::dbExecute(con, sprintf("SET threads TO %d", as.integer(threads))))

  seqs_sql <- as.character(DBI::dbQuoteString(
    con,
    seqs_parquet))

  seq_export <- data.frame(
    SeqID = derep_order,
    Seq_rank = seq_along(derep_order),
    stringsAsFactors = FALSE)
  DBI::dbWriteTable(con, "seq_export", seq_export, temporary = TRUE, overwrite = TRUE)

  input_stats <- DBI::dbGetQuery(con, sprintf("
    SELECT
      COUNT(*) AS n_input_rows,
      COUNT(DISTINCT SeqID) AS n_input_unique
    FROM parquet_scan(%s)
  ", seqs_sql))

  dup_stats <- DBI::dbGetQuery(con, sprintf("
    SELECT COUNT(*) AS n_dup_ids
    FROM (
      SELECT SeqID
      FROM parquet_scan(%s)
      GROUP BY SeqID
      HAVING COUNT(*) > 1
    )
  ", seqs_sql))

  has_duplicate_ids <- dup_stats$n_dup_ids[[1]] > 0

  if(has_duplicate_ids){
    cat("... Duplicated SeqIDs detected in the input parquet - using deduplication export path\n")

    invisible(DBI::dbExecute(con, sprintf("
      CREATE OR REPLACE TEMP VIEW seqs_indexed AS
      SELECT
        ROW_NUMBER() OVER () AS input_row_num,
        SeqID,
        Sequence
      FROM parquet_scan(%s)
    ", seqs_sql)))

    invisible(DBI::dbExecute(con, "
      CREATE OR REPLACE TEMP VIEW seqs_dedup AS
      SELECT
        SeqID,
        Sequence
      FROM seqs_indexed
      QUALIFY ROW_NUMBER() OVER (PARTITION BY SeqID ORDER BY input_row_num) = 1
    "))

    export_stats <- DBI::dbGetQuery(con, "
      SELECT
        (SELECT COUNT(*) FROM seq_export) AS n_requested,
        (SELECT COUNT(*) FROM seq_export e INNER JOIN seqs_dedup s USING (SeqID)) AS n_exported,
        (SELECT COUNT(*) FROM seq_export e LEFT JOIN seqs_dedup s USING (SeqID) WHERE s.SeqID IS NULL) AS n_missing
    ")

    export_query <- "
      SELECT
        '>' || e.SeqID || chr(10) || s.Sequence AS fasta_record
      FROM seq_export e
      INNER JOIN seqs_dedup s USING (SeqID)
      ORDER BY e.Seq_rank
    "
  } else {
    cat("... Input parquet has unique SeqIDs - using direct parquet export path\n")

    export_stats <- DBI::dbGetQuery(con, sprintf("
      SELECT
        (SELECT COUNT(*) FROM seq_export) AS n_requested,
        (SELECT COUNT(*) FROM seq_export e INNER JOIN parquet_scan(%s) s USING (SeqID)) AS n_exported,
        (SELECT COUNT(*) FROM seq_export e LEFT JOIN parquet_scan(%s) s USING (SeqID) WHERE s.SeqID IS NULL) AS n_missing
    ", seqs_sql, seqs_sql))

    export_query <- sprintf("
      SELECT
        '>' || e.SeqID || chr(10) || s.Sequence AS fasta_record
      FROM seq_export e
      INNER JOIN parquet_scan(%s) s USING (SeqID)
      ORDER BY e.Seq_rank
    ", seqs_sql)
  }

  cat("... Total number of sequence rows in input parquet: ", input_stats$n_input_rows[[1]], "\n")
  cat("... Total number unique SeqIDs in input parquet: ", input_stats$n_input_unique[[1]], "\n")
  cat("... Number of dereplicated SeqIDs requested from the abundance table: ", export_stats$n_requested[[1]], "\n")
  cat("... Number of dereplicated SeqIDs to export: ", export_stats$n_exported[[1]], "\n")
  cat("... Number of dereplicated SeqIDs missing from the input parquet: ", export_stats$n_missing[[1]], "\n")

  if(dup_stats$n_dup_ids[[1]] > 0){
    cat("WARNING: duplicated SeqIDs detected in the input parquet - keeping the first row for ",
        dup_stats$n_dup_ids[[1]], " SeqIDs\n", sep = "")
  }

  if(export_stats$n_missing[[1]] > 0){
    cat("WARNING: not all dereplicated SeqIDs from the abundance table were found in the input parquet\n")
  }

  out_con <- pipe(
    sprintf("pigz -p%d > %s",
      as.integer(threads),
      shQuote(file)),
    "wt")

  result <- DBI::dbSendQuery(con, export_query)

  repeat {
    chunk <- DBI::dbFetch(result, n = as.integer(chunk_size))
    if(nrow(chunk) == 0L){
      break
    }
    writeLines(chunk[[1L]], out_con)
  }
}


######################################
###################################### Load and process the data with DuckDB
######################################

cat("\n..Initializing DuckDB\n")
con <- DBI::dbConnect(duckdb::duckdb())
on.exit({
  if(exists("con") && !is.null(con)){
    try(DBI::dbDisconnect(con, shutdown = TRUE), silent = TRUE)
  }
}, add = TRUE)

invisible(DBI::dbExecute(con, sprintf("SET threads TO %d", as.integer(CPUTHREADS))))

seqtab_sql <- as.character(DBI::dbQuoteString(con, SEQTAB))
uc_sql     <- as.character(DBI::dbQuoteString(con, UCF))

cat("\n..Registering Parquet inputs\n")
invisible(DBI::dbExecute(con, sprintf("
  CREATE OR REPLACE TEMP VIEW seqtab AS
  SELECT
    SeqID,
    SampleID,
    Abundance,
    MEEP
  FROM parquet_scan(%s)
", seqtab_sql)))

uc_field_names <- names(DBI::dbGetQuery(con, sprintf("SELECT * FROM parquet_scan(%s) LIMIT 0", uc_sql)))

if(all(c("SeqID", "DerepID") %in% uc_field_names)){
  seqid_col <- "SeqID"
  derep_col <- "DerepID"
} else if(length(uc_field_names) >= 2L){
  seqid_col <- uc_field_names[[1]]
  derep_col <- uc_field_names[[2]]
  cat("..Using first two columns from the dereplication UC parquet as SeqID and DerepID: ",
      seqid_col, ", ", derep_col, "\n", sep = "")
} else {
  stop("UC parquet must contain at least two columns for SeqID and DerepID.", call. = FALSE)
}

seqid_col_sql <- as.character(DBI::dbQuoteIdentifier(con, seqid_col))
derep_col_sql <- as.character(DBI::dbQuoteIdentifier(con, derep_col))

invisible(DBI::dbExecute(con, sprintf("
  CREATE OR REPLACE TEMP VIEW uc AS
  SELECT
    %s AS SeqID,
    %s AS DerepID
  FROM parquet_scan(%s)
", seqid_col_sql, derep_col_sql, uc_sql)))

cat("..Inspecting sequence tables\n")
seqtab_stats <- DBI::dbGetQuery(con, "
  SELECT
    COUNT(*) AS n_records,
    COUNT(DISTINCT SeqID) AS n_seq,
    COUNT(DISTINCT SampleID) AS n_samples
  FROM seqtab
")

cat("... Total number of records: ", seqtab_stats$n_records[[1]], "\n")
cat("... Total number unique sequences: ", seqtab_stats$n_seq[[1]], "\n")
cat("... Total number unique samples (files): ", seqtab_stats$n_samples[[1]], "\n")

cat("..Joining sequence table with dereplicated-sequence membership\n")
invisible(DBI::dbExecute(con, "
  CREATE OR REPLACE TEMP VIEW joined AS
  SELECT
    s.SeqID,
    s.SampleID,
    s.Abundance,
    s.MEEP,
    u.DerepID
  FROM seqtab s
  LEFT JOIN uc u USING (SeqID)
"))

missing_derep <- DBI::dbGetQuery(con, "
  SELECT
    COUNT(*) AS n_missing,
    COALESCE(SUM(Abundance), 0) AS abundance_missing
  FROM joined
  WHERE DerepID IS NULL
")

if(missing_derep$n_missing[[1]] > 0){
  cat("WARNING: not all sequences were assigned to dereplicated-sequence IDs\n")
  cat("..Removing missing/excluded sequences\n")
  cat(".. ", missing_derep$n_missing[[1]], " sequences with total abundance ",
      missing_derep$abundance_missing[[1]], " reads will be excluded\n", sep = "")
}

invisible(DBI::dbExecute(con, "
  CREATE OR REPLACE TEMP VIEW joined_filtered AS
  SELECT *
  FROM joined
  WHERE DerepID IS NOT NULL
"))

cat("\n..Finding singleton dereplicated sequences\n")
invisible(DBI::dbExecute(con, "
  CREATE OR REPLACE TEMP VIEW singleton_derep AS
  SELECT
    DerepID,
    SUM(Abundance) AS Abundance
  FROM joined_filtered
  GROUP BY DerepID
  HAVING SUM(Abundance) < 2
"))

singleton_count <- DBI::dbGetQuery(con, "SELECT COUNT(*) AS n_singletons FROM singleton_derep")
cat("... Number of singleton dereplicated sequences: ", singleton_count$n_singletons[[1]], "\n")

if(RECOV_SINGLET == TRUE){
  invisible(DBI::dbExecute(con, sprintf("
    CREATE OR REPLACE TEMP VIEW excluded_singletons AS
    SELECT DISTINCT
      s.DerepID
    FROM singleton_derep s
    JOIN joined_filtered j
      ON j.SeqID = s.DerepID
    WHERE j.MEEP > %s
  ", sql_num(MAXMEEP))))

  singleton_filtered <- DBI::dbGetQuery(con, "SELECT COUNT(*) AS n_filtered FROM excluded_singletons")
  cat("... Number of singleton dereplicated sequences after filtering by MEEP score: ", singleton_filtered$n_filtered[[1]], "\n")
} else {
  invisible(DBI::dbExecute(con, "
    CREATE OR REPLACE TEMP VIEW excluded_singletons AS
    SELECT DerepID
    FROM singleton_derep
  "))
}

excluded_singletons <- DBI::dbGetQuery(con, "SELECT COUNT(*) AS n_excluded FROM excluded_singletons")

if(excluded_singletons$n_excluded[[1]] > 0){
  cat("..Removing singleton dereplicated sequences\n")
}

invisible(DBI::dbExecute(con, "
  CREATE OR REPLACE TEMP VIEW joined_pruned AS
  SELECT
    j.*
  FROM joined_filtered j
  LEFT JOIN excluded_singletons e USING (DerepID)
  WHERE e.DerepID IS NULL
"))

records_after_singletons <- DBI::dbGetQuery(con, "SELECT COUNT(*) AS n_records FROM joined_pruned")
cat("... Number of records in sequence table after removing singleton dereplicated sequences: ",
    records_after_singletons$n_records[[1]], "\n")

cat("\n..Summarizing dereplicated-sequence abundance\n")
sample_expr_sql <- if(MERGE_SAMPLES == TRUE){
  cat("\n... Merging sample replicates (e.g., re-sequenced samples)\n")
  "CASE WHEN strpos(SampleID, '__') > 0 THEN split_part(SampleID, '__', 2) ELSE NULL END"
} else {
  "SampleID"
}

if(MERGE_SAMPLES == TRUE){
  cat(".... Summarizing abundance by sample and dereplicated sequence\n")
} else {
  cat("... Summarizing abundance by sample and dereplicated sequence\n")
}

RES <- DBI::dbGetQuery(con, sprintf("
  SELECT
    DerepID,
    %s AS SampleID,
    SUM(Abundance) AS Abundance
  FROM joined_pruned
  GROUP BY 1, 2
  ORDER BY 1, 2
", sample_expr_sql))
setDT(RES)

if(nrow(RES) == 0L){
  RES <- data.table(DerepID = character(), SampleID = character(), Abundance = numeric())
} else {
  RES[, Abundance := as.numeric(Abundance)]
}

cat("... Number of rows in summarized long table: ", nrow(RES), "\n")

cat("\n..Disconnecting DuckDB\n")
DBI::dbDisconnect(con, shutdown = TRUE)
con <- NULL


######################################
###################################### Prepare long and wide outputs
######################################

setorder(RES, SampleID, -Abundance, DerepID)
setattr(x = RES, name = "Samples_merged", value = MERGE_SAMPLES)

cat("\nPreparing wide output metadata\n")
n_derep <- uniqueN(RES$DerepID)
n_smp   <- uniqueN(RES$SampleID)
n_cll   <- as.double(n_derep) * as.double(n_smp)

cat("...In total, there are ", n_derep, " dereplicated sequences and ", n_smp, " samples\n")
cat("...The total number of cells in the wide table will be ", format(n_cll, scientific = FALSE), "\n")

if(nrow(RES) > 0L){
  derep_order <- RES[, .(TotalAbundance = sum(Abundance, na.rm = TRUE)), by = "DerepID"][
    order(-TotalAbundance, DerepID), DerepID]
  sample_order <- sort(unique(RES$SampleID), na.last = TRUE)
  setindex(RES, DerepID)
} else {
  derep_order  <- character()
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
  derep_order = derep_order,
  sample_order = sample_order,
  file = "OTU_table_wide.txt.gz")

cat("..Preparing wide table [R]\n")
if(WIDE_RDATA_FORMAT == "dense"){
  WIDE_OBJ <- build_dense_wide(
    res = RES,
    derep_order = derep_order,
    sample_order = sample_order)
} else {
  WIDE_OBJ <- build_sparse_wide(
    res = RES,
    derep_order = derep_order,
    sample_order = sample_order)
}
attr(WIDE_OBJ, "Samples_merged") <- MERGE_SAMPLES

cat("..Exporting wide table [R]\n")
saveRDS.gz(
  object  = WIDE_OBJ,
  file    = "OTU_table_wide.RData",
  threads = CPUTHREADS)


######################################
###################################### Export dereplicated sequences
######################################

cat("\nExporting dereplicated sequences to FASTA\n")
cat("..Preparing sequences\n")
cat("... Exporting dereplicated sequences from parquet with DuckDB\n")
export_derep_fasta(
  seqs_parquet = SEQS_PQ,
  derep_order  = derep_order,
  file         = "OTUs.fa.gz",
  threads      = CPUTHREADS)


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
