#!/usr/bin/env Rscript

## Aim - evenly distribute sequence clusters across a specified number of buckets.
## The goal is to have the total length of sequences in each bucket as equal as possible.

## Number of buckets can be automatically selected
## (e.g., to avoid the DADA2s' error message `long vectors not supported yet`, related with > 2^31 elements)

## Notes - allocating large number of CPUs will not add much (used only for data.table).
#  Currently, exporting FASTA files is single-threaded
#  (with multi-threaded export, RAM usage will high, as for each thread full sequnece list will be duplicated)

## Usage examples:
# bucketize_db.R \
#  --db         Cluster_Membership.txt \
#  --fasta      Input.fa.gz \
#  --summary    bucket_summary.txt \
#  --numbuckets 10 \
#  --threads    10

## Input data:
# - DB: TSV, no header, two columns (ClusterID, QueryID)
# - FASTA file with sequences, header format: >QueryID;size=Abundance


## Check time
start_time <- Sys.time()

cat("\nParsing input options and arguments...\n")

suppressPackageStartupMessages(require(optparse))

## Parse arguments
option_list <- list(
  make_option(c("-d", "--db"),         action="store", default="DB_clu.tsv",         type='character', help="Clustering database"),
  make_option(c("-f", "--fasta"),      action="store", default="Input.fa.gz",        type='character', help="Input sequences in FASTA format"),
  make_option(c("-s", "--summary"),    action="store", default="bucket_summary.txt", type='character', help="Output file summary information"),
  make_option(c("-n", "--numbuckets"), action="store", default=NA,   type='integer', help="Number of output buckets (NA, for automatic selection)"),
  make_option(c("-t", "--threads"),    action="store", default=4,    type='integer', help="Number of CPU threads to use")
)
opt <- parse_args(OptionParser(option_list=option_list))

# Validation of the required arguments
if(is.na(opt$fasta)){
  stop("Input file with sequences is not specified\n")
}
if(is.na(opt$db)){
  stop("Clustering results are not specified\n")
}
if(!is.na(opt$numbuckets) & opt$numbuckets <= 1){
  stop("Number of buckets should be > 1\n")
}


## Assign variables
DATABASE <- opt$db
FASTA    <- opt$fasta
SUMMARY  <- opt$summary
NBUCKETS <- opt$numbuckets
THREADS  <- opt$threads

## Log assigned variables
cat("\nParameters specified:\n")
cat(paste("Clustering database: " ,       DATABASE, "\n", sep = ""))
cat(paste("Input sequences (FASTA): " ,   FASTA,    "\n", sep = ""))
cat(paste("Output with bucket summary: ", SUMMARY,  "\n", sep = ""))
if(is.na(NBUCKETS)){
  cat(paste("Number of buckets: ",        "auto",   "\n", sep = ""))
} else {
  cat(paste("Number of buckets: ",        NBUCKETS, "\n", sep = ""))
}

cat(paste("CPU threads: ",                THREADS,  "\n", sep = ""))
cat("\n")


############################################## Data for debugging

# DATABASE <- "DB_clu.tsv"
# FASTA    <- "Dereplicated.fa.gz"
# SUMMARY  <- "bucket_summary.txt"
# NBUCKETS <- 400
# THREADS  <- 10


############################################## Load packages

cat("Loading R packages...\n")

load_pckg <- function(pkg = "data.table"){
  suppressPackageStartupMessages( library(package = pkg, character.only = TRUE) )
  cat(paste(pkg, packageVersion(pkg), "\n"))
}

load_pckg("data.table")
load_pckg("Biostrings")

if(THREADS < 1){ THREADS <- 1 }
setDTthreads(threads = THREADS)

cat("\n")


############################################## Workflow

## Function to extract sequence size from sequence IDs
extract_seq_size <- function(seq_ids){
  size_field <- tstrsplit(seq_ids, split = ";", fixed = TRUE, keep = 2L)[[1L]]
  size_num <- suppressWarnings(as.numeric(sub(pattern = "^size=", replacement = "", x = size_field)))
  size_num[ is.na(size_num) ] <- 0
  size_num
}

## Helper function for exporting FASTA files
## it opens the output file once, then writes many small DNAStringSet chunks into that same open handle
## The main idea is to avoid the full-bucket allocation, which can be problematic for large data
write_bucket_fasta <- function(seq_index, filepath, seqs, chunk_nseq){
  # seq_index  = indices of the sequences to export
  # filepath   = path to the output file
  # seqs       = DNAStringSet object with all sequences
  # chunk_nseq = number of sequences to write in each chunk
  
  nseq <- length(seq_index)
  
  ## Open a compressed file handle
  filexp_list <- XVector:::open_output_file(
    filepath = filepath,
    append   = FALSE,
    compress = TRUE,
    compression_level = NA)
  on.exit(Biostrings:::.close_filexp_list(filexp_list), add = TRUE)

  ## If there are no sequences to export, write an empty DNAStringSet
  if(nseq == 0L){
    Biostrings:::.write_XStringSet_to_fasta(
      x = seqs[integer()],
      filexp_list = filexp_list,
      width = 9999)

    return(invisible(NULL))
  }

  ## Write chunks of sequences to the compressed file handle
  for(chunk_start in seq.int(1L, nseq, by = chunk_nseq)){
    chunk_end <- min(chunk_start + chunk_nseq - 1L, nseq)

    Biostrings:::.write_XStringSet_to_fasta(
      x = seqs[ seq_index[chunk_start:chunk_end] ],
      filexp_list = filexp_list,
      width = 9999)
  }

  invisible(NULL)
}

## Load seq stats
cat("..Loading input sequences\n")
seqs <- readDNAStringSet(filepath = FASTA)
seq_names     <- names(seqs)
seq_widths    <- width(seqs)
max_seq_width <- max(seq_widths)

## Bound per-write allocations during FASTA export
EXPORT_CHUNK_BASES <- 5e7
EXPORT_CHUNK_NSEQ <- max(1L, as.integer(floor(EXPORT_CHUNK_BASES / max_seq_width)))

cat("..FASTA export chunk size: up to ", format(EXPORT_CHUNK_NSEQ, big.mark = ","), " sequences per write\n", sep = "")

## Load clustering file
cat("..Loading clustering file\n")
DB <- fread(file = DATABASE,
  sep = "\t", header = FALSE,
  col.names = c("Cluster", "Member"))

## Match clustering members to FASTA entries
cat("..Matching clustering members to FASTA entries\n")
DB[ , SeqIndex := match(Member, seq_names) ]

if(anyNA(DB$SeqIndex)){
  missing_ids <- DB[ is.na(SeqIndex), head(Member, 5L) ]
  stop(
    "Could not match ", sum(is.na(DB$SeqIndex)),
    " clustering members to FASTA entries. Examples: ",
    paste(missing_ids, collapse = ", "),
    "\n")
}

## Estimate sequence length and size
cat("..Estimating total length of the sequences\n")
DB[ , Len := seq_widths[SeqIndex] ]
DB[ , Size := extract_seq_size(Member) ]
rm(seq_names, seq_widths)
invisible(gc(FALSE))

## Estimate number of sequences per cluster and the total length of sequences
cat("..Estimating cluster sizes\n")
datt <- DB[ , .(num_seqs = .N, sum_len = sum(Len, na.rm = TRUE)), by = "Cluster" ]

## Sort clusters by the number of sequenes in descending order
cat("..Sorting clusters\n")
setorder(datt, -sum_len, -num_seqs)
cluster_ids      <- datt[[ "Cluster"  ]]
cluster_num_seqs <- datt[[ "num_seqs" ]]
cluster_sum_len  <- datt[[ "sum_len"  ]]


cat("..Bucketizing\n")

if(is.na(NBUCKETS)){
  cat("...Number of buckets is not specified, using automatic selection\n")

  ## For DADA2, a matrix with quality values is required  `as(Biostrings::quality(fq), "matrix")`
  ## It should not exceed 2^31 (2147483648) elemens,
  ## Meaning that `num_seq * len_seq` must be < 2^31

  ## Calculate approximate estimate for the maximum number of sequences per bucket
  maxseqs <- 2^31 / max_seq_width   # quantile(x = DB$Len, probs = 0.99)

  ## Number of buckets
  NBUCKETS <- ceiling(nrow(DB) / maxseqs)

  cat("...The sugested number of buckets is ", NBUCKETS, "\n")
}


## Initializing bucket assignments and bucket sizes
cluster_bucket      <- integer(nrow(datt))
bucket_size_numseqs <- numeric(NBUCKETS)
bucket_size_lenseqs <- numeric(NBUCKETS)

## Distributing files into buckets
## By starting with the largest files and placing each one in the currently smallest bucket, 
## we try to prevent any single bucket from becoming significantly larger than the others
for (i in seq_len(nrow(datt))) {
  
  ## Find the bucket with the minimum total sequence length
  min_bucket_index <- which.min(bucket_size_lenseqs)
  
  ## Assign the cluster to the chosen bucket
  cluster_bucket[i] <- min_bucket_index
  
  # Updating the total sequence length of the chosen bucket
  bucket_size_lenseqs[ min_bucket_index ] <- bucket_size_lenseqs[min_bucket_index] + cluster_sum_len[i]
  bucket_size_numseqs[ min_bucket_index ] <- bucket_size_numseqs[min_bucket_index] + cluster_num_seqs[i]

}

cat("..Bucket summary:\n\n")

## Prepare bucket summary
smr <- data.table(
    BucketID     = seq_len(NBUCKETS),
    Num_clusters = tabulate(cluster_bucket, nbins = NBUCKETS),
    sum_len      = bucket_size_lenseqs,
    num_seqs     = bucket_size_numseqs)

print(smr)

## Add percentages
smr[ , NumClust_Percent := round(Num_clusters / sum(Num_clusters) * 100, 2) ]
smr[ , TotLen_Percent   := round(sum_len   / sum(sum_len)   * 100, 2) ]
smr[ , TotSeqs_Percent  := round(num_seqs  / sum(num_seqs)  * 100, 2) ]

## Bucket summary
cat("..Exporting bucket summary\n")
fwrite(x = smr, file = SUMMARY, sep = "\t", col.names = TRUE)


cat("\n\n..Exporting FASTA file for each bucket\n")

## Prepare export index once instead of scanning DB for every bucket
cat("..Preparing export index\n")
cluster_bucket_dt <- data.table(Cluster = cluster_ids, BucketID = cluster_bucket)
DB[ cluster_bucket_dt, on = "Cluster", BucketID := i.BucketID ]
DB[ , c("Cluster", "Member", "Len") := NULL ]
setcolorder(DB, c("BucketID", "Size", "SeqIndex"))
setorder(DB, BucketID, -Size, SeqIndex)

export_ranges <- DB[ , .(start = .I[1L], end = .I[.N]), by = BucketID ]
export_start <- integer(NBUCKETS)
export_end   <- integer(NBUCKETS)
export_start[ export_ranges$BucketID ] <- export_ranges$start
export_end[ export_ranges$BucketID ]   <- export_ranges$end
export_seq_idx <- DB[["SeqIndex"]]

rm(DB, datt, cluster_bucket, cluster_bucket_dt, export_ranges, cluster_ids, cluster_num_seqs, cluster_sum_len)
invisible(gc(FALSE))

## Exporting function
export_bucket <- function(clustnum = 1){

  cat("...Bucket ", clustnum, "\n")

  ## Cluster ID with leading zero
  cl <- sprintf(paste0("%0", nchar(NBUCKETS), "d"), clustnum)

  ## Extract and export in bounded chunks
  row_start <- export_start[clustnum]
  row_end   <- export_end[clustnum]
  seq_index <- if(row_start > 0L) export_seq_idx[ row_start:row_end ] else integer()

  write_bucket_fasta(
    seq_index  = seq_index,
    filepath   = paste0("bucket_", cl, ".fa.gz"),
    seqs       = seqs,
    chunk_nseq = EXPORT_CHUNK_NSEQ)

  invisible(gc(FALSE))

}

for(clustnum in seq_len(NBUCKETS)){
  export_bucket(clustnum)
}


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
