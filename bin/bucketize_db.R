#!/usr/bin/env Rscript

## Aim - evenly distribute sequence clusters across a specified number of buckets.
## The goal is to have the total length of sequences in each bucket as equal as possible.

## Number of buckets can be automatically selected
## (e.g., to avoid the DADA2s' error message `long vectors not supported yet`, related with > 2^31 elements)

## Usage examples:
# bucketize_db.R \
#  --db         stat_clusters.txt \
#  --fasta      Input.fa.gz \
#  --summary    bucket_summary.txt \
#  --numbuckets 10 \
#  --threads    10


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


############################################## Load packages

cat("Loading R packages...\n")

load_pckg <- function(pkg = "data.table"){
  suppressPackageStartupMessages( library(package = pkg, character.only = TRUE) )
  cat(paste(pkg, packageVersion(pkg), "\n"))
}

load_pckg("data.table")
load_pckg("Biostrings")
load_pckg("plyr")

if(THREADS < 1){ THREADS <- 1 }
if(THREADS > 1){
  cat("Preparing multi-threaded setup\n")

  load_pckg("doFuture")
  registerDoFuture()
  plan(multicore, workers = THREADS)
  options(future.globals.maxSize = 6e10)  # 60GB
  
  setDTthreads(threads = THREADS)         # for data.table

  parall <- TRUE

} else {
  parall <- FALSE
  setDTthreads(threads = 1)
}

cat("\n")


############################################## Workflow

## Load seq stats
cat("..Loading input sequences\n")
seqs <- readDNAStringSet(filepath = FASTA)

## Load clustering file
cat("..Loading clustering file\n")
DB <- fread(file = DATABASE,
  sep = "\t", header = FALSE,
  col.names = c("Cluster", "Member"))

## Estimate sequence length
cat("..Estimating total length of the sequences\n")
seqt <- data.table(Member = names(seqs), Len = width(seqs))
DB <- merge(x = DB, y = seqt, by = "Member", all.x = TRUE)
rm(seqt)

## Estimate number of sequences per cluster and the total length of sequences
cat("..Estimating cluster sizes\n")
datt <- DB[ , .(num_seqs = .N, sum_len = sum(Len, na.rm = TRUE)), by = "Cluster" ]

## Sort clusters by the number of sequenes in descending order
cat("..Sorting clusters\n")
setorder(datt, -sum_len, -num_seqs)


cat("..Bucketizing\n")

if(is.na(NBUCKETS)){
  cat("...Number of buckets is not specified, using automatic selection\n")

  ## For DADA2, a matrix with quality values is required  `as(Biostrings::quality(fq), "matrix")`
  ## It should not exceed 2^31 (2147483648) elemens,
  ## Meaning that `num_seq * len_seq` must be < 2^31

  ## Calculate approximate estimate for the maximum number of sequences per bucket
  maxseqs <- 2^31 / max(DB$Len)   # quantile(x = DB$Len, probs = 0.99)

  ## Number of buckets
  NBUCKETS <- ceiling(nrow(DB) / maxseqs)

  cat("...The sugested number of buckets is ", NBUCKETS, "\n")
}


## Initializing buckets and bucket sizes
buckets <- vector("list", length = NBUCKETS)
bucket_size_numseqs <- numeric(NBUCKETS)
bucket_size_lenseqs <- numeric(NBUCKETS)

## Distributing files into buckets
## By starting with the largest files and placing each one in the currently smallest bucket, 
## we try to prevent any single bucket from becoming significantly larger than the others
for (i in 1:nrow(datt)) {
  
  ## Find the bucket with the minimum total sequence length
  min_bucket_index <- which.min(bucket_size_lenseqs)
  
  ## Add the cluster ID to the chosen bucket
  buckets[[ min_bucket_index ]] <- c(
    buckets[[ min_bucket_index ]],
    datt[i, ]$Cluster
    )
  
  # Updating the total sequence length of the chosen bucket
  bucket_size_lenseqs[ min_bucket_index ] <- bucket_size_lenseqs[min_bucket_index] + datt[i, ]$sum_len
  bucket_size_numseqs[ min_bucket_index ] <- bucket_size_numseqs[min_bucket_index] + datt[i, ]$num_seqs

}

cat("..Bucket summary:\n\n")

## Prepare bucket summary
smr <- data.table(
    BucketID  = 1:length(buckets),
    Num_files = laply(.data = buckets, .fun = function(x){ length(x) }),
    sum_len   = bucket_size_lenseqs,
    num_seqs  = bucket_size_numseqs)

print(smr)

## Add percentages
smr[ , NumFiles_Percent := round(Num_files / sum(Num_files) * 100, 2) ]
smr[ , TotLen_Percent   := round(sum_len   / sum(sum_len)   * 100, 2) ]
smr[ , TotSeqs_Percent  := round(num_seqs  / sum(num_seqs)  * 100, 2) ]


cat("\n\n..Exporting FASTA file for each bucket\n")

## Exporting function
export_bucket <- function(clustnum = 1){

  ## IDs of cluster representatives
  clustids <- buckets[[ clustnum ]]

  ## Find sequence IDs to export
  ids <- data.table(SeqID = DB[ Cluster %in% clustids ]$Member)

  ## Sort sequences by size
  ids[ , Size := tstrsplit(SeqID, split = ";", keep = 2) ]
  ids[ , Size := as.numeric( sub(pattern = "size=", replacement = "", x = Size) ) ]
  setorder(ids, -Size, SeqID)

  ## Cluster ID with leading zero
  cl <- sprintf(paste0("%0", nchar(NBUCKETS), "d"), clustnum)

  ## Extract and export
  writeXStringSet(
    x = seqs[ ids$SeqID ],
    filepath = paste0("bucket_", cl, ".fa.gz"),
    compress = TRUE,
    format = "fasta",
    width = 9999)

}

a_ply(
  .data     = seq_along(buckets),
  .margins  = 1,
  .fun      = export_bucket,
  .parallel = parall)


## Bucket summary
cat("..Exporting bucket summary\n")
fwrite(x = smr, file = SUMMARY, sep = "\t", col.names = TRUE)


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
