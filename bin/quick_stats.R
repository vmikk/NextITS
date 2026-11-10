#!/usr/bin/env Rscript

## Summarise number of reads (demultiplexed and primer-checked)

# quick_stats.R \
#   --raw          Counts_1.RawData.txt \
#   --qc           Counts_2.QC.txt \
#   --demuxed      Counts_3.Demux.txt \
#   --primer       Counts_4.PrimerCheck.txt \
#   --primerartef  Counts_4.PrimerArtefacts.txt \
#   --threads      4



############################################## Parse input parameters

## Check time
start_time <- Sys.time()


cat("Parsing input options and arguments...\n")

suppressPackageStartupMessages(require(optparse))

## Parse arguments
option_list <- list(

  make_option("--raw",        action="store", default=NA, type='character', help="Raw read counts"),
  make_option("--qc",         action="store", default=NA, type='character', help="Counts of reads passed QC"),
  make_option("--demuxed",    action="store", default=NA, type='character', help="Counts of demultiplexed reads"),
  make_option("--primer",     action="store", default=NA, type='character', help="Counts of reads with both primers detected"),
  make_option("--primerartef",action="store", default=NA, type='character', help="Counts of primer artefacts"),
  make_option("--reoriented", action="store", default=NA, type='character', help="Illumina: per-sample input and reoriented read pairs"),
  make_option("--merged",     action="store", default=NA, type='character', help="Illumina: per-sample merged reads"),
  make_option("--joined",     action="store", default=NA, type='character', help="Illumina: counts of joined non-merged read pairs"),
  make_option(c("-t", "--threads"), action="store", default=4L, type='integer', help="Number of CPU threads for arrow, default 4")
)
opt <- parse_args(OptionParser(option_list=option_list))


## Function to convert text "NA"s to NA
to_na <- function(x){
  if(is.null(x) || length(x) == 0 || x %in% c("NA", "null", "Null", "")){ x <- NA }
  return(x)
}

## Assign variables
RAW         <- to_na( opt$raw )
QC          <- to_na( opt$qc )
DEMUXED     <- to_na( opt$demuxed )
PRIMER      <- to_na( opt$primer )
PRIMERARTEF <- to_na( opt$primerartef )
REORIENTED  <- to_na( opt$reoriented )
MERGED      <- to_na( opt$merged )
JOINED      <- to_na( opt$joined )
CPUTHREADS  <- as.numeric( opt$threads )

## Is the input file present and non-empty (and not a `no_*` placeholder)?
is_usable <- function(x){
  if(is.na(x)){ return(FALSE) }
  file.exists(x) && !dir.exists(x) && file.size(x) > 0 && !grepl("^no_", basename(x))
}

## Illumina paired-end data?
ILLUMINA <- is_usable(REORIENTED)

## Log assigned variables
cat(paste("Counts - RawData: " ,     RAW, "\n", sep=""))
cat(paste("Counts - QC: " ,          QC, "\n", sep=""))
cat(paste("Counts - Demux: " ,       DEMUXED, "\n", sep=""))
cat(paste("Counts - PrimerCheck: " , PRIMER, "\n", sep=""))
cat(paste("Counts - Primer Artefacts: " ,   PRIMERARTEF, "\n", sep=""))
cat(paste("Illumina - reoriented pairs: ", REORIENTED, "\n", sep=""))
cat(paste("Illumina - merged reads: ",     MERGED,     "\n", sep=""))
cat(paste("Illumina - joined pairs: ",     JOINED,     "\n", sep=""))
cat(paste("Number of CPU threads to use: ", CPUTHREADS,     "\n", sep=""))

cat("\n")


############################################## data for debuging

# RAW         <- "Counts_1.RawData.txt"
# QC          <- "Counts_2.QC.txt"
# DEMUXED     <- "Counts_3.Demux.txt"
# PRIMER      <- "Counts_4.PrimerCheck.txt"
# PRIMERARTEF <- "Counts_4.PrimerArtefacts.txt"
# REORIENTED  <- NA
# MERGED      <- NA
# JOINED      <- NA
# CPUTHREADS  <- 6


############################################## Load packages and data

cat("Loading R packages...\n")

load_pckg <- function(pkg = "data.table"){
  suppressPackageStartupMessages( library(package = pkg, character.only = TRUE) )
  cat(paste(pkg, packageVersion(pkg), "\n"))
}

load_pckg("data.table")
load_pckg("openxlsx")

cat("\n")

## Set CPU thread number
cat("Setting number of CPU threads to: ", CPUTHREADS, "\n")
setDTthreads(threads = CPUTHREADS)  # for data.table


######################################
###################################### Load the data
######################################

cat("\nLoading input data\n")

## Read a seqkit table (NULL if absent)
read_seqkit <- function(x){
  if(!is_usable(x)){ return(NULL) }
  res <- fread(x, colClasses = list(character = "file"))
  if(nrow(res) == 0){ return(NULL) }
  return(res)
}

## Read a per-sample table (NULL if absent)
read_custom <- function(x){
  if(!is_usable(x)){ return(NULL) }
  res <- fread(x, colClasses = list(character = "SampleID"))
  if(nrow(res) == 0){ return(NULL) }
  return(res)
}


#### Per-dataset stats

cat("..Loading raw counts\n")
RAW <- read_seqkit(RAW)

cat("..Loading QC counts\n")
QC <- read_seqkit(QC)

#### Per-sample stats

SEQKITCOUNTS <- list()
CUSTOMCOUNTS <- list()

cat("..Loading demux counts\n")
if(ILLUMINA){
  ## Demultiplexed read pairs = input of the reorientation step
  REORIENT_data <- read_custom(REORIENTED)
  if(!is.null(REORIENT_data)){
    SEQKITCOUNTS$DEMUXED <- REORIENT_data[ , .(file = SampleID, num_seqs = Input_Pairs) ]
  }
} else {
  SEQKITCOUNTS$DEMUXED <- read_seqkit(DEMUXED)
}

if(ILLUMINA){
  cat("..Loading Illumina-specific counts (reorientation, merging, QC, joining)\n")
  if(!is.null(REORIENT_data)){
    CUSTOMCOUNTS$REORIENTED <- REORIENT_data[ , .(SampleID, NumReads = Reoriented_Pairs) ]
  }

  MERGE_data <- read_custom(MERGED)
  if(!is.null(MERGE_data)){
    CUSTOMCOUNTS$MERGED    <- MERGE_data[ , .(SampleID, NumReads = Merged_Reads) ]
    CUSTOMCOUNTS$NOTMERGED <- MERGE_data[ , .(SampleID, NumReads = NotMerged_Pairs) ]
  }

  ## QC is performed per sample (on merged reads)
  if(!is.null(QC)){
    SEQKITCOUNTS$QC <- copy(QC)
  }

  SEQKITCOUNTS$JOINED <- read_seqkit(JOINED)
}

cat("..Loading primer-checked data counts\n")
SEQKITCOUNTS$PRIMER      <- read_seqkit(PRIMER)
SEQKITCOUNTS$PRIMERARTEF <- read_seqkit(PRIMERARTEF)


## Remove NULL-files
null_seqk <- vapply(SEQKITCOUNTS, function(x){ is.null(x) || nrow(x) == 0 }, logical(1))
null_cust <- vapply(CUSTOMCOUNTS, function(x){ is.null(x) || nrow(x) == 0 }, logical(1))

if(any(null_seqk) || any(null_cust)){
  cat("Some files with counts are missing:\n")
  cat(".. ", paste(c(names(SEQKITCOUNTS)[ null_seqk ], names(CUSTOMCOUNTS)[ null_cust ]), collapse = ", "), "\n")
}
SEQKITCOUNTS <- SEQKITCOUNTS[ ! null_seqk ]
CUSTOMCOUNTS <- CUSTOMCOUNTS[ ! null_cust ]


## Process seqkit counts
seqkit_process <- function(x){
  if(nrow(x) > 0){

    ## Remove reudndant columns
    x <- x[ , .(file, num_seqs) ]

    ## Remove file extensions
    x[ , file := sub(pattern = ".fastq.gz$",        replacement = "", x = file) ]
    x[ , file := sub(pattern = ".fq.gz$",           replacement = "", x = file) ]
    x[ , file := sub(pattern = ".fa.gz$",           replacement = "", x = file) ]
    x[ , file := sub(pattern = "_PrimerChecked$",   replacement = "", x = file) ]
    x[ , file := sub(pattern = "_PrimerArtefacts$", replacement = "", x = file) ]

    ## Joined non-merged Illumina reads are pooled with the merged reads of the same sample
    x[ , file := sub(pattern = "_JoinedPE$", replacement = "", x = file) ]
    x <- x[ , .(num_seqs = sum(num_seqs, na.rm = TRUE)), by = "file" ]

  }
  return(x)
}

## Process custom counts (`SampleID` & `NumReads`)
custom_process <- function(x){
  if(nrow(x) > 0){
    x <- copy(x)
    setnames(x = x, old = "SampleID", new = "file")
    x <- x[ , .(NumReads = sum(NumReads, na.rm = TRUE)), by = "file" ]
  }
  return(x)
}


cat("Processing data\n")
SEQKITCOUNTS <- lapply(SEQKITCOUNTS, seqkit_process)
CUSTOMCOUNTS <- lapply(CUSTOMCOUNTS, custom_process)

## Rename columns
rename_counts <- function(lst, key, old, new){
  if(!is.null(lst[[key]])){
    setnames(x = lst[[key]], old = old, new = new, skip_absent = TRUE)
  }
  return(lst)
}
SEQKITCOUNTS <- rename_counts(SEQKITCOUNTS, "DEMUXED",     "num_seqs", "Demultiplexed_Reads")
SEQKITCOUNTS <- rename_counts(SEQKITCOUNTS, "QC",          "num_seqs", "QC_Reads")
SEQKITCOUNTS <- rename_counts(SEQKITCOUNTS, "JOINED",      "num_seqs", "Joined_Reads")
SEQKITCOUNTS <- rename_counts(SEQKITCOUNTS, "PRIMER",      "num_seqs", "PrimerChecked_Reads")
SEQKITCOUNTS <- rename_counts(SEQKITCOUNTS, "PRIMERARTEF", "num_seqs", "PrimerArtefacts_Reads")
CUSTOMCOUNTS <- rename_counts(CUSTOMCOUNTS, "REORIENTED",  "NumReads", "Reoriented_Pairs")
CUSTOMCOUNTS <- rename_counts(CUSTOMCOUNTS, "MERGED",      "NumReads", "Merged_Reads")
CUSTOMCOUNTS <- rename_counts(CUSTOMCOUNTS, "NOTMERGED",   "NumReads", "NotMerged_Pairs")

## Merge seqkit and custom counts into a single list
cat("Pooling per-sample counts\n")
COUNTS <- c(SEQKITCOUNTS, CUSTOMCOUNTS)

## Pool per-file estimates
merge_dt <- function(x,y){ merge(x, y, by = "file", all = TRUE) }
if(length(COUNTS) > 0){
  PER_SAMPLE_COUNTS_merged <- Reduce(f = merge_dt, x = COUNTS)
} else {
  cat("..No per-sample counts found\n")
  PER_SAMPLE_COUNTS_merged <- data.table(file = character())
}

## Columns that must be present even if the input is missing
for(cc in c("Demultiplexed_Reads", "PrimerChecked_Reads", "PrimerArtefacts_Reads")){
  if(! cc %in% colnames(PER_SAMPLE_COUNTS_merged)){
    PER_SAMPLE_COUNTS_merged[ , (cc) := 0 ]
  }
}
rm(cc)

## Replace NAs with zeros
cat("Replacing NAs with zero\n")
for (j in setdiff(colnames(PER_SAMPLE_COUNTS_merged), "file")){
  set(PER_SAMPLE_COUNTS_merged, which(is.na(PER_SAMPLE_COUNTS_merged[[j]])), j, 0)
}
rm(j)

## Estimate percentage of primer artefacts
## (no reads at all - 0%; no good reads - 100%)
cat("Estimating percentage of primer artefacts\n")
PER_SAMPLE_COUNTS_merged[ , PrimerArtefacts_Percent := 0 ]
PER_SAMPLE_COUNTS_merged[
  PrimerChecked_Reads + PrimerArtefacts_Reads > 0,
  PrimerArtefacts_Percent := round(
    PrimerArtefacts_Reads / (PrimerChecked_Reads + PrimerArtefacts_Reads) * 100,
    2)
  ]

## Percentage of merged read pairs (Illumina)
if("Merged_Reads" %in% colnames(PER_SAMPLE_COUNTS_merged) && "Reoriented_Pairs" %in% colnames(PER_SAMPLE_COUNTS_merged)){
  PER_SAMPLE_COUNTS_merged[ , Merged_Percent := 0 ]
  PER_SAMPLE_COUNTS_merged[ Reoriented_Pairs > 0, Merged_Percent := round( Merged_Reads / Reoriented_Pairs * 100, 2) ]
}

## Reorder columns
setcolorder(PER_SAMPLE_COUNTS_merged,
  skip_absent = TRUE,
  neworder = c(
    "file", "Demultiplexed_Reads",
    "Reoriented_Pairs", "Merged_Reads", "NotMerged_Pairs", "Merged_Percent", "Joined_Reads", "QC_Reads",
    "PrimerChecked_Reads", "PrimerArtefacts_Reads", "PrimerArtefacts_Percent"))

setorder(PER_SAMPLE_COUNTS_merged, file)


## Prepare per-run stats
cat("Preparing per-run stats\n")
sum_or_na <- function(x){ if(is.null(x)){ NA_real_ } else { as.numeric(sum(x$num_seqs, na.rm = TRUE)) } }

PER_RUN_COUNTS_merged <- data.table(
  Total_Number_Of_Reads = sum_or_na(RAW),
  Reads_Demultiplexed   = sum(PER_SAMPLE_COUNTS_merged$Demultiplexed_Reads, na.rm = TRUE),
  Reads_Passed_QC       = sum_or_na(QC),
  Reads_PrimerChecked   = sum(PER_SAMPLE_COUNTS_merged$PrimerChecked_Reads, na.rm = TRUE)
  )

if(ILLUMINA){
  PER_RUN_COUNTS_merged[ , Pairs_Reoriented := sum(PER_SAMPLE_COUNTS_merged$Reoriented_Pairs, na.rm = TRUE) ]
  PER_RUN_COUNTS_merged[ , Reads_Merged     := sum(PER_SAMPLE_COUNTS_merged$Merged_Reads,     na.rm = TRUE) ]
  PER_RUN_COUNTS_merged[ , Percentage_Merged := NA_real_ ]
  PER_RUN_COUNTS_merged[ Pairs_Reoriented > 0, Percentage_Merged := round(Reads_Merged / Pairs_Reoriented * 100, 1) ]
  if("Joined_Reads" %in% colnames(PER_SAMPLE_COUNTS_merged)){
    PER_RUN_COUNTS_merged[ , Reads_Joined := sum(PER_SAMPLE_COUNTS_merged$Joined_Reads, na.rm = TRUE) ]
  }
}

## Estimate percentages (relative to the number of raw reads or read pairs)
cat("..Estimating per-run percentages\n")
pct <- function(x, y, digits = 1){
  ifelse(is.finite(x) & is.finite(y) & y > 0, round(x / y * 100, digits), NA_real_)
}

PER_RUN_COUNTS_merged[ , Percentage_QC :=
  pct(Reads_Passed_QC, Total_Number_Of_Reads) ]

PER_RUN_COUNTS_merged[ , Percentage_Demultiplexed :=
  pct(Reads_Demultiplexed, Total_Number_Of_Reads) ]

PER_RUN_COUNTS_merged[ , Percentage_PrimerChecked :=
  pct(Reads_PrimerChecked, Total_Number_Of_Reads) ]

PER_RUN_COUNTS_merged[ , PrimerArtefacts_Reads := sum(PER_SAMPLE_COUNTS_merged$PrimerArtefacts_Reads, na.rm = TRUE) ]
PER_RUN_COUNTS_merged[ , PrimerArtefacts_Percent := 0 ]
PER_RUN_COUNTS_merged[ Reads_PrimerChecked + PrimerArtefacts_Reads > 0,
  PrimerArtefacts_Percent := round(PrimerArtefacts_Reads / (Reads_PrimerChecked + PrimerArtefacts_Reads) * 100, 2) ]


## Export summary stats
cat("Exporting results\n")
fwrite(PER_SAMPLE_COUNTS_merged, file = "per_sample.tsv", sep = "\t")
fwrite(PER_RUN_COUNTS_merged,    file = "per_run.tsv",    sep = "\t")
write.xlsx(list(
  "per_sample" = PER_SAMPLE_COUNTS_merged,
  "per_run"    = PER_RUN_COUNTS_merged
  ),
  file = "Run_summary.xlsx", colNames = TRUE)


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
