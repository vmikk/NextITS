#!/usr/bin/env Rscript

## Script to pool quality-filtered and trimmed sequences from multiple sequencing runs
## And summarize sequence abundance at OTU level (per sample)
## 

# Input:
#   1. UC file                         (`UC_Pooled.parquet`)
#   2. FASTA file with OTU sequences   (`Clustered.fa.gz`)
#   3. Max MEEP score
#   4. Max de novo chimera score

# Outputs:
#  - OTU table in long format    (`OTU_table_long.txt.gz` & `OTU_table_long.RData`)
#  - OTU table in wide format    (`OTU_table_wide.txt.gz` & `OTU_table_wide.RData`)
#  - FASTA file with sequences   (`OTUs.fa.gz`)

## Usage:
# ./pool_seq_runs.R \
#    --uc   "UC_Pooled.parquet" \
#    --otus "Clustered.fa.gz" \
#    --maxmeep 0.5 \
#    --maxchim 0.6 \
#    --recoverdenovo TRUE \
#    --recoversinglet TRUE \
#    --mergesamples TRUE \
#    --threads 4

## Do-novo chimera recovery:
# if a sequence identified as putative chimera was observed in the other samples,
# where there is no evidence that it is chimeric, it will be recovered

## Singleton recovery
# if a within-sequencing run singleton sequence with relatively low qualty (based on MEEP)
# was observed in the other samples, it will be recovered


## TO DO:
# - load sequence tables (TAB) in parallel


############################################## Parse input parameters

## Check time
start_time <- Sys.time()


cat("Parsing input options and arguments:\n")

suppressPackageStartupMessages(require(optparse))

## Parse arguments
option_list <- list(
  make_option("--seqtab",  action="store", default=NA,  type='character', help="Sequence tables in long format with de novo chimeras removed (Parquet format)"),
  make_option("--uc",      action="store", default=NA,  type='character', help="UC file (Parquet format)"),
  make_option("--otus",    action="store", default=NA,  type='character', help="FASTA file with OTU sequences"),
  make_option("--maxmeep", action="store", default=0.5, type='double', help="Max MEEP score"),
  make_option("--recoversinglet", action="store", default=TRUE, type='logical', help="Recover singletons"),
  make_option(c("-m", "--mergesamples"), action="store", default=FALSE, type='logical', help="Merge sample replicates (default, false)"),
  make_option(c("-t", "--threads"), action="store", default=4L, type='integer', help="Number of CPU threads for arrow, default 4")
  # make_option(c("-s", "--scriptdir"),  action="store", default=getwd(), type='character', help="Directory containing source scripts")
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
if(is.na(opt$seqtab)){
  cat("Input file is not specified: sequence tables in Parquet format.\n", file=stderr())
  stop()
}
if(is.na(opt$uc)){
  cat("Input file is not specified: UC file is required.\n", file=stderr())
  stop()
}
if(is.na(opt$otus)){
  cat("Input file is not specified: FASTA file with OTU sequences.\n", file=stderr())
  stop()
}
if(opt$recoversinglet == TRUE && is.na(opt$maxmeep)){
  cat("For singleton recovery, the max MEEP score must be specified.\n", file=stderr())
  stop()
}

## Assign variables
SEQTAB        <- opt$seqtab
UCF           <- opt$uc
MAXMEEP       <- as.numeric( opt$maxmeep )
RECOV_SINGLET <- as.logical(opt$recoversinglet)
MERGE_SAMPLES <- as.logical(opt$mergesamples)
OTUS          <- opt$otus

CPUTHREADS <- as.numeric( opt$threads )
# SCRIPTDIR  <- opt$scriptdir

## Log assigned variables
cat(paste("Sequence tables (Parquet format): ", SEQTAB,        "\n", sep=""))
cat(paste("UC file (Parquet format): ",         UCF,           "\n", sep=""))
cat(paste("Max MEEP score: ",                   MAXMEEP,       "\n", sep=""))
cat(paste("Low-quality singleton recovery: ",   RECOV_SINGLET, "\n", sep=""))
cat(paste("Merge sample replicates: ",          MERGE_SAMPLES, "\n", sep=""))
cat(paste("OTU sequences: ",                    OTUS,          "\n", sep=""))
cat(paste("Number of CPU threads to use: ",     CPUTHREADS,    "\n", sep=""))
# cat(paste("Directory containing source scripts: ", SCRIPTDIR, "\n", sep=""))

cat("\n")



############################################## Data for debugging

# SEQTAB        <- "Seqs.parquet"
# UCF           <- "UC_Pooled.parquet"
# MAXMEEP       <- 0.5
# RECOV_SINGLET <- TRUE
# MERGE_SAMPLES <- TRUE
# OTUS          <- "Clustered.fa.gz"
# CPUTHREADS    <- 4


############################################## Load packages and data

cat("Loading R packages:\n")

load_pckg <- function(pkg = "data.table"){
  suppressPackageStartupMessages( library(package = pkg, character.only = TRUE) )
  cat(paste(pkg, packageVersion(pkg), "\n"))
}

load_pckg("data.table")
load_pckg("plyr")
load_pckg("metagMisc")
load_pckg("Biostrings")
load_pckg("arrow")


cat("\n")


# cat("Loading additional R funcitons...\n")
# source(file.path(SCRIPTDIR, "R_functions.R"))
# cat("\n")


## Set CPU thread number
cat("Setting number of CPU threads to: ", CPUTHREADS, "\n")
setDTthreads(threads = CPUTHREADS)     # for data.table
set_cpu_count(CPUTHREADS)              # for arrow

######################################
###################################### Load the data
######################################

## Load sequence tables
cat("\n..Loading sequence tables\n")
TAB <- arrow::read_parquet(SEQTAB)
setDT(TAB)
cat("... Total number of records: ",             nrow(TAB), "\n")
cat("... Total number unique sequences: ",       length(unique(TAB$Sequence)), "\n")
cat("... Total number unique samples (files): ", length(unique(TAB$SampleID)), "\n")

## Load UC file for globally dereplicated sequences
cat("..Loading pooled UC file\n")
UC <- open_dataset(UCF) |> dplyr::collect() |> setDT()

## Add OTU IDs to seq table
cat("... Adding OTU IDs to sequence table\n")
cat(".... Number of records in sequence table before merging: ", nrow(TAB), "\n")
TAB <- merge(x = TAB, y = UC, by = "SeqID", all.x = TRUE)
cat(".... Number of records in sequence table after merging: ",  nrow(TAB), "\n")

## Remove NA OTUs -- probably excluded seqs
if(any(is.na(TAB$OTU))){
  cat("WARNING: not all sequences were assigned to OTUs\n")
  cat("..Removing missing/excluded sequences\n")
  cat(".. ", sum(is.na(TAB$OTU)), " sequences with total abundance ", 
      sum(TAB[ is.na(OTU) ]$Abundance, na.rm = TRUE), " reads will be excluded\n")
  TAB <- TAB[ ! is.na(OTU) ]
}


## Find singleton OTUs
cat("\n..Finding singleton OTUs\n")
SINGLETONS <- TAB[ , .(Abundance = sum(Abundance, na.rm = TRUE)), by = .(OTU) ][ Abundance < 2 ]
cat("... Number of singleton OTUs: ", nrow(SINGLETONS), "\n")

## If singleton recovery is reqired
if(RECOV_SINGLET == TRUE && nrow(SINGLETONS) > 0){

  ## Add quality scores
  SINGLETONS <- merge(x = SINGLETONS, y = TAB[ , .(SeqID, MEEP)], by.x = "OTU", by.y = "SeqID", all.x = TRUE)

  ## Filter by MEEP score
  SINGLETONS <- SINGLETONS[ MEEP > MAXMEEP ]
  cat("... Number of singleton OTUs after filtering by MEEP score: ", nrow(SINGLETONS), "\n")

}

if(nrow(SINGLETONS) > 0){
  cat("..Removing singleton OTUs\n")
  TAB <- TAB[ ! OTU %in% SINGLETONS$OTU ]
  cat("... Number of records in sequence table after removing singleton OTUs: ", nrow(TAB), "\n")
}


## Summarize abundance by sample and OTU
cat("\n..Summarizing OTU abundance\n")
if(MERGE_SAMPLES == TRUE){

  cat("\n... Merging sample replicates (e.g., re-sequenced samples)\n")

  ## Extract sample names
  cat(".... Extracting sample names\n")
  TAB[ , SampleName := tstrsplit(x = SampleID, split = "__", keep = 2) ]

  cat(".... Summarizing abundance by sample and OTU\n")
  RES <- TAB[ , 
    .( Abundance  = sum(Abundance,  na.rm = TRUE) ),
    by = c("OTU", "SampleName") ]

  setnames(x = RES, old = "SampleName", new = "SampleID")

} else {

  cat("... Summarizing abundance by sample and OTU\n")
  RES <- TAB[ , 
    .( Abundance  = sum(Abundance,  na.rm = TRUE) ),
    by = c("OTU", "SampleID") ]

}

#### Reshape to wide table
cat("\nReshaping table into wide format\n")

## Check if we can reshape the table in a single pass
n_otu <- length(unique(RES$OTU))
n_smp <- length(unique(RES$SampleID))
n_cll <- as.numeric(n_otu) * as.numeric(n_smp)
cat("...In total, there are ", n_otu, " OTUs and ",  n_smp, " samples\n")
cat("...The total number of cells in the wide table will be ", n_cll, "\n")

## Reshape data in one pass
if(n_cll < 50000000){
  REW <- dcast(data = RES, 
    formula = OTU ~ SampleID, 
    fun.aggregate = sum, fill = 0, value.var = "Abundance")
} else {
## Split data into chunks, reshape, and merge back

  cat("..The input table is too large to reshape in a single pass, reshaping by chunks\n")

  ## Function to split vector into N chunks
  chunk <- function(x, n){
    if(n > 1) { res <- split(x, cut(seq_along(x), n, labels = FALSE)) }
    if(n == 1){ res <- list();  res[[1]] <- x }
    return(res)
  }

  ## Choose the number of chunks
  n_chunks <- data.table::fcase(
      n_cll <  9e7,                2L,
      n_cll >= 9e7 & n_cll < 5e8,  5L,
      n_cll >= 5e8 & n_cll < 5e9,  6L,
      n_cll >= 5e9 & n_cll < 5e10, 7L,
      n_cll >= 5e10,               8L)

  cat("...The number of chunks to process: , ", n_chunks, "\n")

  ch <- chunk(x = sort(unique(RES$SampleID)), n = n_chunks)

  ## Chunk-and-reshape loop
  REWL <- plyr::llply(
    .data = ch,
    .fun = function(x){
      
      ## Reshape to wide
      res <- dcast(
        data = RES[ SampleID %in% x , ],
        formula = OTU ~ SampleID,
        fill = 0, fun.aggregate = sum, value.var = "Abundance")
      
      ## Create key on a data.table (should improve merging speed)
      setkey(res, OTU)
      
      return(res)
      },
    .progress = "text")

  cat("...Chunk reshaping finished\n")
  cat("..Merging data into a single wide table\n")

  ## Merge chunks into a single wide table
  merge_dt <- function(x,y){ data.table::merge.data.table(x, y, by = "OTU", all = TRUE) }
  REW <- Reduce(f = merge_dt, x = REWL)
  cat("...Merging finished\n")

  ## Clean up
  cat("...Cleaning up\n")
  rm(REWL); gc()

  ## Replace NAs with zeros
  cat("...Filling missing values with zeros\n")
  for (j in seq_len(ncol(REW))){
    set(REW, which(is.na(REW[[j]])), j, 0)
  }

} ## end of reshaping

cat("...Reshaping to the wide format done!\n")


## Reorder OTU rows
cat("\n..Reordering OTU rows by total abundance\n")
otu_tots <- rowSums(REW[, -1], na.rm = TRUE)
REW <- REW[ order(otu_tots, decreasing = T), ]

## Add attributes if samples were merged
setattr(x = RES, name = "Samples_merged", value = MERGE_SAMPLES)
setattr(x = REW, name = "Samples_merged", value = MERGE_SAMPLES)


cat("\nExporting results\n")

## Export data
saveRDS.gz <- function(object, file, threads = parallel::detectCores()) {
  con <- pipe(paste0("pigz -p",threads," > ",file),"wb")
  saveRDS(object, file = con)
  close(con)
}

cat("..Exporting long table [R]\n")
saveRDS.gz(object = RES,
  file = "OTU_table_long.RData",
  threads = CPUTHREADS)

cat("..Exporting wide table [R]\n")
saveRDS.gz(object = REW,
  file = "OTU_table_wide.RData",
  threads = CPUTHREADS)

cat("..Exporting long table [tab-delimited]\n")
fwrite(x = RES, file = "OTU_table_long.txt.gz", sep = "\t", compress = "gzip")

cat("..Exporting wide table [tab-delimited]\n")
fwrite(x = REW, file = "OTU_table_wide.txt.gz", sep = "\t", compress = "gzip")


cat("\nExporting OTU sequences to FASTA\n")

cat("..Preparing sequences\n")

## Take sequences from the data (NB! there are a several different sequence per OTU)
# SQS <- unique(RES[, .(OTU) ])
# tmp_OTUs <- unique(TAB[ OTU %in% SQS$OTU & SeqID == OTU , .(OTU, Sequence) ])
# SQS <- merge(x = SQS, y = tmp_OTUs, by = "OTU", all.x = TRUE)
# rm(tmp_OTUs)
# 
# cat("...Preparing XStringSet object\n")
# SQF <- DNAStringSet(x = SQS$Sequence)
# names(SQF) <- SQS$OTU

## Take sequnces from the OTU file
cat("... Loading FASTA file\n")
SQS <- readDNAStringSet(filepath = OTUS, format="fasta")
cat("... Extracting sequence IDs\n")
names(SQS) <- tstrsplit(x = names(SQS), split = ";", keep = 1)[[1]]

if(any(duplicated(names(SQS)))){
  cat("WARNING: duplicated OTU names detected!\n")
}

cat("... Subsetting OTUs\n")
SQF <- SQS[ names(SQS) %in% unique(REW$OTU) ]

cat("....Total number of OTUs in input sequences: ", length(SQS), "\n")
cat("....Number of OTUs to export: ",                length(SQF), "\n")
cat("....Number of OTUs in the OTU table: ",         nrow(REW),   "\n")

cat("... Writing FASTA file\n")
writeXStringSet(x = SQF,
  filepath = "OTUs.fa.gz",
  compress = TRUE, format = "fasta", width = 9999)


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
