#!/usr/bin/env Rscript

## Summarise number of reads per process

# read_count_summary.R \
#   --raw          Counts_1.RawData.txt \
#   --qc           Counts_2.QC.txt \
#   --demuxed      Counts_3.Demux.txt \
#   --primer       Counts_4.PrimerCheck.txt \
#   --primermulti  Counts_4.PrimerMultiArtifacts.txt \
#   --itsx         Counts_5.ITSx_or_PrimTrim.txt \
#   --chimrefn     Counts_6.ChimRef_reads.txt \
#   --chimrefu     Counts_6.ChimRef_uniqs.txt \
#   --chimdenovo   Counts_7.ChimDenov.txt \
#   --chimrecovn   Counts_8.ChimRecov_reads.txt \
#   --chimrecovu   Counts_8.ChimRecov_uniqs.txt \
#   --tj           TagJump_OTUs.RData \
#   --seqtab       Seqs.RData \
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
  make_option("--primermulti",action="store", default=NA, type='character', help="Counts of multi-primer artifacts"),
  make_option("--itsx",       action="store", default=NA, type='character', help="Read counts after ITSx or primer removal"),
  make_option("--chimrefn",   action="store", default=NA, type='character', help="Number of reads for reference-based chimeras"),
  make_option("--chimrefu",   action="store", default=NA, type='character', help="Number of unique sequences detected as reference-based chimeras"),
  make_option("--chimdenovo", action="store", default=NA, type='character', help="Number of de novo chimeras"),
  make_option("--chimrecovn", action="store", default=NA, type='character', help="Number of resued reads for de novo chimeras (false positives)"),
  make_option("--chimrecovu", action="store", default=NA, type='character', help="Number of resued unique sequences detected as de novo chimeras (false positives)"),
  make_option("--tj",         action="store", default=NA, type='character', help="Tag jump removal data"),
  make_option("--seqtab",     action="store", default=NA, type='character', help="Final seq table"),
  make_option(c("-t", "--threads"), action="store", default=4L, type='integer', help="Number of CPU threads for arrow, default 4")
)
opt <- parse_args(OptionParser(option_list=option_list))


## Validation of the required argiments
# if(is.na(opt$raw)){
#   cat("Input file is not specified: ....\n", file=stderr())
#   stop()
# }



## Function to convert text "NA"s to NA
to_na <- function(x){ 
  if(x %in% c("NA", "null", "Null")){ x <- NA }
  return(x)
}

## Assign variables
RAW         <- opt$raw
QC          <- opt$qc
DEMUXED     <- opt$demuxed
PRIMER      <- opt$primer
PRIMERMULTI <- opt$primermulti
ITSX        <- opt$itsx
CHIMREFN    <- opt$chimrefn
CHIMREFU    <- opt$chimrefu
CHIMDENOVO  <- opt$chimdenovo
CHIMRECOVN  <- opt$chimrecovn
CHIMRECOVU  <- opt$chimrecovu
TJ          <- opt$tj
SEQTAB      <- opt$seqtab
CPUTHREADS  <- as.numeric( opt$threads )

## Log assigned variables
cat(paste("Counts - RawData: " ,     RAW, "\n", sep=""))
cat(paste("Counts - QC: " ,          QC, "\n", sep=""))
cat(paste("Counts - Demux: " ,       DEMUXED, "\n", sep=""))
cat(paste("Counts - PrimerCheck: " , PRIMER, "\n", sep=""))
cat(paste("Counts - Primer Multi Artifacts: " ,                        PRIMERMULTI, "\n", sep=""))
cat(paste("Counts - ITSx or Primer Trim: " ,                           ITSX, "\n", sep=""))
cat(paste("Counts - Chimera Ref-based, reads: " ,                      CHIMREFN, "\n", sep=""))
cat(paste("Counts - Chimera Ref-based, unique sequences: " ,           CHIMREFU, "\n", sep=""))
cat(paste("Counts - Chimera de novo: " ,                               CHIMDENOVO, "\n", sep=""))
cat(paste("Counts - Chimera Ref-based recoverd, reads: " ,             CHIMRECOVN, "\n", sep=""))
cat(paste("Counts - Chimera Ref-based recoverd, unique sequences: " ,  CHIMRECOVU, "\n", sep=""))
cat(paste("Tag-jump data: " ,                TJ, "\n", sep=""))
cat(paste("Final sequence table: " ,         SEQTAB, "\n", sep=""))
cat(paste("Number of CPU threads to use: ",  CPUTHREADS,    "\n", sep=""))

cat("\n")


############################################## data for debuging

# RAW         <- "Counts_1.RawData.txt"
# QC          <- "Counts_2.QC.txt"
# DEMUXED     <- "Counts_3.Demux.txt"
# PRIMER      <- "Counts_4.PrimerCheck.txt"
# PRIMERMULTI <- "Counts_4.PrimerMultiArtifacts.txt"
# ITSX        <- "Counts_5.ITSx_or_PrimTrim.txt"
# CHIMREFN    <- "Counts_6.ChimRef_reads.txt"
# CHIMREFU    <- "Counts_6.ChimRef_uniqs.txt"
# CHIMDENOVO  <- "Counts_7.ChimDenov.txt"
# CHIMRECOVN  <- "Counts_8.ChimRecov_reads.txt"
# CHIMRECOVU  <- "Counts_8.ChimRecov_uniqs.txt"
# TJ          <- "Seqs.RData"
# SEQTAB      <- "TagJump_OTUs.RData"
# CPUTHREADS  <- 6


############################################## Load packages and data

cat("Loading R packages...\n")

load_pckg <- function(pkg = "data.table"){
  suppressPackageStartupMessages( library(package = pkg, character.only = TRUE) )
  cat(paste(pkg, packageVersion(pkg), "\n"))
}

load_pckg("data.table")
load_pckg("plyr")
load_pckg("metagMisc")
load_pckg("openxlsx")

cat("\n")

## Set CPU thread number
cat("Setting number of CPU threads to: ", CPUTHREADS, "\n")
setDTthreads(threads = CPUTHREADS)  # for data.table


######################################
###################################### Load the data
######################################

cat("\nLoading input data\n")


#### Per-dataset stats

## Load ASV table
cat("..Loading raw counts\n")
RAW <- fread(RAW)

cat("..Loading QC counts\n")
QC <- fread(QC)

#### Per-sample stats

SEQKITCOUNTS <- list()
CUSTOMCOUNTS <- list()

cat("..Loading demux counts\n")
SEQKITCOUNTS$DEMUXED <- fread(DEMUXED)

cat("..Loading primer-checked data counts\n")
SEQKITCOUNTS$PRIMER      <- fread(PRIMER)
SEQKITCOUNTS$PRIMERMULTI <- fread(PRIMERMULTI)

cat("..Loading ITSx or primer trim counts\n")
CUSTOMCOUNTS$ITSX <- fread(ITSX)

cat("..Loading ref-based chimera counts\n")
CUSTOMCOUNTS$CHIMREFN <- fread(CHIMREFN)
SEQKITCOUNTS$CHIMREFU <- fread(CHIMREFU)

cat("..Loading de novo chimera counts\n")
CHIMDENOVO <- fread(CHIMDENOVO)                 # incorporate to the main table

cat("..Loading rescued ref-based chimera counts\n")
CUSTOMCOUNTS$CHIMRECOVN <- fread(CHIMRECOVN)
SEQKITCOUNTS$CHIMRECOVU <- fread(CHIMRECOVU)

# cat("..Loading tag-jump filtration data\n")
# TJ

# cat("..Loading sequence table\n")
# SEQTAB


## Remove NULL-files
null_cust <- laply(.data = CUSTOMCOUNTS, .fun = nrow)
null_seqk <- laply(.data = SEQKITCOUNTS, .fun = nrow)

if(any(null_cust == 0)){
  cat("Some files with counts are missing:\n")
  to_rm <- which(null_cust == 0)
  cat(".. ", paste(names(CUSTOMCOUNTS)[ to_rm ], collapse = ", "), "\n")
  CUSTOMCOUNTS[ to_rm ] <- NULL
  rm(to_rm)
}

if(any(null_seqk == 0)){
  cat("Some files with counts are missing:\n")
  to_rm <- which(null_seqk == 0)
  cat(".. ", paste(names(SEQKITCOUNTS)[ to_rm ], collapse = ", "), "\n")
  SEQKITCOUNTS[ to_rm ] <- NULL
  rm(to_rm)
}


## Process seqkit counts
seqkit_process <- function(x){
  if(nrow(x) > 0){

    ## Remove reudndant columns
    x <- x[ , .(file, num_seqs) ]

    ## Remove file extensions
    x[ , file := sub(pattern = ".fastq.gz",        replacement = "", x = file) ]
    x[ , file := sub(pattern = ".fq.gz",           replacement = "", x = file) ]
    x[ , file := sub(pattern = ".fa.gz",           replacement = "", x = file) ]
    x[ , file := sub(pattern = ".full.fasta",      replacement = "", x = file) ]
    x[ , file := sub(pattern = ".ITS1.fasta.gz",   replacement = "", x = file) ]
    x[ , file := sub(pattern = ".ITS2.fasta.gz",   replacement = "", x = file) ]
    x[ , file := sub(pattern = "_PrimerChecked$",  replacement = "", x = file) ]
    x[ , file := sub(pattern = "_Mutiprimer$",     replacement = "", x = file) ]
    x[ , file := sub(pattern = "_Chimera$",        replacement = "", x = file) ]
    x[ , file := sub(pattern = "_RescuedChimera$", replacement = "", x = file) ]
    x[ , file := sub(pattern = "^Rescued_Chimeric_sequences.part_", replacement = "", x = file) ]

  }
  return(x)
}

## Process custom counts
custom_process <- function(x){
  if(nrow(x) > 0){

    ## There should be just two columns - `SampleID` & `NumReads`

    ## Rename "SampleID" into "file"
    setnames(x = x, old = "SampleID", new = "file")

    ## Remove file extensions
    x[ , file := sub(pattern = ".full.fasta",         replacement = "", x = file) ]
    x[ , file := sub(pattern = "_Chimera.fa",         replacement = "", x = file) ]
    x[ , file := sub(pattern = "_RescuedChimera.fa$", replacement = "", x = file) ]
    x[ , file := sub(pattern = "^Rescued_Chimeric_sequences.part_", replacement = "", x = file) ]

  }
  return(x)
}


cat("Processing data\n")
SEQKITCOUNTS <- llply(.data = SEQKITCOUNTS, .fun = seqkit_process)
CUSTOMCOUNTS <- llply(.data = CUSTOMCOUNTS, .fun = custom_process)


## Rename columns
if(!is.null(SEQKITCOUNTS$DEMUXED)){
setnames(x = SEQKITCOUNTS$DEMUXED, old = "num_seqs", new = "Demultiplexed_Reads", skip_absent = TRUE)
}
if(!is.null(SEQKITCOUNTS$PRIMER)){
setnames(x = SEQKITCOUNTS$PRIMER,  old = "num_seqs", new = "PrimerChecked_Reads", skip_absent = TRUE)
}
if(!is.null(SEQKITCOUNTS$PRIMERMULTI)){
setnames(x = SEQKITCOUNTS$PRIMERMULTI, old = "num_seqs", new = "MultiprimerArtifacts_Reads", skip_absent = TRUE)
}
if(!is.null(SEQKITCOUNTS$CHIMREFU)){
setnames(x = SEQKITCOUNTS$CHIMREFU, old = "num_seqs", new = "ReferenceBasedChimera_NumUniqSequences", skip_absent = TRUE)
}
if(!is.null(SEQKITCOUNTS$CHIMRECOVU)){
setnames(x = SEQKITCOUNTS$CHIMRECOVU, old = "num_seqs", new = "Recovered_ReferenceBasedChimea_NumUniqSequences", skip_absent = TRUE)
}

if(!is.null(CUSTOMCOUNTS$ITSX)){
setnames(x = CUSTOMCOUNTS$ITSX,       old = "NumReads", new = "ITSx_Extracted_Reads", skip_absent = TRUE)
}
if(!is.null(CUSTOMCOUNTS$CHIMREFN)){
setnames(x = CUSTOMCOUNTS$CHIMREFN,   old = "NumReads", new = "ReferenceBasedChimera_Reads", skip_absent = TRUE)
}
if(!is.null(CUSTOMCOUNTS$CHIMRECOVN)){
setnames(x = CUSTOMCOUNTS$CHIMRECOVN, old = "NumReads", new = "Recovered_ReferenceBasedChimea_Reads", skip_absent = TRUE)
}

## Merge seqkit and custom counts into a single list
cat("Pooling per-sample counts\n")
COUNTS <- c(SEQKITCOUNTS, CUSTOMCOUNTS)

## Pool per-file estimates
merge_dt <- function(x,y){ merge(x, y, by = "file", all = TRUE) }
PER_SAMPLE_COUNTS_merged <- Reduce(f = merge_dt, x = COUNTS)


### ... update



## Export summary stats
write.xlsx(list(
  "per_sample" = SEQKITCOUNTS$DEMUXED
  # , "per_run"
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
