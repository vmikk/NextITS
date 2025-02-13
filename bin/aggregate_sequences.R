#!/usr/bin/env Rscript

## Script to aggregate sequences from multiple runs into a single file (for dereplication and subsequent clustering)
## Also, performs removal of de novo chimeras with high scores (with option to recover sequences that occurred in multiple runs)


## Function to load packages
load_pckg <- function(pkg = "data.table"){
    suppressPackageStartupMessages( library(package = pkg, character.only = TRUE) )
    cat(".. ", paste(pkg, packageVersion(pkg), "\n"))
}

cat("Loading packages:\n")

load_pckg("optparse")
load_pckg("data.table")
load_pckg("Biostrings")
load_pckg("plyr")
load_pckg("arrow")
# load_pckg("dplyr")


cat("Parsing input options and arguments...\n")

option_list <- list(
  make_option("--seqtabs", action="store", default=NA,  type='character', help = "Direcotry containing long tables with quality-filtered sequences (Parquet format)"),
  make_option("--maxchim", action="store", default=0.6, type='numeric',   help = "Maximum de novo chimera score to remove"),
  make_option("--recoverdenovo", action="store", default=TRUE, type='logical', help="Recover de-novo chimeras (logical)"),
  make_option("--output",  action="store", default="Seqs", type='character', help = "Output prefix"),
  make_option("--threads", action="store", default=4,   type='integer',   help = "Number of CPU threads to use")
)

opt <- parse_args(OptionParser(option_list=option_list))

## Function to convert text "NA"s to NA
to_na <- function(x){
  if(x %in% c("NA", "null", "Null")){ x <- NA }
  return(x)
}

## Replaces "null"s from Nextflow with NA
opt <- lapply(X = opt, FUN = to_na)


## Validation of the required arguments
if (is.na(opt$seqtabs)) {
  stop("Input directory with quality-filtered sequences is not specified\n")
}

## Assign variables
SEQTABS      <- opt$seqtabs
MAXCHIM      <- opt$maxchim
RECOV_DENOVO <- opt$recoverdenovo
OUTPUT       <- opt$output
CPUTHREADS   <- as.numeric( opt$threads )

## Log assigned variables
cat(paste("Path to sequence tables: ",   SEQTABS,      "\n", sep=""))
cat(paste("Max de novo chimera score: ", MAXCHIM,      "\n", sep=""))
cat(paste("De novo chimera recovery: ",  RECOV_DENOVO, "\n", sep=""))
cat(paste("Output prefix: ",             OUTPUT,       "\n", sep=""))
cat(paste("CPU threads: ",               CPUTHREADS,   "\n", sep=""))

cat("\n")

## Set CPU thread number
cat("\nSetting number of CPU threads to: ", CPUTHREADS, "\n")
setDTthreads(threads = CPUTHREADS)     # for data.table
set_cpu_count(CPUTHREADS)              # for arrow


######################################
###################################### Process the data
######################################

## Load sequence tables
cat("..Looking for sequence tables\n")
TABS <- list.files(path = SEQTABS, pattern = ".parquet", full.names = TRUE, recursive = TRUE)
cat("... Tables found: ", length(TABS), "\n")

cat("..Loading sequence tables\n")
TAB <- alply(.data = TABS, .margins = 1, .fun = function(x){
  res <- arrow::read_parquet(x)
  setDT(res)
  return(res)
})
TAB <- rbindlist(TAB, use.names = TRUE, fill = TRUE)
cat("... Total number of records: ",             nrow(TAB), "\n")
cat("... Total number unique sequences: ",       length(unique(TAB$Sequence)), "\n")
cat("... Total number unique samples (files): ", length(unique(TAB$SampleID)), "\n")


## Filter sequences by chimeric score (MAXCHIM)
if(!is.na(MAXCHIM)){

  cat("..Filtering data by max de novo chimera score\n")
  nrecs <- nrow(TAB)
  nabun <- sum(TAB$Abundance, na.rm = TRUE)
  
  ## If no chimera recovery is required
  if(RECOV_DENOVO == FALSE){

    TAB <- TAB[ DeNovo_Chimera_Score < MAXCHIM | is.na(DeNovo_Chimera_Score) ]
  
  ## If we need to recover chimeras
  } else {

    ## Find putative chimeras
    CHIMERAS <- TAB[ DeNovo_Chimera_Score >= MAXCHIM, .(SeqID___SampleID, DeNovo_Chimera_Score, Sequence, Abundance) ]
    NONCHIMERAS <- TAB[ ! SeqID___SampleID %in% CHIMERAS$SeqID___SampleID ]

    ## Recover false-positives
    chim_seqs    <- unique(CHIMERAS$Sequence)
    nonchim_seqs <- unique(NONCHIMERAS$Sequence)
    fp_chims <- chim_seqs %in% nonchim_seqs
    if(any(fp_chims)){
      cat(".... Probably there are a few false-positive chimeras\n")
      cat(".... Recovering ", sum(fp_chims), "sequences\n")
      fp_seqs <- chim_seqs[ fp_chims ]
      CHIMERAS <- CHIMERAS[ ! Sequence %in% fp_seqs ]
      rm(fp_seqs)
    }

    TAB <- TAB[ ! Sequence %in% CHIMERAS$Sequence ]
    rm(CHIMERAS, NONCHIMERAS)
  
  } # end of chimera recovery

  ## Data summary after filtering
  nrecs_delta <- nrecs - nrow(TAB)
  nabun_delta <- nabun - sum(TAB$Abundance, na.rm = TRUE)

  cat("... Records removed: ", nrecs_delta, " (", round(nrecs_delta/nrecs * 100, 1), "%)\n")
  cat("... Reads removed: ", nabun_delta, " (", round(nabun_delta/nabun * 100, 1), "%)\n")

  rm(nrecs_delta, nabun_delta)

} # end of MAXCHIM filtering


cat("..Sorting table by abundance, quality score\n")
setorder(x = TAB, -Abundance, -PhredScore, SeqID)

cat("..Preparing FASTA file\n")

SQF <- DNAStringSet(x = TAB$Sequence)
names(SQF) <- paste0(TAB$SeqID, ";size=", TAB$Abundance)   # , ";sample=", TAB$SampleID, ";"

## Export FASTA
cat("..Exporting FASTA file with filtered sequences\n")

writeXStringSet(
  x = SQF,
  filepath = paste0(OUTPUT, ".fa.gz"),
  compress = TRUE, format = "fasta", width = 9999)

## Export FASTA
cat("..Exporting filtered table\n")

write_parquet(
  x    = TAB,
  sink = paste0(OUTPUT, ".parquet"),
  compression       = "zstd",
  compression_level = 10,
  use_dictionary    = TRUE)

cat("All done.\n")
