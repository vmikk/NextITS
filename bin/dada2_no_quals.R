#!/usr/bin/env Rscript

## Perform sequence denoising with DADA2

############################################## Parse input parameters

## Check time
start_time <- Sys.time()


cat("\nParsing input options and arguments...\n")

suppressPackageStartupMessages(require(optparse))

## Parse arguments
option_list <- list(
  make_option(c("-i", "--input"),            action="store", default=NA,    type='character', help=""),
  make_option(c("-n", "--nbases"),           action="store", default=1e6,   type='double', help=""),
  make_option(c("-b", "--bandsize"),         action="store", default=16,    type='double', help=""),
  make_option(c("-s", "--detectsingletons"), action="store", default=TRUE,  type='logical', help=""),
  make_option(c("-A", "--omegaA"),           action="store", default=1e-20, type='double', help=""),
  make_option(c("-C", "--omegaC"),           action="store", default=1e-40, type='double', help=""),
  make_option(c("-P", "--omegaP"),           action="store", default=1e-4,  type='double', help=""),
  make_option(c("-x", "--maxconsist"),       action="store", default=10,    type='integer', help=""),
  make_option("--match",                     action="store", default=4,     type='double', help=""),
  make_option("--mismatch",                  action="store", default=-5,    type='double', help=""),
  make_option("--gappenalty",                action="store", default=-8,    type='double', help=""),
  make_option("--hpgap",                     action="store", default=NULL,  type='double', help=""),
  make_option(c("-t", "--threads"),          action="store", default=4L,    type='integer', help="Number of CPU threads for arrow, default 4")
)
opt <- parse_args(OptionParser(option_list=option_list))

# Validation of the required argiments
if(is.na(opt$input)){
  cat("Input file is not specified: ....\n", file=stderr())
  stop()
}


## Function to convert text "NA"s to NA
# to_na <- function(x){ 
#   if(x %in% c("NA", "null", "Null")){ x <- NA }
#   return(x)
# }

## Assign variables
INPUT              <- opt$input
NBASES             <- opt$nbases
BAND_SIZE          <- opt$bandsize
DETECT_SINGLETONS  <- opt$detectsingletons
OMEGA_A            <- opt$omegaA
OMEGA_C            <- opt$omegaC
OMEGA_P            <- opt$omegaP
MAX_CONSIST        <- opt$maxconsist
MATCH              <- opt$match
MISMATCH           <- opt$mismatch
GAP_PENALTY        <- opt$gappenalty
HOMOPOLYMER_GAP_PENALTY <- opt$hpgap     # PacBio CCS does not make homopolymer errors at a higher rate than normal indels -> NULL
CPUTHREADS         <- opt$threads


## Log assigned variables
cat("\nParameters specified:\n")
cat(paste("Input file: " ,     INPUT, "\n", sep=""))
cat(paste("Number of bases to use for error rate learning: ", NBASES, "\n", sep = ""))
cat(paste("Band size for the Needleman-Wunsch alignment: ", BAND_SIZE, "\n", sep = ""))
cat(paste("Singleton detection: ", DETECT_SINGLETONS, "\n", sep = ""))
cat(paste("OMEGA_A: ", OMEGA_A, "\n", sep = ""))
cat(paste("OMEGA_C: ", OMEGA_C, "\n", sep = ""))
cat(paste("OMEGA_P: ", OMEGA_P, "\n", sep = ""))
cat(paste("Number of iterations of the self-consistency loop: ", MAX_CONSIST, "\n", sep = ""))
cat(paste("Alignment for matches: ", MATCH, "\n", sep = ""))
cat(paste("Alignment for mismatches: ", MISMATCH, "\n", sep = ""))
cat(paste("Gap penalty: ", GAP_PENALTY, "\n", sep = ""))
cat(paste("Homopolymer gap penalty: ", HOMOPOLYMER_GAP_PENALTY, "\n", sep = ""))
cat(paste("Number of CPU threads to use: ",  CPUTHREADS,    "\n", sep=""))

cat("\n")



############################################## Load packages

cat("Loading R packages...\n")

load_pckg <- function(pkg = "data.table"){
  suppressPackageStartupMessages( library(package = pkg, character.only = TRUE) )
  cat(paste(pkg, packageVersion(pkg), "\n"))
}

load_pckg("Biostrings")
load_pckg("ShortRead")
load_pckg("data.table")
load_pckg("dada2")

cat("\n")

## Set CPU thread number
cat("Setting number of CPU threads to: ", CPUTHREADS, "\n")
setDTthreads(threads = CPUTHREADS)  # for data.table

## Set seed
set.seed(111)
