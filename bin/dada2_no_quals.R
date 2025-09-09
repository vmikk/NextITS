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

## Set DADA options
cat("Setting DADA2 options\n")
setDadaOpt(
  BAND_SIZE         = BAND_SIZE,           # dada2 default, 16
  DETECT_SINGLETONS = DETECT_SINGLETONS,   # dada2 default, FALSE
  OMEGA_A           = OMEGA_A,             # dada2 default, 1e-40
  OMEGA_C           = OMEGA_C,             # dada2 default, 1e-40
  OMEGA_P           = OMEGA_P,             # dada2 default, 1e-4
  MAX_CONSIST       = MAX_CONSIST,         # dada2 default, 10
  GAP_PENALTY       = GAP_PENALTY,         # dada2 default, -8
  MATCH             = MATCH,               # dada2 default, 4
  MISMATCH          = MISMATCH,            # dada2 default, -5
  HOMOPOLYMER_GAP_PENALTY = HOMOPOLYMER_GAP_PENALTY  # PacBio CCS does not make homopolymer errors at a higher rate than normal indels
  )

## Get DADA options
# getDadaOpt()

############################################## Workflow


## Load FASTQ file
fq <- readFastq(dirPath = INPUT, qualityType = "FastqQuality")

## Extract sequence headers
sq <- as.data.table(fq@id)
setnames(x = sq, new = "SeqName")
sq[ , c("SeqID", "Abundance") := tstrsplit(x = SeqName, split = ";size=", keep = 1:2) ]
sq[ , Abundance := as.numeric(Abundance) ]
sq[ , Sequence := as.character(sread(fq))]

## Extract sequence qualities
seq_quals <- as(quality(fq), "matrix")
# dada2:::qtables2(fq)

## Summary stats
num_seqs  <- nrow(sq)
num_singl <- nrow(sq[ Abundance < 2 ])
num_reads <- sum(sq$Abundance, na.rm = TRUE)
perc_nonsingleton <- (num_seqs - num_singl) / num_seqs * 100


## Manually create a derep-class object
##   See also https://github.com/benjjneb/dada2/blob/004ce26909268e1318a2f68e0ea26807412c7a2d/R/sequenceIO.R#L240-L242
#             https://github.com/benjjneb/dada2/blob/004ce26909268e1318a2f68e0ea26807412c7a2d/R/sequenceIO.R#L45

## Prepare derep-class object
uniques <- sq$Abundance
names(uniques) <- as.character(sread(fq))   # names = full amplicon sequence
rownames(seq_quals) <- names(uniques)

derep <- list(
  uniques = uniques,
  quals   = seq_quals,
  map     = NULL,
  SeqID   = sq$SeqID         # add allso sequence IDs
  )

derep <- as(derep, "derep")


## Estimate error rates for each type of transition while ignoring quality scores
errors <- try(
  learnErrors(
    fls    = derep,
    nbases = NBASES,
    errorEstimationFunction = noqualErrfun,
    qualityType = "FastqQuality",
    verbose     = 1,
    multithread = CPUTHREADS
    )
  )


## Export results
saveRDS(object = errors,
  file = "DADA2_ErrorRates_noqualErrfun.RData",
  compress = "xz")


## Plot observed and estimated error rates
# plotErrors(errors)


## Run sample inference with DADA2
dadares <- dada(
  derep = derep,
  err   = errors,
  errorEstimationFunction = noqualErrfun,
  selfConsist = FALSE,
  verbose     = 1,
  multithread = CPUTHREADS)
saveRDS(object = dadares,
  file = "DADA2_InferedSeqs_noqualErrfun.RData",
  compress = "xz")


res <- data.table(
  Sequence  = dadares$sequence,
  Abundance = dadares$denoised)

## Add sequence IDs
res[ , SeqNumID := .I ]
res <- merge(
  x = res,
  y = sq[, .(SeqID, Sequence)],
  by = "Sequence", all.x = TRUE)

## Sort by abundance
setorder(res, -Abundance, SeqID, na.last = TRUE)

## Export denoised sequences
ASVS <- DNAStringSet(x = res$Sequence)
names(ASVS) <- paste0(res$SeqID, ";size=", res$Abundance)

writeXStringSet(
  x = ASVS,
  filepath = "DADA2_denoised.fa.gz",
  compress = TRUE,
  format   = "fasta",
  width    = 20000)


UC <- data.table(
  DerepSeqID  = derep$SeqID,
  SeqNumID    = dadares$map,
  Abundance   = derep$uniques)

UC <- merge(
  x = UC,
  y = res[ , .(SeqNumID, SeqID) ],
  by = "SeqNumID", all.x = TRUE)

setorder(UC, SeqNumID, na.last = TRUE)
setnames(x = UC, old = "SeqID", new = "ASV")

## Export pre-UC file
saveRDS(
  object = UC,
  file = "DADA2_UC.RData",
  compress = "xz")


