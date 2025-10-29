#!/usr/bin/env Rscript

## Summarise number of reads per process

# read_count_summary.R \
#   --raw          Counts_1.RawData.txt \
#   --qc           Counts_2.QC.txt \
#   --demuxed      Counts_3.Demux.txt \
#   --primer       Counts_4.PrimerCheck.txt \
#   --primerartef  Counts_4.PrimerArtefacts.txt \
#   --itsx         Counts_5.ITSx_or_PrimTrim.txt \
#   --homopolymer  Counts_5.Homopolymers.txt \
#   --chimrefn     Counts_6.ChimRef_reads.txt \
#   --chimrefu     Counts_6.ChimRef_uniqs.txt \
#   --chimdenovo   Counts_7.ChimDenov.txt \
#   --chimrecovn   Counts_8.ChimRecov_reads.txt \
#   --chimrecovu   Counts_8.ChimRecov_uniqs.txt \
#   --tj           TagJump_scores.qs \
#   --seqtab       Seqs.parquet \
#   --maxchim     0.6 \
#   --threads      4


############################################## Parse input parameters

## Check time
start_time <- Sys.time()


cat("\nParsing input options and arguments...\n")

suppressPackageStartupMessages(require(optparse))

## Parse arguments
option_list <- list(

  make_option("--raw",        action="store", default=NA, type='character', help="Raw read counts"),
  make_option("--qc",         action="store", default=NA, type='character', help="Counts of reads passed QC"),
  make_option("--demuxed",    action="store", default=NA, type='character', help="Counts of demultiplexed reads"),
  make_option("--primer",     action="store", default=NA, type='character', help="Counts of reads with both primers detected"),
  make_option("--primerartef",action="store", default=NA, type='character', help="Counts of primer artefacts"),
  make_option("--itsx",       action="store", default=NA, type='character', help="Read counts after ITSx or primer removal"),
  make_option("--homopolymer",action="store", default=NA, type='character', help="Homopolymer correction results"),
  make_option("--chimrefn",   action="store", default=NA, type='character', help="Number of reads for reference-based chimeras"),
  make_option("--chimrefu",   action="store", default=NA, type='character', help="Number of unique sequences detected as reference-based chimeras"),
  make_option("--chimdenovo", action="store", default=NA, type='character', help="Number of de novo chimeras"),
  make_option("--chimrecovn", action="store", default=NA, type='character', help="Number of resued reads for de novo chimeras (false positives)"),
  make_option("--chimrecovu", action="store", default=NA, type='character', help="Number of resued unique sequences detected as de novo chimeras (false positives)"),
  make_option("--tj",         action="store", default=NA, type='character', help="Tag jump removal data (serialized in qs format)"),
  make_option("--seqtab",     action="store", default=NA, type='character', help="Final seq table (Parquet format)"),
  make_option("--maxchim",    action="store", default=0.6, type='numeric',  help = "Maximum de novo chimera score to remove"),
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
PRIMERARTEF <- opt$primerartef
ITSX        <- opt$itsx
HOMOPOLY    <- opt$homopolymer
CHIMREFN    <- opt$chimrefn
CHIMREFU    <- opt$chimrefu
CHIMDENOVO  <- opt$chimdenovo
CHIMRECOVN  <- opt$chimrecovn
CHIMRECOVU  <- opt$chimrecovu
TJ          <- opt$tj
SEQTAB      <- opt$seqtab
MAXCHIM     <- opt$maxchim
CPUTHREADS  <- as.numeric( opt$threads )

## Log assigned variables
cat(paste("Counts - RawData: " ,     RAW, "\n", sep=""))
cat(paste("Counts - QC: " ,          QC, "\n", sep=""))
cat(paste("Counts - Demux: " ,       DEMUXED, "\n", sep=""))
cat(paste("Counts - PrimerCheck: " , PRIMER, "\n", sep=""))
cat(paste("Counts - Primer Artefacts: " ,                              PRIMERARTEF, "\n", sep=""))
cat(paste("Counts - ITSx or Primer Trim: " ,                           ITSX, "\n", sep=""))
cat(paste("Counts - Homopolymer correction results: " ,                HOMOPOLY, "\n", sep=""))
cat(paste("Counts - Chimera Ref-based, reads: " ,                      CHIMREFN, "\n", sep=""))
cat(paste("Counts - Chimera Ref-based, unique sequences: " ,           CHIMREFU, "\n", sep=""))
cat(paste("Counts - Chimera de novo: " ,                               CHIMDENOVO, "\n", sep=""))
cat(paste("Counts - Chimera Ref-based recoverd, reads: " ,             CHIMRECOVN, "\n", sep=""))
cat(paste("Counts - Chimera Ref-based recoverd, unique sequences: " ,  CHIMRECOVU, "\n", sep=""))
cat(paste("Tag-jump data: " ,                TJ, "\n", sep=""))
cat(paste("Final sequence table: " ,         SEQTAB,  "\n", sep=""))
cat(paste("Maximum de novo chimera score: ", MAXCHIM, "\n", sep=""))
cat(paste("Number of CPU threads to use: ",  CPUTHREADS,    "\n", sep=""))

cat("\n")


############################################## data for debuging

# RAW         <- "Counts_1.RawData.txt"
# QC          <- "Counts_2.QC.txt"
# DEMUXED     <- "Counts_3.Demux.txt"
# PRIMER      <- "Counts_4.PrimerCheck.txt"
# PRIMERARTEF <- "Counts_4.PrimerArtefacts.txt"
# ITSX        <- "Counts_5.ITSx_or_PrimTrim.txt"
# HOMOPOLY    <- "Counts_5.Homopolymers.txt"
# CHIMREFN    <- "Counts_6.ChimRef_reads.txt"
# CHIMREFU    <- "Counts_6.ChimRef_uniqs.txt"
# CHIMDENOVO  <- "Counts_7.ChimDenov.txt"
# CHIMRECOVN  <- "Counts_8.ChimRecov_reads.txt"
# CHIMRECOVU  <- "Counts_8.ChimRecov_uniqs.txt"
# TJ          <- "TagJump_scores.qs"
# SEQTAB      <- "Seqs.parquet"
# MAXCHIM     <- 0.6
# CPUTHREADS  <- 6


############################################## Load packages and data

cat("Loading R packages...\n")

load_pckg <- function(pkg = "data.table"){
  suppressPackageStartupMessages( library(package = pkg, character.only = TRUE) )
  cat(paste(pkg, packageVersion(pkg), "\n"))
}

load_pckg("data.table")
load_pckg("plyr")
load_pckg("arrow")
# load_pckg("dplyr")
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
SEQKITCOUNTS$PRIMERARTEF <- fread(PRIMERARTEF)

cat("..Loading ITSx or primer trim counts\n")
CUSTOMCOUNTS$ITSX <- fread(ITSX)

cat("..Loading homopolymer correction results\n")
HOMOPOLY_data <- fread(HOMOPOLY)

cat("..Loading ref-based chimera counts\n")
CUSTOMCOUNTS$CHIMREFN <- fread(CHIMREFN)
SEQKITCOUNTS$CHIMREFU <- fread(CHIMREFU)

cat("..Loading de novo chimera counts\n")
CHIMDENOVO <- fread(CHIMDENOVO)                 # incorporate to the main table

cat("..Loading rescued ref-based chimera counts\n")
CUSTOMCOUNTS$CHIMRECOVN <- fread(CHIMRECOVN)
SEQKITCOUNTS$CHIMRECOVU <- fread(CHIMRECOVU)

if(!is.na(TJ) && TJ != "no_tj" && file.exists(TJ)){
  cat("..Loading tag-jump filtration data\n")
  TJ <- qs::qread(TJ)
  tjdata <- TRUE
} else {
  cat("..No tag-jump filtration data found\n")
  tjdata <- FALSE
} 

cat("..Loading sequence table\n")
SEQTAB <- arrow::open_dataset(SEQTAB)


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
    x[ , file := sub(pattern = ".fastq.gz$",            replacement = "", x = file) ]
    x[ , file := sub(pattern = ".fq.gz$",               replacement = "", x = file) ]
    x[ , file := sub(pattern = ".fa.gz$",               replacement = "", x = file) ]
    x[ , file := sub(pattern = ".full.fasta$",          replacement = "", x = file) ]
    x[ , file := sub(pattern = ".ITS1.fasta.gz$",       replacement = "", x = file) ]
    x[ , file := sub(pattern = ".ITS2.fasta.gz$",       replacement = "", x = file) ]
    x[ , file := sub(pattern = "_PrimerChecked$",       replacement = "", x = file) ]
    x[ , file := sub(pattern = "_PrimerArtefacts$",     replacement = "", x = file) ]
    x[ , file := sub(pattern = "_Chimera$",             replacement = "", x = file) ]
    x[ , file := sub(pattern = "_RescuedChimera$",      replacement = "", x = file) ]
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
    x[ , file := sub(pattern = ".full.fasta$",          replacement = "", x = file) ]
    x[ , file := sub(pattern = "_ITS1_58S_ITS2.fasta$", replacement = "", x = file) ]
    x[ , file := sub(pattern = "_Chimera.fa$",          replacement = "", x = file) ]
    x[ , file := sub(pattern = "_RescuedChimera.fa$",   replacement = "", x = file) ]
    x[ , file := sub(pattern = "^Rescued_Chimeric_sequences.part_", replacement = "", x = file) ]

  }
  return(x)
}


cat("Processing data\n")
SEQKITCOUNTS <- llply(.data = SEQKITCOUNTS, .fun = seqkit_process)
CUSTOMCOUNTS <- llply(.data = CUSTOMCOUNTS, .fun = custom_process)

cat("Estimating homopolymer stats\n")
if(nrow(HOMOPOLY_data) > 0){
  HOMOPOLY_counts <- HOMOPOLY_data[ , .(
    N_UniqSequences_AfterITSx_or_PrimerTrimming = .N,
    N_UniqSequences_AfterHomopolymerCorrection = length(unique(Target))
    ),
    by = "SampleID" ]
} else {
  cat("..No homopolymer correction data found\n")
}

# HOMOPOLY_counts[, Num_HomopolymerCorrectedSequences := 
#   N_UniqSequences_AfterITSx_or_PrimerTrimming - N_UniqSequences_AfterHomopolymerCorrection ]


## Rename columns
if(!is.null(SEQKITCOUNTS$DEMUXED)){
setnames(x = SEQKITCOUNTS$DEMUXED, old = "num_seqs", new = "Demultiplexed_Reads", skip_absent = TRUE)
}
if(!is.null(SEQKITCOUNTS$PRIMER)){
setnames(x = SEQKITCOUNTS$PRIMER,  old = "num_seqs", new = "PrimerChecked_Reads", skip_absent = TRUE)
}
if(!is.null(SEQKITCOUNTS$PRIMERARTEF)){
setnames(x = SEQKITCOUNTS$PRIMERARTEF, old = "num_seqs", new = "PrimerArtefacts_Reads", skip_absent = TRUE)
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

## If there are no primer artefacts
if(is.null(SEQKITCOUNTS$PRIMERARTEF)){
  PER_SAMPLE_COUNTS_merged[ , PrimerArtefacts_Reads := 0 ]
}

## Estimate percentage of primer artefacts
cat("Estimating percentage of primer artefacts\n")
PER_SAMPLE_COUNTS_merged[ , 
 PrimerArtefacts_Percent := round(
  PrimerArtefacts_Reads / (PrimerChecked_Reads + PrimerArtefacts_Reads) * 100,
  2)
 ]

## Estimate tag-jump stats
if(tjdata == TRUE){
  cat("Estimating tag-jump removal yields\n")
  TJ_stats <- TJ[ TagJump == TRUE, .(
      TagJump_Events = .N,
      TagJump_Reads  = sum(Abundance, na.rm = TRUE)),
    by = "SampleID" ]

  if(nrow(TJ_stats) > 0){
    PER_SAMPLE_COUNTS_merged <- merge(
      x = PER_SAMPLE_COUNTS_merged,
      y = TJ_stats,
      by.x = "file", by.y = "SampleID", all.x = TRUE)
  } else {
    PER_SAMPLE_COUNTS_merged[ , TagJump_Events := 0 ]
    PER_SAMPLE_COUNTS_merged[ , TagJump_Reads  := 0 ]
  }
}

## Add homopolymer stats
if(nrow(HOMOPOLY_data) > 0){
  cat("Adding homopolymer stats\n")
  PER_SAMPLE_COUNTS_merged <- merge(
    x = PER_SAMPLE_COUNTS_merged, 
    y = HOMOPOLY_counts,
    by.x = "file", by.y = "SampleID", all.x = TRUE)
}

## Add de novo chimera stats
cat("Adding de novo chimera stats\n")
denovo_stats <- SEQTAB |>
  dplyr::filter(DeNovo_Chimera_Score >= MAXCHIM) |>
  dplyr::group_by(SampleID) |>
  dplyr::summarize(
    DeNovoChimeras_NumReads    = sum(Abundance, na.rm = TRUE),
    DeNovoChimeras_NumUniqSeqs = n()) |>
  dplyr::collect() |>
  setDT()

PER_SAMPLE_COUNTS_merged <- merge(
  x = PER_SAMPLE_COUNTS_merged,
  y = denovo_stats,
  by.x = "file", by.y = "SampleID", all.x = TRUE)



## Number of reads and unique sequences in the sequence table (per sample)
## Excluding putative de novo chimeras (with score >= MAXCHIM)
cat("Adding sequence table stats\n")
seqtab_stats <- SEQTAB |>
  dplyr::filter(DeNovo_Chimera_Score < MAXCHIM | is.na(DeNovo_Chimera_Score) ) |>
  dplyr::group_by(SampleID) |>
  dplyr::summarize(
    SeqTable_NumReads    = sum(Abundance, na.rm = TRUE),
    SeqTable_NumUniqSeqs = n()) |>
  dplyr::collect() |>
  setDT()

PER_SAMPLE_COUNTS_merged <- merge(
  x = PER_SAMPLE_COUNTS_merged,
  y = seqtab_stats,
  by.x = "file", by.y = "SampleID", all.x = TRUE)


## Replace NAs with zeros
cat("Replacing NAs with zero\n")
for (j in seq_len(ncol(PER_SAMPLE_COUNTS_merged))){
  set(PER_SAMPLE_COUNTS_merged, which(is.na(PER_SAMPLE_COUNTS_merged[[j]])), j, 0)
}
rm(j)

## Estimate percentage of reads retained (starting from demultiplexed reads)
cat("Estimating percentage of reads retained\n")
PER_SAMPLE_COUNTS_merged[ , Percentage_Reads_Retained := round( SeqTable_NumReads / Demultiplexed_Reads * 100, 2) ]

## Estimate percentage of reads retained after ITSx
if("ITSx_Extracted_Reads" %in% colnames(PER_SAMPLE_COUNTS_merged)){
  PER_SAMPLE_COUNTS_merged[ , ITSx_Yield_Percent := round( ITSx_Extracted_Reads / PrimerChecked_Reads * 100, 2) ]
}

### TODO:
# .. estimate percentages
# .. add per-run positive / negative counts (based on default sample names)


## Reorder columns
setcolorder(PER_SAMPLE_COUNTS_merged,
  skip_absent = TRUE,
  neworder = c(
    "file", "Demultiplexed_Reads",
    "PrimerChecked_Reads", "PrimerArtefacts_Reads", "PrimerArtefacts_Percent", 
    "ReferenceBasedChimera_Reads", "ReferenceBasedChimera_NumUniqSequences", 
    "Recovered_ReferenceBasedChimea_Reads", "Recovered_ReferenceBasedChimea_NumUniqSequences", 
    "DeNovoChimeras_NumReads", "DeNovoChimeras_NumUniqSeqs", 
    "ITSx_Extracted_Reads", "ITSx_Yield_Percent",
    "N_UniqSequences_AfterITSx_or_PrimerTrimming", 
    "N_UniqSequences_AfterHomopolymerCorrection", 
    # "Num_HomopolymerCorrectedSequences",
    "TagJump_Reads", "TagJump_Events", 
    "SeqTable_NumReads", "SeqTable_NumUniqSeqs", 
    "Percentage_Reads_Retained"))




## Prepare per-run stats
cat("Preparing per-run stats\n")
PER_RUN_COUNTS_merged <- data.table(
  Total_Number_Of_Reads = sum(RAW$num_seqs, na.rm = TRUE),
  Reads_Demultiplexed   = sum(PER_SAMPLE_COUNTS_merged$Demultiplexed_Reads, na.rm = TRUE),
  Reads_Passed_QC       = sum(QC$num_seqs, na.rm = TRUE),
  Reads_PrimerChecked   = sum(PER_SAMPLE_COUNTS_merged$PrimerChecked_Reads, na.rm = TRUE)
  )

if("ITSx_Extracted_Reads" %in% colnames(PER_SAMPLE_COUNTS_merged)){
  PER_RUN_COUNTS_merged[ , Reads_ITSx_Extracted := sum(PER_SAMPLE_COUNTS_merged$ITSx_Extracted_Reads, na.rm = TRUE) ]
}

## Estimate percentage of reads passed primer checking
cat("..Estimating per-run percentages\n")
PER_RUN_COUNTS_merged[ , Percentage_Demultiplexed := 
  round(Reads_Demultiplexed / Total_Number_Of_Reads * 100, 1) ]

PER_RUN_COUNTS_merged[ , Percentage_PrimerChecked := 
  round(Reads_PrimerChecked / Total_Number_Of_Reads * 100, 1) ]

## Final per-run num reads
PER_RUN_COUNTS_merged[ , SeqTable_NumReads := sum(PER_SAMPLE_COUNTS_merged$SeqTable_NumReads, na.rm = TRUE) ]
PER_RUN_COUNTS_merged[ , Percentage_Reads_Retained := round( SeqTable_NumReads / Total_Number_Of_Reads * 100, 2) ]

NumUniqSeqs <- SEQTAB |> dplyr::select(SeqID) |>   dplyr::summarize(N = n()) |> dplyr::collect()
PER_RUN_COUNTS_merged$SeqTable_NumUniqueSequences <- NumUniqSeqs$N

## Export summary stats
cat("Exporting results\n")
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
