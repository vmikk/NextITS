#!/usr/bin/env Rscript

## Script to pool quality-filtered and trimmed sequences from multiple sequencing runs
## And summarize sequence abundance at OTU level (per sample)
## 

# Input:
#   1. UC file from global dereplication  (`Dereplicated.uc.gz`)
#   2. UC file from clustering            (`Clustered.uc.gz`)
#   3. FASTA file with OTU sequences      (`Clustered.fa.gz`)
#   4. Max MEEP score
#   5. Max de novo chimera score

# Outputs:
#  - OTU table in long format    (`OTU_table_long.txt.gz` & `OTU_table_long.RData`)
#  - OTU table in wide format    (`OTU_table_wide.txt.gz` & `OTU_table_wide.RData`)
#  - FASTA file with sequences   (`OTUs.fa.gz`)

## Usage:
# ./pool_seq_runs.R \
#    --ucderep "Dereplicated.uc.gz" \
#    --ucclust "Clustered.uc.gz" \
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



############################################## Parse input parameters

## Check time
start_time <- Sys.time()


cat("Parsing input options and arguments...\n")

suppressPackageStartupMessages(require(optparse))

## Parse arguments
option_list <- list(
  make_option("--ucderep", action="store", default=NA,  type='character', help="UC file from global dereplication"),
  make_option("--ucpreclust", action="store", default="NoPrecluster",  type='character', help="UC file from pre-clustering (optional; use 'NoPrecluster' to disable this step)"),
  make_option("--ucclust", action="store", default=NA,  type='character', help="UC file from clustering"),
  make_option("--otus",    action="store", default=NA,  type='character', help="FASTA file with OTU sequences"),
  make_option("--maxmeep", action="store", default=0.5, type='double', help="Max MEEP score"),
  make_option("--maxchim", action="store", default=0.6, type='double', help="Max de novo chimera score"),
  
  make_option("--recoverdenovo", action="store", default=TRUE, type='logical', help="Recover de-novo chimeras"),
  make_option("--recoversinglet", action="store", default=TRUE, type='logical', help="Recover singletons"),
  
  make_option(c("-m", "--mergesamples"), action="store", default=FALSE, type='logical', help="Merge sample replicates (default, true)"),
  
  make_option(c("-t", "--threads"), action="store", default=4L, type='integer', help="Number of CPU threads for arrow, default 4")
  # make_option(c("-s", "--scriptdir"),  action="store", default=getwd(), type='character', help="Directory containing source scripts")
)
opt <- parse_args(OptionParser(option_list=option_list))


## Validation of the required argiments
if(is.na(opt$ucderep)){
  cat("Input file is not specified: UC file from global dereplication.\n", file=stderr())
  stop()
}
if(is.na(opt$ucclust)){
  cat("Input file is not specified: UC file from clustering.\n", file=stderr())
  stop()
}
if(is.na(opt$otus)){
  cat("Input file is not specified: FASTA file with OTU sequences.\n", file=stderr())
  stop()
}

## Function to convert text "NA"s to NA
to_na <- function(x){ 
  if(x %in% c("NA", "null", "Null")){ x <- NA }
  return(x)
}

## Assign variables
UCDEREP    <- to_na( opt$ucderep )
UCPRECLUST <- to_na( opt$ucpreclust )
UCCLUST    <- to_na( opt$ucclust )
MAXMEEP    <- as.numeric( to_na( opt$maxmeep ) )
MAXCHIM    <- as.numeric( to_na( opt$maxchim ) )

RECOV_DENOVO  <- as.logical(opt$recoverdenovo)
RECOV_SINGLET <- as.logical(opt$recoversinglet)
MERGE_SAMPLES <- as.logical(opt$mergesamples)
OTUS          <- opt$otus

CPUTHREADS <- as.numeric( opt$threads )
# SCRIPTDIR  <- opt$scriptdir

## Log assigned variables
cat(paste("UC file from global dereplication: ", UCDEREP, "\n", sep=""))
cat(paste("UC file from pre-clustering or denoising: ", UCPRECLUST, "\n", sep=""))
cat(paste("UC file from clustering: ",           UCCLUST, "\n", sep=""))
cat(paste("Max MEEP score: ",                    MAXMEEP, "\n", sep=""))
cat(paste("Max de novo chimera score: ",         MAXCHIM, "\n", sep=""))
cat(paste("De novo chimera recovery: ",          RECOV_DENOVO,  "\n", sep=""))
cat(paste("Low-quality singleton recovery: ",    RECOV_SINGLET, "\n", sep=""))
cat(paste("Merge sample replicates: ",           MERGE_SAMPLES, "\n", sep=""))
cat(paste("OTU sequences: ",                     OTUS,          "\n", sep=""))
cat(paste("Number of CPU threads to use: ",      CPUTHREADS,    "\n", sep=""))
# cat(paste("Directory containing source scripts: ", SCRIPTDIR, "\n", sep=""))

cat("\n")



############################################## Data for debugging

# UCDEREP       <- "Dereplicated.uc.gz"
# UCCLUST       <- "Clustered.uc.gz"
# UCPRECLUST    <- "NoPrecluster"         # UNOISE.uc.gz"
# MAXMEEP       <- 0.5
# MAXCHIM       <- 0.6
# RECOV_DENOVO  <- TRUE
# RECOV_SINGLET <- TRUE
# MERGE_SAMPLES <- TRUE
# OTUS          <- "Clustered.fa.gz"
# CPUTHREADS    <- 4


############################################## Load packages and data

cat("Loading R packages...\n")

load_pckg <- function(pkg = "data.table"){
  suppressPackageStartupMessages( library(package = pkg, character.only = TRUE) )
  cat(paste(pkg, packageVersion(pkg), "\n"))
}

load_pckg("data.table")
load_pckg("plyr")
load_pckg("metagMisc")
load_pckg("Biostrings")

cat("\n")


# cat("Loading additional R funcitons...\n")
# source(file.path(SCRIPTDIR, "R_functions.R"))
# cat("\n")


## Set CPU thread number
cat("Setting number of CPU threads to: ", CPUTHREADS, "\n")
setDTthreads(threads = CPUTHREADS)  # for data.table


######################################
###################################### Load the data
######################################

## Load sequence tables
cat("..Loading sequence tables\n")
TABS <- list.files(path = ".", pattern = "Seqs.RData", full.names = TRUE, recursive = TRUE)
cat("... Tables found: ", length(TABS), "\n")

cat("..Loading sequence tables\n")
TAB <- alply(.data = TABS, .margins = 1, .fun = readRDS)
TAB <- rbindlist(TAB, use.names = TRUE, fill = TRUE)
cat("... Total number of records: ",             nrow(TAB), "\n")
cat("... Total number unique sequences: ",       length(unique(TAB$Sequence)), "\n")
cat("... Total number unique samples (files): ", length(unique(TAB$SampleID)), "\n")

## Load UC file for globally dereplicated sequences
cat("..Loading dereplication UC file\n")
UCA <- metagMisc::parse_uc(x = UCDEREP, map_only = TRUE, package = "data.table")
UCA <- unique(UCA)
colnames(UCA) <- c("SeqID", "DerepID")

## Load UC file for pre-clustering or denoising (e.g., UNOISE) - optional
cat("..Loading pre-clustering or denoising UC file\n")
if(UCPRECLUST %in% "NoPrecluster"){
  cat("...No file provided, skipping this step\n")
} else {
  UCP <- metagMisc::parse_uc(x = UCPRECLUST, map_only = TRUE)
  UCP <- unique(UCP)
  colnames(UCP) <- c("DerepID", "PreclusterID")
}

## Load UC file for clustered sequences (OTUs / SWARM clusters / UNOISE)
cat("..Loading clustering UC file\n")
UCO <- metagMisc::parse_uc(x = UCCLUST, map_only = TRUE, package = "data.table")
UCO <- unique(UCO)
colnames(UCO) <- c("DerepID", "OTU")

## Filter sequences by the max number of expected error per 100bp (MEEP)
if(!is.na(MAXMEEP)){
  cat("..Filtering data by max MEEP score\n")
  nrecs <- nrow(TAB)                           # number of records with all low-quality seqs
  nabun <- sum(TAB$Abundance, na.rm = TRUE)    # total abundance with low-quality seqs
  
  ## If no sequence recovery is reqired
  if(RECOV_SINGLET == FALSE){

    ## Remove low-quality sequences
    TAB <- TAB[ MEEP < MAXMEEP | is.na(MEEP) ]

  ## If we need to recover low-quality sequences with multiple occurrence
  } else {

    ## Get low-quality singleton
    SINGLETONS_lowquality <- TAB[ MEEP >= MAXMEEP, .(SeqID___SampleID, Sequence) ]

    ## Recover non-unique singletons
    SINGLETONS_dups <- duplicated(SINGLETONS_lowquality$Sequence)

    if(any(SINGLETONS_dups)){
      cat(".... There are a few sequencing-run-level singletons with relatively low quality, that occurr in the other sequncing runs.\n")
      cat(".... N = ", sum(SINGLETONS_dups), "\n")

      SINGLETONS_torecover <- unique( SINGLETONS_lowquality$Sequence[ SINGLETONS_dups ] )
      cat(".... Recovering ", length(SINGLETONS_torecover), " unique sequences\n")
      SINGLETONS_lowquality <- SINGLETONS_lowquality[ ! Sequence %in% SINGLETONS_torecover ]

      rm(SINGLETONS_torecover)
    }

    cat(".... In total, ", nrow(SINGLETONS_lowquality), " low-quality singletons will be removed\n")
    
    if(nrow(SINGLETONS_lowquality) > 0){
      TAB <- TAB[ ! SeqID___SampleID %in% SINGLETONS_lowquality$SeqID___SampleID ]
    }

    rm(SINGLETONS_lowquality, SINGLETONS_dups)

  } # end of recovery

  nrecs_delta <- nrecs - nrow(TAB)                          # num records after filtering
  nabun_delta <- nabun - sum(TAB$Abundance, na.rm = TRUE)   # total abundance after filtering

  cat("... Records removed: ", nrecs_delta, " (", round(nrecs_delta/nrecs * 100, 1), "%)\n")
  cat("... Reads removed: ", nabun_delta, " (", round(nabun_delta/nabun * 100, 1), "%)\n")

  rm(nrecs_delta, nabun_delta)

} # end of MAXMEEP filtering



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


## Prepare sequence IDs without pre-clustering or denoising
if(UCPRECLUST %in% "NoPrecluster"){
  
  ## Prepare OTU IDs
  cat("..Preparing OTU IDs [no pre-clustering or denoising ]\n")
  UCA <- merge(x = UCA, y = UCO, by = "DerepID", all.x = TRUE)

  ## Add OTU IDs to seq table
  cat("... Adding OTU IDs to sequence table\n")
  cat(".... Number of records in sequence table before merging: ", nrow(TAB), "\n")
  TAB <- merge(x = TAB, y = UCA, by = "SeqID", all.x = TRUE)
  cat(".... Number of records in sequence table after merging: ",  nrow(TAB), "\n")

} else {

  ## Prepare OTU IDs
  cat("..Preparing OTU IDs [with pre-clustering or denoising ]\n")
  
  cat("... Adding pre-cluster or denoised IDs\n")
  UCA <- merge(x = UCA, y = UCP, by = "DerepID", all.x = TRUE)

  cat("... Adding clustering IDs\n")
  UCA <- merge(x = UCA, y = UCO, by.x = "PreclusterID", by.y = "DerepID", all.x = TRUE)

  cat(".... Number of sequences without OTU ID: ", sum(is.na(UCA$OTU)), "\n")

  cat("... Adding OTU IDs to sequence table\n")
  cat(".... Number of records in sequence table before merging: ", nrow(TAB), "\n")
  TAB <- merge(x = TAB, y = UCA, by = "SeqID", all.x = TRUE)
  cat(".... Number of records in sequence table after merging: ",  nrow(TAB), "\n")

}



## Remove NA OTUs -- probably excluded seqs
if(any(is.na(TAB$OTU))){
  cat("WARNING: not all sequences were assigned to OTUs\n")
  cat("..Removing missing/excluded sequences\n")
  cat(".. ", sum(is.na(TAB$OTU)), " sequences with total abundance ", 
      sum(TAB[ is.na(OTU) ]$Abundance, na.rm = TRUE), " reads will be excluded\n")
  TAB <- TAB[ ! is.na(OTU) ]
}


## Summarize abundance by sample and OTU
if(MERGE_SAMPLES == TRUE){

  cat("Summarizing OTU abundance and merging sample replicates (e.g., re-sequenced samples)\n")

  ## Extract sample names
  cat("..Extracting sample names\n")
  TAB[ , SampleName := tstrsplit(x = SampleID, split = "__", keep = 2) ]

  cat("..Summarizing abundance by sample and OTU\n")
  RES <- TAB[ , 
    .( Abundance  = sum(Abundance,  na.rm = TRUE) ),
    by = c("OTU", "SampleName") ]

  setnames(x = RES, old = "SampleName", new = "SampleID")

} else {

  cat("Summarizing OTU abundance\n")

  cat("..Summarizing abundance by sample and OTU\n")
  RES <- TAB[ , 
    .( Abundance  = sum(Abundance,  na.rm = TRUE) ),
    by = c("OTU", "SampleID") ]

}

#### Reshape to wide table
cat("Reshaping table into wide format\n")

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

cat("..Reshaping to the wide format done!\n")


## Reorder OTU rows
cat("..Reordering OTU rows by total abundance\n")
otu_tots <- rowSums(REW[, -1], na.rm = TRUE)
REW <- REW[ order(otu_tots, decreasing = T), ]

## Add attributes if samples were merged
setattr(x = RES, name = "Samples_merged", value = MERGE_SAMPLES)
setattr(x = REW, name = "Samples_merged", value = MERGE_SAMPLES)


cat("Exporting results\n")

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


cat("..Exporting OTU sequences to FASTA\n")

cat("...Preparing sequences\n")

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
cat("....Loading FASTA file\n")
SQS <- readDNAStringSet(filepath = OTUS, format="fasta")
cat("....Extracting sequence IDs\n")
names(SQS) <- tstrsplit(x = names(SQS), split = ";", keep = 1)[[1]]

if(any(duplicated(names(SQS)))){
  cat("WARNING: duplicated OTU names detected!\n")
}

cat("....Subsetting OTUs\n")
SQF <- SQS[ names(SQS) %in% unique(REW$OTU) ]

cat("....Total number of OTUs in the input file: ", length(SQS), "\n")
cat("....Number of OTUs to export: ",               length(SQF), "\n")
cat("....Number of OTUs in the OTU talbe: ",        nrow(REW),   "\n")

cat("...Writing FASTA file\n")
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
