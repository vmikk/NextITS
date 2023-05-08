#!/usr/bin/env Rscript

## Script to pool quality-filtered and trimmed sequences from multiple sequencing runs
## And summarize sequence abundance at OTU level (per sample)

############################################## Parse input parameters

## Check time
start_time <- Sys.time()


cat("Parsing input options and arguments...\n")

suppressPackageStartupMessages(require(optparse))

## Parse arguments
option_list <- list(
  make_option("--ucderep", action="store", default=NA,  type='character', help="UC file from global dereplication"),
  make_option("--ucclust", action="store", default=NA,  type='character', help="UC file from clustering"),
  make_option("--maxmeep", action="store", default=0.5, type='double', help="Max MEEP score"),
  make_option("--maxchim", action="store", default=0.6, type='double', help="Max de novo chimera score"),
  
  make_option("--recoverdenovo", action="store", default=TRUE, type='logical', help="Recover de-novo chimeras"),
  make_option("--recoversinglet", action="store", default=TRUE, type='logical', help="Recover singletons"),
  
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

## Function to convert text "NA"s to NA
to_na <- function(x){ 
  if(x %in% c("NA", "null", "Null")){ x <- NA }
  return(x)
}

## Assign variables
UCDEREP    <- to_na( opt$ucderep )
UCCLUST    <- to_na( opt$ucclust )
MAXMEEP    <- as.numeric( to_na( opt$maxmeep ) )
MAXCHIM    <- as.numeric( to_na( opt$maxchim ) )

RECOV_DENOVO  <- as.logical(opt$recoverdenovo)
RECOV_SINGLET <- as.logical(opt$recoversinglet)

CPUTHREADS <- as.numeric( opt$threads )
# SCRIPTDIR  <- opt$scriptdir

## Log assigned variables
cat(paste("UC file from global dereplication: ", UCDEREP, "\n", sep=""))
cat(paste("UC file from clustering: ",           UCCLUST, "\n", sep=""))
cat(paste("Max MEEP score: ",                    MAXMEEP, "\n", sep=""))
cat(paste("Max de novo chimera score: ",         MAXCHIM, "\n", sep=""))
cat(paste("De novo chimera recovery: ",          RECOV_DENOVO,  "\n", sep=""))
cat(paste("Low-quality singleton recovery: ",    RECOV_SINGLET, "\n", sep=""))
cat(paste("Number of CPU threads to use: ",      CPUTHREADS,    "\n", sep=""))
# cat(paste("Directory containing source scripts: ", SCRIPTDIR, "\n", sep=""))

cat("\n")


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

## Load ASV table
cat("..Loading sequence tables\n")
TABS <- list.files(path = ".", pattern = "ASVs.RData", full.names = TRUE, recursive = TRUE)
cat("... Tables found: ", length(TABS), "\n")

cat("..Loading sequence tables\n")
TAB <- alply(.data = TABS, .margins = 1, .fun = readRDS)
TAB <- rbindlist(TAB, use.names = TRUE, fill = TRUE)
cat("... Total number of records: ", nrow(TAB), "\n")

## Load UC file for globally dereplicated sequences
cat("..Loading dereplication UC file\n")
UCA <- metagMisc::parse_uc(x = UCDEREP, map_only = TRUE, package = "data.table")
UCA <- unique(UCA)
colnames(UCA) <- c("SeqID", "DerepID")

## Load UC file for clustered sequences (OTUs / SWARM clusters / UNOISE)
cat("..Loading clustering UC file\n")
UCO <- metagMisc::parse_uc(x = UCCLUST, map_only = TRUE, package = "data.table")
UCO <- unique(UCO)
colnames(UCO) <- c("DerepID", "OTU")

## Filter sequences by MEEP
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



## Filter sequences by MEEP
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


## Extract sample names
cat("..Extracting sample names\n")
TAB[ , SampleName := tstrsplit(x = SampleID, split = "__", keep = 2) ]

## Prepare OTU IDs
cat("..Preparing OTU IDs\n")
UCA <- merge(x = UCA, y = UCO, by = "DerepID", all.x = TRUE)

## Add OTU IDs to seq table
cat("..Adding OTU IDs to sequence table\n")
nrow(TAB)
TAB <- merge(x = TAB, y = UCA, by = "SeqID", all.x = TRUE)
nrow(TAB)

## Remove NA OTUs -- probably excluded seqs
if(any(is.na(TAB$OTU))){
  cat("WARNING: not all sequences were assigned to OTUs\n")
  cat("..Removing missing/excluded sequences\n")
  cat(".. ", sum(is.na(TAB$OTU)), " sequences with total abundance ", 
      sum(TAB[ is.na(OTU) ]$Abundance, na.rm = TRUE), " reads will be excluded\n")
  TAB <- TAB[ ! is.na(OTU) ]
}

## Summarize abundance by sample and OTU
cat("..Summarizing abundance by sample and OTU\n")
RES <- TAB[ , .( Abundance  = sum(Abundance,  na.rm = TRUE) ),
  by = c("OTU", "SampleName", "Sequence") ]

## Reshape to wide table
cat("..Reshaping table into wide format\n")
REW <- dcast(data = RES, 
  formula = OTU ~ SampleName, 
  fun.aggregate = sum, fill = 0, value.var = "Abundance")

## Reorder OTU rows
cat("..Reordering OTU rows by total abundance\n")
otu_tots <- rowSums(REW[, -1], na.rm = TRUE)
REW <- REW[ order(otu_tots, decreasing = T), ]



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
