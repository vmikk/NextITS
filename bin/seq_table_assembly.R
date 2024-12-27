#!/usr/bin/env Rscript

## Script to perform tag-jump removal

## To do:
#  - add HMM profile ID if ITSx was used
#  - export data in Excel format? (if table is not too large)
#  - compress output tabs?


# Input is given as positional arguments:
#   1. non-filtered Seq table   (`Seq_tab_not_filtered.txt.gz`)
#   2. Sequences in fasta       (`Seq_not_filtered.fa.gz`)
#   3. sequence mapping to OTUs (`Sample_mapping.uc.gz`)
#   4. tag-jumped OTU list      (`TagJump_OTUs.RData`)
#   5. de novo chimera scores   (`DeNovo_Chimera.txt`)
#   6. sequence qualities       (`SeqQualities.parquet`)

# Outputs:
#  - FASTA with filtered Seqs       `Seqs.fa.gz`
#  - Seq table in long format       `Seqs.txt.gz`  (with additional sequence info)
#  - Seq table in wide format       `Seq_tab.txt.gz`
#  - Data in R-serialization format `Seqs.RData`


args <- commandArgs(trailingOnly = TRUE)

## Debug:
# args <- c(
#   "Seq_tab_not_filtered.txt.gz",
#   "Seq_not_filtered.fa.gz",
#   "Sample_mapping.uc.gz",
#   "TagJump_OTUs.RData",
#   "DeNovo_Chimera.txt",
#   "SeqQualities.parquet"
#   )

suppressMessages(library(data.table))
suppressMessages(library(Biostrings))
suppressMessages(library(plyr))
suppressMessages(library(arrow))
suppressMessages(library(openxlsx))


######################################
###################################### Load the data
######################################

## Load sequnece table
cat("..Loading non-filtered sequence table\n")
TAB <- fread(
  file = args[1],
  sep = "\t", header = FALSE,
  col.names = c("SeqID", "Abundance", "SampleID"))


## Load sequences in fasta format
cat("..Loading sequneces in FASTA format\n")
SQS <- readDNAStringSet(filepath = args[2])


## Load sequence mapping to OTUs
cat("..Loading sequence mapping table\n")
MAP <- fread(
  file = args[3],
  header = FALSE, sep = "\t")

MAP <- MAP[ V1 != "S" ]
MAP[, c("SeqID", "SampleID") := tstrsplit(V9, ";", keep = c(1,3)) ]
MAP[, OTU := tstrsplit(V10, ";", keep = 1) ]
MAP[V1 == "C", OTU := SeqID ]
MAP[, SampleID := gsub(pattern = "sample=", replacement = "", x = SampleID) ]
MAP <- MAP[, .(SeqID, SampleID, OTU) ]


## Load list of tag-jumped OTUs
cat("..Loading list of tag-jumped OTUs\n")
JMP <- readRDS( args[4] )


## Load de novo chimera scores
cat("..Loading de novo chimera scores\n")

CHI <- try(
  fread(
    file = args[5],
    header = FALSE, sep = "\t",
    col.names = c("SeqID", "DeNovo_Chimera_Score", "SampleID"))
  )

if("try-error" %in% class(CHI)){
  cat("\nCould not read the file with de novo chimeric scores\n")
  cat("Most likely, the file file is empty (no de novo chimeras)\n")
  
  ## Initialize empty data table
  CHI <- data.table(SeqID = character(), DeNovo_Chimera_Score = numeric(), SampleID = character())
}


## Load sequence quality scores
cat("..Loading sequence quality scores\n")
QLT <- arrow::read_parquet(file = args[6])
setDT(QLT)
clz <- c("SampleID", "Hash", "Length", "AvgPhredScore", "MaxEE", "MEEP")
QLT <- QLT[ , ..clz ]
setnames(QLT,
  old = c("Hash", "Length", "AvgPhredScore"),
  new = c("SeqID", "SeqLen", "PhredScore"))

# old header: c("SampleID", "SeqID", "SeqLen", "PhredScore", "MaxEE", "MEEP")
# new header: c("SampleID", "Hash", "PacBioID", "PhredScore", "MaxEE", "MEEP", "Sequence", "Quality", "Length")



## Create SeqID___SampleID column
TAB[, SeqID___SampleID := paste0(SeqID, "___", SampleID) ]
QLT[, SeqID___SampleID := paste0(SeqID, "___", SampleID) ]
MAP[, SeqID___SampleID := paste0(SeqID, "___", SampleID) ]

MAP[, c("SeqID", "SampleID") := NULL ]

######################################
###################################### Remove tag-jumps
######################################

cat("..Removing tag-jumped sequences\n")

if(nrow(JMP) > 0){

  # JMP[ , SeqID___SampleID := paste0(OTU, "___", SampleID) ]
  JMP[ , TagJump := TRUE ]

  ## Add OTU ID to sequences
  cat("...Adding OTU IDs to sequences\n")

  TAB <- merge(x = TAB, y = MAP,
    by = "SeqID___SampleID", all.x = TRUE)

  ## Add tag-jump information to the sequence table
  TAB <- merge(x = TAB, y = JMP,
    by = c("OTU", "SampleID"),
    all.x = TRUE)

  cat("... ", sum(TAB$TagJump, na.rm = TRUE), " tag-jump occurrences detected\n")

  ## Remove tag-jumps
  if(any(TAB$TagJump)){
    TAB <- TAB[(!TagJump %in% TRUE)]
    # because of NAs, TAB[ ! TagJump ] does not work properly
  }

  TAB[, TagJump := NULL ]
  TAB[, OTU := NULL ]

# end of `nrow(JMP) > 0`
} else {

  cat("...no tag-jumps found\n")

}


######################################
###################################### Add singleton qualities
######################################

cat("..Looking for singleton sequences\n")

## Find singletons
cat("...Counting sequence occurrence\n")
SNG <- TAB[ , .(Occurrence = .N), by = "SeqID" ]
SNG <- SNG[ Occurrence < 2 ]$SeqID

cat("... ", length(SNG), " sequences with single occurrence found\n")

## Remove non-singleton seqs from the quality table
if(length(SNG) > 0){
  cat("...Subsetting quality table\n")
  
  ## Remove sequences with multiple reads
  dups <- QLT$SeqID___SampleID[ which(duplicated(QLT$SeqID___SampleID)) ]
  if(length(dups) > 0){
    QLT <- QLT[ ! SeqID___SampleID %in% dups ]
  }

  ## Keep only single-occurrence sequence
  QLT <- QLT[ SeqID %in% SNG ]
  QLT[ , SeqLen := NULL ]
}

## Add Phred scores to the main table
if(nrow(QLT) > 0){

cat("...Adding Phred-scores to the main table\n")
TAB <- merge(x = TAB,
  y = QLT[, .(SeqID___SampleID, PhredScore, MaxEE, MEEP) ],
  by = c("SeqID___SampleID"), all.x = TRUE)

## Remove scores for non-singleton sequences
TAB[ Abundance > 1, PhredScore := NA ]
TAB[ Abundance > 1, MaxEE := NA ]
TAB[ Abundance > 1, MEEP := NA ]

} else {
## No singletons
TAB[ , PhredScore := NA ]
TAB[ , MaxEE      := NA ]
TAB[ , MEEP       := NA ]
}

## Convert variables to numberic scores
TAB[ , PhredScore := as.numeric(PhredScore) ]
TAB[ , MaxEE      := as.numeric(MaxEE) ]
TAB[ , MEEP       := as.numeric(MEEP) ]

# with(TAB, plot(Abundance, PhredScore))


######################################
###################################### Add chimera info
######################################

cat("..Adding info about de novo chimeric sequences\n")

if(nrow(CHI) > 0){

  TAB <- merge(x = TAB, y = CHI,
    by = c("SeqID", "SampleID"), all.x = TRUE)

  ## Convert variables to numeric scores
  TAB[ , DeNovo_Chimera_Score := as.numeric(DeNovo_Chimera_Score) ]

  ## Classify sequences into putative chimeras
  TAB[ !is.na(DeNovo_Chimera_Score), DeNovo_Chimera := TRUE  ]
  TAB[  is.na(DeNovo_Chimera_Score), DeNovo_Chimera := FALSE ]

  cat("... ", sum( TAB$DeNovo_Chimera), " putative de novo chimeras found\n")
  cat("... ", sum(!TAB$DeNovo_Chimera), " non-chimeric sequences\n")

} else {

  ## No de novo chimeras

  TAB[ , DeNovo_Chimera_Score := as.numeric(NA) ]
  TAB[ , DeNovo_Chimera       := FALSE  ]

  cat("... ", 0,         " putative de novo chimeras found\n")
  cat("... ", nrow(TAB), " non-chimeric sequences\n")

}


######################################
###################################### Add sequences
######################################

cat("..Processing sequences\n")

SQTAB <- data.table(
  SeqHeader = names(SQS),
  Sequence = as.character(SQS))

## Split the header  (`feb76b9;size=1;sample=ABCD;` )
SQTAB[ , c("SeqID", "SampleID") := tstrsplit(x = SeqHeader, split = ";", keep = c(1,3)) ]
SQTAB[ , SeqHeader := NULL ]
SQTAB[ , SampleID := gsub(pattern = "sample=", replacement = "", x = SampleID) ]


SQTAB[ , SeqID___SampleID := paste0(SeqID, "___", SampleID) ]
SQTAB[ , c("SeqID", "SampleID") := NULL ]

cat("..Adding sequences to the main table\n")

TAB <- merge(x = TAB, y = SQTAB,
  by = c("SeqID___SampleID"), all.x = TRUE)


cat("..Sorting table by SampleID and number of reads per sequence\n")

setorder(x = TAB, SampleID, -Abundance)


cat("..Preparing FASTA file with filtered sequences\n")

SQF <- DNAStringSet(x = TAB$Sequence)
names(SQF) <- paste0(TAB$SeqID, ";size=", TAB$Abundance, ";sample=", TAB$SampleID, ";")

## Export FASTA
cat("..Exporting FASTA file with filtered sequences\n")

writeXStringSet(x = SQF,
  filepath = "Seqs.fa.gz",
  compress = TRUE, format = "fasta", width = 9999)




######################################
###################################### Export results
######################################

cat("..Reshaping sequence table into wide format\n")

TABW <- dcast(data = TAB,
  formula = SeqID ~ SampleID, 
  value.var = "Abundance",
  fill = 0)


cat("..Exporting result\n")

cat("...Exporting RData\n")
saveRDS(object = TAB, file = "Seqs.RData", compress = "xz")

## Long table
cat("...Exporting long table\n")

TAB[ , SeqID___SampleID := NULL ]
setcolorder(
  x = TAB,
  neworder = c("SampleID", "SeqID", "Abundance",
               "PhredScore", "MaxEE", "MEEP",
               "DeNovo_Chimera", "DeNovo_Chimera_Score",
               "Sequence"))

fwrite(x = TAB, file = "Seqs.txt.gz", sep = "\t", compress = "gzip")

## Wide table
cat("...Exporting wide table\n")
fwrite(x = TABW, file = "Seq_tab.txt.gz", sep = "\t", compress = "gzip")


cat("All done.")
