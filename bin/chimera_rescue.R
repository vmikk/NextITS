#!/usr/bin/env Rscript

## Script to rescue sequences that were annotated as chimeric, 
## but have high occurrence within sequenceing run (occurrence > 2) 

# Input is given as positional arguments:
#   1. List of all chimeric sequences (`All_chimeras.txt.gz`)
#   2. Min sequence occurrence to be preserved (e.g., 2)
#   3. Output file name (`Rescued_Chimeric_sequences.fa.gz`)

suppressMessages(library(data.table)); setDTthreads(threads = 1)
suppressMessages(library(Biostrings))

args <- commandArgs(trailingOnly = TRUE)

MINOCC <- as.numeric( args[2] )

## Load sequences
cat("..Loading chimeric sequences\n")
CH <- try(
  fread(file = args[1],
    sep = "\t", header = F,
    col.names = c("SeqID", "Seq"))
  )

if("try-error" %in% class(CH)){
  cat("\nCould not read the file with chimeric sequences\n")
  cat("Most likely, the file file is empty (no chimeras)\n")
   q(save = "no", status = 0, runLast = FALSE)
}

cat("..Total number of chimeric records: ", nrow(CH), "\n")

if(nrow(CH) > 0){

  ## Extract sample name and sequencing run ID
  CH[, SampleID := tstrsplit(x = SeqID, split = ";", keep = 2) ]
  CH[, SampleID := gsub(pattern = "sample=", replacement = "", x = SampleID) ]

  ## Estimate sequence frequency
  cat("..Estimating chimera occurrence\n")
  CF <- CH[, .(Occurrence = .N), by = "Seq"]

  ## Exclude sequences with low occurrence (most probably chimeric)
  ## Sequences with higher occurrence should be "real" sequences
  CF <- CF[ Occurrence > MINOCC ]

  cat("..Total number of unique chimeric sequences: ", length(unique(CH$Seq)), ".\n")

  ## Export sequences
  if(nrow(CF) > 0){
    cat("..There are", nrow(CF), "unique sequence to rescue.\n")
    NCH <- CH[ Seq %in% CF$Seq ]
    SQS <- DNAStringSet(x = NCH$Seq)
    names(SQS) <- NCH$SeqID

    cat("..Exporting rescued sequences\n")
    writeXStringSet(x = SQS, filepath = args[3],
      compress=TRUE, format="fasta", width=9999)
  } else {
    cat("..No sequences were rescued.\n")
  }

}