#!/usr/bin/env Rscript

## Script to combine Phred scores and MaxEE estimates

# Input is given as positional arguments:
#   1. Phred score table   (`tmp_hash_table.txt`)
#   2. MaxEE table         (`tmp_ee.txt`)
#   3. Output file name    (`${sampID}_hash_table.txt`)


args <- commandArgs(trailingOnly = TRUE)

## Debug:
# args <- c(
#   "tmp_hash_table.txt",
#   "tmp_ee.txt",
#   "res_hash_table.txt"
#   )

suppressMessages(library(data.table))


## Load table with Phred scores
cat("..Loading Phred scores\n")
T1 <- fread(
  file = args[1],
  sep = "\t", header = FALSE,
  col.names  = c("SeqID", "SeqHash", "Len", "PhredScore"),
  colClasses = c("character", "character", "numeric", "numeric"))

if(any(is.na(T1$Len))){
  cat("WARNING: non-numeric data detected. Maybe there are some empty sequences\n")
}

## Load table with Phred scores
cat("..Loading MaxEE estimates\n")
T2 <- fread(
  file = args[2],
  sep = "\t", header = FALSE,
  col.names = c("SeqID", "MaxEE"))

## Merge tables
cat("..Merging tables\n")
TAB <- merge(x = T1, y = T2, by = "SeqID", all.x = TRUE)

## Estimate the MEEP score (Koparde et al., DOI:10.1504/IJCBDD.2017.10006006)
## Maximum number of probable incorrect base calls per every 100 bases in the read
cat("..Estimating MEEP score\n")
TAB[ , MEEP := 100 * MaxEE / Len ]

## Export results
cat("..Exporting results\n")
fwrite(x = TAB,
  file = args[3],
  sep = "\t",
  compress = "none")

