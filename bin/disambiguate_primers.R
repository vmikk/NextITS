#!/usr/bin/Rscript

## The script to disambiguate sequences
## (expand ambiguous nucleotides into all combinations)
## Based on IUPAC codes

# Input is given as positional arguments:
#   1. A text string (e.g., "ACTGNK")
#   2. output file name (e.g., "Primer_F.fasta")

# Output:
#  - FASTA with disambiguated sequences

args <- commandArgs(trailingOnly = TRUE)

cat("..Loading packages\n")
suppressMessages(library(DECIPHER))
suppressMessages(library(Biostrings))

## Convert input string into DNAStringSet object
cat("..Preparing DNAStringSet\n")
dna <- DNAStringSet(args[1])

## Disambiguate
cat("..Disambiguating\n")
res <- Disambiguate(dna)[[1]]

## Assign names
names(res) <- paste0("seq", 1:length(res), sep = "")

## Export FASTA
cat("..Exporting FASTA\n")
writeXStringSet(x = res,
  filepath = args[2],
  compress=FALSE, format="fasta", width=9999)

cat("..done\n")
