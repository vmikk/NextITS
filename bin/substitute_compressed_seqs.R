#!/usr/bin/Rscript

## Script to replace homopolymer-comressed sequences with non-compressed seqs
## + update size annotation

# Input is given as positional arguments:
#   1. Uncompressed sequences (`inp_tab.txt`)
#   2. Homopolymer-compressed sequences (`clust_tab.txt`)
#   3. Name of the output FASTA file (`res.fa`)

suppressMessages(library(data.table)); setDTthreads(threads = 1)
suppressMessages(library(Biostrings))

args <- commandArgs(trailingOnly = TRUE)

## Load data - Uncompressed (inp_tab.txt)
cat("..Loading original sequences\n")
d1 <- fread(file = args[1],
  header=FALSE, sep = "\t", quote = F, col.names = c("SeqID", "Seq_OK"), selec = 1:2)

## Load data - Compressed  (clust_tab.txt)
cat("..Loading compressed sequences\n")
d2 <- fread(file = args[2],
  header=FALSE, sep = "\t", quote = F, col.names = c("SeqID", "Seq_Compr"), selec = 1:2)

cat("...Number of raw sequences: ", nrow(d1), "\n")
cat("...Number of compressed sequences: ", nrow(d2), "\n")

cat("..Processing data\n")

## Remove multiple separators
d1[, SeqID := gsub(pattern = ";;", replacement = ";", x = SeqID)]
d2[, SeqID := gsub(pattern = ";;", replacement = ";", x = SeqID)]

## Split seq ID
d1[, c("Hash", "Size") := tstrsplit(SeqID, ";", keep=1:2)]
d2[, c("Hash", "Size") := tstrsplit(SeqID, ";", keep=1:2)]

## Drop seq ID
d1[, SeqID := NULL ]
d2[, SeqID := NULL ]

## Replace seqs
res <- merge(
  x = d2[, .(Hash, Size)],
  y = d1[, .(Seq_OK, Hash)],
  by = "Hash", all.x = TRUE)

res[, SeqID := do.call(paste, c(.SD, sep = ";")), .SDcols = c("Hash", "Size")]

## Verify the number of reads - should be the same
# sum(as.numeric(gsub(pattern = "size=", replacement = "", x = d1$Size)))
# sum(as.numeric(gsub(pattern = "size=", replacement = "", x = d2$Size)))
# sum(as.numeric(gsub(pattern = "size=", replacement = "", x = res$Size)))

cat("...Total number of reads: ",
  sum(as.numeric(gsub(pattern = "size=", replacement = "", x = res$Size))),
  "\n")

## Prepare sequences
cat("..Exporting results\n")
sqs <- DNAStringSet(x = res$Seq_OK)
names(sqs) <- res$SeqID

## Export FASTA
writeXStringSet(x = sqs, filepath = args[3],
  compress=FALSE, format="fasta", width=9999)
