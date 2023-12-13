#!/usr/bin/env Rscript

## Perform sequence denoising with DADA2

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
