#!/usr/bin/env Rscript

## Script to document the Step-1 workflow of the NextITS pipeline.

## Function to load packages
load_pckg <- function(pkg = "data.table"){
  suppressPackageStartupMessages( library(package = pkg, character.only = TRUE) )
  cat(paste(pkg, packageVersion(pkg), "\n"))
}

load_pckg("glue")
load_pckg("data.table")
load_pckg("yaml")

## Parse arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: document_s1.R <software_versions.yml> <params_table.csv/tsv> [output_path]\n")
  stop()
}

versions_path <- args[[1]]
params_path   <- args[[2]]
output_path   <- ifelse(length(args) >= 3, args[[3]], "README_Step1_Methods.txt")

## Citation registry
citation_db <- list(
  nextits   = "Mikryukov V, Anslan S, Tedersoo L (2025) NextITS - A pipeline for metabarcoding fungi and other eukaryotes with full-length ITS sequenced with PacBio. DOI:10.5281/zenodo.15074882",
  nextflow  = "Di Tommaso P, et al. (2017) Nextflow enables reproducible computational workflows. Nat Biotechnol 35, 316-319, DOI:10.1038/nbt.3820",
  lima      = "Pacific Biosciences (2025) LIMA - The PacBio barcode demultiplexer and primer remover. URL: https://lima.how/",
  seqkit    = "Shen W, Sipos B, Zhao L (2024) SeqKit2: A Swiss Army Knife for Sequence and Alignment Processing. iMeta e191. DOI:10.1002/imt2.191",
  csvtk     = "Shen W (2025) csvtk - a cross-platform, efficient and practical CSV/TSV toolkit. URL: https://github.com/shenwei356/csvtk",
# brename   = "Shen W (2025) brename - batch renaming safely, URL: https://github.com/shenwei356/brename",
  cutadapt  = "Martin M (2011) Cutadapt removes adapter sequences. EMBnet.journal 17(1):10-12, DOI:10.14806/ej.17.1.200",
  itsx      = "Bengtsson-Palme J, et al (2013) Improved software detection and extraction of ITS1 and ITS2 from ribosomal ITS sequences of fungi and other eukaryotes for analysis of environmental sequencing data. Methods Ecol Evol 4:914-919, DOI:10.1111/2041-210X.12073",
  vsearch   = "Rognes T, Flouri T, Nichols B, Quince C, Mahé F (2016) VSEARCH: a versatile open source tool for metagenomics. PeerJ 4:e2584. DOI:10.7717/peerj.2584",
  uchime2   = "Edgar RC (2016) UCHIME2: improved chimera prediction for amplicon sequencing. bioRxiv 074252. DOI:10.1101/074252",
  uncross2  = "Edgar RC (2018) UNCROSS2: identification of cross-talk in 16S rRNA OTU tables. bioRxiv 400762. DOI:10.1101/400762",
  chimscore = "Nilsson RH, et al. (2015) A Comprehensive, Automatically Updated Fungal ITS Sequence Dataset for Reference-Based Chimera Control in Environmental Sequencing Efforts. Microbes Environ. 30(2), 145-50. DOI:10.1264/jsme2.ME14121",
  bedtools  = "Quinlan AR, Hall IM (2010) BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics 26:841-842. DOI:10.1093/bioinformatics/btq033",
  duckdb    = "Raasveldt M, Mühleisen H (2019) DuckDB: an Embeddable Analytical Database. SIGMOD '19: Proceedings of the 2019 International Conference on Management of Data, 1981-1984. DOI:10.1145/3299869.332021",
  parallel  = "Tange O (2011) GNU Parallel: The command-line power tool. Usenix Mag 36 (1), 42",
  eukaryome = "Tedersoo L, et al. (2024). EUKARYOME: the rRNA gene reference database for identification of all eukaryotes. Database (Oxford) 12:baae043. DOI:10.1093/database/baae043",
  R         = "R Core Team (2025) R: A Language and Environment for Statistical Computing. R Foundation for Statistical Computing, Vienna, Austria. URL: https://www.R-project.org/",
  arrow     = "Richardson N, Cook I, Crane N, Dunnington D, François R, Keane J, Moldovan-Grünfeld D, Ooms J, Wujciak-Jens J, and Apache Arrow (2025) arrow: Integration to Apache Arrow. URL: https://github.com/apache/arrow/",
  ggplot2   = "Wickham H (2016) ggplot2: Elegant Graphics for Data Analysis. Springer. DOI:10.1007/978-3-319-24277-4",
  biostrings= "Pagès H, Aboyoun P, Gentleman R, DebRoy S (2025) Biostrings: Efficient manipulation of biological strings. DOI:10.18129/B9.bioc.Biostrings",
  datatable = "Barrett T, Dowle M, Srinivasan A, Gorecki J, Chirico M, Hocking T, Schwendinger B, Krylov I (2025) data.table: Extension of data.frame. URL: <https://r-datatable.com>"
)


##################################
################################## Helpers
##################################

## Get version number
getv <- function(v, process, tool){
  # v       = list (from YAML file)
  # process = process name
  # tool    = tool name

  if(is.null(v[[process]]) || is.null(v[[process]][[tool]])){ return(NA_character_) }
  as.character( v[[process]][[tool]] )
}
# E.g., getv(versions, "demux", "lima")


## Get parameter
getp <- function(p, name, default = NA){
  # p       =  table with parameters (two columns: name and value)
  # name    = parameter name
  # default = default value if parameter is not found

  if(is.null(p[[name]]) || is.na(p[[name]])){ return(default) }
  p[[name]]
}
# E.g., getp(params, "lima_minscore", 93)


## Remove NAs and empty strings (to curate the citations)
trim_na <- function(x){
  x[ !is.na(x) & nzchar(x) ]
}

##################################
################################## Body builders
##################################


emit_demux_pacbio <- function(p, v) {
  ms <- getp(p, "lima_minscore", 93)
  mb <- getp(p, "lima_barcodetype", "dual_symmetric")
  vs <- getv(v, "demux", "lima")
  switch(mb,
    "single"          = {barcode_type <- "single-end barcodes"},
    "dual_symmetric"  = {barcode_type <- "symmetric dual-end barcodes"},
    "dual_asymmetric" = {barcode_type <- "asymmetric dual-end barcodes"},
    "dual"            = {barcode_type <- "combination of symmetric and asymmetric dual-end barcodes"})

  glue("Demultiplexed PacBio reads using LIMA v.{vs} (Pacific Biosciences) with min score {ms} and {barcode_type}.")
}


##################################
################################## Workflow-dependent method descriptions
##################################

## Function to assembly the workflow description and references
build_docs <- function(versions, params){
  body <- character()
  tools_used <- character()

      body <- c(body, emit_demux_pacbio(params, versions))
      tools_used <- c(tools_used, c("lima"))
  tools_used <- unique(tools_used)
  citations <- trim_na( unlist(citation_db[tools_used]) )
  citations <- sort(unique(citations))

  res <- list(
    body = body,
    citations = citations)
  
  return(res)
}


##################################
################################## Assemble body and citations
##################################


## Load inputs
cat("Loading versions YAML...\n")
versions <- yaml::read_yaml(versions_path)

cat("Loading params table...\n")
params <- data.table::fread(params_path, sep = "\t", header = TRUE, na.strings = c("", "NA"))
setnames(params, new = c("name", "value"))

