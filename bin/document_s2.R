#!/usr/bin/env Rscript

## Script to document the Step-2 workflow of the NextITS pipeline.

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
  cat("Usage: document_s2.R <software_versions.yml> <params_table.csv/tsv> [output_path]\n")
  stop()
}

versions_path <- args[[1]]
params_path   <- args[[2]]
output_path   <- ifelse(length(args) >= 3, args[[3]], "README_Step2_Methods.txt")


## Validation
if(is.null(versions_path) || versions_path == ""){
  stop("Versions YAML not specified")
}
if(is.null(params_path) || params_path == ""){
  stop("Params table not specified")
}

if(!file.exists(versions_path)){
  stop(glue("Versions YAML not found: {versions_path}"))
}
if(!file.exists(params_path)){
  stop(glue("Params table not found: {params_path}"))
}



##################################
################################## References
##################################

## Citation registry
citation_db <- list(
  nextits   = "Mikryukov V, Anslan S, Tedersoo L (2025) NextITS - A pipeline for metabarcoding fungi and other eukaryotes with full-length ITS sequenced with PacBio. DOI:10.5281/zenodo.15074882",
  nextflow  = "Di Tommaso P, et al. (2017) Nextflow enables reproducible computational workflows. Nat Biotechnol 35, 316-319, DOI:10.1038/nbt.3820",
  vsearch   = "Rognes T, Flouri T, Nichols B, Quince C, Mahé F (2016) VSEARCH: a versatile open source tool for metagenomics. PeerJ 4:e2584. DOI:10.7717/peerj.2584",
  dada2     = "Callahan BJ, et al. (2016) DADA2: High-resolution sample inference from Illumina amplicon data. Nat Methods 13:581-583. DOI:10.1038/nmeth.3869",
  unoise    = "Edgar RC (2016) UNOISE2: improved error-correction for Illumina 16S and ITS amplicon sequencing. bioRxiv 081257. DOI:10.1101/081257",
  swarm     = "Mahé F, Czech L, Stamatakis A, Quince C, de Vargas C, Dunthorn M, Rognes T. (2021) Swarm v3: towards tera-scale amplicon clustering. Bioinformatics 38(1), 267-269. DOI:10.1093/bioinformatics/btab493",
  lulu      = "Frøslev TG, et al. (2017) Algorithm for post-clustering curation of DNA amplicon data yields reliable biodiversity estimates. Nat Commun 8:1188. DOI:10.1038/s41467-017-01312-x",
  mumu      = "Mahé F (2025) MUMU: C++ implementation of LULU, a R package for post-clustering curation of metabarcoding data. URL: https://github.com/frederic-mahe/mumu",
  ucs       = "Mikryukov V (2025) ucs - USEARCH cluster file parser. URL: https://github.com/vmikk/ucs",
  duckdb    = "Raasveldt M, Mühleisen H (2019) DuckDB: an Embeddable Analytical Database. SIGMOD '19: Proceedings of the 2019 International Conference on Management of Data, 1981-1984. DOI:10.1145/3299869.332021",
  R         = "R Core Team (2025) R: A Language and Environment for Statistical Computing. R Foundation for Statistical Computing, Vienna, Austria. URL: https://www.R-project.org/",
  arrow     = "Richardson N, Cook I, Crane N, Dunnington D, François R, Keane J, Moldovan-Grünfeld D, Ooms J, Wujciak-Jens J, and Apache Arrow (2025) arrow: Integration to Apache Arrow. URL: https://github.com/apache/arrow/",
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

  if(is.null(v[[process]]) || is.null(v[[process]][[tool]])){ return("") }
  as.character( v[[process]][[tool]] )
}
# E.g., getv(versions, "dereplication", "vsearch")


## Get parameter
getp <- function(p, pname, default = NA){
  # p       =  table with parameters (two columns: name and value)
  # pname   = parameter name
  # default = default value if parameter is not found

  pp <- p[ name == pname ]$value
  if(is.null(pp) || is.na(pp)){ return(default) }
  return(pp)
}
# E.g., getp(params, "otu_id", 0.98)


## Remove NAs and empty strings (to curate the citations)
trim_na <- function(x){
  x[ !is.na(x) & nzchar(x) ]
}

##################################
################################## Body builders
##################################

emit_nextits <- function(v) {
  nextits_v <- if(!is.null(v$NextITS$version)){ as.character(v$NextITS$version) } else { "" }
  glue("Bioinformatic processing was performed using the \\
    NextITS pipeline v.{nextits_v} (Mikryukov et al., 2025).")
}

emit_nextflow <- function(v) {
  nextflow_v <- if(!is.null(v$Nextflow$version)){ as.character(v$Nextflow$version) } else { "" }
  glue("Workflow management was performed using \\
    Nextflow v.{nextflow_v} (Di Tommaso et al., 2017).")
}

emit_aggregation <- function(p, v) {
  glue("Sequences from all sequencing runs were aggregated and \\
    de novo chimeric sequences with chimera score >= {getp(p,'max_ChimeraScore',0.6)} were removed.")
}

emit_dereplication <- function(p, v) {
  minlen <- getp(p, "ampliconlen_min", NA)
  maxlen <- getp(p, "ampliconlen_max", NA)
  
  length_filter <- ""
  if(!is.na(minlen) && !is.na(maxlen)){
    length_filter <- glue(" Sequences shorter than {minlen} nt or longer than {maxlen} nt were excluded.")
  } else if(!is.na(minlen)){
    length_filter <- glue(" Sequences shorter than {minlen} nt were excluded.")
  } else if(!is.na(maxlen)){
    length_filter <- glue(" Sequences longer than {maxlen} nt were excluded.")
  }
  
  glue("Global sequence dereplication was performed using \\
    VSEARCH v.{getv(v,'dereplication','vsearch')} (Rognes et al., 2016).\\
    {length_filter}")
}

emit_preclustering <- function(p, v) {
  preclustering_method <- getp(p, "preclustering", "none")
  
  res <- switch(preclustering_method,
    
    "none" = "", # No pre-clustering or denoising was performed
    
    "homopolymer" = glue(
        "Global homopolymer correction was performed using an algorithm implemented in NextITS \\
        with support of VSEARCH v.{getv(v,'homopolymer','vsearch')} (Rognes et al., 2016)."),

    "unoise" = glue(
        "Sequence denoising was performed using the UNOISE3 algorithm (Edgar, 2016) \\
        implemented in VSEARCH v.{getv(v,'unoise','vsearch')} (Rognes et al., 2016) \\
        with alpha parameter {getp(p,'unoise_alpha',6.0)} and minimum size {getp(p,'unoise_minsize',1)}."),

    "dada2" = glue(
        "Sequence denoising was performed using \\
        DADA2 v.{getv(v,'dada2','dada2')} (Callahan et al., 2016)"),
        # using {getp(p,'dada2_pooling','global')} pooling strategy."

    "swarm_d1" = glue(
        "Pre-clustering was performed using \\
        SWARM v.{getv(v,'precluster_swarm','swarm')} (Mahé et al., 2021) \\
        with d=1 and fastidious option enabled.")
    )
    
  return(res)
}



##################################
################################## Workflow-dependent method descriptions
##################################

## Function to assembly the workflow description and references
build_docs <- function(versions, params){
  body <- character()
  tools_used <- character()

  ## Pipeline version
  body <- c(body, emit_nextits(versions))
  tools_used <- c(tools_used, "nextits")

  ## Nextflow version
  body <- c(body, emit_nextflow(versions))
  tools_used <- c(tools_used, "nextflow")

  ## Sequence aggregation
  body <- c(body, emit_aggregation(params, versions))

  ## Sequence dereplication and amplicon length filtering
  body <- c(body, emit_dereplication(params, versions))
  tools_used <- c(tools_used, "vsearch")

  ## Conditional: pre-clustering/denoising
  preclustering_method <- getp(params, "preclustering", "none")
  if(preclustering_method != "none" && !is.na(preclustering_method)){
    body <- c(body, emit_preclustering(params, versions))
    
    switch(preclustering_method,
      "homopolymer" = {
        tools_used <- c(tools_used, "vsearch")
      },
      "unoise" = {
        tools_used <- c(tools_used, c("vsearch", "unoise"))
      },
      "dada2" = {
        tools_used <- c(tools_used, "dada2")
      },
      "swarm_d1" = {
        tools_used <- c(tools_used, "swarm")
      }
    )
  }


  ## Generate citations
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

## Build body and citations
res <- build_docs(versions, params)

## Write output
con <- file(output_path, open = "wt")

writeLines("Methods:", con)
writeLines(res$body, con)

writeLines("", con)

writeLines("References:", con)
writeLines(paste0("- ", res$citations), con)

close(con)

cat("All done.\n")

