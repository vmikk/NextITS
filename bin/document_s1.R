#!/usr/bin/env Rscript

## Script to document the Step-1 workflow of the NextITS pipeline.

## Usage:
##   Rscript document_s1.R <software_versions.yml> <params.csv_or_tsv> [output_path]

## Input:
##  - software_versions.yml
##  - pipeline_params.tsv

## Output:
##  - README_Step1_Methods.txt
##    with two sections: methods and references


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
################################## Data for debugging
##################################

# versions_path <- "software_versions.yml"
# params_path   <- "pipeline_params.tsv"
# output_path   <- "README_Step1_Methods.txt"


##################################
################################## References
##################################

## Citation registry
citation_db <- list(
  nextits   = "Mikryukov V, Anslan S, Tedersoo L (2025) NextITS - A pipeline for metabarcoding fungi and other eukaryotes with full-length ITS sequenced with PacBio. DOI:10.5281/zenodo.15074882",
  nextflow  = "Di Tommaso P, et al. (2017) Nextflow enables reproducible computational workflows. Nat Biotechnol 35, 316-319, DOI:10.1038/nbt.3820",
  lima      = "Pacific Biosciences (2025) LIMA - The PacBio barcode demultiplexer and primer remover. URL: https://lima.how/",
  seqkit    = "Shen W, Sipos B, Zhao L (2024) SeqKit2: A Swiss Army Knife for Sequence and Alignment Processing. iMeta e191. DOI:10.1002/imt2.191",
  csvtk     = "Shen W (2025) csvtk - a cross-platform, efficient and practical CSV/TSV toolkit. URL: https://github.com/shenwei356/csvtk",
# brename   = "Shen W (2025) brename - batch renaming safely, URL: https://github.com/shenwei356/brename",
  cutadapt  = "Martin M (2011) Cutadapt removes adapter sequences. EMBnet.journal 17(1):10-12, DOI:10.14806/ej.17.1.200",
  usearch   = "Edgar RC (2010) Search and clustering orders of magnitude faster than BLAST. Bioinformatics 26(19):2460-2461. DOI:10.1093/bioinformatics/btq461",
  pemerge   = "Edgar RC, Flyvbjerg H (2015) Error filtering, pair assembly and error correction for next-generation sequencing reads. Bioinformatics 31(21):3476-3482. DOI:10.1093/bioinformatics/btv401",
  fastp     = "Chen S, Zhou Y, Chen Y, Gu J (2018) fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics 34(17):i884-i890. DOI:10.1093/bioinformatics/bty560",
  itsx      = "Bengtsson-Palme J, et al (2013) Improved software detection and extraction of ITS1 and ITS2 from ribosomal ITS sequences of fungi and other eukaryotes for analysis of environmental sequencing data. Methods Ecol Evol 4:914-919, DOI:10.1111/2041-210X.12073",
  itsx2     = "ITSx2: structure-aware delimitation of the nuclear rRNA cistron. URL: https://github.com/Changwanseo/pre-itsx2",
  infernal  = "Nawrocki EP, Eddy SR (2013) Infernal 1.1: 100-fold faster RNA homology searches. Bioinformatics 29(22):2933-2935. DOI:10.1093/bioinformatics/btt509",
  rfam      = "Ontiveros-Palacios N, et al. (2025) Rfam 15: RNA families database in 2025. Nucleic Acids Res 53(D1):D258-D267. DOI:10.1093/nar/gkae1023",
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

  if(is.null(v[[process]]) || is.null(v[[process]][[tool]])){ return("") }
  as.character( v[[process]][[tool]] )
}
# E.g., getv(versions, "demux", "lima")


## Get parameter
getp <- function(p, pname, default = NA){
  # p       =  table with parameters (two columns: name and value)
  # pname   = parameter name
  # default = default value if parameter is not found

  pp <- p[ name == pname ]$value
  if(is.null(pp) || is.na(pp)){ return(default) }
  return(pp)
}
# E.g., getp(params, "lima_minscore", 93)


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

emit_qc_pacbio <- function(p, v) {
  glue("Quality control was performed using \\
    VSEARCH v.{getv(v,'qc_se','vsearch')} (Rognes et al., 2016) and \\
    seqkit v.{getv(v,'qc_se','seqkit')} (Shen et al., 2024). \\
    Reads with the number of ambiguous bases >= {getp(p,'qc_maxn',4)}, \\
    expected error rate >= {getp(p,'qc_maxeerate',0.01)}, \\
    or homopolymer stretches longer than {getp(p,'qc_maxhomopolymerlen',25)} nt were removed.")
}

is_true <- function(x){ tolower(as.character(x)) %in% c("true", "t", "1") }
is_set  <- function(x){ !is.na(x) && nzchar(x) && !tolower(as.character(x)) %in% c("null", "false") }

emit_demux_illumina <- function(p, v) {
  bt <- getp(p, "illumina_barcodetype", "dual_symmetric")
  vs <- getv(v, "demux_illumina", "cutadapt")
  common <- glue("allowing up to {getp(p,'barcode_errors',2)} mismatches (no indels) \\
    within the first {getp(p,'barcode_window',30)} nt of the reads")
  switch(bt,
    "single" = glue("Paired-end reads were demultiplexed using cutadapt v.{vs} (Martin, 2011) \\
      based on a sample tag at the 5' end of either mate, {common}."),
    "dual_symmetric" = glue("Paired-end reads were demultiplexed using cutadapt v.{vs} (Martin, 2011), \\
      {common}. Read pairs were assigned to samples only if both mates carried the tag of the same sample; \\
      pairs in which the mates carried tags of different samples were discarded as tag jumps.\\
      {if(is_true(getp(p,'illumina_demux_rescue',TRUE))) ' Pairs with a single readable tag were assigned based on that tag.' else ''}"),
    "dual_asymmetric" = glue("Paired-end reads were demultiplexed using cutadapt v.{vs} (Martin, 2011) \\
      based on combinations of forward and reverse tags (in both mate orientations), {common}. \\
      Pairs carrying an unexpected combination of known tags were discarded as tag jumps.\\
      {if(is_true(getp(p,'illumina_demux_rescue',TRUE))) ' Pairs with a single readable tag were assigned based on that tag, if the tag was used by a single sample.' else ''}"),
    glue("Paired-end reads were demultiplexed using cutadapt v.{vs} (Martin, 2011)."))
}

emit_merge_illumina <- function(p, v) {
  merger <- getp(p, "pe_merger", "usearch")
  mv <- getv(v, "merge_pe", merger)
  res <- glue("Read pairs were oriented based on primer sequences using \\
    cutadapt v.{getv(v,'reorient_pe','cutadapt')}.")
  if(is_set(getp(p, "qc_polyglen", NA))){
    res <- c(res, glue("Poly-G tails (>= {getp(p,'qc_polyglen',10)} nt) were trimmed using \\
      fastp v.{getv(v,'merge_pe','fastp')} (Chen et al., 2018)."))
  }
  if(merger %in% "usearch"){
    res <- c(res, glue("Read pairs were merged using USEARCH v.{mv} (Edgar, 2010; Edgar & Flyvbjerg, 2015) \\
      with a minimum overlap of {getp(p,'pe_minoverlap',16)} nt, \\
      at most {getp(p,'pe_maxdiffs',10)} mismatches and >= {getp(p,'pe_pctid',80)}% identity in the overlap, \\
      and a minimum merged length of {getp(p,'pe_minlen',50)} nt."))
  } else {
    res <- c(res, glue("Read pairs were merged using VSEARCH v.{mv} (Rognes et al., 2016) \\
      with a minimum overlap of {getp(p,'pe_minoverlap',16)} nt, \\
      at most {getp(p,'pe_maxdiffs',10)} mismatches and >= {getp(p,'pe_pctid',80)}% identity in the overlap, \\
      and a minimum merged length of {getp(p,'pe_minlen',50)} nt."))
  }
  if(is_true(getp(p, "illumina_keep_notmerged", FALSE))){
    res <- c(res, glue("Non-merged read pairs were quality-filtered and joined with a padding of \\
      '{getp(p,'illumina_joinpadgap','NNNNNNNNNN')}' using VSEARCH v.{getv(v,'join_pe','vsearch')}."))
  }
  paste0(res, collapse = " ")
}

emit_qc_illumina <- function(p, v) {
  maxee <- getp(p, "qc_maxee", NA)
  ee <- if(is_set(maxee)) glue("more than {maxee} expected errors, ") else ""
  glue("Quality control of merged reads was performed using \\
    VSEARCH v.{getv(v,'qc_se','vsearch')} (Rognes et al., 2016) and \\
    seqkit v.{getv(v,'qc_se','seqkit')} (Shen et al., 2024). \\
    Reads with {ee}expected error rate >= {getp(p,'qc_maxeerate',0.01)}, \\
    the number of ambiguous bases >= {getp(p,'qc_maxn',4)}, \\
    or homopolymer stretches longer than {getp(p,'qc_maxhomopolymerlen',25)} nt were removed.")
}

emit_primer_check <- function(p, v) {
  glue("Primers were trimmed using \\
    cutadapt v.{getv(v, 'primer_check', 'cutadapt')} (Martin, 2011) \\
    with <= {getp(p, 'primer_mismatches', 2)} mismatches. \\
    Reads without both primers were discarded.")
}

emit_long_chimera <- function(p, v) {
  minbits <- getp(p, 'long_chimera_minbits', 50)
  mincov  <- suppressWarnings(as.numeric(getp(p, 'long_chimera_mincov', 0.4)))
  mincov  <- if(is.na(mincov)) 40 else round(mincov * 100, 1)
  glue("Long chimeras were removed from the dereplicated sequences of each sample. \\
    Sequences were searched on both strands using cmsearch from \\
    Infernal v.{getv(v,'long_chimera','Infernal')} (Nawrocki & Eddy, 2013) in HMM-only mode \\
    against the covariance models of eukaryotic 5.8S rRNA (Rfam RF00002), \\
    the 3' end of SSU rRNA (RF01960) and the 5' end of LSU rRNA (RF02543) \\
    (Ontiveros-Palacios et al., 2025), as provided with ITSx2. \\
    Hits with bit score >= {minbits} that covered >= {mincov}% of the model length were considered. \\
    Sequences with such hits on both strands (inverted duplicates and chimeras of amplicons \\
    joined in opposite orientation) were discarded, as were sequences with two non-overlapping hits \\
    of the same model on the same strand that shared >= {mincov}% of the model (repeated rRNA regions).")
}

emit_itsx <- function(p, v) {
  switch(getp(p,'its_region','full'),
    "full"            = {its_region <- "full-length ITS"},
    "SSU"             = {its_region <- "SSU"},
    "ITS1"            = {its_region <- "ITS1"},
    # "5_8S"            = {its_region <- "5.8S"},   # not-yet-implemented
    "ITS2"            = {its_region <- "ITS2"},
    "LSU"             = {its_region <- "LSU"},
    "ITS1_5.8S_ITS2"  = {its_region <- "near-full-length ITS"})


  ## Either ITSx v1.x (HMMER profiles) or ITSx2 (Infernal covariance models)
  if(identical(getp(p, 'itsx_tool', 'ITSx'), 'ITSx2')){
    glue("Extraction of rRNA regions ({its_region}) was performed using \\
      ITSx2 v.{getv(v,'itsx2','ITSx2')}, which delimits the rRNA cistron \\
      pan-eukaryotically with covariance models using \\
      Infernal v.{getv(v,'itsx2','Infernal')} (Nawrocki & Eddy, 2013).")
  } else {
    glue("Extraction of rRNA regions ({its_region}) was performed using \\
      ITSx v.{getv(v,'itsx','ITSx')} (Bengtsson-Palme et al., 2013) \\
      with the '{getp(p, 'ITSx_tax', 'all')}' taxonomic profile(s).")
  }
}

emit_assemble_its <- function(p, v) {
  glue("To assemble near-full-length ITS sequences, we ... (TODO)")
}

emit_hp_and_chimeras <- function(p, v, did_hp) {
  res <- character()
  if(isTRUE(did_hp)){
    res <- c(res, glue(
      "Homopolymer correction of sequences was performed using an algorithm implemented in NextITS \\
      with support of VSEARCH v.{getv(v,'homopolymer','vsearch')} and seqkit v.{getv(v,'homopolymer','seqkit')}.") )
  } else {
    res <- c(res, "Homopolymer correction of sequences was not performed.")
  }
  res <- c(res,
    glue(
      "Two-step chimera detection was done using VSEARCH v.{getv(v,'chimera_denovo','vsearch')}: 
       - de novo using UCHIME2 algorithm (Edgar, 2016) with max score {getp(p,'max_ChimeraScore',0.6)} (Nilsson et al., 2015), 
       - then reference-based against the EUKARYOME database (Tedersoo et al., 2024).")
  )
  res <- paste0(res, collapse = "\n")
  return(res)
}

emit_tj <- function(p, v) {
  glue("Tag-jump detection and removal was performed using \\
    UNCROSS2 algorithm (Edgar, 2018) with the parameter f = {getp(p,'tj_f',0.01)}.")
}

emit_seqtab <- function(p, v) {
  glue("Sequence counts table was generated using \\
    R v.{getv(v,'prep_seqtab','R')} (R Core Team, 2025), \\
    data.table v.{getv(v,'prep_seqtab','data.table')} (Barrett et al., 2025), \\
    and Apache Arrow v.{getv(v,'prep_seqtab','arrow')} (Richardson et al., 2025) \\
    packages.")
}



##################################
################################## Workflow-dependent method descriptions
##################################

## Function to assembly the workflow description and references
build_docs <- function(versions, params){
  body <- character()
  tools_used <- character()

  body <- c(body, emit_nextits(versions))
  tools_used <- c(tools_used, "nextits")

  body <- c(body, emit_nextflow(versions))
  tools_used <- c(tools_used, "nextflow")

  demuxed <- tolower(as.character(getp(params, "demultiplexed", FALSE))) %in% c("true", "t", "1")
  platform <- getp(params, "seqplatform", "PacBio")

  if(platform %in% "Illumina"){

    ## Illumina paired-end reads: demultiplexing (if required), reorientation, merging, QC
    if(!demuxed){
      body <- c(body, emit_demux_illumina(params, versions))
      tools_used <- c(tools_used, c("cutadapt"))
    }

    body <- c(body, emit_merge_illumina(params, versions))
    merger <- getp(params, "pe_merger", "usearch")
    tools_used <- c(tools_used, c("cutadapt", if(merger %in% "usearch") c("usearch", "pemerge") else "vsearch"))
    if(is_set(getp(params, "qc_polyglen", NA))){
      tools_used <- c(tools_used, "fastp")
    }

    body <- c(body, emit_qc_illumina(params, versions))
    tools_used <- c(tools_used, c("vsearch", "seqkit"))

  } else if(!demuxed){

    body <- c(body, emit_demux_pacbio(params, versions))
    tools_used <- c(tools_used, c("lima"))

    body <- c(body, emit_qc_pacbio(params, versions))
    tools_used <- c(tools_used, c("vsearch", "seqkit"))

  } else {

    body <- c(body, emit_qc_pacbio(params, versions))
    tools_used <- c(tools_used, c("vsearch", "seqkit"))

  }

  ## Primer trimming
  body <- c(body, emit_primer_check(params, versions))
  tools_used <- c(tools_used, c("cutadapt"))

  ## Long-chimera removal (optional)
  if(is_true(getp(params, "long_chimera", FALSE))){
    body <- c(body, emit_long_chimera(params, versions))
    tools_used <- c(tools_used, c("infernal", "rfam", "itsx2"))
  }

  ## ITS extraction
  its_region <- getp(params, "its_region", "full")
  extractor  <- if(identical(getp(params, "itsx_tool", "ITSx"), "ITSx2")) c("itsx2", "infernal") else "itsx"
  if(its_region %in% c("full", "ITS1", "ITS2", "SSU", "LSU", "ITS1_5.8S_ITS2")){
    body <- c(body, emit_itsx(params, versions))
    tools_used <- c(tools_used, extractor, c("vsearch", "duckdb", "seqkit", "cutadapt"))
  }

  did_hp <- tolower(as.character(getp(params, "hp", TRUE))) %in% c("true", "t", "1")
  body <- c(body, emit_hp_and_chimeras(params, versions, did_hp))
  tools_used <- c(tools_used, c("vsearch", "uchime2", "eukaryome"))

  body <- c(body, emit_tj(params, versions))
  tools_used <- c(tools_used, c("uncross2"))

  body <- c(body, emit_seqtab(params, versions))
  tools_used <- c(tools_used, c("arrow", "datatable", "R"))

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
