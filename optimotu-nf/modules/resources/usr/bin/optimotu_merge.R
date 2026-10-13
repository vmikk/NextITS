#!/usr/bin/env Rscript

## Script: optimotu_merge.R
## Purpose: Concatenate per-component OptimOTU results into one set of tables
##
##          OTUs are identified by their representative sequence IDs, which are
##          unique across components, so only the UC cluster numbers need to be
##          shifted into one global series
##
##          This is only sound when every input was a complete connected
##          component at a distance cutoff no stricter than the thresholds'
##          maximum clustering distance -- see optimotu_split_components.R --
##          or a part of a taxonomic split (optimotu_cluster.R --split-rank)
##          Otherwise the parts must be clustered together instead
##
##          A part can be empty (a split whose first pass finished nothing);
##          its UC and FASTA files are then empty and skipped
##
## Usage:
##   Rscript optimotu_merge.R \
##      --membership compA.membership.tsv,compB.membership.tsv \
##      --otu-taxonomy compA.otu_taxonomy.tsv,compB.otu_taxonomy.tsv \
##      --otus-fasta compA.otus.fasta,compB.otus.fasta \
##      --uc compA.uc,compB.uc \
##      --output-prefix all
##
## Output:
##   - <prefix>.membership.tsv
##   - <prefix>.otu_taxonomy.tsv    (only with --otu-taxonomy)
##   - <prefix>.uc                  (only with --uc)
##   - <prefix>.otus.fasta          (only with --otus-fasta)

options(stringsAsFactors = FALSE, warn = 1)

## Check time
start_time <- Sys.time()

script_dir <- dirname(normalizePath(sub(
  "^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]]
)))
source(file.path(script_dir, "optimotu_common.R"))
SCRIPT_START <- start_time

banner("optimotu_merge.R", "Merge per-component OptimOTU results")

cat("Parsing input options and arguments...\n")

suppressPackageStartupMessages(require(optparse))

## Parse arguments
option_list <- list(
  make_option(c("--membership"),   default=NULL,  type='character', help="Comma-separated membership TSV files [required]"),
  make_option(c("--output-prefix"), default=NULL, type='character', help="Prefix for the merged outputs [required]"),
  make_option(c("--otu-taxonomy"), default=NULL,  type='character', help="Matching otu_taxonomy TSV files, comma-separated"),
  make_option(c("--otus-fasta"),   default=NULL,  type='character', help="Matching representative FASTA files, comma-separated"),
  make_option(c("--uc"),           default=NULL,  type='character', help="Matching .uc files, comma-separated"),
  make_option(c("--overwrite"),    action='store_true', default=FALSE, help="Replace existing output files")
)

opt <- parse_args(OptionParser(
  option_list = option_list,
  usage = "%prog --membership FILES --output-prefix PATH [options]",
  description = "\nMerge per-component OptimOTU results.\n",
  epilogue = paste0(
    "Every input must have been a complete connected component: the merge only\n",
    "concatenates and renumbers, it never re-clusters across component borders.\n",
    "Overlapping inputs are rejected rather than silently deduplicated.\n")
))

## Validation of required arguments
require_options(opt, c("membership", "output-prefix"))

############################################## Load packages

cat("\nLoading R packages...\n")

load_pckg("data.table")
require_package("Biostrings", "2.70.0")

############################################## Assign and validate parameters

prefix <- path.expand(opt$`output-prefix`)

read_all <- function(value, label) {
  paths <- as_list_arg(value)
  if (!length(paths)) return(character())
  missing_paths <- paths[!file.exists(paths)]
  if (length(missing_paths)) {
    fail(label, " file(s) not readable: ", paste(missing_paths, collapse = ", "))
  }
  normalizePath(paths)
}

membership_paths <- read_all(opt$membership, "Membership")
if (!length(membership_paths)) fail("--membership listed no files")
otu_taxonomy_paths <- read_all(opt$`otu-taxonomy`, "OTU taxonomy")
uc_paths <- read_all(opt$uc, "UC")
fasta_paths <- read_all(opt$`otus-fasta`, "Representative FASTA")

## Log assigned variables
cat("\nParameters specified:\n")
report_params(list(
  "Membership files"   = length(membership_paths),
  "OTU taxonomy files" = length(otu_taxonomy_paths),
  "UC files"           = length(uc_paths),
  "Representative FASTA files" = length(fasta_paths),
  "Output prefix"      = prefix,
  "Overwrite"          = isTRUE(opt$overwrite)
))
cat("\n")

############################################## Helpers

## Read and stack a set of TSV files with identical columns
##
## `fread` + `rbindlist` instead of `lapply(read.delim)` + `do.call(rbind)`:
## rbind() reallocates the whole result once per file
## na.strings = "NA" reproduces read.delim's default exactly
bind_tsv <- function(paths, label) {
  parts <- lapply(paths, function(p) {
    data.table::fread(p, sep = "\t", header = TRUE, quote = "",
                      colClasses = "character", na.strings = "NA",
                      showProgress = FALSE, data.table = FALSE)
  })
  columns <- lapply(parts, names)
  if (length(unique(lapply(columns, sort))) != 1L) {
    fail(label, " files do not share the same columns")
  }
  out <- data.table::rbindlist(parts, use.names = TRUE)
  data.table::setcolorder(out, columns[[1L]])
  data.table::setDF(out)
  out
}

############################################## Main workflow

cat("\n--- Main Workflow ---\n\n")

step("Reading ", length(membership_paths), " membership table(s) ...")
membership <- bind_tsv(membership_paths, "Membership")
if (!all(c("sequence_id", "representative_id") %in% names(membership))) {
  fail("Membership files must have sequence_id and representative_id columns")
}
note("  ", nrow(membership), " membership rows")

duplicates <- unique(membership$sequence_id[duplicated(membership$sequence_id)])
if (length(duplicates)) {
  fail("The same sequence appears in more than one component (e.g. ",
       paste(utils::head(duplicates, 5L), collapse = ", "),
       "); the inputs were not disjoint")
}
representative_ids <- unique(membership$representative_id)
unknown <- setdiff(representative_ids, membership$sequence_id)
if (length(unknown)) {
  fail("Representatives that are not member sequences (e.g. ",
       paste(utils::head(unknown, 5L), collapse = ", "), ")")
}
n_otus <- length(representative_ids)
note("Merged ", length(membership_paths), " components: ",
     nrow(membership), " sequences in ", n_otus, " OTUs")

outputs <- paste0(prefix, ".membership.tsv")
payload <- list(membership)

if (length(otu_taxonomy_paths)) {
  step("Merging ", length(otu_taxonomy_paths), " OTU taxonomy table(s) ...")
  otu_taxonomy <- bind_tsv(otu_taxonomy_paths, "OTU taxonomy")
  if (!"representative_id" %in% names(otu_taxonomy)) {
    fail("OTU taxonomy files must have a representative_id column")
  }
  if (anyDuplicated(otu_taxonomy$representative_id)) {
    fail("OTU taxonomy tables name the same representative more than once")
  }
  if (!setequal(otu_taxonomy$representative_id, representative_ids)) {
    fail("OTU taxonomy tables do not cover exactly the OTUs in the membership tables")
  }
  outputs <- c(outputs, paste0(prefix, ".otu_taxonomy.tsv"))
  payload <- c(payload, list(otu_taxonomy))
  note("  ", nrow(otu_taxonomy), " OTU taxonomy rows")
}

merged_uc <- NULL
if (length(uc_paths)) {
  step("Merging ", length(uc_paths), " UC table(s) ...")
  # Cluster numbers are 0-based per component, so each file's are shifted past
  # the clusters of the files before it
  offset <- 0L
  parts <- vector("list", length(uc_paths))
  for (i in seq_along(uc_paths)) {
    if (!file.size(uc_paths[[i]])) next
    part <- data.table::fread(uc_paths[[i]], sep = "\t", header = FALSE, quote = "",
                              colClasses = "character", na.strings = "NA",
                              showProgress = FALSE, data.table = FALSE)
    cluster <- as.integer(part[[2L]])
    if (anyNA(cluster)) fail("Non-numeric cluster number in ", uc_paths[[i]])
    part[[2L]] <- as.character(cluster + offset)
    offset <- offset + max(cluster) + 1L
    parts[[i]] <- part
  }
  merged_uc <- as.data.frame(data.table::rbindlist(parts, use.names = FALSE))
  hit_rows <- merged_uc[[1L]] != "C"
  if (!setequal(merged_uc[[9L]][hit_rows], membership$sequence_id) ||
      anyDuplicated(merged_uc[[9L]][hit_rows])) {
    fail("UC files do not list every member sequence exactly once")
  }
  # S/H records first, then the C summary records, as VSEARCH writes them
  merged_uc <- merged_uc[order(!hit_rows), , drop = FALSE]
  outputs <- c(outputs, paste0(prefix, ".uc"))
  note("  ", nrow(merged_uc), " UC records")
}

merged_fasta <- NULL
if (length(fasta_paths)) {
  step("Merging ", length(fasta_paths), " representative FASTA file(s) ...")
  merged_fasta <- do.call(c, lapply(fasta_paths[file.size(fasta_paths) > 0], function(p) {
    Biostrings::readDNAStringSet(p, format = "fasta")
  }))
  ids <- strip_seq_ids(names(merged_fasta))
  if (anyDuplicated(ids) || !setequal(ids, representative_ids)) {
    fail("Representative FASTA files do not hold exactly one sequence per OTU")
  }
  outputs <- c(outputs, paste0(prefix, ".otus.fasta"))
  note("  ", length(merged_fasta), " representative sequences")
}

step("Writing ", length(outputs), " output file(s) ...")

with_atomic_outputs(outputs, overwrite = isTRUE(opt$overwrite), writer = function(staged) {
  i <- 0L
  for (table in payload) {
    i <- i + 1L
    data.table::fwrite(table, staged[[i]], sep = "\t", quote = FALSE, na = "")
  }
  if (!is.null(merged_uc)) {
    i <- i + 1L
    data.table::fwrite(merged_uc, staged[[i]], sep = "\t", quote = FALSE,
                       col.names = FALSE, na = "*")
  }
  if (!is.null(merged_fasta)) {
    i <- i + 1L
    Biostrings::writeXStringSet(merged_fasta, staged[[i]], format = "fasta",
                                width = 80L)
  }
})

note("Wrote ", paste(basename(outputs), collapse = ", "), " under ", dirname(prefix))

############################################## Session info

finish(start_time)
