#!/usr/bin/env Rscript

## Script: optimotu_split_components.R
## Purpose: Split a FASTA into connected components of a sequence-similarity
##          graph, so that OptimOTU can be run on the parts in parallel without
##          changing the result
##
##          The invariant: two sequences may only end up in different files if no
##          chain of pairwise distances below the largest threshold OptimOTU might
##          use connects them.  The graph is therefore built at a cutoff at least
##          as permissive as the thresholds' `max_clustering_distance`, and the
##          script refuses to run otherwise
##
## Usage:
##   Rscript optimotu_split_components.R \
##      --input-fasta query.fasta \
##      --output-dir components \
##      --max-dist 0.4 \
##      --threads 4
##
## Output:
##   - <output-dir>/<prefix>_NNNN.fasta   one file per connected component
##   - <output-dir>/components.tsv        manifest: component_id, n_sequences

options(stringsAsFactors = FALSE, warn = 1)

## Check time
start_time <- Sys.time()

script_dir <- dirname(normalizePath(sub(
  "^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]]
)))
source(file.path(script_dir, "optimotu_common.R"))
SCRIPT_START <- start_time

banner("optimotu_split_components.R",
       "Split a FASTA into connected components of a similarity graph")

cat("Parsing input options and arguments...\n")

suppressPackageStartupMessages(require(optparse))

## Parse arguments
option_list <- list(
  make_option(c("--input-fasta"),  default=NULL,  type='character', help="Sequences to split [required]"),
  make_option(c("--output-dir"),   default=NULL,  type='character', help="Destination for the per-component FASTA files [required]"),
  make_option(c("--thresholds"),   default=NULL,  type='character', help="Thresholds TSV; its maximum clustering distance is checked against --max-dist"),
  make_option(c("--max-dist"),     default=0.4,   type='double',    help="Graph distance cutoff [default %default]"),
  make_option(c("--min-identity"), default=NULL,  type='double',    help="Equivalent VSEARCH identity [default: 1 - max-dist]"),
  make_option(c("--prefix"),       default="component", type='character', help="Component file prefix [default %default]"),
  make_option(c("--threads"),      default=1L,    type='integer',   help="VSEARCH threads [default %default]"),
  make_option(c("--vsearch"),      default=NULL,  type='character', help="VSEARCH executable [default: autodetected from $VSEARCH, then PATH]"),
  make_option(c("--force"),        action='store_true', default=FALSE, help="Split even if the cutoff is too strict"),
  make_option(c("--overwrite"),    action='store_true', default=FALSE, help="Replace an existing output directory")
)

opt <- parse_args(OptionParser(
  option_list = option_list,
  usage = "%prog --input-fasta FILE --output-dir DIR [options]",
  description = "\nSplit a FASTA into connected components of a similarity graph.\n",
  epilogue = paste0(
    "The graph is built with --allpairs_global, which aligns every pair, so no\n",
    "near-threshold edge is missed -- but that also makes it quadratic. At the\n",
    "identities OptimOTU's coarser thresholds imply (often 60-70%) most sequences\n",
    "tend to land in one component anyway, so leave splitting off unless a trial\n",
    "run shows the input really does break apart.\n")
))

## Validation of required arguments
require_options(opt, c("input-fasta", "output-dir"))

############################################## Load packages

cat("\nLoading R packages...\n")

load_pckg("data.table")
require_package("optimotu", OPTIMOTU_MIN_VERSION)
require_package("Biostrings", "2.70.0")
require_package("processx", "3.8.0")

############################################## Assign and validate parameters

max_dist <- as_probability(opt$`max-dist`, "max-dist", 0.4)
min_identity <- as_probability(opt$`min-identity`, "min-identity", 1 - max_dist)
prefix <- sanitize_id(opt$prefix, "--prefix")
threads <- as_count(opt$threads, "threads", 1L)
vsearch <- find_vsearch(opt$vsearch)

input_path <- normalizePath(opt$`input-fasta`, mustWork = TRUE)
output_dir <- path.expand(opt$`output-dir`)
thresholds_path <- if (is.null(opt$thresholds)) NULL else {
  normalizePath(opt$thresholds, mustWork = TRUE)
}

## Log assigned variables
cat("\nParameters specified:\n")
report_params(list(
  "Input FASTA"    = input_path,
  "Output dir"     = output_dir,
  "Thresholds"     = thresholds_path,
  "Max distance"   = max_dist,
  "Min identity"   = min_identity,
  "Component prefix" = prefix,
  "Threads"        = threads,
  "Force"          = isTRUE(opt$force),
  "Overwrite"      = isTRUE(opt$overwrite)
))
report_tool("vsearch", vsearch)
cat("\n")

############################################## Helpers

## Connected components by vectorized label propagation
##
## Each vertex starts labelled with its own index; every round pushes the
## smallest label across each edge and then jumps every pointer to its current
## target, which converges in O(log n) rounds without a per-edge R loop
##
## The per-vertex minimum is a data.table grouped aggregation rather than
## `tapply()`: tapply built a factor over the vertex ids and returned the result
## keyed by their *character* representation, so every round paid for a
## character conversion and an as.integer() round-trip back
connected_components <- function(n, from, to) {
  label <- seq_len(n)
  if (!length(from)) return(label)
  vertex <- c(from, to)
  repeat {
    edge_min <- pmin(label[from], label[to])
    proposed <- data.table::data.table(
      vertex = vertex, value = c(edge_min, edge_min)
    )[, list(value = min(value)), by = "vertex"]
    candidate <- label
    candidate[proposed$vertex] <- pmin(candidate[proposed$vertex], proposed$value)
    candidate <- candidate[candidate]
    if (identical(candidate, label)) break
    label <- candidate
  }
  label
}

############################################## Main workflow

cat("\n--- Main Workflow ---\n\n")

if (!is.null(thresholds_path)) {
  step("Checking --max-dist against the thresholds ...")
  thresholds <- data.table::fread(thresholds_path, sep = "\t", header = TRUE,
                                  quote = "", na.strings = "NA",
                                  showProgress = FALSE, data.table = FALSE)
  if (!"threshold" %in% names(thresholds)) {
    fail("Thresholds file has no 'threshold' column")
  }
  needed <- max(optimotu::threshold_as_dist(
    suppressWarnings(as.numeric(thresholds$threshold))
  ), na.rm = TRUE)
  note("  thresholds' maximum clustering distance: ", format(needed, digits = 6))
  if (max_dist < needed) {
    msg <- paste0("--max-dist ", max_dist, " is stricter than the thresholds' ",
                  "maximum clustering distance ", format(needed, digits = 6),
                  "; components would not be independent")
    if (isTRUE(opt$force)) note("WARNING: ", msg, " (--force given)") else {
      fail(msg, ". Raise --max-dist or pass --force.")
    }
  }
}

step("Reading ", input_path, " ...")
seqs <- read_fasta(input_path, "Input FASTA")
ids <- names(seqs)
n <- length(ids)
note(n, " sequences; building the similarity graph at identity >= ",
     format(min_identity, digits = 4), " ...")

## --- edges ------------------------------------------------------------------

query_file <- tempfile(fileext = ".fasta")
hits_file <- tempfile(fileext = ".b6")
on.exit(unlink(c(query_file, hits_file), force = TRUE), add = TRUE)
Biostrings::writeXStringSet(seqs, query_file, format = "fasta", width = 80L)

step("Running VSEARCH --allpairs_global on ", threads, " thread(s) ...")
# --allpairs_global aligns every pair rather than using a k-mer prefilter,
# so --id only filters the output and no near-threshold edge can be missed.
# --maxaccepts/--maxrejects 0 additionally disables the per-query hit caps.
result <- processx::run(
  vsearch,
  c("--allpairs_global", query_file,
    "--id", format(min_identity, scientific = FALSE),
    "--blast6out", hits_file,
    "--maxaccepts", "0", "--maxrejects", "0",
    "--threads", as.character(threads), "--quiet"),
  error_on_status = FALSE,
  echo = FALSE
)
if (result$status != 0L) fail("VSEARCH all-pairs search failed: ", result$stderr)

## --- connected components ---------------------------------------------------

step("Reading the all-pairs hit table ...")
n_edges <- 0L
from <- integer()
to <- integer()
if (file.exists(hits_file) && file.size(hits_file) > 0L) {
  # the all-pairs table is the quadratic-sized intermediate here, and only its
  # first two columns are edges, so read just those
  hits <- data.table::fread(hits_file, sep = "\t", header = FALSE, quote = "",
                            select = 1:2, colClasses = "character",
                            na.strings = "NA", showProgress = FALSE,
                            data.table = FALSE)
  from <- match(strip_seq_ids(hits[[1L]]), ids)
  to <- match(strip_seq_ids(hits[[2L]]), ids)
  keep <- !is.na(from) & !is.na(to) & from != to
  from <- from[keep]
  to <- to[keep]
  n_edges <- length(from)
}
note("  ", n_edges, " edges")

step("Finding connected components ...")
component <- connected_components(n, from, to)
component <- match(component, unique(component))
n_components <- max(component)
sizes <- tabulate(component, nbins = n_components)
note(n_edges, " edges; ", n_components, " components (largest ", max(sizes),
     " sequences, ", sum(sizes == 1L), " singletons)")

## --- write ------------------------------------------------------------------

step("Writing ", n_components, " component FASTA file(s) ...")

if (dir.exists(output_dir)) {
  if (!isTRUE(opt$overwrite)) {
    fail("Output directory already exists; use --overwrite: ", output_dir)
  }
  unlink(list.files(output_dir, pattern = paste0("^", prefix, "_"),
                    full.names = TRUE), force = TRUE)
} else {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(output_dir)) fail("Cannot create output directory: ", output_dir)
}

width <- max(4L, nchar(n_components))
manifest <- data.frame(
  component_id = sprintf("%s_%0*d", prefix, width, seq_len(n_components)),
  n_sequences = sizes
)
for (k in seq_len(n_components)) {
  Biostrings::writeXStringSet(
    seqs[component == k],
    file.path(output_dir, paste0(manifest$component_id[[k]], ".fasta")),
    format = "fasta", width = 80L
  )
}
data.table::fwrite(manifest, file.path(output_dir, "components.tsv"), sep = "\t",
                   quote = FALSE)

note("Wrote ", n_components, " component FASTA files to ", output_dir)

############################################## Session info

finish(start_time)
