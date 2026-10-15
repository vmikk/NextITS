#!/usr/bin/env Rscript

## Script: optimotu_merge_thresholds.R
## Purpose: Reduce the per-chunk outputs of a partitioned training run back into
##          the single thresholds table (plus JSON sidecar) that
##          optimotu_cluster.R expects
##
##          Because optimotu_partition_reference.R assigns every
##          (superrank, supertaxon) pair to exactly one chunk, merging is a
##          concatenation -- there is no precedence rule to get wrong.  What this
##          script does instead is check that the assumption held: keys must be
##          unique, the chunks must agree on the distance configuration the
##          thresholds are only valid for, and the merged table must satisfy the
##          root-fallback requirement that each chunk on its own could not
##
## Usage:
##   Rscript optimotu_merge_thresholds.R \
##      --thresholds chunk_001.thresholds.tsv,chunk_002.thresholds.tsv \
##      --sidecars chunk_001.thresholds.json,chunk_002.thresholds.json \
##      --expected-rows partition/expected_rows.tsv \
##      --output-prefix reference
##
## Output:
##   - <prefix>.thresholds.tsv
##   - <prefix>.thresholds.json
##   - <prefix>.thresholds.provenance.tsv

options(stringsAsFactors = FALSE, warn = 1)

## Check time
start_time <- Sys.time()

script_dir <- dirname(normalizePath(sub(
  "^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]]
)))
source(file.path(script_dir, "optimotu_common.R"))
SCRIPT_START <- start_time

banner("optimotu_merge_thresholds.R", "Merge per-chunk OptimOTU threshold tables")

cat("Parsing input options and arguments...\n")

suppressPackageStartupMessages(require(optparse))

## Parse arguments
option_list <- list(
  make_option(c("--thresholds"),    default=NULL, type='character', help="Per-chunk .thresholds.tsv files, comma-separated [required]"),
  make_option(c("--sidecars"),      default=NULL, type='character', help="Matching .thresholds.json files, comma-separated [required]"),
  make_option(c("--output-prefix"), default=NULL, type='character', help="Prefix for <prefix>.thresholds.tsv / .json [required]"),
  make_option(c("--expected-rows"), default=NULL, type='character', help="expected_rows.tsv from the partitioner; enables the completeness report and the underpowered chunk-root filter"),
  make_option(c("--label"),         default=NULL, type='character', help="Free-text label recorded in the JSON sidecar"),
  make_option(c("--keep-underpowered"), action='store_true', default=FALSE, help="Keep chunk-root rows whose subset is below the partitioner's --min-taxa/--min-refseq (a monolithic run would not have produced them)"),
  make_option(c("--overwrite"),     action='store_true', default=FALSE, help="Replace existing output files")
)

opt <- parse_args(OptionParser(
  option_list = option_list,
  usage = paste0("%prog --thresholds FILE[,FILE...] --sidecars FILE[,FILE...] ",
                 "--output-prefix PATH [options]"),
  description = "\nMerge per-chunk OptimOTU threshold tables.\n",
  epilogue = paste0(
    "--thresholds and --sidecars are matched by position; the file stems are\n",
    "checked so a caller that collects the two lists separately cannot silently\n",
    "pair a table with another chunk's configuration.\n")
))

## Validation of required arguments
require_options(opt, c("thresholds", "sidecars", "output-prefix"))

############################################## Load packages

cat("\nLoading R packages...\n")

load_pckg("data.table")
require_package("optimotu", OPTIMOTU_MIN_VERSION)
require_package("jsonlite", "1.8.0")

############################################## Assign and validate parameters

threshold_files <- as_list_arg(opt$thresholds)
sidecar_files <- as_list_arg(opt$sidecars)
if (!length(threshold_files)) fail("--thresholds is empty")
if (length(threshold_files) != length(sidecar_files)) {
  fail("--thresholds has ", length(threshold_files), " file(s) but --sidecars ",
       "has ", length(sidecar_files))
}
threshold_files <- vapply(threshold_files, normalizePath, character(1),
                          mustWork = TRUE, USE.NAMES = FALSE)
sidecar_files <- vapply(sidecar_files, normalizePath, character(1),
                        mustWork = TRUE, USE.NAMES = FALSE)

# The two lists are matched by position, so a caller that collects them
# separately can silently pair each table with another chunk's configuration.
# The stems always agree, so check them.
stems <- list(sub("\\.thresholds\\.tsv$", "", basename(threshold_files)),
              sub("\\.thresholds\\.json$", "", basename(sidecar_files)))
mismatched <- which(stems[[1L]] != stems[[2L]])
if (length(mismatched)) {
  fail("--thresholds and --sidecars are not in the same order: position ",
       mismatched[[1L]], " pairs ", basename(threshold_files[[mismatched[[1L]]]]),
       " with ", basename(sidecar_files[[mismatched[[1L]]]]))
}
prefix <- path.expand(opt$`output-prefix`)
expected_path <- if (is.null(opt$`expected-rows`)) NULL else {
  normalizePath(opt$`expected-rows`, mustWork = TRUE)
}

## Log assigned variables
cat("\nParameters specified:\n")
report_params(list(
  "Chunk tables"      = length(threshold_files),
  "Chunk sidecars"    = length(sidecar_files),
  "Expected rows"     = expected_path,
  "Output prefix"     = prefix,
  "Label"             = opt$label,
  "Keep underpowered" = isTRUE(opt$`keep-underpowered`),
  "Overwrite"         = isTRUE(opt$overwrite)
))
cat("\n")

############################################## Main workflow

cat("\n--- Main Workflow ---\n\n")

## --- read the chunks --------------------------------------------------------

KEY <- c("rank", "superrank", "supertaxon", "measure")
SUBSET_KEY <- c("rank", "superrank", "supertaxon")

read_chunk <- function(path) {
  # na.strings = "NA" reproduces read.delim's default exactly
  tab <- data.table::fread(path, sep = "\t", header = TRUE, quote = "",
                           colClasses = "character", check.names = FALSE,
                           na.strings = "NA", showProgress = FALSE,
                           data.table = FALSE)
  missing_columns <- setdiff(c(KEY, "threshold", "value"), names(tab))
  if (length(missing_columns)) {
    fail(basename(path), " lacks required columns: ",
         paste(missing_columns, collapse = ", "))
  }
  tab$threshold <- suppressWarnings(as.numeric(tab$threshold))
  tab$value <- suppressWarnings(as.numeric(tab$value))
  tab
}

step("Reading ", length(threshold_files), " chunk threshold table(s) ...")
chunks <- lapply(threshold_files, read_chunk)

step("Reading ", length(sidecar_files), " chunk sidecar(s) ...")
sidecars <- lapply(sidecar_files, function(p) jsonlite::read_json(p, simplifyVector = TRUE))
chunk_ids <- vapply(seq_along(sidecars), function(i) {
  id <- sidecars[[i]]$chunk$chunk_id
  if (is.null(id) || !length(id) || !nzchar(id)) {
    sub("\\.thresholds\\.tsv$", "", basename(threshold_files[[i]]))
  } else {
    as.character(id)[[1L]]
  }
}, character(1))
# rep(): a chunk left with nothing to optimize (--min-species-seqs) has no rows
for (i in seq_along(chunks)) chunks[[i]]$chunk_id <- rep(chunk_ids[[i]], nrow(chunks[[i]]))

# rbindlist instead of do.call(rbind, ...): rbind() reallocates the whole result
# once per chunk, which is quadratic in the number of chunks
thresholds <- data.table::rbindlist(chunks, use.names = TRUE)
note("Merged ", nrow(thresholds), " rows from ", length(chunks), " chunk(s)")

## --- keys must be disjoint --------------------------------------------------

step("Checking that the chunk key sets are disjoint ...")

# duplicated(by = KEY) hashes the key columns directly, instead of materializing
# a pasted key string for every row
offending <- duplicated(thresholds, by = KEY)
if (any(offending)) {
  dup <- unique(thresholds[offending, KEY, with = FALSE])
  detail <- vapply(seq_len(min(5L, nrow(dup))), function(i) {
    row <- dup[i]
    hit <- thresholds[row, on = KEY, nomatch = NULL]
    paste0(paste(unlist(row, use.names = FALSE), collapse = "/"), " (in ",
           paste(unique(hit$chunk_id), collapse = ", "), ")")
  }, character(1))
  fail(nrow(dup), " (rank, superrank, supertaxon, measure) key(s) were ",
       "produced by more than one chunk, so the partition is not disjoint and ",
       "the merge would silently pick one at random: ",
       paste(detail, collapse = "; "))
}

## --- the chunks must describe the same experiment ---------------------------

step("Checking that the chunks share one training configuration ...")

# Ranks legitimately differ between chunks: a chunk rooted at an order is run
# with `order,family,genus,species` so that rank becomes its own fallback level.
# Every chunk's rank vector must therefore be a suffix of the deepest one, which
# is the rank vector of the whole reference.
chunk_ranks <- lapply(sidecars, function(sc) as.character(sc$ranks))
ranks <- chunk_ranks[[which.max(lengths(chunk_ranks))]]
for (i in seq_along(chunk_ranks)) {
  tail_of <- utils::tail(ranks, length(chunk_ranks[[i]]))
  if (!identical(chunk_ranks[[i]], tail_of)) {
    fail("Chunk ", chunk_ids[[i]], " was trained with ranks ",
         paste(chunk_ranks[[i]], collapse = ","),
         ", which is not a tail of the reference ranks ",
         paste(ranks, collapse = ","), "; the chunks do not describe one taxonomy")
  }
}
root_rank <- ranks[[1L]]

# The root chunk (the one carrying the full rank vector and owning the root
# rank) describes the whole reference, so its sidecar is the base of the merged
# one.  Picking it by content rather than by position keeps the merge
# independent of how the chunks happen to be named or ordered
is_root_chunk <- vapply(seq_along(sidecars), function(i) {
  identical(chunk_ranks[[i]], ranks) &&
    (is.null(sidecars[[i]]$chunk) || root_rank %in% sidecars[[i]]$chunk$own_superranks)
}, logical(1))
if (sum(is_root_chunk) != 1L) {
  fail("Expected exactly one chunk rooted at ", root_rank, " (full rank vector ",
       "and ", root_rank, " among its own superranks), found ", sum(is_root_chunk),
       if (any(is_root_chunk)) {
         paste0(": ", paste(chunk_ids[is_root_chunk], collapse = ", "))
       })
}
base_index <- which(is_root_chunk)
base <- sidecars[[base_index]]
note("  root chunk: ", chunk_ids[[base_index]])

# optimotu_cluster.R refuses a thresholds file whose distance configuration does
# not match the one it is about to cluster with, so a merged table assembled
# from chunks that disagreed would be meaningless.
config_of <- function(sc) list(distance = sc$distance, clustering = sc$clustering,
                               optimization = sc$optimization[c(
                                 "measures", "dist_min", "dist_max", "dist_step",
                                 "n_thresholds", "min_taxa", "min_refseq",
                                 "min_species_seqs", "drop_unpaired",
                                 "reference_root_rank")])
reference_config <- config_of(base)
for (i in seq_along(sidecars)) {
  if (!identical(config_of(sidecars[[i]]), reference_config)) {
    fail("Chunk ", chunk_ids[[i]], " was trained with a different configuration ",
         "than ", chunk_ids[[base_index]], "; thresholds from different distance, ",
         "clustering, measure or threshold-grid settings cannot be merged")
  }
}

# The letters map SINTAX-style prefixes to ranks for every downstream reader,
# so all chunks must agree on them
letters_of <- function(sc) {
  x <- unlist(sc$rank_letters)
  if (is.null(x)) x else x[order(names(x))]
}
for (i in seq_along(sidecars)) {
  if (!identical(letters_of(sidecars[[i]]), letters_of(base))) {
    fail("Chunk ", chunk_ids[[i]], " has different rank letters than the root ",
         "chunk ", chunk_ids[[base_index]])
  }
}

measures <- as.character(base$optimization$measures)
note("  reference ranks: ", paste(ranks, collapse = ","))
note("  measures: ", paste(measures, collapse = ", "))

## --- compare against the plan -----------------------------------------------

underpowered <- 0L
if (!is.null(expected_path)) {
  step("Comparing the merged table against the partitioner's plan ...")
  expected <- data.table::fread(expected_path, sep = "\t", header = TRUE,
                                quote = "", colClasses = "character",
                                check.names = FALSE, na.strings = "NA",
                                showProgress = FALSE)
  missing_columns <- setdiff(c(SUBSET_KEY, "optimized"), names(expected))
  if (length(missing_columns)) {
    fail("expected_rows.tsv lacks required columns: ",
         paste(missing_columns, collapse = ", "))
  }
  optimized <- toupper(expected$optimized) %in% c("TRUE", "T", "1")
  want <- unique(expected[optimized, SUBSET_KEY, with = FALSE])
  got <- unique(thresholds[, SUBSET_KEY, with = FALSE])

  # optimize_thresholds() always keeps its ranks[1] subsets, however small, so a
  # chunk emits chunk-root rows that a monolithic run would have filtered out.
  # Dropping them here is what makes chunked and unchunked runs comparable.
  extra <- data.table::fsetdiff(got, want)
  if (nrow(extra) && !isTRUE(opt$`keep-underpowered`)) {
    drop <- thresholds[extra, on = SUBSET_KEY, which = TRUE, nomatch = NULL]
    underpowered <- length(drop)
    thresholds <- thresholds[!drop]
    note("Dropped ", underpowered, " row(s) in ", nrow(extra),
         " subset(s) below the partitioner's --min-taxa/--min-refseq ",
         "(pass --keep-underpowered to retain them)")
  }
  absent <- data.table::fsetdiff(
    want, unique(thresholds[, SUBSET_KEY, with = FALSE])
  )
  if (nrow(absent)) {
    note("NOTE: ", nrow(absent), " planned subset(s) produced no threshold ",
         "(expected where a subsampled chunk's sample, or --min-species-seqs, ",
         "left a subset below the size filters), e.g. ",
         paste(apply(utils::head(absent, 3L), 1L, paste, collapse = "/"),
               collapse = "; "))
  }
}

if (!nrow(thresholds)) fail("Merging left no thresholds")

## --- the check each chunk could not make ------------------------------------

step("Checking the ", root_rank, "-level fallback coverage ...")

# optimotu::calc_taxon_thresholds() falls back up the taxonomy and needs a
# root-level row for every rank; optimotu_cluster.R rejects a table without one.
coverage <- root_fallback_gaps(thresholds, ranks, measures)
complete <- measures[lengths(coverage) == 0L]
if (!length(complete)) {
  fail("The merged table has no ", root_rank, "-level threshold for every rank ",
       "under any measure (missing at ",
       paste(sort(unique(unlist(coverage))), collapse = ", "),
       "). The chunk rooted at ", root_rank, " must cover every rank; lower ",
       "--min-taxa/--min-refseq or raise the partitioner's --subsample-max.")
}
for (m in measures[lengths(coverage) > 0L]) {
  note("Measure ", m, " has no ", root_rank, "-level threshold for ",
       paste(coverage[[m]], collapse = ", "),
       "; do not select --metric ", m, " at clustering time")
}
note("Measures usable as --metric: ", paste(complete, collapse = ", "))

## --- outputs ----------------------------------------------------------------

step("Sorting and writing the merged thresholds ...")

# order() rather than setorderv(): the sort is on `supertaxon`/`measure`, and
# data.table sorts character columns in the C locale while order() uses the
# session's collation.  The table is at most a few thousand rows, so keeping the
# original row order is worth more here than the faster sort
data.table::setDF(thresholds)
thresholds <- thresholds[
  order(match(thresholds$rank, ranks), match(thresholds$superrank, ranks),
        thresholds$supertaxon, thresholds$measure),
  , drop = FALSE
]
rownames(thresholds) <- NULL

max_clustering_distance <- max(
  optimotu::threshold_as_dist(thresholds$threshold), na.rm = TRUE
)

# The chunk_id column is provenance, not part of the contract optimotu_cluster.R
# reads, so it is published beside the table rather than inside it.
provenance <- thresholds[, c(KEY, "chunk_id")]
out_table <- thresholds[, setdiff(names(thresholds), "chunk_id"), drop = FALSE]

rows_per_chunk <- table(thresholds$chunk_id)
outputs <- paste0(prefix, c(".thresholds.tsv", ".thresholds.json",
                            ".thresholds.provenance.tsv"))

# Sidecars written before file_record() carry an absolute `path`, and their
# training_fasta `sequences` counts the chunk's subset rather than the FASTA
merged_input <- function(rec, chunked) {
  if (is.null(rec) || is.null(rec$path)) return(rec)
  out <- list(file = basename(rec$path), md5 = rec$md5)
  if (!chunked && !is.null(rec$sequences)) out$sequences <- rec$sequences
  out
}
base_chunked <- !is.null(base$chunk$seq_ids)
inputs <- list(
  training_fasta = merged_input(base$inputs$training_fasta, base_chunked),
  training_taxonomy = merged_input(base$inputs$training_taxonomy, base_chunked)
)
inputs <- inputs[!vapply(inputs, is.null, logical(1))]

# threads (or its older forms) describes one chunk's job, not the merged run
optimization <- base$optimization[setdiff(
  names(base$optimization),
  c("threads", "requested_threads", "effective_dist_threads", "elapsed_seconds")
)]
optimization$elapsed_seconds <- sum(vapply(sidecars, function(sc) {
  x <- sc$optimization$elapsed_seconds
  if (is.null(x)) NA_real_ else as.numeric(x)
}, numeric(1)))

counts <- list(
  thresholds = nrow(out_table),
  chunks = length(chunks),
  complete_measures = complete,
  dropped_underpowered = underpowered
)
# counted over the whole reference, so every chunk records the same number
# (absent in chunks trained without --min-species-seqs; NULL leaves it out)
counts$masked_species <- base$counts$masked_species
# rows each chunk dropped with --drop-unpaired (absent in older chunks)
if (isTRUE(base$optimization$drop_unpaired)) {
  counts$unpaired_dropped <- sum(vapply(sidecars, function(sc) {
    x <- sc$counts$unpaired_dropped
    if (is.null(x)) 0L else as.integer(x)
  }, integer(1)))
}

with_atomic_outputs(outputs, overwrite = isTRUE(opt$overwrite), writer = function(staged) {
  data.table::fwrite(out_table, staged[[1L]], sep = "\t", quote = FALSE, na = "")
  write_json_file(
    list(
      schema_version = base$schema_version,
      created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
      label = opt$label,
      ranks = ranks,
      root_rank = root_rank,
      root_taxon = base$root_taxon,
      rank_letters = base$rank_letters,
      inputs = inputs,
      distance = base$distance,
      clustering = base$clustering,
      optimization = optimization,
      packages = base$packages,
      counts = counts,
      max_clustering_distance = max_clustering_distance,
      merged_from = lapply(seq_along(chunks), function(i) list(
        chunk_id = chunk_ids[[i]],
        own_superranks = sidecars[[i]]$chunk$own_superranks,
        sequences = sidecars[[i]]$chunk$seq_ids$sequences,
        rows = if (chunk_ids[[i]] %in% names(rows_per_chunk)) {
          as.integer(rows_per_chunk[[chunk_ids[[i]]]])
        } else 0L
      ))
    ),
    staged[[2L]]
  )
  data.table::fwrite(provenance, staged[[3L]], sep = "\t", quote = FALSE, na = "")
})

note("Wrote ", nrow(out_table), " thresholds (maximum clustering distance ",
     format(max_clustering_distance, digits = 6), ") to ",
     paste(basename(outputs), collapse = ", "), " under ", dirname(prefix))

############################################## Session info

finish(start_time)
