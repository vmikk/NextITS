#!/usr/bin/env Rscript

## Script: optimotu_train.R
## Purpose: Optimize taxon-specific OptimOTU clustering thresholds from a
##          reference training set
##
##          Writes a thresholds TSV plus a sidecar JSON recording the distance
##          configuration the thresholds are only valid for
##
## Usage:
##   Rscript optimotu_train.R \
##      --training-fasta reference.training.fasta \
##      --training-taxonomy reference.training_taxonomy.tsv \
##      --output-prefix reference \
##      --metric FM --dist-max 0.4 --dist-step 0.001 \
##      --threads 4
##
## Output:
##   - <prefix>.thresholds.tsv
##   - <prefix>.thresholds.json

options(stringsAsFactors = FALSE, warn = 1)

## Check time
start_time <- Sys.time()

script_dir <- dirname(normalizePath(sub(
  "^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]]
)))
source(file.path(script_dir, "optimotu_common.R"))
SCRIPT_START <- start_time

banner("optimotu_train.R", "Train taxon-specific OptimOTU clustering thresholds")

cat("Parsing input options and arguments...\n")

suppressPackageStartupMessages(require(optparse))

## Parse arguments
option_list <- list(
  make_option(c("--training-fasta"),    default=NULL, type='character', help="Training sequences [required]"),
  make_option(c("--output-prefix"),     default=NULL, type='character', help="Prefix for <prefix>.thresholds.tsv / .json [required]"),
  make_option(c("--training-taxonomy"), default=NULL, type='character', help="Taxonomy table [default: parse the FASTA headers]"),
  make_option(c("--taxonomy-format"),   default="auto", type='character', help="One of: auto, eukaryome, unite, sintax, bold, lineage, tsv [default %default]"),
  make_option(c("--ranks"),             default=NULL, type='character', help="Rank names, most to least inclusive, comma-separated [default: root,kingdom,phylum,class,order,family,genus,species]"),
  make_option(c("--root-taxon"),        default=NULL, type='character', help="Value for the synthetic first rank [default: Eukaryota]"),
  make_option(c("--rank-letters"),      default=NULL, type='character', help="Letter=rank map for prefixed/SINTAX headers, comma-separated"),
  make_option(c("--seq-ids"),           default=NULL, type='character', help="Chunked training: train on just these sequence IDs"),
  make_option(c("--own-superranks"),    default=NULL, type='character', help="Chunked training: emit only rows with these superranks; also relaxes the root-fallback check, which only the merged table can satisfy"),
  make_option(c("--chunk-id"),          default=NULL, type='character', help="Chunked training: chunk name recorded in the JSON sidecar"),
  make_option(c("--metric"),            default=NULL, type='character', help="Measure(s) to optimize, comma-separated, or 'all' [default: all]"),
  make_option(c("--dist-max"),          default=0.4,  type='double',  help="Largest tested distance [default %default]"),
  make_option(c("--dist-step"),         default=0.001, type='double', help="Threshold grid spacing [default %default]"),
  make_option(c("--min-taxa"),          default=5L,   type='integer', help="Minimum subtaxa per optimized subset [default %default]"),
  make_option(c("--min-refseq"),        default=10L,  type='integer', help="Minimum sequences per optimized subset [default %default]"),
  make_option(c("--dist-method"),       default="wfa2", type='character', help="One of: wfa2, ksw2, edlib, hybrid, hamming, usearch [default %default]"),
  make_option(c("--dist-penalties"),    default=NULL, type='character', help="WFA2/KSW2 penalties, e.g. gap_open=2,gap_extend=1 [default: edit distance]"),
  make_option(c("--clust-method"),      default="tree", type='character', help="One of: tree, matrix, index, slink [default %default]"),
  make_option(c("--threads"),           default=1L,   type='integer', help="Worker threads for alignment and clustering [default %default]"),
  make_option(c("--usearch"),           default=NULL, type='character', help="usearch executable for --dist-method usearch [default: autodetected from $USEARCH, then PATH]"),
  make_option(c("--label"),             default=NULL, type='character', help="Free-text label recorded in the JSON sidecar"),
  make_option(c("--no-qualify-homonyms"), action='store_true', default=FALSE, help="Keep taxon names that occur under several parents"),
  make_option(c("--keep-orphan-ranks"), action='store_true', default=FALSE, help="Bridge Incertae sedis gaps with a synthetic taxon instead of truncating the lineage there"),
  make_option(c("--verbose"),           action='store_true', default=FALSE, help="Print OptimOTU progress"),
  make_option(c("--overwrite"),         action='store_true', default=FALSE, help="Replace existing output files")
)

opt <- parse_args(OptionParser(
  option_list = option_list,
  usage = "%prog --training-fasta FILE --output-prefix PATH [options]",
  description = "\nTrain taxon-specific OptimOTU clustering thresholds.\n",
  epilogue = paste0(
    "Thresholds are only valid for the distance and clustering configuration\n",
    "they were trained with; the JSON sidecar records it and optimotu_cluster.R\n",
    "refuses a mismatch.  For chunked training see optimotu_partition_reference.R\n",
    "and optimotu_merge_thresholds.R.\n")
))

## Validation of required arguments
require_options(opt, c("training-fasta", "output-prefix"))

############################################## Load packages

cat("\nLoading R packages...\n")

load_pckg("data.table")
require_package("optimotu", OPTIMOTU_MIN_VERSION)
require_package("Biostrings", "2.70.0")
require_package("jsonlite", "1.8.0")

############################################## Assign and validate parameters


ranks <- resolve_ranks(opt$ranks)
root_rank <- ranks[[1L]]
root_taxon <- if (is.null(opt$`root-taxon`)) DEFAULT_ROOT_TAXON else opt$`root-taxon`
rank_letters <- as_kv_arg(opt$`rank-letters`, "rank-letters", DEFAULT_RANK_LETTERS)
taxonomy_format <- as_choice(opt$`taxonomy-format`, "taxonomy-format",
                             TAXONOMY_FORMATS, "auto")

measures <- if (is.null(opt$metric) || identical(tolower(opt$metric), "all")) {
  MEASURES
} else {
  requested <- toupper(as_list_arg(opt$metric))
  unknown <- setdiff(requested, MEASURES)
  if (length(unknown)) {
    fail("Unsupported --metric: ", paste(unknown, collapse = ", "),
         " (choose from ", paste(MEASURES, collapse = ", "), " or 'all')")
  }
  requested
}

dist_max <- as_number(opt$`dist-max`, "dist-max", 0.4, 0, inclusive = FALSE)
dist_step <- as_number(opt$`dist-step`, "dist-step", 0.001, 0, inclusive = FALSE)
if (dist_step > dist_max) fail("--dist-step cannot exceed --dist-max")
min_taxa <- as_count(opt$`min-taxa`, "min-taxa", 5L)
min_refseq <- as_count(opt$`min-refseq`, "min-refseq", 10L)
if (min_refseq < min_taxa) fail("--min-refseq cannot be smaller than --min-taxa")

dist_method <- as_choice(opt$`dist-method`, "dist-method", DIST_METHODS, "wfa2")
clust_method <- as_choice(opt$`clust-method`, "clust-method", CLUST_METHODS, "tree")
penalties <- as_kv_arg(opt$`dist-penalties`, "dist-penalties")
threads <- as_count(opt$threads, "threads", 1L)

own_superranks <- as_list_arg(opt$`own-superranks`, NULL)
if (!is.null(own_superranks)) {
  unknown <- setdiff(own_superranks, ranks)
  if (length(unknown)) {
    fail("--own-superranks names ranks that are not in --ranks: ",
         paste(unknown, collapse = ", "))
  }
}
chunked <- !is.null(own_superranks) || !is.null(opt$`seq-ids`)
chunk_id <- opt$`chunk-id`

training_path <- normalizePath(opt$`training-fasta`, mustWork = TRUE)
seq_ids_path <- if (is.null(opt$`seq-ids`)) NULL else {
  normalizePath(opt$`seq-ids`, mustWork = TRUE)
}
taxonomy_path <- if (is.null(opt$`training-taxonomy`)) NULL else {
  normalizePath(opt$`training-taxonomy`, mustWork = TRUE)
}
prefix <- path.expand(opt$`output-prefix`)

usearch <- if (identical(dist_method, "usearch")) {
  find_executable("usearch", opt$usearch)
} else NULL

## Log assigned variables
cat("\nParameters specified:\n")
report_params(list(
  "Training FASTA"    = training_path,
  "Training taxonomy" = taxonomy_path,
  "Taxonomy format"   = taxonomy_format,
  "Output prefix"     = prefix,
  "Ranks"             = ranks,
  "Root taxon"        = root_taxon,
  "Rank letters"      = rank_letters,
  "Sequence IDs"      = seq_ids_path,
  "Chunk ID"          = chunk_id,
  "Own superranks"    = own_superranks,
  "Measures"          = measures,
  "Max distance"      = dist_max,
  "Distance step"     = dist_step,
  "Min taxa"          = min_taxa,
  "Min refseq"        = min_refseq,
  "Distance method"   = dist_method,
  "Distance penalties" = penalties,
  "Clustering method" = clust_method,
  "Threads"           = threads,
  "Label"             = opt$label,
  "Qualify homonyms"  = !isTRUE(opt$`no-qualify-homonyms`),
  "Keep orphan ranks" = isTRUE(opt$`keep-orphan-ranks`),
  "Verbose"           = isTRUE(opt$verbose),
  "Overwrite"         = isTRUE(opt$overwrite)
))
if (!is.null(usearch)) report_tool("usearch", usearch)
cat("\n")

############################################## Configuration

dist_config <- build_dist_config(dist_method, penalties, threads, usearch)
clust_config <- build_clust_config(clust_method)

############################################## Main workflow

cat("\n--- Main Workflow ---\n\n")

step("Reading training sequences from ", training_path, " ...")
training <- read_fasta(training_path, "Training FASTA")
n_training_total <- length(training)
note("  ", length(training), " sequences")

# A chunk trains on one taxonomic subtree out of a shared reference, so it is
# handed the whole FASTA plus the ID list it owns.  Subsetting happens before
# the taxonomy is cleaned so that require_single_root, and the lineage
# truncation, see only this chunk's records.
if (!is.null(seq_ids_path)) {
  wanted <- unique(trimws(readLines(seq_ids_path, warn = FALSE)))
  wanted <- wanted[nzchar(wanted)]
  if (!length(wanted)) fail("--seq-ids is empty: ", seq_ids_path)
  absent <- setdiff(wanted, names(training))
  if (length(absent)) {
    fail("--seq-ids names ", length(absent), " sequence(s) missing from the ",
         "training FASTA (e.g. ",
         paste(utils::head(absent, 5L), collapse = ", "), ")")
  }
  training <- training[names(training) %in% wanted]
  note("  ", length(training), " sequences selected by --seq-ids")
}

step("Parsing and cleaning the training taxonomy ...")
taxonomy <- if (is.null(taxonomy_path)) {
  read_reference_taxonomy(training_path, ranks, taxonomy_format, rank_letters,
                          headers = NULL)
} else {
  read_reference_taxonomy(taxonomy_path, ranks, taxonomy_format, rank_letters)
}
if (!is.null(seq_ids_path)) {
  taxonomy <- taxonomy[taxonomy$seq_id %in% names(training), , drop = FALSE]
}
taxonomy <- add_root_rank(taxonomy, root_rank, root_taxon)
taxonomy <- clean_reference_taxonomy(
  taxonomy, ranks, "Training taxonomy",
  require_single_root = TRUE,
  qualify = !isTRUE(opt$`no-qualify-homonyms`),
  keep_orphan_ranks = isTRUE(opt$`keep-orphan-ranks`)
)

missing_tax <- setdiff(names(training), taxonomy$seq_id)
extra_tax <- setdiff(taxonomy$seq_id, names(training))
if (length(missing_tax) || length(extra_tax)) {
  fail("Training taxonomy IDs do not match the FASTA (missing: ",
       paste(utils::head(missing_tax, 5L), collapse = ", "), "; extra: ",
       paste(utils::head(extra_tax, 5L), collapse = ", "), ")")
}
# optimize_thresholds() reorders internally, but SLINK is order-sensitive, so
# align the two up front.
taxonomy <- taxonomy[match(names(training), taxonomy$seq_id), , drop = FALSE]
rownames(taxonomy) <- NULL

# The effective root taxon is whatever the cleaned taxonomy actually carries at
# ranks[1]; --root-taxon only supplies a value where the column was absent or
# empty, so a reference already rooted at, say, kingdom Fungi keeps its own name.
root_taxon <- unique(taxonomy[[root_rank]])

# one pass per rank column, rather than materializing a full logical matrix
assigned <- vapply(ranks, function(r) sum(!is.na(taxonomy[[r]])), integer(1))
note("Assigned per rank: ",
     paste(sprintf("%s=%d", ranks, assigned), collapse = " "))

## --- optimize ---------------------------------------------------------------

step("Building the threshold grid ...")

# The threshold grid is named explicitly. optimotu's WFA2 path sets numeric row
# names on the cluster matrix itself, but the usearch path routes through
# distmx_cluster(), which sets none -- and find_best_threshold() asserts
# `row.names = "unique"` and then reads them back with as.numeric(). Supplying
# thresh_names makes both backends work, so long as the names round-trip as
# numbers.
threshold_grid <- round(seq(0, dist_max, by = dist_step), 10L)
threshold_names <- format(threshold_grid, trim = TRUE, scientific = FALSE,
                          drop0trailing = TRUE)
stopifnot(!anyDuplicated(threshold_names),
          isTRUE(all.equal(as.numeric(threshold_names), threshold_grid)))
n_thresholds <- length(threshold_grid)
note("Optimizing ", paste(measures, collapse = "/"), " over ", n_thresholds,
     " thresholds (0-", dist_max, " step ", dist_step, ") with ", dist_method,
     " distances on ", length(training), " sequences ...")

step("Running optimotu::optimize_thresholds() ...")
started <- Sys.time()
# MCC (and, on tiny training sets, other confusion-matrix measures) is undefined
# for a subset containing only one true taxon: every between-taxon pair count is
# zero, so its denominator vanishes.  optimotu signals that with a bare
# `max(NA, na.rm = TRUE)` warning per affected subset, which is noise; collect
# them and report the outcome as coverage instead.
thresholds <- withCallingHandlers(
  optimotu::optimize_thresholds(
    taxonomy = taxonomy,
    refseq = training,
    ranks = ranks,
    dist_config = dist_config,
    threshold_config = optimotu::threshold_uniform(0, dist_max, dist_step,
                                                   threshold_names),
    clust_config = clust_config,
    parallel_config = optimotu::parallel_concurrent(threads),
    min_taxa = min_taxa,
    min_refseq = min_refseq,
    measures = measures,
    verbose = isTRUE(opt$verbose)
  ),
  warning = function(w) {
    if (grepl("no non-missing arguments to max|number of items to replace",
              conditionMessage(w))) {
      invokeRestart("muffleWarning")
    }
  }
)
elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))

thresholds <- as.data.frame(thresholds, stringsAsFactors = FALSE)
if (!nrow(thresholds)) fail("Threshold optimization produced no thresholds")
thresholds <- thresholds[!is.na(thresholds$threshold), , drop = FALSE]
if (!nrow(thresholds)) fail("Every optimized threshold was NA")

# Chunks overlap in what they *can* compute -- a chunk rooted at an order also
# sees its families -- but each owns a disjoint slice, so that merging is a
# concatenation with no precedence rules.  Drop everything outside this chunk's
# slice before the coverage check reads the table.
if (!is.null(own_superranks)) {
  before <- nrow(thresholds)
  thresholds <- thresholds[thresholds$superrank %in% own_superranks, , drop = FALSE]
  note("Keeping ", nrow(thresholds), " of ", before, " rows for superrank ",
       paste(own_superranks, collapse = "/"))
  if (!nrow(thresholds)) {
    fail("No optimized threshold has superrank ",
         paste(own_superranks, collapse = "/"), "; this chunk owns nothing. ",
         "Check --own-superranks against --ranks.")
  }
}

# Every rank needs a root-level fallback row, otherwise
# optimotu::calc_taxon_thresholds() has nothing to fall back to for sequences
# whose parent taxon was never optimized.  Coverage is per measure: a measure
# that is undefined on a single-taxon subset simply has no row there, which only
# matters if that measure is the one selected at clustering time.
root_rows <- thresholds[thresholds$superrank == root_rank, , drop = FALSE]
coverage <- lapply(measures, function(m) {
  setdiff(ranks[-1L], unique(root_rows$rank[root_rows$measure == m]))
})
names(coverage) <- measures
complete <- measures[lengths(coverage) == 0L]
if (!length(complete) && !chunked) {
  fail("No measure has a ", root_taxon, "-level threshold for every rank; the ",
       "training set has too few sequences assigned at ",
       paste(unique(unlist(coverage)), collapse = ", "),
       ". Lower --min-taxa/--min-refseq or use a broader training set.")
}
if (chunked) {
  # Only the chunk rooted at the reference's own root rank can supply the
  # fallback row every rank needs; for the others an incomplete measure is
  # expected.  optimotu_merge_thresholds.R makes the check that matters, on the
  # merged table.
  note("Chunk mode: deferring the root-fallback check to the merge step")
} else {
  for (m in measures[lengths(coverage) > 0L]) {
    note("Measure ", m, " has no ", root_taxon, "-level threshold for ",
         paste(coverage[[m]], collapse = ", "),
         " (undefined where the training set holds a single taxon); ",
         "do not select --metric ", m, " at clustering time")
  }
  note("Measures usable as --metric: ", paste(complete, collapse = ", "))
}

# order() rather than setorderv(): the sort is on `supertaxon`/`measure`, and
# data.table sorts character columns in the C locale while order() uses the
# session's collation
thresholds <- thresholds[
  order(match(thresholds$rank, ranks), match(thresholds$superrank, ranks),
        thresholds$supertaxon, thresholds$measure),
  , drop = FALSE
]
rownames(thresholds) <- NULL

max_clustering_distance <- max(
  optimotu::threshold_as_dist(thresholds$threshold), na.rm = TRUE
)
note("Optimized ", nrow(thresholds), " thresholds in ", round(elapsed, 1),
     " s; maximum clustering distance ", format(max_clustering_distance, digits = 6))

## --- outputs ----------------------------------------------------------------

outputs <- paste0(prefix, c(".thresholds.tsv", ".thresholds.json"))

step("Writing ", length(outputs), " output file(s) ...")

with_atomic_outputs(outputs, overwrite = isTRUE(opt$overwrite), writer = function(staged) {
  data.table::fwrite(thresholds, staged[[1L]], sep = "\t", quote = FALSE, na = "")
  write_json_file(
    list(
      schema_version = 1L,
      created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
      label = opt$label,
      chunk = if (!chunked) NULL else list(
        chunk_id = chunk_id,
        own_superranks = own_superranks,
        seq_ids = if (is.null(seq_ids_path)) NULL else {
          file_record(seq_ids_path, sequences = length(training))
        }
      ),
      ranks = ranks,
      root_rank = root_rank,
      root_taxon = root_taxon,
      rank_letters = as.list(rank_letters),
      inputs = list(
        # `sequences` is the whole FASTA; a --seq-ids subset is in chunk$seq_ids
        training_fasta = file_record(training_path, sequences = n_training_total),
        training_taxonomy = if (is.null(taxonomy_path)) NULL else {
          file_record(taxonomy_path)
        }
      ),
      distance = c(list(method = dist_method),
                   dist_config[setdiff(names(dist_config), c("method", "call"))]),
      clustering = list(method = clust_method),
      optimization = list(
        measures = measures,
        dist_min = 0,
        dist_max = dist_max,
        dist_step = dist_step,
        n_thresholds = n_thresholds,
        min_taxa = min_taxa,
        min_refseq = min_refseq,
        threads = threads,
        elapsed_seconds = elapsed
      ),
      packages = list(R = as.character(getRversion()),
                      optimotu = as.character(utils::packageVersion("optimotu"))),
      counts = list(
        thresholds = nrow(thresholds),
        complete_measures = complete,
        assigned = as.list(stats::setNames(as.integer(assigned), ranks))
      ),
      max_clustering_distance = max_clustering_distance
    ),
    staged[[2L]]
  )
})

note("Wrote ", paste(basename(outputs), collapse = ", "), " under ", dirname(prefix))

############################################## Session info

finish(start_time)
