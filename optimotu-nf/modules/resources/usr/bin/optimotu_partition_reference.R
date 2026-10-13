#!/usr/bin/env Rscript

## Script: optimotu_partition_reference.R
## Purpose: Plan an embarrassingly parallel threshold-training run
##
##          `optimotu::optimize_thresholds()` optimizes one threshold per
##          (superrank, rank, supertaxon) subset, and every such subset is a
##          taxonomic subtree -- so any two subsets are either nested or
##          disjoint.  Running the optimizer on one subtree therefore reproduces,
##          bit for bit, every row whose `superrank` lies at or below that
##          subtree's own rank.  Only rows whose `superrank` sits *above* the
##          split need sequences from across the database
##
##          This script turns that observation into a work plan.  It walks the
##          taxonomy top-down and, for each node, either claims the whole subtree
##          as an EXACT job (when its all-pairs cost fits the budget) or emits a
##          SUBSAMPLED job that owns only that node's own superrank and recurses
##          into the children.  Ownership is disjoint by construction, so merging
##          the per-job tables is a concatenation
##
##          It reads taxonomy only -- never sequences -- so it runs in seconds
##
## Usage:
##   Rscript optimotu_partition_reference.R \
##      --training-taxonomy reference.training_taxonomy.tsv \
##      --output-dir partition \
##      --budget-pairs 2.5e7 \
##      --subsample-max 6000
##
## Output:
##   - <output-dir>/partition.tsv       one row per training job
##   - <output-dir>/expected_rows.tsv   every threshold row a monolithic run
##                                      would produce, mapped to its owning chunk
##   - <output-dir>/partition.json      provenance sidecar
##   - <output-dir>/<chunk>.ids.txt     sequence IDs for each job

options(stringsAsFactors = FALSE, warn = 1)

## Check time
start_time <- Sys.time()

script_dir <- dirname(normalizePath(sub(
  "^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]]
)))
source(file.path(script_dir, "optimotu_common.R"))
SCRIPT_START <- start_time

banner("optimotu_partition_reference.R",
       "Plan a chunked, parallel OptimOTU threshold-training run")

cat("Parsing input options and arguments...\n")

suppressPackageStartupMessages(require(optparse))

## Parse arguments
option_list <- list(
  make_option(c("--training-taxonomy"), default=NULL, type='character', help="Cleaned, qualified taxonomy TSV from optimotu_reference_prep.R [required]"),
  make_option(c("--output-dir"),        default=NULL, type='character', help="Destination for the plan and per-chunk ID lists [required]"),
  make_option(c("--ranks"),             default=NULL, type='character', help="Rank names, most to least inclusive, comma-separated [default: root,kingdom,phylum,class,order,family,genus,species]"),
  make_option(c("--budget-pairs"),      default=2.5e7, type='double',  help="Largest all-pairs cost an exact chunk may have [default %default, about 7000 sequences]"),
  make_option(c("--subsample-max"),     default=6000L, type='integer', help="Sequence cap for a subsampled chunk [default %default]"),
  make_option(c("--subsample-min-per-taxon"), default=2L, type='integer', help="When a subsampled chunk must drop taxa, keep this many sequences from each it keeps [default %default]"),
  make_option(c("--min-taxa"),          default=5L,   type='integer',  help="Minimum subtaxa per optimized subset [default %default]"),
  make_option(c("--min-refseq"),        default=10L,  type='integer',  help="Minimum sequences per optimized subset [default %default]"),
  make_option(c("--seed"),              default=42L,  type='integer',  help="Subsampling seed [default %default]"),
  make_option(c("--prefix"),            default="chunk", type='character', help="Chunk ID prefix [default %default]"),
  make_option(c("--overwrite"),         action='store_true', default=FALSE, help="Replace existing plan files")
)

opt <- parse_args(OptionParser(
  option_list = option_list,
  usage = "%prog --training-taxonomy FILE --output-dir DIR [options]",
  description = "\nPlan a chunked, parallel OptimOTU threshold-training run.\n",
  epilogue = paste0(
    "The plan is disjoint by construction: every (superrank, supertaxon) pair is\n",
    "claimed by exactly one chunk, so optimotu_merge_thresholds.R only has to\n",
    "concatenate the per-chunk tables -- and it verifies that this held.\n")
))

## Validation of required arguments
require_options(opt, c("training-taxonomy", "output-dir"))

############################################## Load packages

cat("\nLoading R packages...\n")

load_pckg("data.table")
require_package("optimotu", OPTIMOTU_MIN_VERSION)
require_package("jsonlite", "1.8.0")

############################################## Assign and validate parameters

ranks <- resolve_ranks(opt$ranks)
budget <- as_number(opt$`budget-pairs`, "budget-pairs", 2.5e7, 1)
subsample_max <- as_count(opt$`subsample-max`, "subsample-max", 6000L)
subsample_min_per_taxon <- as_count(opt$`subsample-min-per-taxon`,
                                    "subsample-min-per-taxon", 2L)
min_taxa <- as_count(opt$`min-taxa`, "min-taxa", 5L)
min_refseq <- as_count(opt$`min-refseq`, "min-refseq", 10L)
seed <- as_count(opt$seed, "seed", 42L)
prefix <- sanitize_id(opt$prefix, "--prefix")
if (min_refseq < min_taxa) fail("--min-refseq cannot be smaller than --min-taxa")

taxonomy_path <- normalizePath(opt$`training-taxonomy`, mustWork = TRUE)
output_dir <- path.expand(opt$`output-dir`)

## Log assigned variables
cat("\nParameters specified:\n")
report_params(list(
  "Training taxonomy" = taxonomy_path,
  "Output dir"        = output_dir,
  "Ranks"             = ranks,
  "Budget (pairs)"    = format(budget, big.mark = ",", scientific = FALSE),
  "Subsample max"     = subsample_max,
  "Subsample min per taxon" = subsample_min_per_taxon,
  "Min taxa"          = min_taxa,
  "Min refseq"        = min_refseq,
  "Seed"              = seed,
  "Chunk prefix"      = prefix,
  "Overwrite"         = isTRUE(opt$overwrite)
))
cat("\n")

############################################## Main workflow

cat("\n--- Main Workflow ---\n\n")

## --- taxonomy ---------------------------------------------------------------

# The table is already cleaned and homonym-qualified by optimotu_reference_prep.R;
# re-running either step per chunk would change the true partition, so this
# script deliberately only reads it.
step("Reading the training taxonomy ...")
# na.strings = "NA" reproduces read.delim's default exactly
taxonomy <- data.table::fread(taxonomy_path, sep = "\t", header = TRUE,
                              quote = "", colClasses = "character",
                              check.names = FALSE, na.strings = "NA",
                              showProgress = FALSE)
missing_columns <- setdiff(c("seq_id", ranks), names(taxonomy))
if (length(missing_columns)) {
  fail("Training taxonomy lacks required columns: ",
       paste(missing_columns, collapse = ", "))
}
taxonomy <- taxonomy[, c("seq_id", ranks), with = FALSE]
# set() rewrites each rank column in place rather than copying it
for (rank in ranks) {
  value <- taxonomy[[rank]]
  value[!nzchar(trimws(value))] <- NA_character_
  data.table::set(taxonomy, j = rank, value = value)
}
# optimotu::summarize_by_rank() asserts a data.frame, and the partition walk
# indexes columns with data.frame semantics throughout
data.table::setDF(taxonomy)
n_seq_total <- nrow(taxonomy)
if (!n_seq_total) fail("Training taxonomy is empty")
# Chunks are addressed by sequence ID, so a duplicate would put the same record
# in two chunks and break the disjointness the merge relies on.
duplicates <- unique(taxonomy$seq_id[duplicated(taxonomy$seq_id)])
if (length(duplicates)) {
  fail("Training taxonomy has duplicate sequence IDs: ",
       paste(utils::head(duplicates, 10L), collapse = ", "))
}
if (anyNA(taxonomy[[ranks[[1L]]]])) {
  fail("Training taxonomy has records with no ", ranks[[1L]], " assignment")
}
note(n_seq_total, " sequences; ranks: ", paste(ranks, collapse = ", "))

## --- the rows a monolithic run would produce --------------------------------

# summarize_by_rank() is the same C++ enumeration optimize_thresholds() uses, so
# taking the plan's inventory from it guarantees the two agree about which
# subsets exist.  Only the counts are kept; the seq_id/true_partition list
# columns would be hundreds of megabytes on a full database.
step("Enumerating subsets ...")
subsets <- optimotu::summarize_by_rank(taxonomy, ranks, "seq_id")
subsets <- data.frame(
  rank       = as.character(subsets$rank),
  superrank  = as.character(subsets$superrank),
  supertaxon = as.character(subsets$supertaxon),
  n_seq      = as.integer(subsets$n_seq),
  n_taxa     = as.integer(subsets$n_taxa),
  stringsAsFactors = FALSE
)
# Exactly the filter optimize_thresholds() applies, including its exemption for
# subsets rooted at ranks[1]: those rows are the fallbacks
# optimotu::calc_taxon_thresholds() walks up to, so they are kept however small.
# A chunk rooted deeper down gets the same exemption for *its* root, which is
# what the merge step filters back out -- hence this column has to describe the
# monolithic run, not the chunked one.
subsets$optimized <-
  (subsets$n_seq >= min_refseq & subsets$n_taxa >= min_taxa) |
  subsets$superrank == ranks[[1L]]
note("  ", nrow(subsets), " subsets, ", sum(subsets$optimized),
     " above --min-taxa/--min-refseq")

## --- subsampling ------------------------------------------------------------

## Stratified subsample of a set of records.
##
## Strata are whole lineages (the deepest assigned path), not a single rank:
## a high-rank job needs the sequences that stop at order just as much as the
## ones identified to species, and stratifying on one `require_rank` -- as
## optimotu_reference_prep.R does -- would throw them away.  Within a stratum the
## per-taxon cap is chosen as large as the sequence budget allows, so the sample
## keeps as many within-taxon pairs as it can afford before it starts dropping
## whole taxa.
stratified_sample <- function(idx, cap) {
  if (length(idx) <= cap) return(idx)
  key <- do.call(paste, c(lapply(ranks, function(r) taxonomy[[r]][idx]),
                          list(sep = "\r")))
  groups <- unname(split(idx, key))
  sizes <- lengths(groups)
  draw <- function(g, m) if (length(g) <= m) g else g[sample.int(length(g), m)]

  # Largest per-stratum cap whose total still fits, by binary search over m.
  lo <- 0L
  hi <- max(sizes)
  while (lo < hi) {
    mid <- as.integer((lo + hi + 1L) %/% 2L)
    if (sum(pmin(sizes, mid)) <= cap) lo <- mid else hi <- mid - 1L
  }
  if (lo >= 1L) {
    return(sort(unlist(lapply(groups, draw, m = lo), use.names = FALSE)))
  }

  # Even one sequence per stratum overflows, so strata have to be dropped.  Keep
  # `subsample_min_per_taxon` sequences from as many strata as fit rather than
  # one sequence from more strata: a taxon represented once contributes no
  # within-taxon pair, and a sample made entirely of such taxa presents the
  # optimizer with a partition of singletons, whose best threshold is a
  # degenerate zero.
  shuffled <- sample.int(length(groups))
  kept <- pmin(sizes[shuffled], subsample_min_per_taxon)
  shuffled <- shuffled[cumsum(kept) <= cap]
  if (!length(shuffled)) {
    fail("--subsample-max ", cap, " is too small to hold even one taxon at ",
         subsample_min_per_taxon, " sequence(s)")
  }
  sort(unlist(lapply(groups[shuffled], draw, m = subsample_min_per_taxon),
              use.names = FALSE))
}

## --- the partition ----------------------------------------------------------

pairs_of <- function(n) n * (n - 1) / 2

jobs <- list()
emit <- function(kind, rank_index, supertaxon, label, own_ranks, idx) {
  jobs[[length(jobs) + 1L]] <<- list(
    kind = kind, rank_index = rank_index, supertaxon = supertaxon,
    label = label, own_ranks = own_ranks, idx = idx
  )
}

# The tip rank has no subtaxa, so it is never a superrank and never roots a job.
last_superrank <- length(ranks) - 1L

## Claim a taxonomic subtree, splitting it only when it is too expensive.
##
## `i` is the index in `ranks` of the node's own rank, `idx` the rows under it
## (all of which have that rank assigned).
partition <- function(i, supertaxon, idx) {
  n <- length(idx)
  # a subset needs at least two sequences before any pair, and hence any
  # threshold, exists
  if (n < 2L || i > last_superrank) return(invisible(NULL))

  own <- ranks[[i]]
  # Collapse a chain of single children that covers the whole node -- root ->
  # kingdom Fungi is the usual case.  Both superranks see exactly the same
  # sequences, so one job can own both instead of clustering them twice.
  # `supertaxon` stays the taxon at the job's own rank; `label` follows the chain
  # down, because the deepest name is the one that says what the chunk holds.
  label <- supertaxon
  j <- i
  while (j < last_superrank) {
    child <- taxonomy[[ranks[[j + 1L]]]][idx]
    values <- unique(child[!is.na(child)])
    if (length(values) != 1L || anyNA(child)) break
    j <- j + 1L
    own <- c(own, ranks[[j]])
    label <- values[[1L]]
  }

  if (pairs_of(n) <= budget) {
    # the whole subtree fits, so this job owns every superrank inside it
    emit("exact", i, supertaxon, label, ranks[i:last_superrank], idx)
    return(invisible(NULL))
  }

  emit("subsampled", i, supertaxon, label,
       own, stratified_sample(idx, subsample_max))
  if (j >= last_superrank) return(invisible(NULL))
  child <- taxonomy[[ranks[[j + 1L]]]][idx]
  for (value in sort(unique(child[!is.na(child)]))) {
    partition(j + 1L, value, idx[!is.na(child) & child == value])
  }
  invisible(NULL)
}

step("Walking the taxonomy and claiming subtrees ...")
set.seed(seed)
partition(1L, taxonomy[[ranks[[1L]]]][[1L]], seq_len(n_seq_total))
if (!length(jobs)) fail("Partitioning produced no training jobs")
note("  ", length(jobs), " candidate job(s)")

## --- ownership --------------------------------------------------------------

# Map every enumerated subset onto the job that owns it: the job whose own
# superrank set contains the subset's superrank and whose sequences contain the
# subset's supertaxon.  Each (superrank, supertaxon) pair is owned exactly once,
# which is what makes the merge a plain concatenation.
step("Mapping each subset onto its owning chunk ...")

chunk_ids <- character(length(jobs))
width <- max(3L, nchar(length(jobs)))
for (k in seq_along(jobs)) {
  chunk_ids[[k]] <- sprintf("%s_%0*d_%s", prefix, width, k,
                            sanitize_id(jobs[[k]]$label, "chunk label"))
}

# One keyed table instead of an environment used as a hash map: the previous
# form assigned into `owner` once per (job, rank, taxon) with a pasted string
# key, then probed it with a vapply(exists) and a vapply(get) over every subset
# -- three R-level loops over the whole taxonomy
own_map <- data.table::rbindlist(lapply(seq_along(jobs), function(k) {
  job <- jobs[[k]]
  data.table::rbindlist(lapply(job$own_ranks, function(r) {
    taxon <- unique(taxonomy[[r]][job$idx])
    taxon <- taxon[!is.na(taxon)]
    if (!length(taxon)) return(NULL)
    data.table::data.table(superrank = r, supertaxon = taxon, job = k)
  }))
}))
# the assign() loop let a later job overwrite an earlier claim on the same key,
# so keep the last one; the partition is disjoint, so this should never bind
own_map <- unique(own_map, by = c("superrank", "supertaxon"), fromLast = TRUE)

job_of <- own_map[
  list(subsets$superrank, subsets$supertaxon),
  on = c("superrank", "supertaxon"),
  x.job
]
found <- !is.na(job_of)
# A supertaxon holding a single sequence gets no job -- there is no pair to
# cluster -- so its subsets are legitimately unowned.  Every subset that
# optimize_thresholds() would actually optimize must have an owner, though; a
# gap there means the partition has lost data and is a hard error.
if (!all(found | !subsets$optimized)) {
  orphan <- unique(subsets[!found & subsets$optimized,
                           c("superrank", "supertaxon")])
  fail("Partition does not cover ", nrow(orphan), " optimizable ",
       "(superrank, supertaxon) combination(s), e.g. ",
       paste(utils::head(sprintf("%s=%s", orphan$superrank, orphan$supertaxon), 5L),
             collapse = ", "))
}
subsets$chunk_id <- NA_character_
subsets$chunk_id[found] <- chunk_ids[job_of[found]]

## --- prune jobs that cannot produce a row -----------------------------------

# A job whose every subset is below --min-taxa/--min-refseq would still emit its
# own chunk-root rows, which the merge step drops again.  Skipping it here saves
# the task entirely; the long tail of tiny orders is most of the job list.
productive <- unique(subsets$chunk_id[subsets$optimized])
keep <- chunk_ids %in% productive
if (!any(keep)) fail("No training job can produce an optimized threshold; ",
                     "lower --min-taxa/--min-refseq")
dropped <- sum(!keep)
jobs <- jobs[keep]
chunk_ids <- chunk_ids[keep]
subsets <- subsets[is.na(subsets$chunk_id) | subsets$chunk_id %in% chunk_ids,
                   , drop = FALSE]

## --- outputs ----------------------------------------------------------------

step("Assembling the plan ...")

plan <- data.frame(
  chunk_id   = chunk_ids,
  job_kind   = vapply(jobs, `[[`, character(1), "kind"),
  superrank  = vapply(jobs, function(j) ranks[[j$rank_index]], character(1)),
  supertaxon = vapply(jobs, `[[`, character(1), "supertaxon"),
  # ranks[1] of the chunk must be the chunk's own root, so that
  # clean_reference_taxonomy(require_single_root = TRUE) holds and
  # optimize_thresholds() treats it as the fallback level
  chunk_ranks = vapply(jobs, function(j) {
    paste(ranks[j$rank_index:length(ranks)], collapse = ",")
  }, character(1)),
  own_superranks = vapply(jobs, function(j) paste(j$own_ranks, collapse = ","),
                          character(1)),
  n_seq     = vapply(jobs, function(j) length(j$idx), integer(1)),
  est_pairs = vapply(jobs, function(j) pairs_of(length(j$idx)), numeric(1)),
  ids_file  = paste0(chunk_ids, ".ids.txt"),
  stringsAsFactors = FALSE
)
plan <- plan[order(-plan$est_pairs), , drop = FALSE]
rownames(plan) <- NULL

expected <- subsets[order(match(subsets$rank, ranks),
                          match(subsets$superrank, ranks),
                          subsets$supertaxon), , drop = FALSE]
rownames(expected) <- NULL

id_paths <- file.path(output_dir, plan$ids_file)
outputs <- c(file.path(output_dir, c("partition.tsv", "expected_rows.tsv",
                                     "partition.json")), id_paths)

step("Writing the plan and ", nrow(plan), " ID list(s) ...")

with_atomic_outputs(outputs, overwrite = isTRUE(opt$overwrite), writer = function(staged) {
  data.table::fwrite(plan, staged[[1L]], sep = "\t", quote = FALSE, na = "")
  data.table::fwrite(expected, staged[[2L]], sep = "\t", quote = FALSE, na = "")
  write_json_file(
    list(
      schema_version = 1L,
      created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
      ranks = ranks,
      inputs = list(training_taxonomy = file_record(taxonomy_path,
                                                    sequences = n_seq_total)),
      partitioning = list(budget_pairs = budget, subsample_max = subsample_max,
                          subsample_min_per_taxon = subsample_min_per_taxon,
                          min_taxa = min_taxa, min_refseq = min_refseq,
                          seed = seed),
      counts = list(
        chunks = nrow(plan),
        exact = sum(plan$job_kind == "exact"),
        subsampled = sum(plan$job_kind == "subsampled"),
        pruned = dropped,
        expected_rows = nrow(expected),
        optimized_rows = sum(expected$optimized),
        total_pairs = sum(plan$est_pairs),
        monolithic_pairs = pairs_of(n_seq_total)
      )
    ),
    staged[[3L]]
  )
  for (k in seq_len(nrow(plan))) {
    job <- jobs[[match(plan$chunk_id[[k]], chunk_ids)]]
    writeLines(taxonomy$seq_id[job$idx], staged[[3L + k]])
  }
})

note("Planned ", nrow(plan), " training jobs (",
     sum(plan$job_kind == "exact"), " exact, ",
     sum(plan$job_kind == "subsampled"), " subsampled; ", dropped,
     " pruned as unproductive)")
note("  total ", format(sum(plan$est_pairs), big.mark = ",", scientific = FALSE),
     " pairs vs ", format(pairs_of(n_seq_total), big.mark = ",", scientific = FALSE),
     " for a single job (", round(pairs_of(n_seq_total) / sum(plan$est_pairs), 1),
     "x less)")
note("  largest job: ", plan$chunk_id[[1L]], " (", plan$n_seq[[1L]],
     " sequences)")
note("Wrote partition.tsv, expected_rows.tsv, partition.json and ",
     nrow(plan), " ID lists under ", output_dir)

############################################## Session info

finish(start_time)
