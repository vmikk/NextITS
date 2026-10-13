#!/usr/bin/env Rscript

## Script: optimotu_cluster.R
## Purpose: Taxonomically guided OTU clustering with OptimOTU: sequences plus a classification plus trained per-taxon thresholds in,
##          sequence-to-OTU membership out
##
## Usage:
##   Rscript optimotu_cluster.R \
##      --input-fasta query.fasta \
##      --thresholds reference.thresholds.tsv \
##      --tax-prob query.taxprob.tsv \
##      --output-prefix query \
##      --metric FM --confidence 0.5 --threads 4
##
## Output:
##   - <prefix>.membership.tsv
##   - <prefix>.otus.fasta
##   - <prefix>.otu_taxonomy.tsv
##   - <prefix>.sequence_taxonomy.tsv
##   - <prefix>.run.json
##   - <prefix>.uc                  (unless --no-uc)
##
##   With --split-rank R the clustering stops after rank R. Sequences still
##   classified or attached at R are written as parts of at most --split-target
##   sequences (whole rank-R taxa, never split), to be clustered separately:
##   - <split-dir>/pNNN.fasta
##   - <split-dir>/pNNN.taxprob.tsv  (ranks down to R as certain calls)
##   The outputs above then cover only the sequences finished here, those
##   clustered de novo at or above R. See optimotu_loop.R for why this is exact
##
##   An OTU is identified by its representative, the ID of its most abundant
##   member sequence, as VSEARCH identifies a cluster by its centroid. The IDs
##   stay unique when parts of a dataset are clustered separately and merged

options(stringsAsFactors = FALSE, warn = 1)

## Check time
start_time <- Sys.time()

script_dir <- dirname(normalizePath(sub(
  "^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]]
)))
source(file.path(script_dir, "optimotu_common.R"))
SCRIPT_START <- start_time

banner("optimotu_cluster.R", "Cluster sequences into OTUs with OptimOTU")

cat("Parsing input options and arguments...\n")

suppressPackageStartupMessages(require(optparse))

## Parse arguments
option_list <- list(
  make_option(c("--input-fasta"),    default=NULL, type='character', help="Query sequences [required]"),
  make_option(c("--thresholds"),     default=NULL, type='character', help="Thresholds TSV from optimotu_train.R [required]"),
  make_option(c("--tax-prob"),       default=NULL, type='character', help="Classification TSV from optimotu_classify.R [required]"),
  make_option(c("--output-prefix"),  default=NULL, type='character', help="Prefix for the output files [required]"),
  make_option(c("--metric"),         default=NULL, type='character', help="Measure to select from the thresholds file [default: FM]"),
  make_option(c("--confidence"),     default=0.5,  type='double',    help="Classifier probability cutoff [default %default]"),
  make_option(c("--force-denovo"),   default=NULL, type='character', help="Ranks whose classifier calls are ignored, comma-separated, or 'none' [default: species]"),
  make_option(c("--ranks"),          default=NULL, type='character', help="Rank names, most to least inclusive, comma-separated [default: from the thresholds sidecar, else root,kingdom,...,species]"),
  make_option(c("--dist-method"),    default="wfa2", type='character', help="One of: wfa2, ksw2, edlib, hybrid, hamming, usearch [default %default]"),
  make_option(c("--dist-penalties"), default=NULL, type='character', help="WFA2/KSW2 penalties, e.g. gap_open=2,gap_extend=1 [default: edit distance]"),
  make_option(c("--clust-method"),   default="tree", type='character', help="One of: tree, matrix, index, slink [default %default]"),
  make_option(c("--threads"),        default=1L,   type='integer',   help="Worker threads for alignment and clustering [default %default]"),
  make_option(c("--usearch"),        default=NULL, type='character', help="usearch executable for --dist-method usearch [default: autodetected from $USEARCH, then PATH]"),
  make_option(c("--chunk-id"),       default=NULL, type='character', help="Label for the chunk_id column [default: the input basename]"),
  make_option(c("--seed"),           default=42L,  type='integer',   help="RNG seed; OptimOTU breaks distance ties at random [default %default]"),
  make_option(c("--split-rank"),     default=NULL, type='character', help="Stop after this rank and write the sequences resolved there as per-taxon parts [default: cluster everything]"),
  make_option(c("--split-target"),   default=20000L, type='integer', help="With --split-rank: pack small taxa into parts of up to this many sequences [default %default]"),
  make_option(c("--split-dir"),      default=NULL, type='character', help="With --split-rank: directory for the parts [default: <prefix>.parts]"),
  make_option(c("--closed-ref-from"), default=NULL, type='character', help="Attach unclassified sequences to classified ones only at this rank and below; above it they are clustered de novo [default: every rank]"),
  make_option(c("--no-uc"),          action='store_true', default=FALSE, help="Skip the VSEARCH-style .uc file"),
  make_option(c("--dereplicate"),    action='store_true', default=FALSE, help="Collapse identical sequences instead of failing"),
  make_option(c("--force"),          action='store_true', default=FALSE, help="Proceed despite a thresholds/distance mismatch"),
  make_option(c("--verbose"),        action='store_true', default=FALSE, help="Print OptimOTU progress"),
  make_option(c("--overwrite"),      action='store_true', default=FALSE, help="Replace existing output files")
)

opt <- parse_args(OptionParser(
  option_list = option_list,
  usage = paste0("%prog --input-fasta FILE --thresholds FILE --tax-prob FILE ",
                 "--output-prefix PATH [options]"),
  description = "\nCluster sequences into OTUs with OptimOTU.\n",
  epilogue = paste0(
    "Thresholds are only valid for the distance configuration they were trained\n",
    "with, so a mismatch against the thresholds' JSON sidecar is an error unless\n",
    "--force is given.\n")
))

## Validation of required arguments
require_options(opt, c("input-fasta", "thresholds", "tax-prob", "output-prefix"))

############################################## Load packages

cat("\nLoading R packages...\n")

load_pckg("data.table")
require_package("optimotu", OPTIMOTU_MIN_VERSION)
require_package("Biostrings", "2.70.0")
require_package("jsonlite", "1.8.0")
# needs the optimotu namespace, so only after the version check
source(file.path(script_dir, "optimotu_loop.R"))

############################################## Assign and validate parameters


metric <- toupper(if (is.null(opt$metric)) "FM" else opt$metric)
if (!metric %in% MEASURES) fail("Unsupported --metric: ", metric)
confidence <- as_probability(opt$confidence, "confidence", 0.5)
dist_method <- as_choice(opt$`dist-method`, "dist-method", DIST_METHODS, "wfa2")
clust_method <- as_choice(opt$`clust-method`, "clust-method", CLUST_METHODS, "tree")
penalties <- as_kv_arg(opt$`dist-penalties`, "dist-penalties")
threads <- as_count(opt$threads, "threads", 1L)
split_rank <- opt$`split-rank`
if (!is.null(split_rank) && tolower(split_rank) %in% c("", "none")) split_rank <- NULL
split_target <- as_count(opt$`split-target`, "split-target", 20000L)
closed_ref_from <- opt$`closed-ref-from`
if (!is.null(closed_ref_from) && tolower(closed_ref_from) %in% c("", "none")) closed_ref_from <- NULL
seed <- as_count(opt$seed, "seed", 42L)
# optimotu::closed_ref_cluster() calls sample.int() to pick between equally
# distant references, so without a fixed seed the same input can give different
# OTUs from run to run.
set.seed(seed)

input_path <- normalizePath(opt$`input-fasta`, mustWork = TRUE)
thresholds_path <- normalizePath(opt$thresholds, mustWork = TRUE)
tax_prob_path <- normalizePath(opt$`tax-prob`, mustWork = TRUE)
prefix <- path.expand(opt$`output-prefix`)
split_dir <- path.expand(if (is.null(opt$`split-dir`)) paste0(prefix, ".parts") else opt$`split-dir`)

usearch <- if (identical(dist_method, "usearch")) {
  find_executable("usearch", opt$usearch)
} else NULL

## Log assigned variables
cat("\nParameters specified:\n")
report_params(list(
  "Input FASTA"        = input_path,
  "Thresholds"         = thresholds_path,
  "Classification"     = tax_prob_path,
  "Output prefix"      = prefix,
  "Metric"             = metric,
  "Confidence"         = confidence,
  "Distance method"    = dist_method,
  "Distance penalties" = penalties,
  "Clustering method"  = clust_method,
  "Threads"            = threads,
  "Seed"               = seed,
  "Closed-ref from"    = closed_ref_from,
  "Split rank"         = split_rank,
  "Split target"       = if (is.null(split_rank)) NULL else split_target,
  "Emit UC"            = !isTRUE(opt$`no-uc`),
  "Dereplicate"        = isTRUE(opt$dereplicate),
  "Force"              = isTRUE(opt$force),
  "Verbose"            = isTRUE(opt$verbose),
  "Overwrite"          = isTRUE(opt$overwrite)
))
if (!is.null(usearch)) report_tool("usearch", usearch)
cat("\n")

############################################## Main workflow

cat("\n--- Main Workflow ---\n\n")

## --- thresholds and their sidecar -------------------------------------------

sidecar_path <- sub("\\.tsv$", ".json", thresholds_path)
sidecar <- if (identical(sidecar_path, thresholds_path) ||
               !file.exists(sidecar_path)) {
  NULL
} else {
  tryCatch(
    jsonlite::read_json(sidecar_path, simplifyVector = TRUE),
    error = function(e) {
      note("Ignoring unreadable thresholds sidecar ", sidecar_path, ": ",
           conditionMessage(e))
      NULL
    }
  )
}

ranks <- if (!is.null(opt$ranks)) {
  resolve_ranks(opt$ranks)
} else if (!is.null(sidecar$ranks)) {
  as.character(sidecar$ranks)
} else {
  DEFAULT_RANKS
}
root_rank <- ranks[[1L]]
# the ranks can come from three places, so say which one won
note("  ranks (", if (!is.null(opt$ranks)) "--ranks" else if (!is.null(sidecar$ranks)) {
  "thresholds sidecar"
} else "built-in default", "): ", paste(ranks, collapse = ","))

dist_config <- build_dist_config(dist_method, penalties, threads, usearch)

# Optimized thresholds are only meaningful for the distance algorithm they were
# optimized with, so a mismatch is an error rather than a warning.
if (!is.null(sidecar$distance$method)) {
  trained_method <- tolower(as.character(sidecar$distance$method))
  if (!identical(trained_method, dist_method)) {
    msg <- paste0("Thresholds were trained with '", trained_method,
                  "' distances but --dist-method is '", dist_method, "'")
    if (isTRUE(opt$force)) note("WARNING: ", msg, " (--force given)") else {
      fail(msg, ". Retrain, match --dist-method, or pass --force.")
    }
  }
  # Compare the effective penalties (defaults included), not only the ones
  # given on the command line
  if (dist_method %in% PENALIZED_DIST_METHODS && identical(trained_method, dist_method)) {
    shared <- intersect(names(sidecar$distance),
                        setdiff(names(dist_config), c("method", "call")))
    differing <- shared[
      vapply(shared,
             function(k) !identical(as.integer(sidecar$distance[[k]]),
                                    as.integer(dist_config[[k]])),
             logical(1))
    ]
    if (length(differing)) {
      msg <- paste0(toupper(dist_method), " penalties differ from training: ",
                    paste(differing, collapse = ", "))
      if (isTRUE(opt$force)) note("WARNING: ", msg, " (--force given)") else {
        fail(msg, ". Match the trained penalties or pass --force.")
      }
    }
  }
}

step("Reading the thresholds table ...")
# na.strings = "NA" reproduces read.delim's default exactly
thresholds <- data.table::fread(thresholds_path, sep = "\t", header = TRUE,
                                quote = "", colClasses = "character",
                                na.strings = "NA", showProgress = FALSE,
                                data.table = FALSE)
required_columns <- c("rank", "superrank", "supertaxon", "threshold")
missing_columns <- setdiff(required_columns, names(thresholds))
if (length(missing_columns) || !nrow(thresholds)) {
  fail("Thresholds file has an invalid schema or is empty; missing: ",
       paste(missing_columns, collapse = ", "))
}
thresholds$threshold <- suppressWarnings(as.numeric(thresholds$threshold))
if ("value" %in% names(thresholds)) {
  thresholds$value <- suppressWarnings(as.numeric(thresholds$value))
}
# `measure` is the current column name; `metric` is the name older OptimOTU
# threshold tables use, and calc_taxon_thresholds() accepts either.
measure_column <- intersect(c("measure", "metric"), names(thresholds))
if (length(measure_column)) {
  measure_column <- measure_column[[1L]]
  if (!metric %in% thresholds[[measure_column]]) {
    fail("Thresholds file has no '", metric, "' rows (available: ",
         paste(sort(unique(thresholds[[measure_column]])), collapse = ", "), ")")
  }
} else {
  measure_column <- NULL
  metric <- NULL
}
thresholds <- thresholds[!is.na(thresholds$threshold), , drop = FALSE]
if (!nrow(thresholds)) fail("Every threshold in the thresholds file is NA")

# The root-level rows are the universal fallback for sequences whose parent taxon
# was never optimized, so the selected measure must have one for every rank.
# Catching it here beats a confusing failure inside calc_taxon_thresholds().
selected <- if (is.null(measure_column)) thresholds else {
  thresholds[thresholds[[measure_column]] == metric, , drop = FALSE]
}
root_rows <- selected[selected$superrank == root_rank, , drop = FALSE]
missing_fallback <- setdiff(ranks[-1L], unique(root_rows$rank))
if (length(missing_fallback)) {
  fail("The thresholds file has no '", root_rank, "'-level ",
       if (is.null(metric)) "threshold" else paste0(metric, " threshold"),
       " for rank(s): ", paste(missing_fallback, collapse = ", "),
       ". Pick another --metric or retrain with a broader training set.")
}

max_clustering_distance <- max(
  optimotu::threshold_as_dist(thresholds$threshold), na.rm = TRUE
)

## --- input sequences --------------------------------------------------------

read_queries <- function(path, dereplicate) {
  seqs <- read_fasta(path, "Input FASTA")
  headers <- names(Biostrings::fasta.seqlengths(path, seqtype = "DNA"))
  ids <- names(seqs)
  text <- toupper(as.character(seqs))
  invalid <- grepl("[^ACGTRYSWKMBDHVN]", text)
  if (any(invalid)) {
    fail("Input FASTA contains gaps or non-IUPAC DNA symbols (e.g. ",
         ids[which(invalid)[1L]], ")")
  }
  abundance <- parse_sizes(headers)

  duplicated_seq <- duplicated(text)
  if (any(duplicated_seq) && !isTRUE(dereplicate)) {
    fail("Input FASTA is not dereplicated: identical sequences occur more than ",
         "once (e.g. ", ids[which(duplicated_seq)[1L]],
         "). Pass --dereplicate to collapse them.")
  }
  # Keep the first occurrence as the cluster-visible sequence and remember the
  # mapping, so membership can be expanded back to every input ID afterwards.
  representative_of <- match(text, text)
  list(
    ids = ids,
    seq = text,
    abundance = abundance,
    input_order = seq_along(ids),
    unique_index = unique(representative_of),
    derep_parent = representative_of
  )
}

step("Reading the query sequences ...")
input <- read_queries(input_path, opt$dereplicate)
n_input <- length(input$ids)
uniq <- input$unique_index
if (length(uniq) < n_input) {
  note("Collapsed ", n_input, " sequences into ", length(uniq), " unique sequences")
  # parts are written per unique sequence, and duplicates would need their
  # membership carried along; dereplicate before splitting instead
  if (!is.null(split_rank)) fail("--split-rank needs dereplicated input")
}

default_chunk <- sub("\\.(fa|fna|fas|fasta)(\\.gz)?$", "", basename(input_path),
                     ignore.case = TRUE)
chunk_id <- sanitize_id(
  if (is.null(opt$`chunk-id`)) default_chunk else opt$`chunk-id`,
  "--chunk-id"
)

## --- classification ---------------------------------------------------------

step("Reading the classification table ...")
tax_prob <- data.table::fread(tax_prob_path, sep = "\t", header = TRUE,
                              quote = "", colClasses = "character",
                              na.strings = "NA", showProgress = FALSE,
                              data.table = FALSE)
if (!all(c("seq_id", "rank", "taxon") %in% names(tax_prob))) {
  fail("Classification TSV must have columns seq_id, rank, taxon [, prob]")
}
if (!"prob" %in% names(tax_prob)) tax_prob$prob <- "1"
tax_prob$prob <- suppressWarnings(as.numeric(tax_prob$prob))
tax_prob <- tax_prob[, c("seq_id", "rank", "taxon", "prob"), drop = FALSE]
tax_prob$taxon[!nzchar(tax_prob$taxon)] <- NA_character_

classifier_calls <- tax_prob
force_denovo <- as_list_arg(opt$`force-denovo`, "species")
if (length(force_denovo) == 1L && tolower(force_denovo) %in% c("none", "")) {
  force_denovo <- character()
}
unknown_ranks <- setdiff(force_denovo, ranks)
if (length(unknown_ranks)) {
  fail("--force-denovo names rank(s) that are not in --ranks: ",
       paste(unknown_ranks, collapse = ", "))
}
if (root_rank %in% force_denovo) {
  fail("--force-denovo cannot include the first rank '", root_rank,
       "'; every sequence must be assigned there")
}
if (length(force_denovo)) {
  note("Ignoring classifier assignments at: ",
       paste(force_denovo, collapse = ", "))
  tax_prob <- tax_prob[!tax_prob$rank %in% force_denovo, , drop = FALSE]
}

# optimotu() requires tax_prob to cover exactly the sequences it is given
tax_prob <- tax_prob[tax_prob$seq_id %in% input$ids[uniq] &
                       tax_prob$rank %in% ranks, , drop = FALSE]
covered <- unique(tax_prob$seq_id[tax_prob$rank == root_rank])
if (!setequal(covered, input$ids[uniq])) {
  fail("Classification TSV does not assign every input sequence at rank '",
       root_rank, "' (", length(setdiff(input$ids[uniq], covered)), " missing)")
}
rownames(tax_prob) <- NULL

# OptimOTU groups sequences by the taxon name at each rank, and its threshold
# lookup (calc_subtaxon_thresholds()) fails on a name with two parents. That
# happens when the classifier database reuses a name under two lineages, e.g.
# decoys or placeholder kingdoms with real lower ranks; say so up front
accepted <- tax_prob[!is.na(tax_prob$taxon) & tax_prob$prob >= confidence, , drop = FALSE]
accepted <- data.table::as.data.table(accepted)
homonyms <- data.table::rbindlist(lapply(seq_along(ranks)[-1L], function(i) {
  child <- accepted[rank == ranks[[i]], list(seq_id, taxon)]
  parent <- accepted[rank == ranks[[i - 1L]], list(seq_id, parent = taxon)]
  pairs <- unique(merge(child, parent, by = "seq_id")[, list(taxon, parent)])
  pairs[, list(n = .N, parents = paste(parent, collapse = " | ")), by = "taxon"][
    n > 1L][, rank := ranks[[i]]]
}))
if (nrow(homonyms)) {
  fail("The classification assigns ", nrow(homonyms), " taxon name(s) under more ",
       "than one parent, which OptimOTU cannot cluster (e.g. ",
       paste0(utils::head(homonyms$rank, 3L), " ", utils::head(homonyms$taxon, 3L),
              " in ", utils::head(homonyms$parents, 3L), collapse = "; "),
       "). Each name must have one lineage in the classifier database.")
}

# Thresholds are looked up by taxon name, so the classifier must use the names
# the thresholds were trained on.  A called taxon without its own rows falls back
# to an ancestor's thresholds, which is normal for small taxa; but a bare name X
# where training has the qualified homonym "Parent X" means the classifier
# database was not built from the training run's reference set
supertaxa <- unique(selected[selected$superrank != root_rank,
                             c("superrank", "supertaxon"), drop = FALSE])
for (r in intersect(ranks, supertaxa$superrank)) {
  called <- unique(accepted$taxon[accepted$rank == r])
  if (!length(called)) next
  known <- unique(supertaxa$supertaxon[supertaxa$superrank == r])
  own <- called %in% known
  note("  ", r, ": ", sum(own), " of ", length(called), " called taxa have ",
       "their own thresholds", if (!all(own)) "; the rest use an ancestor's")
  qualified <- known[grepl(" ", known, fixed = TRUE)]
  suffix <- sub("^.* ", "", qualified)
  bare <- called[!own & called %in% suffix]
  if (length(bare)) {
    examples <- utils::head(bare, 3L)
    note("WARNING: ", length(bare), " ", r, " name(s) called by the classifier ",
         "exist in the thresholds only as qualified homonyms (e.g. ",
         paste0(examples, " -> ", qualified[match(examples, suffix)],
                collapse = "; "),
         "), so they fall back to an ancestor's thresholds. The classifier ",
         "database does not use the training names: use the training run's ",
         "reference set (--reference_set) or rebuild the classifier database ",
         "from it")
  }
}

## --- cluster ----------------------------------------------------------------

# clust_slink() is order-sensitive and cannot accept the randomly ordered
# distances a k-mer indexed backend produces, so it is not the default here even
# though optimotu() itself defaults to it.
clust_config <- build_clust_config(clust_method)
if (identical(clust_method, "slink") && identical(dist_method, "usearch")) {
  fail("--clust-method slink cannot be combined with --dist-method usearch: ",
       "SLINK requires distances in canonical order")
}

seq_table <- data.frame(seq_id = input$ids[uniq], seq = input$seq[uniq])
note("Clustering ", nrow(seq_table), " sequences with ", dist_method,
     " distances (measure ", if (is.null(metric)) "unspecified" else metric,
     ", confidence ", confidence, ") ...")

step("Running the OptimOTU rank loop",
     if (is.null(split_rank)) "" else paste0(" (stopping after ", split_rank, ")"), " ...")
started <- Sys.time()
result <- cluster_by_rank(
  seqs = seq_table,
  tax_prob = tax_prob,
  threshold_optima = thresholds,
  ranks = ranks,
  measure = metric,
  prob_thresh = confidence,
  dist_config = dist_config,
  clust_config = clust_config,
  parallel_config = optimotu::parallel_concurrent(threads),
  verbose = isTRUE(opt$verbose),
  stop_rank = split_rank,
  closed_ref_from = closed_ref_from
)
elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))

clustered <- as.data.frame(result$final, stringsAsFactors = FALSE)
pending <- result$pending

## --- split: keep the finished sequences here, the rest goes to the parts -----

all_input <- input
if (!is.null(split_rank)) {
  note(nrow(clustered), " sequences finished (clustered de novo at or above ",
       split_rank, "), ", nrow(pending), " left for per-", split_rank, " parts")
  keep <- input$ids %in% clustered$seq_id
  n_keep <- sum(keep)
  input <- list(
    ids = input$ids[keep],
    seq = input$seq[keep],
    abundance = input$abundance[keep],
    input_order = seq_len(n_keep),
    unique_index = seq_len(n_keep),
    derep_parent = seq_len(n_keep)
  )
  n_input <- n_keep
  uniq <- input$unique_index
  seq_table <- seq_table[keep, , drop = FALSE]

  # Whole rank-R taxa, largest first, packed next-fit into parts of up to
  # --split-target sequences; a taxon larger than that is a part on its own.
  # A taxon is keyed by its name alone, as the rank loop groups by name
  groups <- split(pending$seq_id, pending[[split_rank]])
  sizes <- lengths(groups)
  part_of_group <- integer(length(groups))
  part <- 0L
  fill <- Inf
  for (g in order(-sizes, names(groups))) {
    if (fill + sizes[[g]] > split_target) {
      part <- part + 1L
      fill <- 0
    }
    part_of_group[[g]] <- part
    fill <- fill + sizes[[g]]
  }
  n_parts <- part
  note("Packed ", length(groups), " ", split_rank, " taxa into ", n_parts, " part(s)",
       if (n_parts) paste0("; largest part ", max(tapply(sizes, part_of_group, sum)),
                           " sequences") else "")
}
if (!all(c("seq_id", ranks) %in% names(clustered))) {
  fail("OptimOTU returned an unexpected result schema")
}
clustered <- clustered[match(seq_table$seq_id, clustered$seq_id),
                       c("seq_id", ranks), drop = FALSE]
if (anyNA(clustered$seq_id) || anyNA(clustered[, ranks, drop = FALSE])) {
  fail("OptimOTU did not assign every sequence at every rank")
}

## --- OTUs -------------------------------------------------------------------

step("Forming OTUs from the clustered rank paths ...")

# Each distinct full rank path -- real taxa plus the pseudotaxon names OptimOTU
# invents for de novo groups -- is one OTU.
path_key <- do.call(paste, c(clustered[, ranks, drop = FALSE], sep = "\037"))
group_levels <- unique(path_key)
unique_group <- match(path_key, group_levels)

# expand from unique sequences back to every input sequence
group_number <- unique_group[match(input$derep_parent, uniq)]
n_otus <- length(group_levels)

# One grouped pass instead of `which(group_number == g)` per OTU, which scanned
# the whole membership vector n_otus times -- O(n_otus * n_input) on inputs where
# both grow together
#
# The representative is the most abundant member, ties broken by input order.
# setorder puts the groups in ascending order, and every group 1..n_otus is
# present, so the result comes back in group order without a reindex
rep_order <- data.table::data.table(
  row = input$input_order,
  group = group_number,
  abundance = input$abundance
)
data.table::setorder(rep_order, group, -abundance, row)
# with --split-rank the first pass can finish no sequence at all
representatives <- if (n_input) as.integer(rep_order[, row[[1L]], by = "group"]$V1) else integer(0)
representative_ids <- input$ids[representatives]
representative_id <- representative_ids[group_number]
is_representative <- input$input_order %in% representatives

membership <- data.frame(
  sequence_id = input$ids,
  representative_id = representative_id,
  is_representative = is_representative,
  chunk_id = rep(chunk_id, n_input)
)
membership <- cbind(membership,
                    clustered[match(input$derep_parent, uniq), ranks, drop = FALSE])
rownames(membership) <- NULL

otu_taxonomy <- data.frame(
  representative_id = representative_ids,
  n_sequences = tabulate(group_number, nbins = n_otus),
  # keyby sorts the integer groups ascending, so this lines up with 1..n_otus;
  # the vapply it replaces re-scanned group_number once per OTU
  total_size = as.numeric(
    data.table::data.table(group = group_number,
                           abundance = input$abundance)[
      , list(total = sum(abundance)), keyby = "group"]$total
  ),
  clustered[match(input$derep_parent[representatives], uniq), ranks, drop = FALSE]
)
rownames(otu_taxonomy) <- NULL

# Split the classifier calls by rank once, rather than re-scanning the whole
# table for every (rank, column) pair -- 2 * length(ranks) full scans before.
# The factor levels keep a rank that the classifier never assigned as an empty
# but correctly typed slice, so its column still comes out all-NA
calls_by_rank <- split(
  classifier_calls[, c("seq_id", "taxon", "prob"), drop = FALSE],
  factor(classifier_calls$rank, levels = ranks)
)
classifier_wide <- function(values_column, prefix) {
  columns <- lapply(ranks, function(rank) {
    rows <- calls_by_rank[[rank]]
    rows[[values_column]][match(input$ids, rows$seq_id)]
  })
  names(columns) <- paste0(prefix, ranks)
  as.data.frame(columns, stringsAsFactors = FALSE, check.names = FALSE)
}
sequence_taxonomy <- data.frame(sequence_id = input$ids,
                                representative_id = representative_id)
sequence_taxonomy <- cbind(
  sequence_taxonomy,
  clustered[match(input$derep_parent, uniq), ranks, drop = FALSE],
  classifier_wide("taxon", "classified_"),
  classifier_wide("prob", "prob_")
)
rownames(sequence_taxonomy) <- NULL

note("Assigned ", n_input, " sequences to ", n_otus, " OTUs in ",
     round(elapsed, 1), " s")

## --- UC ---------------------------------------------------------------------

## Render membership as a VSEARCH-style UC table.
##
## Field 4 carries each member's percent identity to its OTU representative, so
## UC parsers that expect a number there keep working; `S` and `C` records carry
## `*`, as VSEARCH writes them.  The identity is a reported value only -- OptimOTU
## clusters on taxon-specific thresholds, not on a single global identity -- and
## is measured as WFA2 global edit distance, which `align_wfa2_global()` already
## returns normalized by alignment length.
build_uc <- function() {
  if (!n_input) return(data.frame())
  identity <- rep("*", n_input)
  members <- which(!is_representative)
  if (length(members)) {
    rep_seq <- input$seq[representatives][group_number[members]]
    dist <- mapply(optimotu::align_wfa2_global, input$seq[members], rep_seq,
                   USE.NAMES = FALSE)
    identity[members] <- sprintf("%.1f", 100 * (1 - dist))
  }
  # 0-based cluster numbers ordered by first appearance, as VSEARCH does
  hit <- data.frame(
    type = ifelse(is_representative, "S", "H"),
    cluster = group_number - 1L,
    size = nchar(input$seq),
    identity = identity,
    strand = ifelse(is_representative, "*", "+"),
    f6 = "*", f7 = "*", cigar = "*",
    query = input$ids,
    target = ifelse(is_representative, "*", representative_id)
  )
  centroid <- data.frame(
    type = "C",
    cluster = seq_len(n_otus) - 1L,
    size = otu_taxonomy$n_sequences,
    identity = "*", strand = "*", f6 = "*", f7 = "*", cigar = "*",
    query = representative_ids,
    target = "*"
  )
  as.data.frame(data.table::rbindlist(list(hit, centroid), use.names = TRUE))
}

## --- outputs ----------------------------------------------------------------

suffixes <- c(".membership.tsv", ".otus.fasta", ".otu_taxonomy.tsv",
              ".sequence_taxonomy.tsv", ".run.json")
if (!isTRUE(opt$`no-uc`)) suffixes <- c(suffixes, ".uc")
outputs <- paste0(prefix, suffixes)
named <- stats::setNames(as.list(basename(outputs)), sub("^\\.", "", suffixes))

step("Writing ", length(outputs), " output file(s) ...")

with_atomic_outputs(outputs, overwrite = isTRUE(opt$overwrite), writer = function(staged) {
  write_tsv <- function(x, path) {
    data.table::fwrite(x, path, sep = "\t", quote = FALSE, na = "")
  }
  write_tsv(membership, staged[[1L]])

  rep_seqs <- Biostrings::DNAStringSet(input$seq[representatives])
  # `<id>;size=<total>`, as VSEARCH writes centroids with --sizeout
  if (n_otus) {
    names(rep_seqs) <- paste0(
      representative_ids, ";size=",
      format(otu_taxonomy$total_size, scientific = FALSE, trim = TRUE)
    )
  }
  Biostrings::writeXStringSet(rep_seqs, staged[[2L]], format = "fasta", width = 80L)

  write_tsv(otu_taxonomy, staged[[3L]])
  write_tsv(sequence_taxonomy, staged[[4L]])

  write_json_file(
    list(
      schema_version = 1L,
      created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
      input = file_record(input_path, sequences = n_input,
                          unique_sequences = length(uniq),
                          total_size = sum(input$abundance)),
      thresholds = file_record(thresholds_path,
                               sidecar = if (is.null(sidecar)) NULL else basename(sidecar_path),
                               rows = nrow(thresholds),
                               max_clustering_distance = max_clustering_distance),
      classification = file_record(tax_prob_path,
                                   assignments = nrow(classifier_calls)),
      parameters = list(chunk_id = chunk_id, ranks = ranks, metric = metric, confidence = confidence,
                        force_denovo = if (length(force_denovo)) force_denovo else NULL,
                        closed_ref_from = closed_ref_from,
                        seed = seed,
                        distance = dist_method, clustering = clust_method,
                        threads = threads),
      packages = list(R = as.character(getRversion()),
                      optimotu = as.character(utils::packageVersion("optimotu"))),
      split = if (is.null(split_rank)) NULL else list(
        rank = split_rank, target = split_target, pending = nrow(pending),
        parts = n_parts, dir = basename(split_dir)
      ),
      results = list(otus = n_otus, elapsed_seconds = elapsed),
      outputs = named
    ),
    staged[[5L]]
  )

  if (!isTRUE(opt$`no-uc`)) {
    data.table::fwrite(build_uc(), staged[[6L]], sep = "\t", quote = FALSE,
                       col.names = FALSE, na = "*")
  }
})

note("Wrote ", n_input, " memberships in ", n_otus, " OTUs with prefix ", prefix)

## --- parts ------------------------------------------------------------------

## One FASTA and one classification table per part. Ranks down to the split rank
## carry the lineage resolved here, with probability 1, so the per-part run
## treats them as known and only works below; lower ranks keep the classifier's
## own calls. Sequences keep their input order and size annotations
if (!is.null(split_rank)) {
  step("Writing ", n_parts, " part(s) to ", split_dir, " ...")
  if (dir.exists(split_dir)) {
    if (!isTRUE(opt$overwrite)) fail("Split directory already exists; use --overwrite: ", split_dir)
    unlink(split_dir, recursive = TRUE)
  }
  dir.create(split_dir, recursive = TRUE, showWarnings = FALSE)
  resolved_ranks <- ranks[seq_len(match(split_rank, ranks))]
  # p001, p002, ...; p000 is this first pass
  part_names <- sprintf(paste0("p%0", max(3L, nchar(n_parts)), "d"), seq_len(n_parts))
  for (p in seq_len(n_parts)) {
    ids <- unlist(groups[part_of_group == p], use.names = FALSE)
    idx <- sort(match(ids, all_input$ids))
    ids <- all_input$ids[idx]

    part_seqs <- Biostrings::DNAStringSet(all_input$seq[idx])
    names(part_seqs) <- paste0(ids, ";size=",
                               format(all_input$abundance[idx], scientific = FALSE, trim = TRUE))
    Biostrings::writeXStringSet(part_seqs, file.path(split_dir, paste0(part_names[[p]], ".fasta")),
                                format = "fasta", width = 80L)

    lineage <- pending[match(ids, pending$seq_id), , drop = FALSE]
    certain <- data.table::rbindlist(lapply(resolved_ranks, function(rank) {
      data.table::data.table(seq_id = ids, rank = rank, taxon = lineage[[rank]], prob = 1)
    }))
    calls <- classifier_calls[classifier_calls$seq_id %in% ids &
                                !classifier_calls$rank %in% resolved_ranks, , drop = FALSE]
    data.table::fwrite(
      data.table::rbindlist(list(certain, calls), use.names = TRUE),
      file.path(split_dir, paste0(part_names[[p]], ".taxprob.tsv")),
      sep = "\t", quote = FALSE, na = "NA"
    )
  }
}

############################################## Session info

finish(start_time)
