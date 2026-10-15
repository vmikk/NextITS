#!/usr/bin/env Rscript

## Script: optimotu_filter_thresholds.R
## Purpose: Drop or adjust trained threshold rows that the training subset
##          cannot support, after training and merging
##
##          A threshold row (rank, superrank, supertaxon, measure) can pass
##          --min-taxa/--min-refseq and still say little: a genus whose species
##          each have one sequence scores best at 0 by construction, and a row
##          trained on a few named species at a high superrank is inherited by
##          the whole unnamed part of its lineage.  Dropping such a row makes
##          optimotu::calc_taxon_thresholds() fall back to the nearest ancestor
##          that has one, so filtering needs no change upstream
##
##          Every rule is off by default, and with no rule given the output table
##          is byte-identical to the input.  Rows whose superrank is the root rank
##          are the final fallback and are never dropped.  They are changed only by
##          --root-thresholds median / --root-threshold-values, which replace them
##          with the median of the kept rows one rank down (e.g. species|genus|*)
##          or with given values.  This gives a usable fallback when the training
##          set covers one taxon only (Fungi-only UNITE), whose root rows are
##          degenerate copies of that taxon's rows
##
##          Order: subset composition -> drop rules -> value-changing rules
##          (--measure-gap-action min, --plateau-action lower, --min-threshold)
##          -> root rows (--root-thresholds, --root-threshold-values)
##          -> --monotonic clip
##
## Usage:
##   Rscript optimotu_filter_thresholds.R \
##      --thresholds reference.thresholds.tsv \
##      --training-taxonomy reference.training_taxonomy.tsv \
##      --min-taxa-ge2 2 --min-score 0.5 \
##      --output-prefix filtered/reference
##
## Output:
##   - <prefix>.thresholds.tsv              kept and modified rows, input order
##   - <prefix>.thresholds.json             input sidecar plus a "filtering" block
##   - <prefix>.thresholds.diagnostics.tsv  every input row, its subset
##                                          composition, status and reasons

options(stringsAsFactors = FALSE, warn = 1)

## Check time
start_time <- Sys.time()

script_dir <- dirname(normalizePath(sub(
  "^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]]
)))
source(file.path(script_dir, "optimotu_common.R"))
SCRIPT_START <- start_time

banner("optimotu_filter_thresholds.R", "Filter trained OptimOTU threshold rows")

cat("Parsing input options and arguments...\n")

suppressPackageStartupMessages(require(optparse))

## Parse arguments
option_list <- list(
  make_option(c("--thresholds"),        default=NULL, type='character', help="Thresholds TSV; its .json sidecar must sit next to it [required]"),
  make_option(c("--training-taxonomy"), default=NULL, type='character', help="Training taxonomy TSV (seq_id + one column per rank) the thresholds were trained on [required]"),
  make_option(c("--output-prefix"),     default=NULL, type='character', help="Prefix for <prefix>.thresholds.tsv / .json / .diagnostics.tsv [required]"),
  make_option(c("--min-taxa-ge2"),      default=NULL, type='integer',   help="Drop rows whose subset has fewer lower-rank taxa with >= 2 sequences"),
  make_option(c("--min-frac-taxa-ge2"), default=NULL, type='double',    help="Drop rows whose subset has a smaller share of lower-rank taxa with >= 2 sequences (singleton taxa pull the optimum towards 0)"),
  make_option(c("--min-score"),         default=NULL, type='double',    help="Drop rows whose best score (value) is below this"),
  make_option(c("--max-dominance"),     default=NULL, type='character', help="Drop rows whose largest lower taxon holds more than this share of the subset, per rank, e.g. species=0.9,default=0.75"),
  make_option(c("--max-measure-gap"),   default=NULL, type='double',    help="Act on rows whose threshold differs from the --measure-gap-ref threshold of the same subset by more than this"),
  make_option(c("--measure-gap-ref"),   default="MCC", type='character', help="Reference measure for --max-measure-gap [default %default]"),
  make_option(c("--measure-gap-action"), default="drop", type='character', help="drop, or min: lower the threshold to the reference's [default %default]"),
  make_option(c("--max-threshold"),     default=NULL, type='character', help="Drop rows whose threshold exceeds a per-rank cap, e.g. species=0.05,genus=0.2"),
  make_option(c("--min-threshold"),     default=NULL, type='character', help="Raise thresholds below a per-rank floor to the floor, e.g. species=0.005 (a species optimum at 0 only says the marker does not separate the named species)"),
  make_option(c("--drop-grid-edge"),    action='store_true', default=FALSE, help="Drop rows whose threshold is the largest grid value (the sidecar's dist_max)"),
  make_option(c("--drop-grid-zero"),    action='store_true', default=FALSE, help="Drop rows whose threshold is exactly 0"),
  make_option(c("--above-rank"),        default="genus", type='character', help="--min-refseq-above, --min-taxa-above and --min-coverage apply to superranks strictly above this rank [default %default]"),
  make_option(c("--min-refseq-above"),  default=NULL, type='integer',   help="Drop rows above --above-rank whose subset has fewer sequences"),
  make_option(c("--min-taxa-above"),    default=NULL, type='integer',   help="Drop rows above --above-rank whose subset has fewer lower-rank taxa"),
  make_option(c("--min-coverage"),      default=NULL, type='double',    help="Drop rows above --above-rank whose subset holds less than this share of the supertaxon's sequences"),
  make_option(c("--plateau"),           default=NULL, type='character', help="Plateau table (rank, superrank, supertaxon, measure, lo_eps01, hi_eps01, lo_eps02, hi_eps02, ...)"),
  make_option(c("--plateau-eps"),       default="0.01", type='character', help="Plateau tolerance to read from --plateau: 0.01 or 0.02 [default %default]"),
  make_option(c("--plateau-max-width"), default=NULL, type='double',    help="Act on rows whose plateau (hi - lo) is wider than this"),
  make_option(c("--plateau-action"),    default="drop", type='character', help="drop, or lower: move the threshold to the plateau's lower end [default %default]"),
  make_option(c("--monotonic"),         default="none", type='character', help="none, or clip: lower a row above the effective threshold of the next coarser rank in the same supertaxon [default %default]"),
  make_option(c("--root-thresholds"),   default="keep", type='character', help="keep, or median: replace each root-level row by the median threshold of the kept rows of its rank and measure whose superrank is the next coarser rank (species from species|genus|*) [default %default]"),
  make_option(c("--root-threshold-values"), default=NULL, type='character', help="Per-rank values for the root-level rows, applied after --root-thresholds to every measure, e.g. kingdom=0.3,species=0.015"),
  make_option(c("--verbose"),           action='store_true', default=FALSE, help="Print per-rank counts and example rows for each rule"),
  make_option(c("--overwrite"),         action='store_true', default=FALSE, help="Replace existing output files")
)

opt <- parse_args(OptionParser(
  option_list = option_list,
  usage = paste0("%prog --thresholds FILE --training-taxonomy FILE ",
                 "--output-prefix PATH [rules]"),
  description = "\nFilter trained OptimOTU threshold rows.\n",
  epilogue = paste0(
    "Every rule is off unless given; with none, the table is copied unchanged.\n",
    "Rows at the root rank are never dropped, and are changed only by\n",
    "--root-thresholds median and --root-threshold-values. A dropped row makes\n",
    "clustering fall back to the nearest ancestor that has a row. Ranks, the root\n",
    "rank and the threshold grid come from the JSON sidecar.\n")
))

## Validation of required arguments
require_options(opt, c("thresholds", "training-taxonomy", "output-prefix"))

############################################## Load packages

cat("\nLoading R packages...\n")

load_pckg("data.table")
require_package("optimotu", OPTIMOTU_MIN_VERSION)
require_package("jsonlite", "1.8.0")

############################################## Assign and validate parameters

thresholds_path <- normalizePath(opt$thresholds, mustWork = TRUE)
sidecar_path <- sub("\\.tsv$", ".json", thresholds_path)
if (identical(sidecar_path, thresholds_path) || !file.exists(sidecar_path)) {
  fail("No JSON sidecar next to the thresholds table: ", sidecar_path)
}
taxonomy_path <- normalizePath(opt$`training-taxonomy`, mustWork = TRUE)
plateau_path <- if (is.null(opt$plateau)) NULL else {
  normalizePath(opt$plateau, mustWork = TRUE)
}
prefix <- path.expand(opt$`output-prefix`)

# The sidecar is read twice: simplified for the values used here, and as plain
# lists for the copy written back out, so that its arrays stay arrays (a
# simplified length-1 vector would be written as a scalar)
sidecar <- jsonlite::read_json(sidecar_path, simplifyVector = TRUE)
sidecar_raw <- jsonlite::read_json(sidecar_path, simplifyVector = FALSE)
ranks <- as.character(sidecar$ranks)
if (length(ranks) < 2L) fail("The thresholds sidecar lists no ranks: ", sidecar_path)
root_rank <- sidecar$root_rank %||% ranks[[1L]]
if (!identical(root_rank, ranks[[1L]])) {
  fail("The sidecar's root_rank (", root_rank, ") is not its first rank")
}
below_root <- ranks[-1L]
dist_max <- sidecar$optimization$dist_max
dist_step <- sidecar$optimization$dist_step %||% 0.001

## Per-rank limits: `species=0.9,default=0.75` -> c(species = 0.9, default = 0.75)
as_rank_limits <- function(x, name) {
  if (is.null(x)) return(NULL)
  spec <- as_kv_arg(x, name)
  if (!length(spec)) return(NULL)
  unknown <- setdiff(names(spec), c(below_root, "default"))
  if (length(unknown)) {
    fail("--", name, " names unknown rank(s): ", paste(unknown, collapse = ", "),
         " (use ", paste(below_root, collapse = ", "), " or default)")
  }
  value <- suppressWarnings(as.numeric(spec))
  if (anyNA(value) || any(!is.finite(value))) {
    fail("--", name, " values must be finite numbers")
  }
  stats::setNames(value, names(spec))
}

## The limit that applies to each row's rank; NA where the spec has neither the
## rank nor a default, i.e. the rule does not apply there
limit_for <- function(limits, rank) {
  out <- unname(limits[rank])
  if ("default" %in% names(limits)) out[is.na(out)] <- limits[["default"]]
  out
}

rules <- list(
  min_taxa_ge2     = if (is.null(opt$`min-taxa-ge2`)) NULL else {
    as_count(opt$`min-taxa-ge2`, "min-taxa-ge2", 1L)
  },
  min_frac_taxa_ge2 = if (is.null(opt$`min-frac-taxa-ge2`)) NULL else {
    as_probability(opt$`min-frac-taxa-ge2`, "min-frac-taxa-ge2", NULL)
  },
  min_score        = if (is.null(opt$`min-score`)) NULL else {
    as_number(opt$`min-score`, "min-score", NULL, -Inf)
  },
  max_dominance    = as_rank_limits(opt$`max-dominance`, "max-dominance"),
  max_measure_gap  = if (is.null(opt$`max-measure-gap`)) NULL else {
    as_number(opt$`max-measure-gap`, "max-measure-gap", NULL, 0)
  },
  measure_gap_ref  = toupper(as_choice(opt$`measure-gap-ref`, "measure-gap-ref",
                                       tolower(MEASURES), "mcc")),
  measure_gap_action = as_choice(opt$`measure-gap-action`, "measure-gap-action",
                                 c("drop", "min"), "drop"),
  max_threshold    = as_rank_limits(opt$`max-threshold`, "max-threshold"),
  min_threshold    = as_rank_limits(opt$`min-threshold`, "min-threshold"),
  drop_grid_edge   = isTRUE(opt$`drop-grid-edge`),
  drop_grid_zero   = isTRUE(opt$`drop-grid-zero`),
  above_rank       = as_choice(opt$`above-rank`, "above-rank", below_root, "genus"),
  min_refseq_above = if (is.null(opt$`min-refseq-above`)) NULL else {
    as_count(opt$`min-refseq-above`, "min-refseq-above", 1L)
  },
  min_taxa_above   = if (is.null(opt$`min-taxa-above`)) NULL else {
    as_count(opt$`min-taxa-above`, "min-taxa-above", 1L)
  },
  min_coverage     = if (is.null(opt$`min-coverage`)) NULL else {
    as_probability(opt$`min-coverage`, "min-coverage", NULL)
  },
  plateau_eps      = as_choice(opt$`plateau-eps`, "plateau-eps", c("0.01", "0.02"), "0.01"),
  plateau_max_width = if (is.null(opt$`plateau-max-width`)) NULL else {
    as_number(opt$`plateau-max-width`, "plateau-max-width", NULL, 0)
  },
  plateau_action   = as_choice(opt$`plateau-action`, "plateau-action",
                               c("drop", "lower"), "drop"),
  monotonic        = as_choice(opt$monotonic, "monotonic", c("none", "clip"), "none"),
  root_thresholds  = as_choice(opt$`root-thresholds`, "root-thresholds", c("keep", "median"), "keep"),
  root_values      = as_rank_limits(opt$`root-threshold-values`, "root-threshold-values")
)
if ("default" %in% names(rules$root_values)) {
  fail("--root-threshold-values takes ranks only, not default")
}
if (!is.null(rules$root_values) && any(rules$root_values < 0)) {
  fail("--root-threshold-values must not be negative")
}
if (!is.null(rules$plateau_max_width) && is.null(plateau_path)) {
  fail("--plateau-max-width needs --plateau")
}
if (rules$drop_grid_edge && is.null(dist_max)) {
  fail("--drop-grid-edge needs optimization.dist_max in the sidecar")
}

## Log assigned variables
cat("\nParameters specified:\n")
report_params(list(
  "Thresholds"        = thresholds_path,
  "Sidecar"           = sidecar_path,
  "Training taxonomy" = taxonomy_path,
  "Output prefix"     = prefix,
  "Ranks"             = ranks,
  "Min taxa >= 2 seq" = rules$min_taxa_ge2,
  "Min share >= 2 seq" = rules$min_frac_taxa_ge2,
  "Min score"         = rules$min_score,
  "Max dominance"     = rules$max_dominance,
  "Max measure gap"   = rules$max_measure_gap,
  "Measure gap ref"   = if (!is.null(rules$max_measure_gap)) rules$measure_gap_ref,
  "Measure gap action" = if (!is.null(rules$max_measure_gap)) rules$measure_gap_action,
  "Max threshold"     = rules$max_threshold,
  "Min threshold"     = rules$min_threshold,
  "Drop grid edge"    = rules$drop_grid_edge,
  "Drop grid zero"    = rules$drop_grid_zero,
  "Above rank"        = rules$above_rank,
  "Min refseq above"  = rules$min_refseq_above,
  "Min taxa above"    = rules$min_taxa_above,
  "Min coverage"      = rules$min_coverage,
  "Plateau"           = plateau_path,
  "Plateau eps"       = if (!is.null(plateau_path)) rules$plateau_eps,
  "Plateau max width" = rules$plateau_max_width,
  "Plateau action"    = if (!is.null(rules$plateau_max_width)) rules$plateau_action,
  "Monotonic"         = rules$monotonic,
  "Root thresholds"   = rules$root_thresholds,
  "Root values"       = rules$root_values,
  "Verbose"           = isTRUE(opt$verbose),
  "Overwrite"         = isTRUE(opt$overwrite)
))
cat("\n")

############################################## Main workflow

cat("\n--- Main Workflow ---\n\n")

KEY <- c("rank", "superrank", "supertaxon", "measure")
SUBSET_KEY <- c("rank", "superrank", "supertaxon")

# Two thresholds within this distance are the same grid point
TOL <- 1e-9

## --- thresholds -------------------------------------------------------------

step("Reading the thresholds table ...")
# Every column is read as text and written back as it was, so a row that no
# rule changes keeps its exact bytes; na.strings = "" round-trips empty cells
table_in <- data.table::fread(thresholds_path, sep = "\t", header = TRUE,
                              quote = "", colClasses = "character",
                              check.names = FALSE, na.strings = "",
                              showProgress = FALSE)
missing_columns <- setdiff(c(KEY, "threshold", "value"), names(table_in))
if (length(missing_columns)) {
  fail("The thresholds table lacks required columns: ",
       paste(missing_columns, collapse = ", "))
}
n_rows <- nrow(table_in)
if (!n_rows) fail("The thresholds table is empty")
unknown <- setdiff(unique(c(table_in$rank, table_in$superrank)), ranks)
if (length(unknown)) {
  fail("The thresholds table uses ranks missing from its sidecar: ",
       paste(unknown, collapse = ", "))
}
if (anyDuplicated(table_in, by = KEY)) {
  fail("The thresholds table has duplicate (rank, superrank, supertaxon, ",
       "measure) keys")
}
threshold_in <- suppressWarnings(as.numeric(table_in$threshold))
value <- suppressWarnings(as.numeric(table_in$value))
measures <- unique(table_in$measure)
is_root <- table_in$superrank == root_rank
note("  ", n_rows, " rows (", sum(is_root), " at ", root_rank, ", never ",
     "filtered); measures: ", paste(measures, collapse = ", "))

gaps_in <- root_fallback_gaps(table_in, ranks, measures)
complete_in <- measures[lengths(gaps_in) == 0L]

## --- subset composition -----------------------------------------------------

# Read exactly as optimotu_partition_reference.R does: the table is already
# cleaned and homonym-qualified, and cleaning it again could change the names
# the thresholds are keyed by
step("Reading the training taxonomy ...")
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
for (rank in ranks) {
  v <- taxonomy[[rank]]
  v[!nzchar(trimws(v))] <- NA_character_
  data.table::set(taxonomy, j = rank, value = v)
}
n_seq_total <- nrow(taxonomy)
if (!n_seq_total) fail("Training taxonomy is empty")
root_taxon <- unique(taxonomy[[root_rank]])
if (length(root_taxon) != 1L || anyNA(root_taxon)) {
  fail("Training taxonomy must have exactly one ", root_rank, " value")
}
note("  ", n_seq_total, " sequences")

# summarize_by_rank() is the enumeration optimize_thresholds() itself uses: every
# (superrank, rank) pair, skipping sequences unassigned at either rank.  Its
# true_partition is 0-based, so tabulate(x + 1) gives the size of each lower
# taxon.  The counts describe the whole training taxonomy; a row from a
# subsampled training chunk was optimized on fewer sequences than this
step("Summarizing the composition of every subset ...")
data.table::setDF(taxonomy)
subsets <- optimotu::summarize_by_rank(taxonomy, ranks, "seq_id")
sizes <- lapply(subsets$true_partition, function(x) tabulate(x + 1L))
composition <- data.table::data.table(
  rank           = as.character(subsets$rank),
  superrank      = as.character(subsets$superrank),
  supertaxon     = as.character(subsets$supertaxon),
  n_seq          = as.integer(subsets$n_seq),
  n_taxa         = as.integer(subsets$n_taxa),
  n_taxa_ge2     = vapply(sizes, function(k) sum(k >= 2L), integer(1)),
  # as.numeric: k * (k - 1) overflows an integer for a root-level subset
  n_pairs_within = vapply(sizes, function(k) {
    k <- as.numeric(k)
    sum(k * (k - 1) / 2)
  }, numeric(1)),
  largest        = vapply(sizes, max, integer(1))
)
if (!identical(composition$n_seq, vapply(sizes, sum, integer(1))) ||
    !identical(composition$n_taxa, vapply(sizes, function(k) sum(k > 0L), integer(1)))) {
  fail("optimotu::summarize_by_rank() partitions disagree with its n_seq/n_taxa")
}
composition[, dominance := largest / n_seq]
composition[, largest := NULL]
rm(subsets, sizes)
invisible(gc())

# Every sequence under the supertaxon, assigned at the lower rank or not
totals <- data.table::rbindlist(lapply(ranks[-length(ranks)], function(r) {
  v <- taxonomy[[r]]
  data.table::data.table(supertaxon = v[!is.na(v)])[
    , list(n_total = .N), by = "supertaxon"][, superrank := r]
}))
composition <- totals[composition, on = c("superrank", "supertaxon")]
composition[, coverage := n_seq / n_total]
note("  ", nrow(composition), " subsets")

stats <- composition[table_in[, SUBSET_KEY, with = FALSE], on = SUBSET_KEY,
                     list(n_seq, n_taxa, n_taxa_ge2, n_pairs_within, dominance,
                          n_total, coverage)]
composition_missing <- !is_root & is.na(stats$n_seq)
if (any(composition_missing)) {
  note("NOTE: ", sum(composition_missing), " row(s) have a subset the training ",
       "taxonomy does not contain; composition rules leave them alone (e.g. ",
       paste(utils::head(unique(do.call(paste, c(
         table_in[composition_missing, SUBSET_KEY, with = FALSE],
         list(sep = "|"))
       )), 3L), collapse = "; "), ")")
}

## --- plateau widths ---------------------------------------------------------

plateau <- NULL
plateau_missing <- 0L
plateau_na <- 0L
if (!is.null(plateau_path)) {
  step("Reading the plateau table ...")
  suffix <- sub("^0\\.", "", rules$plateau_eps)
  lo_col <- paste0("lo_eps", suffix)
  hi_col <- paste0("hi_eps", suffix)
  plateau_tab <- data.table::fread(plateau_path, sep = "\t", header = TRUE,
                                   quote = "", colClasses = "character",
                                   check.names = FALSE, na.strings = c("", "NA"),
                                   showProgress = FALSE)
  missing_columns <- setdiff(c(KEY, lo_col, hi_col), names(plateau_tab))
  if (length(missing_columns)) {
    fail("The plateau table lacks required columns: ",
         paste(missing_columns, collapse = ", "))
  }
  if (anyDuplicated(plateau_tab, by = KEY)) {
    fail("The plateau table has duplicate (rank, superrank, supertaxon, measure) keys")
  }
  plateau_tab <- plateau_tab[, list(rank, superrank, supertaxon, measure,
                                    lo = suppressWarnings(as.numeric(get(lo_col))),
                                    hi = suppressWarnings(as.numeric(get(hi_col))),
                                    present = TRUE)]
  plateau <- plateau_tab[table_in[, KEY, with = FALSE], on = KEY,
                         list(lo, hi, width = hi - lo, present = !is.na(present))]
  plateau_missing <- sum(!is_root & !plateau$present)
  plateau_na <- sum(!is_root & plateau$present & is.na(plateau$width))
  note("  ", sum(plateau$present), " of ", n_rows, " rows have a plateau (",
       plateau_missing, " non-", root_rank, " rows absent, left alone; ",
       plateau_na, " with no width)")
}

## --- drop rules -------------------------------------------------------------

step("Applying the drop rules ...")

# Every rule is evaluated on the input, independently of the others, so a row
# records each rule it fails and not only the first
DROP_RULES <- c("min_taxa_ge2", "min_frac_taxa_ge2", "min_score", "max_dominance", "max_threshold",
                "grid_edge", "grid_zero", "min_refseq_above", "min_taxa_above",
                "min_coverage", "measure_gap", "plateau_width")
MODIFY_RULES <- c("measure_gap_min", "plateau_lower", "min_threshold", "root_median",
                  "root_value", "monotonic_clip")
fails <- matrix(FALSE, n_rows, length(DROP_RULES),
                dimnames = list(NULL, DROP_RULES))
# NA (no composition, no reference row, no plateau) never fails a rule, and the
# root-level fallback rows are exempt from all of them
hit_of <- function(x) !is.na(x) & x & !is_root

above <- match(table_in$superrank, ranks) < match(rules$above_rank, ranks)

if (!is.null(rules$min_taxa_ge2)) {
  fails[, "min_taxa_ge2"] <- hit_of(stats$n_taxa_ge2 < rules$min_taxa_ge2)
}
if (!is.null(rules$min_frac_taxa_ge2)) {
  fails[, "min_frac_taxa_ge2"] <- hit_of(
    stats$n_taxa_ge2 / stats$n_taxa < rules$min_frac_taxa_ge2 - TOL
  )
}
if (!is.null(rules$min_score)) {
  fails[, "min_score"] <- hit_of(value < rules$min_score)
}
if (!is.null(rules$max_dominance)) {
  fails[, "max_dominance"] <- hit_of(
    stats$dominance > limit_for(rules$max_dominance, table_in$rank) + TOL
  )
}
if (!is.null(rules$max_threshold)) {
  fails[, "max_threshold"] <- hit_of(
    threshold_in > limit_for(rules$max_threshold, table_in$rank) + TOL
  )
}
if (rules$drop_grid_edge) {
  fails[, "grid_edge"] <- hit_of(abs(threshold_in - dist_max) < TOL)
}
if (rules$drop_grid_zero) {
  fails[, "grid_zero"] <- hit_of(abs(threshold_in) < TOL)
}
if (!is.null(rules$min_refseq_above)) {
  fails[, "min_refseq_above"] <- hit_of(above & stats$n_seq < rules$min_refseq_above)
}
if (!is.null(rules$min_taxa_above)) {
  fails[, "min_taxa_above"] <- hit_of(above & stats$n_taxa < rules$min_taxa_above)
}
if (!is.null(rules$min_coverage)) {
  fails[, "min_coverage"] <- hit_of(above & stats$coverage < rules$min_coverage - TOL)
}

# The reference measure's threshold for the same subset.  Where it has none (MCC
# is undefined on a single-taxon subset) this rule cannot judge the row, and the
# composition rules are the ones that should
ref_threshold <- NULL
ref_missing <- NULL
if (!is.null(rules$max_measure_gap)) {
  ref <- data.table::data.table(table_in[, SUBSET_KEY, with = FALSE],
                                ref = threshold_in)[table_in$measure == rules$measure_gap_ref]
  ref_threshold <- ref[table_in[, SUBSET_KEY, with = FALSE], on = SUBSET_KEY, x.ref]
  ref_threshold[table_in$measure == rules$measure_gap_ref] <- NA_real_
  ref_missing <- !is_root & table_in$measure != rules$measure_gap_ref & is.na(ref_threshold)
  gap <- abs(threshold_in - ref_threshold)
  if (identical(rules$measure_gap_action, "drop")) {
    fails[, "measure_gap"] <- hit_of(gap > rules$max_measure_gap + TOL)
  }
}
if (!is.null(rules$plateau_max_width) && identical(rules$plateau_action, "drop")) {
  fails[, "plateau_width"] <- hit_of(plateau$width > rules$plateau_max_width + TOL)
}

dropped <- rowSums(fails) > 0L
first_reason <- rep(NA_character_, n_rows)
first_reason[dropped] <- DROP_RULES[max.col(fails[dropped, , drop = FALSE],
                                            ties.method = "first")]
note("  dropped ", sum(dropped), " of ", n_rows, " rows")

## --- value-changing rules ---------------------------------------------------

# These act only on rows that survived the drop rules; all but --min-threshold
# only ever lower a threshold
threshold_out <- threshold_in
threshold_out[dropped] <- NA_real_
modified_by <- matrix(FALSE, n_rows, length(MODIFY_RULES),
                      dimnames = list(NULL, MODIFY_RULES))

if (!is.null(rules$max_measure_gap) && identical(rules$measure_gap_action, "min")) {
  step("Lowering thresholds that exceed the ", rules$measure_gap_ref,
       " threshold by more than ", rules$max_measure_gap, " ...")
  hit <- hit_of(!dropped & gap > rules$max_measure_gap + TOL &
                  threshold_out > ref_threshold + TOL)
  threshold_out[hit] <- ref_threshold[hit]
  modified_by[, "measure_gap_min"] <- hit
}

if (!is.null(rules$plateau_max_width) && identical(rules$plateau_action, "lower")) {
  step("Moving thresholds on plateaus wider than ", rules$plateau_max_width,
       " to the plateau's lower end ...")
  hit <- hit_of(!dropped & plateau$width > rules$plateau_max_width + TOL &
                  threshold_out > plateau$lo + TOL)
  threshold_out[hit] <- plateau$lo[hit]
  modified_by[, "plateau_lower"] <- hit
}

# The floor is the one rule that raises a threshold.  It runs before the clip,
# so a floored row still cannot end above the next coarser rank
if (!is.null(rules$min_threshold)) {
  step("Raising thresholds below the per-rank floor ...")
  floor_of <- limit_for(rules$min_threshold, table_in$rank)
  hit <- hit_of(!dropped & !is.na(floor_of) & threshold_out < floor_of - TOL)
  threshold_out[hit] <- floor_of[hit]
  modified_by[, "min_threshold"] <- hit
}

## --- root-level rows --------------------------------------------------------

# A root row (rank R, superrank root, measure M) is the threshold for any
# sequence whose lineage has no row of its own, e.g. a non-fungal sequence when
# only Fungi were trained.  "median" sets it to the median of the kept rows at
# rank R and measure M whose superrank is the rank just above R: species from
# every species|genus|* row, genus from genus|family|*, and so on.  For the first
# rank below the root that parent is the root itself, so there is nothing to
# take the median of and the trained row stays (with a warning) unless
# --root-threshold-values gives it a value.  It runs after the floor and before
# the clip, so the clip compares against the new values
root_log <- NULL
if (identical(rules$root_thresholds, "median") || !is.null(rules$root_values)) {
  step("Replacing the ", root_rank, "-level rows ...")
  root_idx <- which(is_root)
  parent_of <- ranks[match(table_in$rank[root_idx], ranks) - 1L]
  n_source <- rep(NA_integer_, length(root_idx))
  method <- rep("kept", length(root_idx))
  if (identical(rules$root_thresholds, "median")) {
    source_ok <- !dropped & !is_root & !is.na(threshold_out)
    src <- data.table::data.table(rank = table_in$rank, superrank = table_in$superrank,
                                  measure = table_in$measure,
                                  thr = threshold_out)[source_ok]
    med <- src[, list(med = stats::median(thr), n = .N), by = c("rank", "superrank", "measure")]
    hit_med <- med[data.table::data.table(rank = table_in$rank[root_idx],
                                          superrank = parent_of,
                                          measure = table_in$measure[root_idx]),
                   on = c("rank", "superrank", "measure"), list(med, n)]
    n_source <- hit_med$n
    ok <- !is.na(hit_med$med) & parent_of != root_rank
    threshold_out[root_idx[ok]] <- hit_med$med[ok]
    modified_by[root_idx[ok], "root_median"] <- TRUE
    method[ok] <- "median"
    no_source <- unique(table_in$rank[root_idx[!ok]])
    if (length(no_source)) {
      note("WARNING: no rows to take the median of for ", root_rank, "-level ",
           paste(no_source, collapse = ", "), "; ",
           if (all(no_source %in% names(rules$root_values))) "--root-threshold-values sets them"
           else "the trained rows are kept")
    }
  }
  if (!is.null(rules$root_values)) {
    given <- unname(rules$root_values[table_in$rank[root_idx]])
    hv <- !is.na(given)
    threshold_out[root_idx[hv]] <- given[hv]
    modified_by[root_idx[hv], "root_value"] <- TRUE
    modified_by[root_idx[hv], "root_median"] <- FALSE
    method[hv] <- "value"
  }
  root_log <- data.table::data.table(
    rank = table_in$rank[root_idx], measure = table_in$measure[root_idx],
    method = method, n_source_rows = n_source,
    threshold_in = threshold_in[root_idx], threshold_out = threshold_out[root_idx]
  )
  for (r in unique(root_log$rank)) {
    x <- root_log[rank == r]
    note("  ", formatC(r, width = -10), " ", paste(sprintf("%s %s->%s (%s, n=%s)",
      x$measure, format(x$threshold_in, digits = 4), format(x$threshold_out, digits = 4),
      x$method, ifelse(is.na(x$n_source_rows), "-", x$n_source_rows)), collapse = "; "))
  }
}

## --- monotonic clip ---------------------------------------------------------

# Within one supertaxon S at superrank R, a rank's threshold should not exceed
# that of the next coarser rank: calc_subtaxon_thresholds() would otherwise
# raise the coarser one with cummax.  The coarser rank's *effective* threshold is
# the one clustering would use: S's own row, else the nearest ancestor's, walking
# up S's lineage exactly as calc_taxon_thresholds() does.  Superranks and ranks
# are processed coarse to fine, so a clipped value is what the finer rows below
# it are compared with
if (identical(rules$monotonic, "clip")) {
  step("Clipping thresholds to the next coarser rank ...")
  work <- data.table::data.table(table_in[, KEY, with = FALSE], thr = threshold_out)
  tax_dt <- data.table::as.data.table(taxonomy)
  lineage <- lapply(seq_along(ranks), function(i) {
    if (i < 2L) return(NULL)
    unique(tax_dt[!is.na(tax_dt[[ranks[[i]]]]), ranks[seq_len(i)], with = FALSE],
           by = ranks[[i]])
  })
  rm(tax_dt)
  ancestor_of <- function(taxa, from, to) {
    if (to == from) return(taxa)
    if (to == 1L) return(rep(root_taxon, length(taxa)))
    lin <- lineage[[from]]
    lin[[ranks[[to]]]][match(taxa, lin[[ranks[[from]]]])]
  }
  # the root-level rows are not clipped, so the walk starts one rank below
  for (i in seq_len(length(ranks) - 1L)[-1L]) {
    if (i + 2L > length(ranks)) next
    for (j in seq.int(i + 2L, length(ranks))) {
      target <- which(work$superrank == ranks[[i]] & work$rank == ranks[[j]] &
                        !is.na(work$thr))
      if (!length(target)) next
      parent <- ranks[[j - 1L]]
      eff <- rep(NA_real_, length(target))
      for (k in rev(seq_len(i))) {
        need <- which(is.na(eff))
        if (!length(need)) break
        query <- data.table::data.table(
          rank = parent, superrank = ranks[[k]],
          supertaxon = ancestor_of(work$supertaxon[target[need]], i, k),
          measure = work$measure[target[need]]
        )
        eff[need] <- work[query, on = KEY, x.thr]
      }
      hit <- !is.na(eff) & work$thr[target] > eff + TOL
      if (!any(hit)) next
      data.table::set(work, i = target[hit], j = "thr", value = eff[hit])
      modified_by[target[hit], "monotonic_clip"] <- TRUE
    }
  }
  threshold_out <- work$thr
  rm(work, lineage)
}

modified <- !dropped & abs(threshold_out - threshold_in) > TOL
modified[is.na(modified)] <- FALSE
status <- ifelse(dropped, "dropped", ifelse(modified, "modified", "kept"))
first_reason[!dropped & modified] <-
  MODIFY_RULES[max.col(modified_by[!dropped & modified, , drop = FALSE],
                       ties.method = "first")]
all_flags <- cbind(fails, modified_by)
reasons <- apply(all_flags, 1L, function(x) paste(colnames(all_flags)[x], collapse = ","))

note("  kept ", sum(status == "kept"), ", modified ", sum(status == "modified"),
     ", dropped ", sum(status == "dropped"))

## Count rows per rule and rank, as a nested rule -> rank -> n list
by_rule_rank <- function(flags) {
  out <- lapply(colnames(flags), function(rule) {
    hit <- flags[, rule]
    if (!any(hit)) return(NULL)
    n <- table(factor(table_in$rank[hit], levels = below_root))
    as.list(n[n > 0L])
  })
  names(out) <- colnames(flags)
  out[!vapply(out, is.null, logical(1))]
}
first_flags <- vapply(DROP_RULES, function(r) dropped & first_reason %in% r,
                      logical(n_rows))
dim(first_flags) <- c(n_rows, length(DROP_RULES))
colnames(first_flags) <- DROP_RULES
dropped_counts <- by_rule_rank(first_flags)
failed_counts <- by_rule_rank(fails)
modified_counts <- by_rule_rank(modified_by & !dropped)

for (rule in DROP_RULES[colSums(fails) > 0L]) {
  note("  ", formatC(rule, width = -17), " fails ", sum(fails[, rule]),
       " row(s), first reason for ", sum(first_flags[, rule]))
}
for (rule in MODIFY_RULES[colSums(modified_by) > 0L]) {
  note("  ", formatC(rule, width = -17), " changes ", sum(modified_by[, rule]),
       " row(s)")
}
if (!is.null(ref_missing)) {
  note("  ", sum(ref_missing), " row(s) have no ", rules$measure_gap_ref,
       " row for their subset; --max-measure-gap leaves them alone")
}
if (isTRUE(opt$verbose)) {
  # the example rows carry every composition column; keep each on one line
  old_width <- options(width = 200L)
  flags <- cbind(fails, modified_by)
  for (rule in colnames(flags)[colSums(flags) > 0L]) {
    hit <- which(flags[, rule])
    note("\n  ", rule, " by rank and measure:")
    print(table(rank = factor(table_in$rank[hit], levels = below_root),
                measure = table_in$measure[hit]))
    note("  e.g.")
    print(utils::head(data.frame(table_in[hit, KEY, with = FALSE],
                                 threshold = threshold_in[hit], value = value[hit],
                                 stats[hit], check.names = FALSE), 10L),
          row.names = FALSE)
  }
  options(old_width)
}

## --- root-level fallback ----------------------------------------------------

step("Checking the ", root_rank, "-level fallback coverage ...")

# Root-level rows are never dropped, so this cannot fail on account of the
# rules; it re-asserts what optimotu_cluster.R requires of the table it is given
kept <- !dropped
gaps_out <- root_fallback_gaps(table_in[kept], ranks, measures)
complete_out <- measures[lengths(gaps_out) == 0L]
lost <- setdiff(complete_in, complete_out)
if (length(lost)) {
  fail("Filtering removed the ", root_rank, "-level fallback of measure(s) ",
       paste(lost, collapse = ", "), "; this is a bug, root-level rows must ",
       "never be dropped")
}
if (!length(complete_out)) {
  fail("The thresholds table has no ", root_rank, "-level threshold for every ",
       "rank under any measure (missing at ",
       paste(sort(unique(unlist(gaps_out))), collapse = ", "), ")")
}
note("Measures usable as --metric: ", paste(complete_out, collapse = ", "))

## --- outputs ----------------------------------------------------------------

step("Writing the filtered thresholds ...")

## Write a changed threshold the way fwrite() writes a double: plain decimal, up
## to 15 significant digits, no trailing zeros
format_threshold <- function(x) {
  vapply(x, function(v) format(v, digits = 15L, scientific = FALSE,
                               drop0trailing = TRUE, trim = TRUE),
         character(1), USE.NAMES = FALSE)
}

threshold_text <- table_in$threshold
threshold_text[modified] <- format_threshold(threshold_out[modified])
threshold_text[dropped] <- NA_character_

# A replaced root row's score belonged to the trained threshold, not the new one
root_replaced <- is_root & modified
value_text <- table_in$value
value_text[root_replaced] <- NA_character_

out_table <- data.table::copy(table_in)
data.table::set(out_table, j = "threshold", value = threshold_text)
data.table::set(out_table, j = "value", value = value_text)
out_table <- out_table[!dropped]

diagnostics <- data.table::data.table(
  table_in[, KEY, with = FALSE],
  threshold_in = table_in$threshold,
  threshold = threshold_text,
  value = table_in$value,
  stats
)
if (!is.null(ref_threshold)) {
  diagnostics[, ref_threshold := ref_threshold]
  diagnostics[, ref_missing := ref_missing]
}
if (!is.null(plateau)) {
  diagnostics[, plateau_lo := plateau$lo]
  diagnostics[, plateau_hi := plateau$hi]
  diagnostics[, plateau_width := plateau$width]
}
diagnostics[, status := status]
diagnostics[, first_reason := first_reason]
diagnostics[, reasons := reasons]

max_clustering_distance <- max(
  optimotu::threshold_as_dist(threshold_out[!dropped]), na.rm = TRUE
)

outputs <- paste0(prefix, c(".thresholds.tsv", ".thresholds.json",
                            ".thresholds.diagnostics.tsv"))

# The input sidecar describes the training run, which filtering does not change;
# only the counts and the clustering distance follow the table
sidecar_out <- sidecar_raw
sidecar_out$counts <- sidecar_out$counts %||% list()
sidecar_out$counts$thresholds <- nrow(out_table)
sidecar_out$counts$complete_measures <- complete_out
sidecar_out$max_clustering_distance <- max_clustering_distance
sidecar_out$filtering <- list(
  created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  inputs = list(
    thresholds = file_record(thresholds_path, rows = n_rows),
    sidecar = file_record(sidecar_path),
    plateau = if (is.null(plateau_path)) NULL else file_record(plateau_path)
  ),
  rules = list(
    min_taxa_ge2 = rules$min_taxa_ge2,
    min_frac_taxa_ge2 = rules$min_frac_taxa_ge2,
    min_score = rules$min_score,
    max_dominance = if (is.null(rules$max_dominance)) NULL else as.list(rules$max_dominance),
    max_measure_gap = rules$max_measure_gap,
    measure_gap_ref = if (is.null(rules$max_measure_gap)) NULL else rules$measure_gap_ref,
    measure_gap_action = if (is.null(rules$max_measure_gap)) NULL else rules$measure_gap_action,
    max_threshold = if (is.null(rules$max_threshold)) NULL else as.list(rules$max_threshold),
    min_threshold = if (is.null(rules$min_threshold)) NULL else as.list(rules$min_threshold),
    drop_grid_edge = rules$drop_grid_edge,
    drop_grid_zero = rules$drop_grid_zero,
    above_rank = rules$above_rank,
    min_refseq_above = rules$min_refseq_above,
    min_taxa_above = rules$min_taxa_above,
    min_coverage = rules$min_coverage,
    plateau_eps = if (is.null(plateau_path)) NULL else as.numeric(rules$plateau_eps),
    plateau_max_width = rules$plateau_max_width,
    plateau_action = if (is.null(rules$plateau_max_width)) NULL else rules$plateau_action,
    monotonic = rules$monotonic,
    root_thresholds = rules$root_thresholds,
    root_threshold_values = if (is.null(rules$root_values)) NULL else as.list(rules$root_values),
    root_rows_exempt = !any(root_replaced)
  ),
  root = if (is.null(root_log)) NULL else list(
    method = rules$root_thresholds,
    source = paste0("median of the kept rows at the same rank and measure whose ",
                    "superrank is the next coarser rank"),
    rows = root_log
  ),
  composition_source = list(
    training_taxonomy = file_record(taxonomy_path, sequences = n_seq_total),
    method = "optimotu::summarize_by_rank()",
    note = paste0("Counts describe the full training taxonomy. Rows from ",
                  "subsampled training chunks were optimized on fewer sequences.")
  ),
  counts = list(
    rows_in = n_rows,
    rows_out = nrow(out_table),
    rows_dropped = sum(dropped),
    rows_modified = sum(modified),
    root_rows_replaced = sum(root_replaced),
    # first failing rule per dropped row, so these add up to rows_dropped
    dropped = dropped_counts,
    # every rule a row fails, so a row can count under several
    failed = failed_counts,
    modified = modified_counts,
    composition_missing = sum(composition_missing),
    measure_gap_ref_missing = if (is.null(ref_missing)) NULL else sum(ref_missing),
    plateau_keys_missing = if (is.null(plateau_path)) NULL else plateau_missing,
    plateau_width_na = if (is.null(plateau_path)) NULL else plateau_na
  ),
  diagnostics = basename(outputs[[3L]])
)

with_atomic_outputs(outputs, overwrite = isTRUE(opt$overwrite), writer = function(staged) {
  data.table::fwrite(out_table, staged[[1L]], sep = "\t", quote = FALSE, na = "")
  write_json_file(sidecar_out, staged[[2L]])
  data.table::fwrite(diagnostics, staged[[3L]], sep = "\t", quote = FALSE, na = "")
})

note("Wrote ", nrow(out_table), " of ", n_rows, " thresholds (maximum ",
     "clustering distance ", format(max_clustering_distance, digits = 6), ") to ",
     paste(basename(outputs), collapse = ", "), " under ", dirname(prefix))

############################################## Session info

finish(start_time)
