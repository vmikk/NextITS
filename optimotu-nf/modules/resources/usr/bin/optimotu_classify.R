#!/usr/bin/env Rscript

## Script: optimotu_classify.R
## Purpose: Classify query sequences and emit the long taxonomic-probability
##          table that optimotu_cluster.R consumes:
##
##            seq_id <tab> rank <tab> taxon <tab> prob
##
##          This is the pipeline's classifier boundary.  Any tool that can
##          produce this table can be plugged in without touching the clustering
##          step
##
## Usage:
##   Rscript optimotu_classify.R \
##      --input-fasta query.fasta \
##      --output query.taxprob.tsv \
##      --classifier sintaxer \
##      --db reference.stx \
##      --threads 4
##
## Classifiers:
##   - sintaxer   SINTAX-style bootstrap classification with sintaxer (default);
##                --db is a sintaxer index (.stx) or a SINTAX-format FASTA,
##                which is then indexed on the fly
##   - sintax     the same algorithm in VSEARCH (`vsearch --sintax`); --db is a
##                SINTAX-format FASTA
##   - vsearch    best-hit transfer with `vsearch --usearch_global`
##   - none       root rank only, i.e. fully de novo clustering
##
## Output:
##   - <output>            long taxonomic-probability TSV
##   - <output>.json       sidecar recording the database and tool versions

options(stringsAsFactors = FALSE, warn = 1)

## Check time
start_time <- Sys.time()

script_dir <- dirname(normalizePath(sub(
  "^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]]
)))
source(file.path(script_dir, "optimotu_common.R"))
SCRIPT_START <- start_time

CLASSIFIERS <- c("sintaxer", "sintax", "vsearch", "none")

banner("optimotu_classify.R",
       "Classify sequences into a long taxonomic-probability table")

cat("Parsing input options and arguments...\n")

suppressPackageStartupMessages(require(optparse))

## Parse arguments
option_list <- list(
  make_option(c("--input-fasta"),     default=NULL,  type='character', help="Query sequences [required]"),
  make_option(c("--output"),          default=NULL,  type='character', help="Destination TSV (seq_id, rank, taxon, prob); a .json sidecar is written beside it [required]"),
  make_option(c("--classifier"),      default="sintaxer", type='character', help="One of: sintaxer, sintax (VSEARCH SINTAX), vsearch, none [default %default]"),
  make_option(c("--db"),              default=NULL,  type='character', help="Reference database (required unless --classifier none): SINTAX-format FASTA, or a sintaxer index (.stx) for --classifier sintaxer"),
  make_option(c("--min-identity"),    default=0.8,   type='double',    help="Minimum identity for --classifier vsearch [default %default]"),
  make_option(c("--threads"),         default=1L,    type='integer',   help="Classifier threads [default %default]"),
  make_option(c("--seed"),            default=42L,   type='integer',   help="PRNG seed; for sintax, 0 uses a random source [default %default]"),
  make_option(c("--vsearch"),         default=NULL,  type='character', help="VSEARCH executable [default: autodetected from $VSEARCH, then PATH]"),
  make_option(c("--sintaxer"),        default=NULL,  type='character', help="sintaxer executable [default: autodetected from $SINTAXER, then PATH]"),
  make_option(c("--ranks"),           default=NULL,  type='character', help="Rank names, most to least inclusive, comma-separated [default: root,kingdom,phylum,class,order,family,genus,species]"),
  make_option(c("--root-taxon"),      default=NULL,  type='character', help="Taxon assigned at the first rank [default: Eukaryota]"),
  make_option(c("--rank-letters"),    default=NULL,  type='character', help="Letter=rank map for SINTAX output, comma-separated"),
  make_option(c("--taxonomy-format"), default="auto", type='character', help="Reference header format for --classifier vsearch [default %default]"),
  make_option(c("--keep-placeholders"), action='store_true', default=FALSE, help="Keep calls to placeholder taxa (e.g. Agaricales_fam_Incertae_sedis, Fungi_sp); by default they are dropped, leaving the sequence unknown at that rank"),
  make_option(c("--overwrite"),       action='store_true', default=FALSE, help="Replace an existing output file")
)

opt <- parse_args(OptionParser(
  option_list = option_list,
  usage = "%prog --input-fasta FILE --output FILE [options]",
  description = "\nClassify sequences into a long taxonomic-probability table.\n",
  epilogue = paste0(
    "'none' assigns only the first rank, which makes the clustering step run\n",
    "fully de novo -- the natural baseline against a plain VSEARCH clustering.\n")
))

## Validation of required arguments
require_options(opt, c("input-fasta", "output"))

############################################## Load packages

cat("\nLoading R packages...\n")

load_pckg("data.table")
require_package("Biostrings", "2.70.0")
require_package("processx", "3.8.0")
require_package("optimotu", OPTIMOTU_MIN_VERSION)

############################################## Assign and validate parameters

classifier <- as_choice(opt$classifier, "classifier", CLASSIFIERS, "sintaxer")
if (!identical(classifier, "none") && is.null(opt$db)) {
  fail("--db is required for --classifier ", classifier)
}
ranks <- resolve_ranks(opt$ranks)
root_rank <- ranks[[1L]]
root_taxon <- if (is.null(opt$`root-taxon`)) DEFAULT_ROOT_TAXON else opt$`root-taxon`
rank_letters <- as_kv_arg(opt$`rank-letters`, "rank-letters", DEFAULT_RANK_LETTERS)
taxonomy_format <- as_choice(opt$`taxonomy-format`, "taxonomy-format",
                             TAXONOMY_FORMATS, "auto")
threads <- as_count(opt$threads, "threads", 1L)
min_identity <- as_probability(opt$`min-identity`, "min-identity", 0.8)
vsearch <- if (classifier %in% c("sintax", "vsearch")) find_vsearch(opt$vsearch) else NULL
sintaxer <- if (identical(classifier, "sintaxer")) find_sintaxer(opt$sintaxer) else NULL
seed <- {
  value <- suppressWarnings(as.integer(opt$seed))
  if (is.na(value) || value < 0L) fail("--seed must be a non-negative integer")
  value
}
# SINTAX bootstraps over randomly sampled k-mers. VSEARCH is only reproducible
# when a seed is combined with a single thread -- it warns about this itself, and
# in practice runs of the same input can otherwise disagree on the taxon name,
# not just the bootstrap value, which would move sequences between OTUs.
# Classification is per sequence, so the way to recover parallelism without
# losing reproducibility is to split the query FASTA across tasks.
# sintaxer needs no such pin: it derives every bootstrap from the seed and the
# sequence itself, so its output is identical at any thread count.
sintax_threads <- threads
if (identical(classifier, "sintax") && seed > 0L && threads > 1L) {
  note("SINTAX is pinned to one thread so that --seed ", seed,
       " makes it reproducible (requested ", threads,
       " threads); pass --seed 0 to trade reproducibility for speed")
  sintax_threads <- 1L
}

input_path <- normalizePath(opt$`input-fasta`, mustWork = TRUE)
db_path <- if (is.null(opt$db)) NULL else normalizePath(opt$db, mustWork = TRUE)
output_path <- path.expand(opt$output)

## Log assigned variables
cat("\nParameters specified:\n")
report_params(list(
  "Input FASTA"       = input_path,
  "Output TSV"        = output_path,
  "Classifier"        = classifier,
  "Database"          = db_path,
  "Min identity"      = if (identical(classifier, "vsearch")) min_identity else NULL,
  "Ranks"             = ranks,
  "Root taxon"        = root_taxon,
  "Rank letters"      = rank_letters,
  "Taxonomy format"   = taxonomy_format,
  "Threads requested" = threads,
  "Threads effective" = if (identical(classifier, "sintax")) sintax_threads else threads,
  "Seed"              = seed,
  "Keep placeholders" = isTRUE(opt$`keep-placeholders`),
  "Overwrite"         = isTRUE(opt$overwrite)
))
report_tool("vsearch", vsearch)
report_tool("sintaxer", sintaxer)
cat("\n")

############################################## Classifier backends

EMPTY_ASSIGNMENTS <- data.frame(seq_id = character(), rank = character(),
                                taxon = character(), prob = numeric())

## Parse SINTAX tabbed output, as written by both `vsearch --sintax --tabbedout`
## and `sintaxer classify`:
##
##   query <tab> d:Fungi(1.00),p:Basidiomycota(0.98),... <tab> strand <tab> ...
##
## Only the first two columns are read -- the cutoff-filtered lineage is ignored,
## because the confidence cutoff is applied later by optimotu()
parse_sintax_tabbed <- function(out_file, ranks, rank_letters) {
  if (!file.exists(out_file)) return(EMPTY_ASSIGNMENTS)
  lines <- readLines(out_file, warn = FALSE)
  if (!length(lines)) return(EMPTY_ASSIGNMENTS)

  ## The whole file is parsed in one pass: tstrsplit() transposes the tab split,
  ## and the lineage terms are exploded into a long vector so the capture regex
  ## runs once over every term rather than once per query sequence
  fields <- data.table::tstrsplit(lines, "\t", fixed = TRUE)
  seq_id <- fields[[1L]]
  lineage <- if (length(fields) >= 2L) fields[[2L]] else NA_character_
  lineage <- rep_len(lineage, length(seq_id))
  lineage[is.na(lineage)] <- ""

  terms <- strsplit(lineage, ",", fixed = TRUE)
  record <- rep.int(seq_along(seq_id), lengths(terms))
  term <- unlist(terms, use.names = FALSE)
  keep <- nzchar(term)
  record <- record[keep]
  term <- term[keep]
  if (!length(term)) return(EMPTY_ASSIGNMENTS)

  # greedy '.+' so the trailing bootstrap is the last parenthesized number
  captures <- regmatches(
    term, regexec("^([A-Za-z]):(.+)\\(([0-9.]+)\\)$", term, perl = TRUE)
  )
  matched <- lengths(captures) == 4L
  if (!any(matched)) return(EMPTY_ASSIGNMENTS)
  # rows are (full match, letter, taxon, bootstrap)
  groups <- matrix(unlist(captures[matched], use.names = FALSE), nrow = 4L)

  out <- data.frame(
    seq_id = seq_id[record[matched]],
    rank = unname(rank_letters[tolower(groups[2L, ])]),
    taxon = trimws(groups[3L, ]),
    prob = as.numeric(groups[4L, ]),
    stringsAsFactors = FALSE
  )
  out[!is.na(out$rank) & out$rank %in% ranks, , drop = FALSE]
}

## Run VSEARCH SINTAX and parse its tabbed output
##
## `--sintax_random` is deliberately omitted: with it, equally good matches are
## broken randomly, which would make OTU membership non-reproducible
classify_sintax <- function(query_file, db, threads, ranks, rank_letters,
                            vsearch, seed) {
  out_file <- tempfile(fileext = ".sintax.tsv")
  on.exit(unlink(out_file, force = TRUE), add = TRUE)
  result <- processx::run(
    vsearch,
    c("--sintax", query_file, "--db", db, "--tabbedout", out_file,
      "--threads", as.character(threads),
      if (seed > 0L) c("--randseed", as.character(seed)),
      "--quiet"),
    error_on_status = FALSE,
    echo = FALSE
  )
  if (result$status != 0L) fail("SINTAX classification failed: ", result$stderr)

  parse_sintax_tabbed(out_file, ranks, rank_letters)
}

## Is `path` a sintaxer index? Indexes start with the 8-byte magic `SINTAXER`
is_sintaxer_index <- function(path) {
  identical(readBin(path, "raw", 8L), charToRaw("SINTAXER"))
}

## Classify with sintaxer and parse its tabbed output
##
## sintaxer classifies against a prebuilt index; a FASTA database is indexed into
## a temporary file first (the Nextflow pipeline does this once per run, in
## SINTAXER_INDEX, and passes the index here instead)
## Its output is deterministic for a given seed at any thread count, and it
## reports ranks down to genus only -- species annotations are not used
classify_sintaxer <- function(query_file, db, threads, ranks, rank_letters,
                              sintaxer, seed) {
  index <- db
  if (!is_sintaxer_index(db)) {
    index <- tempfile(fileext = ".stx")
    on.exit(unlink(index, force = TRUE), add = TRUE)
    note("  indexing ", basename(db), " for sintaxer (a prebuilt .stx index skips this)")
    result <- processx::run(
      sintaxer, c("index", db, "--output", index),
      error_on_status = FALSE, echo = FALSE
    )
    if (result$status != 0L) fail("sintaxer index failed: ", result$stderr)
    note("  ", trimws(result$stderr))
  }

  out_file <- tempfile(fileext = ".sintaxer.tsv")
  on.exit(unlink(out_file, force = TRUE), add = TRUE)
  result <- processx::run(
    sintaxer,
    c("classify", query_file, "--db", index, "--output", out_file,
      "--threads", as.character(threads), "--seed", as.character(seed)),
    error_on_status = FALSE,
    echo = FALSE
  )
  if (result$status != 0L) fail("sintaxer classification failed: ", result$stderr)
  # sintaxer reports a one-line throughput summary on stderr
  note("  ", trimws(result$stderr))

  parse_sintax_tabbed(out_file, ranks, rank_letters)
}

## Best-hit classification with `vsearch --usearch_global`
##
## The hit's whole lineage is transferred to the query, with the alignment
## identity as the confidence at every rank
classify_vsearch <- function(query_file, db, threads, ranks, root_rank,
                             rank_letters, taxonomy_format, min_identity,
                             vsearch) {
  out_file <- tempfile(fileext = ".b6.tsv")
  on.exit(unlink(out_file, force = TRUE), add = TRUE)
  result <- processx::run(
    vsearch,
    c("--usearch_global", query_file, "--db", db,
      "--id", format(min_identity, scientific = FALSE),
      "--blast6out", out_file, "--maxaccepts", "100", "--top_hits_only",
      "--strand", "plus", "--threads", as.character(threads)),
    error_on_status = FALSE,
    echo = FALSE
  )
  if (result$status != 0L) fail("VSEARCH search failed: ", result$stderr)

  if (!file.exists(out_file) || file.size(out_file) == 0L) return(EMPTY_ASSIGNMENTS)
  # only the first three blast6 columns are used, so do not read the other nine
  hits <- data.table::fread(out_file, sep = "\t", header = FALSE, quote = "",
                            select = 1:3, colClasses = "character",
                            na.strings = "NA", showProgress = FALSE,
                            data.table = FALSE)
  if (!nrow(hits)) return(EMPTY_ASSIGNMENTS)
  names(hits)[1:3] <- c("seq_id", "ref_header", "identity")
  hits$identity <- suppressWarnings(as.numeric(hits$identity))
  hits$ref_id <- strip_seq_ids(hits$ref_header)

  # keep one deterministic best hit per query
  hits <- hits[order(hits$seq_id, -hits$identity, hits$ref_id), , drop = FALSE]
  hits <- hits[!duplicated(hits$seq_id), , drop = FALSE]

  ref_tax <- read_reference_taxonomy(db, ranks, taxonomy_format, rank_letters)
  ref_tax <- ref_tax[match(hits$ref_id, ref_tax$seq_id), , drop = FALSE]

  rows <- lapply(setdiff(ranks, root_rank), function(rank) {
    taxon <- ref_tax[[rank]]
    keep <- !is.na(taxon)
    if (!any(keep)) return(NULL)
    data.frame(seq_id = hits$seq_id[keep], rank = rank, taxon = taxon[keep],
               prob = hits$identity[keep] / 100)
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!length(rows)) return(EMPTY_ASSIGNMENTS)
  as.data.frame(data.table::rbindlist(rows, use.names = TRUE))
}

############################################## Main workflow

cat("\n--- Main Workflow ---\n\n")

step("Reading ", input_path, " ...")
seqs <- read_fasta(input_path, "Input FASTA")
ids <- names(seqs)
note(length(ids), " query sequences; classifier: ", classifier)

query_file <- tempfile(fileext = ".fasta")
on.exit(unlink(query_file, force = TRUE), add = TRUE)
Biostrings::writeXStringSet(seqs, query_file, format = "fasta", width = 80L)

step("Classifying with ", classifier, " ...")
raw <- switch(
  classifier,
  sintaxer = classify_sintaxer(query_file, db_path, threads, ranks,
                               rank_letters, sintaxer, seed),
  sintax = classify_sintax(query_file, db_path, sintax_threads, ranks,
                           rank_letters, vsearch, seed),
  vsearch = classify_vsearch(query_file, db_path, threads, ranks, root_rank,
                             rank_letters, taxonomy_format, min_identity,
                             vsearch),
  none = EMPTY_ASSIGNMENTS
)
note("  ", nrow(raw), " raw assignment(s) returned")

step("Normalizing the assignments ...")
raw <- as.data.frame(raw, stringsAsFactors = FALSE)
if (nrow(raw)) {
  raw <- raw[raw$seq_id %in% ids & raw$rank %in% ranks & raw$rank != root_rank &
               !is.na(raw$taxon) & nzchar(raw$taxon), , drop = FALSE]
  raw$prob <- suppressWarnings(as.numeric(raw$prob))
  # taxon names must not carry the delimiters our own tables rely on
  raw$taxon <- gsub("[,;\t]", "_", raw$taxon)
}

# A placeholder (Incertae_sedis, *_sp, unidentified ...) is not a taxon: as a
# call it would group unrelated sequences under one name and pick up that
# name's thresholds.  Only the call at that rank is dropped; deeper calls stay.
placeholders <- 0L
if (nrow(raw) && !isTRUE(opt$`keep-placeholders`)) {
  is_ph <- optimotu::is_placeholder(raw$taxon)
  placeholders <- sum(is_ph)
  note("  dropped ", placeholders, " placeholder call(s)",
       if (placeholders) {
         paste0(", e.g. ", paste(utils::head(unique(raw$taxon[is_ph]), 3L),
                                 collapse = ", "))
       })
  raw <- raw[!is_ph, , drop = FALSE]
}

# Every sequence is assigned at the first rank: optimotu() requires it, and it is
# what makes the root-level threshold the universal fallback.
root <- data.frame(seq_id = ids, rank = root_rank, taxon = root_taxon, prob = 1)
tax_prob <- rbind(root, raw[, c("seq_id", "rank", "taxon", "prob")])

# one deterministic assignment per (sequence, rank): highest probability wins,
# ties broken by taxon name
#
# order() rather than setorderv(): the tie-break is on `taxon`, and data.table
# sorts character columns in the C locale while order() uses the session's
# collation -- switching would silently change which of two tied taxon names wins
tax_prob$sequence_order <- match(tax_prob$seq_id, ids)
tax_prob$rank_order <- match(tax_prob$rank, ranks)
tax_prob <- tax_prob[
  order(tax_prob$sequence_order, tax_prob$rank_order, -tax_prob$prob,
        tax_prob$taxon, na.last = TRUE),
  , drop = FALSE
]
tax_prob <- tax_prob[!duplicated(tax_prob[, c("seq_id", "rank")]), , drop = FALSE]
tax_prob <- tax_prob[, c("seq_id", "rank", "taxon", "prob"), drop = FALSE]
rownames(tax_prob) <- NULL

if (!setequal(unique(tax_prob$seq_id), ids)) {
  fail("Classification did not cover every input sequence")
}

# one pass over the rank column, instead of one pass per rank
per_rank <- stats::setNames(
  tabulate(match(tax_prob$rank, ranks), nbins = length(ranks)), ranks
)
note("Assignments per rank: ",
     paste(sprintf("%s=%d", ranks, per_rank), collapse = " "))

# The sidecar closes the provenance chain: the clustering run record stores this
# table's checksum, and this says which database and tool produced it.
sidecar_path <- if (grepl("\\.tsv$", output_path)) {
  sub("\\.tsv$", ".json", output_path)
} else {
  paste0(output_path, ".json")
}

step("Writing the assignments and the provenance sidecar ...")

with_atomic_outputs(c(output_path, sidecar_path),
                    overwrite = isTRUE(opt$overwrite),
                    writer = function(staged) {
  data.table::fwrite(tax_prob, staged[[1L]], sep = "\t", quote = FALSE, na = "")
  write_json_file(
    list(
      schema_version = 1L,
      created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
      classifier = classifier,
      input = file_record(input_path, sequences = length(ids)),
      database = if (is.null(db_path)) NULL else file_record(db_path),
      parameters = list(ranks = ranks, root_rank = root_rank,
                        root_taxon = root_taxon,
                        rank_letters = as.list(rank_letters),
                        min_identity = if (identical(classifier, "vsearch")) {
                          min_identity
                        } else NULL,
                        keep_placeholders = isTRUE(opt$`keep-placeholders`),
                        seed = seed, requested_threads = threads,
                        effective_threads = if (identical(classifier, "sintax")) {
                          sintax_threads
                        } else threads),
      versions = list(R = as.character(getRversion()),
                      vsearch = if (is.null(vsearch)) NULL else {
                        tool_version(vsearch)
                      },
                      sintaxer = if (is.null(sintaxer)) NULL else {
                        tool_version(sintaxer)
                      }),
      counts = list(assignments = nrow(tax_prob),
                    placeholders_dropped = placeholders,
                    per_rank = as.list(stats::setNames(as.integer(per_rank),
                                                       ranks)))
    ),
    staged[[2L]]
  )
})

note("Wrote ", nrow(tax_prob), " assignments to ", output_path)

############################################## Session info

finish(start_time)
