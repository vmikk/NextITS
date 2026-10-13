#!/usr/bin/env Rscript

## Script: optimotu_reference_prep.R
## Purpose: Turn a large reference database into a training set for OptimOTU
##          threshold optimization, and optionally into a SINTAX classification
##          database
##
## Usage:
##   Rscript optimotu_reference_prep.R \
##      --reference-fasta reference.fasta \
##      --output-prefix reference \
##      --ranks kingdom,phylum,class,order,family,genus,species \
##      --emit-sintax
##
## Output:
##   - <prefix>.training.fasta
##   - <prefix>.training_taxonomy.tsv
##   - <prefix>.training.json
##   - <prefix>.sintax.fasta          (only with --emit-sintax)

options(stringsAsFactors = FALSE, warn = 1)

## Check time
start_time <- Sys.time()

script_dir <- dirname(normalizePath(sub(
  "^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]]
)))
source(file.path(script_dir, "optimotu_common.R"))
SCRIPT_START <- start_time

banner("optimotu_reference_prep.R",
       "Prepare an OptimOTU threshold-training set")

cat("Parsing input options and arguments...\n")

suppressPackageStartupMessages(require(optparse))

## Parse arguments
option_list <- list(
  make_option(c("--reference-fasta"),     default=NULL,  type='character', help="Reference sequences (FASTA, optionally gzipped) [required]"),
  make_option(c("--output-prefix"),       default=NULL,  type='character', help="Prefix for the output files [required]"),
  make_option(c("--taxonomy-tsv"),        default=NULL,  type='character', help="Taxonomy table [default: parse the FASTA headers]"),
  make_option(c("--taxonomy-format"),     default="auto", type='character', help="One of: auto, eukaryome, unite, sintax, bold, lineage, tsv [default %default]"),
  make_option(c("--ranks"),               default=NULL,  type='character', help="Rank names, most to least inclusive, comma-separated [default: root,kingdom,phylum,class,order,family,genus,species]"),
  make_option(c("--root-taxon"),          default=NULL,  type='character', help="Value for the synthetic first rank [default: Eukaryota]"),
  make_option(c("--rank-letters"),        default=NULL,  type='character', help="Letter=rank map for prefixed/SINTAX headers, comma-separated [default: d=kingdom,k=kingdom,p=phylum,...]"),
  make_option(c("--restrict"),            default=NULL,  type='character', help="Keep only these lineages, e.g. kingdom=Fungi"),
  make_option(c("--require-rank"),        default=NULL,  type='character', help="Drop records unassigned at this rank [default: the tip rank]"),
  make_option(c("--min-per-taxon"),       default=2L,    type='integer',   help="Drop --require-rank taxa with fewer records [default %default]"),
  make_option(c("--max-per-taxon"),       default=2L,    type='integer',   help="Keep at most N records per --require-rank taxon [default %default]"),
  make_option(c("--max-seqs"),            default=5000L, type='integer',   help="Overall cap, applied by dropping whole taxa [default %default]"),
  make_option(c("--seed"),                default=42L,   type='integer',   help="Subsampling seed [default %default]"),
  make_option(c("--emit-sintax"),         action='store_true', default=FALSE, help="Also write <prefix>.sintax.fasta"),
  make_option(c("--no-qualify-homonyms"), action='store_true', default=FALSE, help="Keep taxon names that occur under several parents"),
  make_option(c("--keep-orphan-ranks"),   action='store_true', default=FALSE, help="Bridge Incertae sedis gaps with a synthetic taxon instead of truncating the lineage there"),
  make_option(c("--overwrite"),           action='store_true', default=FALSE, help="Replace existing output files")
)

opt <- parse_args(OptionParser(
  option_list = option_list,
  usage = "%prog --reference-fasta FILE --output-prefix PATH [options]",
  description = paste0(
    "\nPrepare an OptimOTU threshold-training set from a reference database.\n"),
  epilogue = paste0(
    "Subsetting is not optional: threshold training compares all pairs of\n",
    "sequences, so its cost grows quadratically and a full database will not\n",
    "finish.  --max-per-taxon keeps the within-taxon pairs that carry the signal\n",
    "while --max-seqs bounds the total by dropping whole taxa.\n")
))

## Validation of required arguments
require_options(opt, c("reference-fasta", "output-prefix"))

############################################## Load packages

cat("\nLoading R packages...\n")

load_pckg("data.table")
require_package("optimotu", OPTIMOTU_MIN_VERSION)
require_package("Biostrings", "2.70.0")
require_package("jsonlite", "1.8.0")

############################################## Assign and validate parameters

ranks <- resolve_ranks(opt$ranks)
root_rank <- ranks[[1L]]
tip_rank <- ranks[[length(ranks)]]
root_taxon <- if (is.null(opt$`root-taxon`)) DEFAULT_ROOT_TAXON else opt$`root-taxon`
rank_letters <- as_kv_arg(opt$`rank-letters`, "rank-letters", DEFAULT_RANK_LETTERS)
taxonomy_format <- as_choice(opt$`taxonomy-format`, "taxonomy-format",
                             TAXONOMY_FORMATS, "auto")
require_rank <- if (is.null(opt$`require-rank`)) tip_rank else opt$`require-rank`
if (!require_rank %in% ranks) {
  fail("--require-rank must be one of --ranks: ", paste(ranks, collapse = ", "))
}
min_per_taxon <- as_count(opt$`min-per-taxon`, "min-per-taxon", 2L)
max_per_taxon <- as_count(opt$`max-per-taxon`, "max-per-taxon", 2L)
if (max_per_taxon < min_per_taxon) {
  fail("--max-per-taxon cannot be smaller than --min-per-taxon")
}
max_seqs <- as_count(opt$`max-seqs`, "max-seqs", 5000L)
seed <- as_count(opt$seed, "seed", 42L)

restrict <- as_kv_arg(opt$restrict, "restrict")
if (length(restrict) && !all(names(restrict) %in% ranks)) {
  fail("--restrict names a rank that is not in --ranks: ",
       paste(setdiff(names(restrict), ranks), collapse = ", "))
}

reference_path <- normalizePath(opt$`reference-fasta`, mustWork = TRUE)
taxonomy_path <- if (is.null(opt$`taxonomy-tsv`)) NULL else {
  normalizePath(opt$`taxonomy-tsv`, mustWork = TRUE)
}
prefix <- path.expand(opt$`output-prefix`)

## Log assigned variables
cat("\nParameters specified:\n")
report_params(list(
  "Reference FASTA"   = reference_path,
  "Taxonomy table"    = taxonomy_path,
  "Taxonomy format"   = taxonomy_format,
  "Output prefix"     = prefix,
  "Ranks"             = ranks,
  "Root taxon"        = root_taxon,
  "Rank letters"      = rank_letters,
  "Restrict"          = restrict,
  "Require rank"      = require_rank,
  "Min per taxon"     = min_per_taxon,
  "Max per taxon"     = max_per_taxon,
  "Max sequences"     = max_seqs,
  "Seed"              = seed,
  "Emit SINTAX DB"    = isTRUE(opt$`emit-sintax`),
  "Qualify homonyms"  = !isTRUE(opt$`no-qualify-homonyms`),
  "Keep orphan ranks" = isTRUE(opt$`keep-orphan-ranks`),
  "Overwrite"         = isTRUE(opt$overwrite)
))
cat("\n")

############################################## Main workflow

cat("\n--- Main Workflow ---\n\n")

## --- taxonomy ---------------------------------------------------------------

step("Reading reference headers from ", reference_path, " ...")
headers <- names(Biostrings::fasta.seqlengths(reference_path, seqtype = "DNA"))
if (!length(headers)) fail("Reference FASTA contains no sequences")
note("  ", length(headers), " records")

step("Parsing the reference taxonomy (format: ", taxonomy_format, ") ...")
taxonomy <- if (is.null(taxonomy_path)) {
  read_reference_taxonomy(reference_path, ranks, taxonomy_format, rank_letters,
                          headers = headers)
} else {
  read_reference_taxonomy(taxonomy_path, ranks, taxonomy_format, rank_letters)
}
note("  ", nrow(taxonomy), " records parsed")

step("Cleaning and validating the taxonomy ...")
taxonomy <- add_root_rank(taxonomy, root_rank, root_taxon)
taxonomy <- clean_reference_taxonomy(
  taxonomy, ranks, "Reference taxonomy",
  qualify = !isTRUE(opt$`no-qualify-homonyms`),
  keep_orphan_ranks = isTRUE(opt$`keep-orphan-ranks`)
)
note("  records assigned per rank: ",
     paste(sprintf("%s=%d", ranks,
                   colSums(!is.na(taxonomy[, ranks, drop = FALSE]))),
           collapse = " "))

if (!is.null(taxonomy_path)) {
  step("Matching the taxonomy table against the FASTA IDs ...")
  fasta_ids <- parse_seq_ids(headers, "Reference FASTA")
  missing_tax <- setdiff(fasta_ids, taxonomy$seq_id)
  if (length(missing_tax)) {
    fail("Taxonomy table is missing ", length(missing_tax),
         " FASTA IDs (e.g. ", paste(utils::head(missing_tax, 5L), collapse = ", "), ")")
  }
  taxonomy <- taxonomy[match(fasta_ids, taxonomy$seq_id), , drop = FALSE]
}

n_parsed <- nrow(taxonomy)

## --- subsetting -------------------------------------------------------------

step("Subsetting to a trainable size ...")

for (rank in names(restrict)) {
  wanted <- as_list_arg(restrict[[rank]])
  keep <- !is.na(taxonomy[[rank]]) & taxonomy[[rank]] %in% wanted
  taxonomy <- taxonomy[keep, , drop = FALSE]
  note("  restricted to ", rank, " in {", paste(wanted, collapse = ", "), "}: ",
       nrow(taxonomy), " records")
  if (!nrow(taxonomy)) fail("--restrict left no records")
}

taxonomy <- taxonomy[!is.na(taxonomy[[require_rank]]), , drop = FALSE]
note("  records assigned at ", require_rank, ": ", nrow(taxonomy))
if (!nrow(taxonomy)) {
  fail("No records are assigned at --require-rank ", require_rank)
}

set.seed(seed)

## Keep whole taxa together: threshold optimization scores within- versus between-taxon pairs,
## so a taxon represented by a single sequence contributes no within-taxon pairs and only inflates the all-pairs cost
##
## Deliberately still base R: `split()` is already linear here, and any other
## grouping would visit the taxa in a different order and so draw a different
## subsample from the same --seed
by_taxon <- split(seq_len(nrow(taxonomy)), taxonomy[[require_rank]])
by_taxon <- by_taxon[lengths(by_taxon) >= min_per_taxon]
if (!length(by_taxon)) {
  fail("No ", require_rank, " taxon has at least --min-per-taxon (",
       min_per_taxon, ") records")
}
note("  distinct ", require_rank, " taxa with >= ", min_per_taxon, " records: ",
     length(by_taxon))

by_taxon <- lapply(by_taxon, function(idx) {
  if (length(idx) <= max_per_taxon) idx else sort(sample(idx, max_per_taxon))
})

if (sum(lengths(by_taxon)) > max_seqs) {
  # drop whole taxa, largest-first ordering unchanged, until the cap is met
  order_taxa <- sample(seq_along(by_taxon))
  cumulative <- cumsum(lengths(by_taxon)[order_taxa])
  keep_taxa <- order_taxa[cumulative <= max_seqs]
  if (!length(keep_taxa)) fail("--max-seqs is smaller than a single taxon")
  note("  capping at --max-seqs ", max_seqs, ": keeping ", length(keep_taxa),
       " of ", length(by_taxon), " ", require_rank, " taxa")
  by_taxon <- by_taxon[sort(keep_taxa)]
}

selected <- sort(unlist(by_taxon, use.names = FALSE))
taxonomy <- taxonomy[selected, , drop = FALSE]
rownames(taxonomy) <- NULL
note("  training set: ", nrow(taxonomy), " sequences, ",
     data.table::uniqueN(taxonomy[[require_rank]]), " ", require_rank, " taxa")

## --- sequences --------------------------------------------------------------

step("Reading reference sequences ...")
reference <- read_fasta(reference_path, "Reference FASTA")
missing_seq <- setdiff(taxonomy$seq_id, names(reference))
if (length(missing_seq)) {
  fail("Selected IDs absent from the FASTA (e.g. ",
       paste(utils::head(missing_seq, 5L), collapse = ", "), ")")
}
training <- reference[taxonomy$seq_id]
rm(reference)
invisible(gc(verbose = FALSE))
note("  kept ", length(training), " of ", length(headers), " sequences")

## --- outputs ----------------------------------------------------------------

suffixes <- c(".training.fasta", ".training_taxonomy.tsv", ".training.json")
if (isTRUE(opt$`emit-sintax`)) suffixes <- c(suffixes, ".sintax.fasta")
outputs <- paste0(prefix, suffixes)

step("Writing ", length(outputs), " output file(s) ...")

with_atomic_outputs(outputs, overwrite = isTRUE(opt$overwrite), writer = function(staged) {
  Biostrings::writeXStringSet(training, staged[[1L]], format = "fasta", width = 80L)
  data.table::fwrite(taxonomy, staged[[2L]], sep = "\t", quote = FALSE, na = "")
  write_json_file(
    list(
      schema_version = 1L,
      created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
      reference = file_record(reference_path, records = length(headers)),
      taxonomy = c(
        if (is.null(taxonomy_path)) list(file = NULL) else file_record(taxonomy_path),
        list(format = taxonomy_format, ranks = ranks, root_rank = root_rank,
             root_taxon = root_taxon, parsed_records = n_parsed)
      ),
      subsetting = list(
        restrict = if (length(restrict)) as.list(restrict) else NULL,
        require_rank = require_rank,
        min_per_taxon = min_per_taxon,
        max_per_taxon = max_per_taxon,
        max_seqs = max_seqs,
        seed = seed
      ),
      counts = list(
        sequences = nrow(taxonomy),
        taxa = data.table::uniqueN(taxonomy[[require_rank]]),
        assigned = as.list(stats::setNames(
          as.integer(colSums(!is.na(taxonomy[, ranks, drop = FALSE]))), ranks
        ))
      )
    ),
    staged[[3L]]
  )
  if (isTRUE(opt$`emit-sintax`)) {
    sintax_db <- training
    names(sintax_db) <- sintax_headers(taxonomy, ranks, rank_letters,
                                       skip_ranks = root_rank)
    # A record with no named rank from kingdom to genus is no use to a
    # classifier, and sintaxer refuses to index a database that has one
    annotated <- grepl(";tax=([^;]*,)?[dkpcofg]:[^,;]*[^,;[:space:]]",
                       names(sintax_db))
    if (!all(annotated)) {
      note("Leaving ", sum(!annotated), " record(s) with no rank between ",
           "kingdom and genus out of the SINTAX database")
      sintax_db <- sintax_db[annotated]
    }
    Biostrings::writeXStringSet(sintax_db, staged[[4L]], format = "fasta",
                                width = 80L)
  }
})

note("Wrote ", paste(basename(outputs), collapse = ", "), " under ", dirname(prefix))

############################################## Session info

finish(start_time)
