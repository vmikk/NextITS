#!/usr/bin/env Rscript

## Script: optimotu_reference_set.R
## Purpose: Bundle trained thresholds with a SINTAX classifier database whose
##          taxon names are the ones the thresholds were trained on
##
##          Training cleans placeholder names and qualifies homonyms with their
##          parent ("Lyophyllaceae Clitocybe", see qualify_homonyms()). A
##          classifier DB built from the raw reference keeps the raw names, so
##          threshold rows keyed by a cleaned or qualified supertaxon are never
##          selected and those sequences fall back to an ancestor's threshold.
##          Each classifier record is rewritten as:
##            - decoy (top-rank name starting with "_", e.g. _pseudogene): that
##              name only; its lower ranks reuse real names
##            - in the training taxonomy: the training lineage
##            - absent from training: its own lineage, cleaned as in training and
##              qualified wherever training qualified the same (parent, name)
##            - nothing left below the root after either: the raw top-rank name
##              only (e.g. Eukaryota_reg_incertae_sedis), which classification
##              treats as unknown; without a top-rank name the record is dropped
##
## Usage:
##   Rscript optimotu_reference_set.R \
##      --fasta 13.Eukaryome_2.2_ITS_full.sintax.fasta.gz \
##      --training-taxonomy reference.training_taxonomy.tsv \
##      --thresholds reference.thresholds.tsv \
##      --label eukaryome_its_v2.2 \
##      --output-dir reference_set
##
## Output (in --output-dir):
##   - thresholds.tsv, thresholds.json  copies of --thresholds and its sidecar
##   - classifier.sintax.fasta.gz       the classifier database
##   - unreachable_thresholds.tsv       threshold rows the classifier cannot select
##   - manifest.json

options(stringsAsFactors = FALSE, warn = 1)

## Check time
start_time <- Sys.time()

script_dir <- dirname(normalizePath(sub(
  "^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]]
)))
source(file.path(script_dir, "optimotu_common.R"))
SCRIPT_START <- start_time

banner("optimotu_reference_set.R",
       "Build a reference set (thresholds + matching classifier database)")

cat("Parsing input options and arguments...\n")

suppressPackageStartupMessages(require(optparse))

## Parse arguments
option_list <- list(
  make_option(c("--fasta"),             default=NULL, type='character', help="Classifier source FASTA, may be gzipped [required]"),
  make_option(c("--training-taxonomy"), default=NULL, type='character', help="reference.training_taxonomy.tsv from the training run [required]"),
  make_option(c("--thresholds"),        default=NULL, type='character', help="Thresholds TSV; its .json sidecar must sit next to it [required]"),
  make_option(c("--output-dir"),        default=NULL, type='character', help="Directory for the reference set [required]"),
  make_option(c("--taxonomy-format"),   default="auto", type='character', help="Header format of --fasta: auto, eukaryome, unite, sintax, bold, lineage [default %default]"),
  make_option(c("--rank-letters"),      default=NULL, type='character', help="Letter=rank map for prefixed/SINTAX headers [default: the sidecar's, else d=kingdom,k=kingdom,p=phylum,...]"),
  make_option(c("--keep-orphan-ranks"), action='store_true', default=FALSE, help="Bridge Incertae sedis gaps in records absent from training, as the training run did with the same flag"),
  make_option(c("--label"),             default=NULL, type='character', help="Label recorded in the manifest"),
  make_option(c("--id"),                default=NULL, type='character', help="Manifest id [default: --label]"),
  make_option(c("--title"),             default=NULL, type='character', help="Manifest title"),
  make_option(c("--database"),          default=NULL, type='character', help="Manifest database name, e.g. EUKARYOME"),
  make_option(c("--database-version"),  default=NULL, type='character', help="Manifest database version"),
  make_option(c("--region"),            default=NULL, type='character', help="Manifest region, e.g. ITS_full"),
  make_option(c("--citation"),          default=NULL, type='character', help="Manifest citation"),
  make_option(c("--index-name"),        default=NULL, type='character', help="Name of the sintaxer index built next to the database, recorded in the manifest"),
  make_option(c("--overwrite"),         action='store_true', default=FALSE, help="Replace existing output files")
)

opt <- parse_args(OptionParser(
  option_list = option_list,
  usage = paste0("%prog --fasta FILE --training-taxonomy FILE --thresholds FILE ",
                 "--output-dir DIR [options]"),
  description = paste0(
    "\nBundle trained thresholds with a classifier database whose taxon names\n",
    "match the threshold names.\n"),
  epilogue = paste0(
    "Ranks, root rank and rank letters come from the thresholds' JSON sidecar.\n",
    "The classifier FASTA may be a superset of the training reference (e.g. with\n",
    "decoys); records absent from training are cleaned and qualified the same way.\n")
))

## Validation of required arguments
require_options(opt, c("fasta", "training-taxonomy", "thresholds", "output-dir"))

############################################## Load packages

cat("\nLoading R packages...\n")

load_pckg("data.table")
require_package("optimotu", OPTIMOTU_MIN_VERSION)
require_package("Biostrings", "2.70.0")
require_package("jsonlite", "1.8.0")

############################################## Assign and validate parameters

fasta_path <- normalizePath(opt$fasta, mustWork = TRUE)
taxonomy_path <- normalizePath(opt$`training-taxonomy`, mustWork = TRUE)
thresholds_path <- normalizePath(opt$thresholds, mustWork = TRUE)
sidecar_path <- sub("\\.tsv$", ".json", thresholds_path)
if (identical(sidecar_path, thresholds_path) || !file.exists(sidecar_path)) {
  fail("No JSON sidecar next to the thresholds table: ", sidecar_path)
}
taxonomy_format <- as_choice(opt$`taxonomy-format`, "taxonomy-format",
                             setdiff(TAXONOMY_FORMATS, "tsv"), "auto")
output_dir <- path.expand(opt$`output-dir`)

sidecar <- jsonlite::read_json(sidecar_path, simplifyVector = TRUE)
ranks <- as.character(sidecar$ranks)
if (length(ranks) < 2L) fail("The thresholds sidecar lists no ranks: ", sidecar_path)
root_rank <- sidecar$root_rank %||% ranks[[1L]]
if (!identical(root_rank, ranks[[1L]])) {
  fail("The sidecar's root_rank (", root_rank, ") is not its first rank")
}
root_taxon <- sidecar$root_taxon %||% DEFAULT_ROOT_TAXON
below_root <- ranks[-1L]
top_rank <- below_root[[1L]]
rank_letters <- if (!is.null(opt$`rank-letters`)) {
  as_kv_arg(opt$`rank-letters`, "rank-letters")
} else if (!is.null(sidecar$rank_letters)) {
  unlist(sidecar$rank_letters)
} else {
  DEFAULT_RANK_LETTERS
}

label <- opt$label
index_name <- opt$`index-name`
if (!is.null(index_name) && basename(index_name) != index_name) {
  fail("--index-name must be a file name, not a path")
}

## Log assigned variables
cat("\nParameters specified:\n")
report_params(list(
  "Classifier FASTA"  = fasta_path,
  "Taxonomy format"   = taxonomy_format,
  "Training taxonomy" = taxonomy_path,
  "Thresholds"        = thresholds_path,
  "Ranks"             = ranks,
  "Root taxon"        = root_taxon,
  "Rank letters"      = rank_letters,
  "Keep orphan ranks" = isTRUE(opt$`keep-orphan-ranks`),
  "Label"             = label,
  "Index name"        = index_name,
  "Output directory"  = output_dir,
  "Overwrite"         = isTRUE(opt$overwrite)
))
cat("\n")

############################################## Main workflow

cat("\n--- Main Workflow ---\n\n")

## --- thresholds and training taxonomy ---------------------------------------

step("Reading the thresholds table ...")
thresholds <- fread(thresholds_path, sep = "\t", colClasses = "character",
                    na.strings = "", quote = "")
missing_columns <- setdiff(c("rank", "superrank", "supertaxon"), names(thresholds))
if (length(missing_columns)) {
  fail("The thresholds table lacks columns: ", paste(missing_columns, collapse = ", "))
}
used_ranks <- unique(c(thresholds$rank, thresholds$superrank))
if (length(setdiff(used_ranks, ranks))) {
  fail("The thresholds table uses ranks missing from its sidecar: ",
       paste(setdiff(used_ranks, ranks), collapse = ", "))
}
note("  ", nrow(thresholds), " threshold rows")

step("Reading the training taxonomy ...")
taxonomy <- fread(taxonomy_path, sep = "\t", colClasses = "character",
                  na.strings = "", quote = "")
if (!identical(names(taxonomy), c("seq_id", ranks))) {
  fail("Training taxonomy columns (", paste(names(taxonomy), collapse = ", "),
       ") do not match seq_id + the sidecar ranks")
}
note("  ", nrow(taxonomy), " records")

# (rank, qualified parent, bare child) for every name training qualified; a
# qualified name is "<parent> <child>"
qualify_map <- rbindlist(lapply(seq_along(ranks)[-1L], function(i) {
  child <- taxonomy[[ranks[[i]]]]
  parent <- taxonomy[[ranks[[i - 1L]]]]
  hit <- !is.na(child) & !is.na(parent) & startsWith(child, paste0(parent, " "))
  unique(data.table(rank = ranks[[i]], parent = parent[hit],
                    bare = substring(child[hit], nchar(parent[hit]) + 2L)))
}))
note("  qualified names per rank: ",
     if (nrow(qualify_map)) {
       paste(qualify_map[, .N, by = "rank"][, sprintf("%s=%d", rank, N)], collapse = " ")
     } else "none")

## --- classifier FASTA -------------------------------------------------------

step("Reading the classifier FASTA ...")
seqs <- tryCatch(
  Biostrings::readDNAStringSet(fasta_path, format = "fasta"),
  error = function(e) fail("Cannot parse the classifier FASTA: ", conditionMessage(e))
)
if (!length(seqs)) fail("The classifier FASTA contains no sequences")
headers <- names(seqs)
ids <- parse_seq_ids(headers, "Classifier FASTA")
note("  ", length(seqs), " records")

if (identical(taxonomy_format, "auto")) {
  taxonomy_format <- detect_taxonomy_format(headers, ranks)
  note("  detected header format: ", taxonomy_format)
}
raw <- as.data.table(read_reference_taxonomy(fasta_path, ranks, taxonomy_format,
                                             rank_letters, headers = headers))
# ',' and ';' are structural in SINTAX headers
raw_top <- gsub("[,;]", "_", raw[[top_rank]])

## --- classify the records ---------------------------------------------------

step("Rewriting the lineages ...")
n <- length(ids)
out <- data.table(seq_id = ids)
for (r in ranks) set(out, j = r, value = NA_character_)
kind <- rep(NA_character_, n)

decoy <- !is.na(raw_top) & startsWith(raw_top, "_")
kind[decoy] <- "decoy"

tax_row <- match(ids, taxonomy$seq_id)
in_training <- !is.na(tax_row) & !decoy
rows <- tax_row[in_training]
has_lineage <- rowSums(!is.na(as.matrix(taxonomy[rows, below_root, with = FALSE]))) > 0L
idx <- which(in_training)[has_lineage]
if (length(idx)) {
  for (r in below_root) set(out, i = idx, j = r, value = taxonomy[[r]][rows[has_lineage]])
}
kind[idx] <- "training"
kind[which(in_training)[!has_lineage]] <- "placeholder_only"
note("  decoys: ", sum(decoy), "; in the training taxonomy: ", sum(in_training),
     " (", sum(!has_lineage), " with no lineage below ", root_rank, ")")
missing_from_db <- sum(!taxonomy$seq_id %chin% ids)
if (missing_from_db) {
  note("  NOTE: ", missing_from_db, " training records are not in the classifier FASTA")
}

absent <- which(is.na(kind))
qualified <- stats::setNames(integer(length(below_root)), below_root)
if (length(absent)) {
  note("  records absent from training: ", length(absent), "; cleaning them as training did")
  ab <- add_root_rank(as.data.frame(raw[absent]), root_rank, root_taxon)
  ab <- clean_reference_taxonomy(ab, ranks, "Records absent from training",
                                 require_single_root = FALSE, qualify = FALSE,
                                 keep_orphan_ranks = isTRUE(opt$`keep-orphan-ranks`))
  setDT(ab)
  # top-down, so that a parent is already qualified when its children are looked up
  for (i in seq_along(ranks)[-1L]) {
    child <- ranks[[i]]
    parent <- ranks[[i - 1L]]
    m <- qualify_map[rank == child]
    if (!nrow(m)) next
    hit <- which(!is.na(ab[[child]]) &
                 paste(ab[[parent]], ab[[child]], sep = "\t") %chin%
                 paste(m$parent, m$bare, sep = "\t"))
    if (length(hit)) {
      set(ab, i = hit, j = child, value = paste(ab[[parent]][hit], ab[[child]][hit]))
      qualified[[child]] <- length(hit)
    }
  }
  if (any(qualified > 0L)) {
    note("  qualified like training (records): ",
         paste(sprintf("%s=%d", names(qualified)[qualified > 0L],
                       qualified[qualified > 0L]), collapse = " "))
  }
  has_lineage <- rowSums(!is.na(as.matrix(ab[, below_root, with = FALSE]))) > 0L
  idx <- absent[has_lineage]
  if (length(idx)) {
    for (r in below_root) set(out, i = idx, j = r, value = ab[[r]][has_lineage])
  }
  kind[idx] <- "absent_cleaned"
  kind[absent[!has_lineage]] <- "placeholder_only"
  rm(ab)
}

# Records left without a lineage keep their raw top-rank name, if they have one
top_only <- which(kind %chin% c("decoy", "placeholder_only"))
if (length(top_only)) set(out, i = top_only, j = top_rank, value = raw_top[top_only])
dropped <- top_only[is.na(raw_top[top_only])]
kind[dropped] <- "dropped"

# Names still under more than one parent: training qualifies homonyms only among
# its own records, so a name can get a second lineage from the records absent
# from training (e.g. genus Chaetopsis, Reticulascaceae in training and one
# Ulidiidae record outside it). OptimOTU cannot cluster such calls, so qualify
# them with their parent, top-down; the lineage training knows keeps the bare
# name, and with it its thresholds
db_qualified <- stats::setNames(integer(length(below_root)), below_root)
db_examples <- character()
for (i in seq_along(below_root)[-1L]) {
  child <- below_root[[i]]
  parent <- below_root[[i - 1L]]
  known <- !is.na(out[[child]]) & !is.na(out[[parent]])
  if (!any(known)) next
  pairs <- unique(data.table(ch = out[[child]][known], pa = out[[parent]][known]))
  homonyms <- pairs[, list(n_parents = .N), by = "ch"][n_parents > 1L, ch]
  if (!length(homonyms)) next
  # (parent, name) of homonyms that training has under a single parent
  trained <- unique(data.table(ch = taxonomy[[child]], pa = taxonomy[[parent]]))
  trained <- trained[ch %chin% homonyms & !is.na(pa)]
  trained <- trained[, list(pa = if (.N == 1L) pa else NA_character_), by = "ch"][!is.na(pa)]
  hit <- which(known & out[[child]] %chin% homonyms)
  bare <- paste(out[[parent]][hit], out[[child]][hit], sep = "\t") %chin%
    paste(trained$pa, trained$ch, sep = "\t")
  hit <- hit[!bare]
  set(out, i = hit, j = child, value = paste(out[[parent]][hit], out[[child]][hit]))
  db_qualified[[child]] <- length(homonyms)
  db_examples <- c(db_examples, paste0(child, ":", utils::head(homonyms, 3L)))
}
if (any(db_qualified > 0L)) {
  note("  qualified names under more than one parent (names): ",
       paste(sprintf("%s=%d", names(db_qualified)[db_qualified > 0L],
                     db_qualified[db_qualified > 0L]), collapse = " "),
       ", e.g. ", paste(utils::head(db_examples, 5L), collapse = ", "))
}

counts <- vapply(c("training", "decoy", "placeholder_only", "absent_cleaned", "dropped"),
                 function(k) sum(kind == k), integer(1))
note("  records per class: ", paste(sprintf("%s=%d", names(counts), counts), collapse = " "))
if (counts[["dropped"]]) {
  note("  dropped ", counts[["dropped"]], " record(s) with no ", top_rank, " name")
}

keep <- kind != "dropped"
out <- out[keep]
seqs <- seqs[keep]
names(seqs) <- sintax_headers(out, ranks, rank_letters, skip_ranks = root_rank)

## --- checks -----------------------------------------------------------------

step("Checking the classifier database against the thresholds ...")
checks <- list()

# threshold rows whose supertaxon the classifier can never return; the root
# taxon is added by the pipeline, not by the classifier
names_by_rank <- function(tab) rbindlist(lapply(below_root, function(r) {
  taxa <- unique(stats::na.omit(tab[[r]]))
  data.table(superrank = rep(r, length(taxa)), supertaxon = taxa)
}))
db_names <- names_by_rank(out)
below <- thresholds[superrank != root_rank]
unreachable <- below[!db_names, on = c("superrank", "supertaxon")]
supertaxa <- unique(below[, c("superrank", "supertaxon")])
supertaxa[, `:=`(in_raw = FALSE, in_db = FALSE)]
supertaxa[names_by_rank(raw), in_raw := TRUE, on = c("superrank", "supertaxon")]
supertaxa[db_names, in_db := TRUE, on = c("superrank", "supertaxon")]
reach <- supertaxa[, list(supertaxa = .N, reachable_raw = sum(in_raw),
                          reachable = sum(in_db)), by = "superrank"]
print(reach)
checks$unreachable_thresholds <- list(
  file = "unreachable_thresholds.tsv",
  rows = nrow(unreachable),
  supertaxa = nrow(unique(unreachable[, c("superrank", "supertaxon")])),
  rows_checked = nrow(below)
)
checks$reachable_supertaxa <- lapply(split(reach, by = "superrank", keep.by = FALSE),
                                     as.list)
if (nrow(unreachable)) {
  warning(nrow(unreachable), " threshold row(s) (",
          checks$unreachable_thresholds$supertaxa, " supertaxa) have a supertaxon ",
          "the classifier database never returns, e.g. ",
          paste(utils::head(unique(unreachable$supertaxon), 5L), collapse = ", "),
          "; see unreachable_thresholds.tsv", call. = FALSE)
}

# bare names of homonyms training qualified: such a name in the database is a
# taxon the thresholds know only under a qualified name
bare <- rbindlist(lapply(below_root, function(r) {
  m <- qualify_map[rank == r]
  if (!nrow(m)) return(NULL)
  hit <- out[[r]] %chin% m$bare
  if (!any(hit)) return(NULL)
  data.table(rank = r, taxon = out[[r]][hit])
}))
bare_names <- if (nrow(bare)) unique(bare) else bare
checks$bare_homonyms <- list(
  names = nrow(bare_names),
  records = nrow(bare),
  examples = I(utils::head(if (nrow(bare_names)) paste0(bare_names$rank, ":", bare_names$taxon)
                           else character(), 10L))
)
if (nrow(bare)) {
  warning(nrow(bare_names), " name(s) in ", nrow(bare), " record(s) equal the bare ",
          "part of a name training qualified with its parent, e.g. ",
          paste(checks$bare_homonyms$examples[seq_len(min(5L, nrow(bare_names)))],
                collapse = ", "),
          "; their parent differs from every parent training saw", call. = FALSE)
}

# names under more than one parent; OptimOTU's threshold lookup assumes a tree
multi <- vapply(seq_along(below_root)[-1L], function(i) {
  pairs <- unique(data.table(ch = out[[below_root[[i]]]], pa = out[[below_root[[i - 1L]]]]))
  pairs <- pairs[!is.na(ch)]
  sum(duplicated(pairs$ch))
}, integer(1))
names(multi) <- below_root[-1L]
checks$multi_parent_names <- as.list(multi)
checks$qualified_db_homonyms <- list(
  names = as.list(db_qualified),
  examples = I(utils::head(db_examples, 10L))
)
# with both names known, every homonym was qualified above
multi_known <- vapply(seq_along(below_root)[-1L], function(i) {
  pairs <- unique(data.table(ch = out[[below_root[[i]]]], pa = out[[below_root[[i - 1L]]]]))
  pairs <- pairs[!is.na(ch) & !is.na(pa)]
  sum(duplicated(pairs$ch))
}, integer(1))
if (any(multi_known > 0L)) {
  fail("Names under more than one parent remain after qualifying them: ",
       paste(sprintf("%s=%d", below_root[-1L][multi_known > 0L], multi_known[multi_known > 0L]),
             collapse = " "))
}
if (any(multi > 0L)) {
  warning("Names under more than one parent: ",
          paste(sprintf("%s=%d", names(multi)[multi > 0L], multi[multi > 0L]),
                collapse = " "), call. = FALSE)
}

## --- outputs ----------------------------------------------------------------

out_names <- c(thresholds = "thresholds.tsv", sidecar = "thresholds.json",
               classifier_db = "classifier.sintax.fasta.gz",
               unreachable = "unreachable_thresholds.tsv",
               manifest = "manifest.json")
outputs <- file.path(output_dir, out_names)

step("Writing ", length(outputs), " file(s) to ", output_dir, " ...")

with_atomic_outputs(outputs, overwrite = isTRUE(opt$overwrite), writer = function(staged) {
  names(staged) <- names(out_names)
  if (!file.copy(thresholds_path, staged[["thresholds"]]) ||
      !file.copy(sidecar_path, staged[["sidecar"]])) {
    fail("Could not copy the thresholds into the reference set")
  }
  Biostrings::writeXStringSet(seqs, staged[["classifier_db"]], compress = TRUE,
                              format = "fasta", width = 80L)
  fwrite(unreachable, staged[["unreachable"]], sep = "\t", quote = FALSE, na = "")

  files <- lapply(c("thresholds", "sidecar", "classifier_db"), function(k) {
    list(name = out_names[[k]], role = k, md5 = file_md5(staged[[k]]))
  })
  if (!is.null(index_name)) {
    files <- c(files, list(list(name = index_name, role = "classifier_index")))
  }
  write_json_file(
    list(
      schema_version = 1L,
      id = opt$id %||% label,
      title = opt$title,
      database = opt$database,
      database_version = opt$`database-version`,
      region = opt$region,
      classifier = "sintaxer",
      classifier_db = out_names[["classifier_db"]],
      classifier_index = index_name,
      thresholds = out_names[["thresholds"]],
      citation = opt$citation,
      label = label,
      created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
      files = files,
      ranks = I(ranks),
      root_rank = root_rank,
      root_taxon = root_taxon,
      rank_letters = as.list(rank_letters),
      distance = sidecar$distance,
      sources = list(
        classifier_fasta = file_record(fasta_path, format = taxonomy_format, records = n),
        training_taxonomy = file_record(taxonomy_path, records = nrow(taxonomy)),
        thresholds = file_record(thresholds_path, rows = nrow(thresholds))
      ),
      counts = c(
        list(records = n, written = length(seqs)),
        as.list(counts),
        list(training_not_in_fasta = missing_from_db,
             qualified_absent = as.list(qualified))
      ),
      checks = checks
    ),
    staged[["manifest"]]
  )
})

note("Wrote ", length(seqs), " classifier records and ", nrow(thresholds),
     " thresholds to ", output_dir)

############################################## Session info

finish(start_time)
