## Shared helpers for the OptimOTU Nextflow pipeline scripts
##
## Sourced by the optimotu_*.R CLIs in the same directory; contains no
## top-level side effects other than defining functions, constants and
## `SCRIPT_START`

DEFAULT_RANKS <- c("root", "kingdom", "phylum", "class", "order", "family",
                   "genus", "species")

DEFAULT_ROOT_TAXON <- "Eukaryota"

# SINTAX / UNITE style single-letter rank codes.  Eukaryome's SINTAX database
# labels its top rank "d" even though the names it carries ("Fungi",
# "Viridiplantae", ...) are kingdoms, so the mapping is a parameter rather than a
# constant.
DEFAULT_RANK_LETTERS <- c(d = "kingdom", k = "kingdom", p = "phylum",
                          c = "class", o = "order", f = "family",
                          g = "genus", s = "species")

MEASURES <- c("MCC", "RI", "ARI", "FMI", "MI", "AMI", "FM")

DIST_METHODS <- c("wfa2", "ksw2", "edlib", "hybrid", "hamming", "usearch")

# Alignment backends that take gap/mismatch penalties (--dist-penalties)
PENALIZED_DIST_METHODS <- c("wfa2", "ksw2")

# optimotu() and optimize_thresholds() exist from 0.9.6, but multithreaded
# WFA2 clustering only gives correct results from 0.10.0
OPTIMOTU_MIN_VERSION <- "0.10.0"

CLUST_METHODS <- c("tree", "matrix", "index", "slink")

TAXONOMY_FORMATS <- c("auto", "eukaryome", "unite", "sintax", "bold", "lineage",
                      "tsv")

## ---------------------------------------------------------------------------
## logging, errors and run bookkeeping
## ---------------------------------------------------------------------------
##
## Everything informational goes to stdout through note()/step(); only fail() and
## R's own warnings go to stderr
## Nextflow keeps the two streams in separate files, so splitting progress
## messages across both is what makes a failed task hard to read

## Wall-clock start of the run, used by step() and finish()
## Each CLI overwrites it as its first statement, so sourcing order is irrelevant
SCRIPT_START <- Sys.time()

fail <- function(...) stop(paste0(...), call. = FALSE)

## Unconditional log line, flushed:
## the last line in the log before a crash is then really the last thing that happened
note <- function(...) {
  cat(paste0(..., "\n"))
  flush(stdout())
}

## Log line prefixed with elapsed time since SCRIPT_START, for phase boundaries
## In a task that hangs, the last [hh:mm:ss] says where it hung
step <- function(...) {
  secs <- as.integer(difftime(Sys.time(), SCRIPT_START, units = "secs"))
  note(sprintf("[%02d:%02d:%02d] ", secs %/% 3600L, (secs %% 3600L) %/% 60L,
               secs %% 60L), ...)
}

## Script banner, printed before anything else happens
banner <- function(script, title) {
  cat("\n========================================\n")
  cat(script, ": ", title, "\n", sep = "")
  cat("========================================\n\n")
}

## Assert a package is installed at >= a version, and log what was found
## The CLIs call everything through `pkg::fn()`, so nothing is attached here
require_package <- function(package, minimum = NULL) {
  if (!requireNamespace(package, quietly = TRUE)) {
    fail("R package '", package, "' is required",
         if (is.null(minimum)) "" else paste0(" (minimum version ", minimum, ")"))
  }
  installed <- utils::packageVersion(package)
  if (!is.null(minimum) && installed < package_version(minimum)) {
    fail("R package '", package, "' is ", installed,
         "; version >= ", minimum, " is required")
  }
  cat(paste0("  ", package, " ", installed, "\n"))
  invisible(installed)
}

## Attach a package and log its version, as extra/rules/00_init_run.R does
## Only for packages whose infix/NSE syntax needs them on the search path
## (optparse, data.table)
load_pckg <- function(pkg, minimum = NULL) {
  require_package(pkg, minimum)
  suppressPackageStartupMessages(library(package = pkg, character.only = TRUE))
}

## Echo the resolved parameters -- the "what am I actually running" record
## `params` is a named list; NULL/empty prints as <none>, logicals as yes/no
report_params <- function(params) {
  if (!length(params)) return(invisible(NULL))
  width <- max(nchar(names(params)))
  for (name in names(params)) {
    value <- params[[name]]
    keys <- names(value)
    label <- if (is.null(value) || !length(value)) {
      "<none>"
    } else if (is.logical(value)) {
      ifelse(value, "yes", "no")
    } else if (!is.null(keys) && any(nzchar(keys))) {
      # key=value options (--rank-letters, --restrict) read back as they were given
      paste0(keys, "=", as.character(value))
    } else {
      as.character(value)
    }
    cat(paste0("  ", formatC(paste0(name, ":"), width = width + 1L, flag = "-"),
               "  ", paste(label, collapse = ", "), "\n"))
  }
  invisible(NULL)
}

## Log an external tool's resolved path together with the version it reports,
## so the log answers "which vsearch did this task actually use?"
report_tool <- function(name, path) {
  if (is.null(path)) {
    cat(paste0("  ", name, ":  <not used>\n"))
  } else {
    cat(paste0("  ", name, ":  ", path, " (",
               tool_version(path, lines = 1L), ")\n"))
  }
  invisible(path)
}

## Report missing required options the way 00_init_run.R does, then exit 2
require_options <- function(opt, required) {
  missing <- required[vapply(
    required,
    function(x) {
      value <- opt[[x]]
      is.null(value) || !length(value) || !nzchar(as.character(value)[[1L]])
    },
    logical(1)
  )]
  if (length(missing)) {
    cat("ERROR: missing required option(s): ",
        paste0("--", missing, collapse = ", "), "\n", sep = "", file = stderr())
    quit(save = "no", status = 2L)
  }
  invisible(NULL)
}

## Closing block: elapsed time and session info, as 00_init_run.R does
finish <- function(start_time = SCRIPT_START) {
  cat("\nAll done.\n")
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "min"))
  cat("\n========================================\n")
  cat("Elapsed time: ", round(elapsed, 3), " minutes\n", sep = "")
  cat("========================================\n")
  cat("\nSession info:\n")
  print(utils::sessionInfo())
  cat("\n")
  invisible(NULL)
}

## ---------------------------------------------------------------------------
## option validation
## ---------------------------------------------------------------------------
##
## optparse applies the type and the default; these enforce range and enum, and
## split the composite (CSV / key=value) options it hands back as plain strings

as_count <- function(x, name, default) {
  if (is.null(x)) return(as.integer(default))
  value <- suppressWarnings(as.numeric(x))
  if (length(value) != 1L || !is.finite(value) || value < 1 || value != floor(value)) {
    fail("--", name, " must be a positive integer")
  }
  as.integer(value)
}

as_number <- function(x, name, default, lower = 0, inclusive = TRUE) {
  if (is.null(x)) return(default)
  value <- suppressWarnings(as.numeric(x))
  ok <- length(value) == 1L && is.finite(value) &&
    if (inclusive) value >= lower else value > lower
  if (!ok) {
    fail("--", name, " must be a finite number ",
         if (inclusive) ">= " else "> ", lower)
  }
  value
}

as_probability <- function(x, name, default) {
  if (is.null(x)) return(default)
  value <- suppressWarnings(as.numeric(x))
  if (length(value) != 1L || !is.finite(value) || value < 0 || value > 1) {
    fail("--", name, " must be between 0 and 1")
  }
  value
}

as_choice <- function(x, name, choices, default = choices[[1L]]) {
  if (is.null(x)) return(default)
  value <- tolower(trimws(x))
  if (!value %in% choices) {
    fail("--", name, " must be one of: ", paste(choices, collapse = ", "))
  }
  value
}

## Split a comma-separated option value into a trimmed character vector.
as_list_arg <- function(x, default = character()) {
  if (is.null(x)) return(default)
  parts <- trimws(strsplit(as.character(x), ",", fixed = TRUE)[[1L]])
  parts[nzchar(parts)]
}

## Parse a `key=value,key=value` option into a named character vector.
as_kv_arg <- function(x, name, default = character()) {
  if (is.null(x)) return(default)
  parts <- as_list_arg(x)
  if (!length(parts)) return(default)
  if (!all(grepl("^[^=]+=[^=]*$", parts))) {
    fail("--", name, " must be a comma-separated list of key=value pairs")
  }
  out <- sub("^[^=]+=", "", parts)
  names(out) <- sub("=.*$", "", parts)
  out
}

## ---------------------------------------------------------------------------
## external tools
## ---------------------------------------------------------------------------
##
## The pipeline no longer passes binary paths: vsearch, usearch and sintaxer ship
## inside the container and are found on PATH
## The options below stay so that a user can point a run at a different build

## Locate an external executable, autodetecting it when no path is given
##
## Precedence: explicit path (--vsearch / --usearch / --sintaxer) -> the upper-case
## environment variable of the same name ($VSEARCH / $USEARCH / $SINTAXER) -> PATH
## The environment variable tier matters in practice:
## an `.Rprofile` that sanitizes `PATH` (e.g., stripping a conda prefix) can hide a perfectly good binary from `Sys.which()`
find_executable <- function(name, explicit = NULL) {
  if (!is.null(explicit) && nzchar(explicit)) {
    if (file.access(explicit, 1L) != 0L) {
      fail(name, " is not executable: ", explicit)
    }
    return(normalizePath(explicit))
  }
  from_env <- Sys.getenv(toupper(name), "")
  if (nzchar(from_env)) {
    if (file.access(from_env, 1L) != 0L) {
      fail(name, " from $", toupper(name), " is not executable: ", from_env)
    }
    return(normalizePath(from_env))
  }
  found <- Sys.which(name)
  if (!nzchar(found)) {
    fail(name, " is required but was not found. Put it on PATH, set $",
         toupper(name), ", or pass --", name)
  }
  unname(normalizePath(found))
}

find_vsearch <- function(explicit = NULL) find_executable("vsearch", explicit)
find_sintaxer <- function(explicit = NULL) find_executable("sintaxer", explicit)

## Capture a tool's self-reported version
##
## `lines` caps how much is kept: the JSON sidecars record the whole banner for
## provenance, while the startup log wants only the first line -- vsearch prints
## its citation, RAM and zlib flags after it
tool_version <- function(executable, flag = "--version", lines = NULL) {
  out <- tryCatch(
    system2(executable, flag, stdout = TRUE, stderr = TRUE),
    error = function(e) {
      fail("Could not execute ", executable, ": ", conditionMessage(e))
    }
  )
  # only the capped form is trimmed; the uncapped string goes into the JSON
  # sidecars verbatim, so it must stay byte-stable
  if (!is.null(lines)) out <- trimws(utils::head(out, lines))
  paste(out, collapse = " ")
}

## ---------------------------------------------------------------------------
## distance and clustering configuration
## ---------------------------------------------------------------------------

## Distance configuration for --dist-method / --dist-penalties
##
## WFA2 and KSW2 take the same penalties (match, mismatch, gap_open,
## gap_extend, gap_open2, gap_extend2); the other methods take none.
## The defaults of both are edit distance.
build_dist_config <- function(method, penalties, threads = 1L, usearch = NULL) {
  if (method %in% PENALIZED_DIST_METHODS) {
    known <- c("match", "mismatch", "gap_open", "gap_extend", "gap_open2",
               "gap_extend2")
    unknown <- setdiff(names(penalties), known)
    if (length(unknown)) {
      fail("Unknown ", toupper(method), " penalty: ", paste(unknown, collapse = ", "),
           " (known: ", paste(known, collapse = ", "), ")")
    }
    values <- lapply(penalties, function(x) {
      v <- suppressWarnings(as.integer(x))
      if (is.na(v)) fail("--dist-penalties must be integers")
      v
    })
    builder <- switch(method, wfa2 = optimotu::dist_wfa2, ksw2 = optimotu::dist_ksw2)
    return(do.call(builder, values))
  }
  if (length(penalties)) {
    fail("--dist-penalties only applies to --dist-method ",
         paste(PENALIZED_DIST_METHODS, collapse = " or "))
  }
  switch(
    method,
    edlib = optimotu::dist_edlib(),
    hybrid = optimotu::dist_hybrid(),
    hamming = optimotu::dist_hamming(),
    # usearch_ncpu must be given explicitly: dist_usearch() stores
    # as.integer(NULL), i.e. integer(0), which trips a length assertion deeper
    # in optimotu. usearch does the alignment here, so give it the threads.
    usearch = optimotu::dist_usearch(usearch = usearch, usearch_ncpu = threads)
  )
}

build_clust_config <- function(method) {
  switch(
    method,
    tree = optimotu::clust_tree(),
    matrix = optimotu::clust_matrix(),
    index = optimotu::clust_index(),
    slink = optimotu::clust_slink()
  )
}

## ---------------------------------------------------------------------------
## ranks
## ---------------------------------------------------------------------------

## Resolve the `--ranks` option, validating shape.
resolve_ranks <- function(x, default = DEFAULT_RANKS) {
  ranks <- as_list_arg(x, default)
  if (length(ranks) < 2L) fail("--ranks needs at least two ranks")
  if (anyDuplicated(ranks)) fail("--ranks contains duplicate rank names")
  ranks
}

## Add a constant top rank to a wide taxonomy table
##
## `optimotu::calc_taxon_thresholds()` asserts that `ranks[1]` holds exactly one
## distinct value, so pan-eukaryotic references need a synthetic root.  Records
## that already carry a value at `root_rank` keep it; missing ones are filled, so
## this is a no-op on a reference whose own top rank is already constant
add_root_rank <- function(tab, root_rank, root_taxon) {
  if (!root_rank %in% names(tab)) {
    tab[[root_rank]] <- NA_character_
    tab <- tab[, c("seq_id", root_rank,
                   setdiff(names(tab), c("seq_id", root_rank))), drop = FALSE]
  }
  value <- trimws(as.character(tab[[root_rank]]))
  value[!nzchar(value) | is.na(value)] <- root_taxon
  tab[[root_rank]] <- value
  tab
}

## ---------------------------------------------------------------------------
## thresholds tables
## ---------------------------------------------------------------------------

## Treat names of the last rank with fewer than `min_seqs` records as unknown
##
## A species with one sequence has no within-species pairs: at every threshold
## above 0 it can only be merged with another species, so a subset made mostly
## of such species scores best at (or near) 0 whatever the real within-species
## variation is.  Names are counted per (parent, name) over the table given, so
## optimotu_partition_reference.R and every optimotu_train.R chunk, which all
## read the same training taxonomy, mask the same names.
##
## @param taxonomy data.frame with one column per rank
## @return list(taxonomy = the table with rare names set to NA,
##              names = number of names masked, records = records affected)
mask_rare_tip_names <- function(taxonomy, ranks, min_seqs) {
  if (min_seqs <= 1L) return(list(taxonomy = taxonomy, names = 0L, records = 0L))
  tip <- ranks[[length(ranks)]]
  parent <- if (length(ranks) > 1L) ranks[[length(ranks) - 1L]] else NULL
  key <- if (is.null(parent)) taxonomy[[tip]] else {
    paste(taxonomy[[parent]], taxonomy[[tip]], sep = "|")
  }
  # read_reference_taxonomy() can leave unknown names as "" (cleaning turns them
  # into NA later); counting them would add one "Genus|" name per genus
  key[is.na(taxonomy[[tip]]) | !nzchar(trimws(taxonomy[[tip]]))] <- NA_character_
  size <- table(key)
  rare <- !is.na(key) & as.vector(size[key]) < min_seqs
  taxonomy[[tip]][rare] <- NA_character_
  list(taxonomy = taxonomy, names = sum(size < min_seqs), records = sum(rare))
}

## Ranks below the root that have no root-level threshold row, per measure
##
## optimotu::calc_taxon_thresholds() falls back up the taxonomy and needs a
## root-level row for every rank; optimotu_cluster.R rejects a table without one.
## Coverage is per measure: a measure that is undefined on a single-taxon subset
## (MCC) simply has no row there, which only matters if that measure is the one
## selected at clustering time
##
## @param thresholds table with `rank`, `superrank` and `measure` columns
## @return named list, measure -> the ranks it lacks a root-level row for
root_fallback_gaps <- function(thresholds, ranks, measures) {
  root <- thresholds$superrank == ranks[[1L]]
  gaps <- lapply(measures, function(m) {
    setdiff(ranks[-1L], unique(thresholds$rank[root & thresholds$measure == m]))
  })
  names(gaps) <- measures
  gaps
}

## ---------------------------------------------------------------------------
## FASTA input
## ---------------------------------------------------------------------------

## Read a FASTA file, returning a DNAStringSet whose names are bare sequence IDs.
##
## IDs are the header up to the first whitespace, `;` or `|`; a trailing
## `;size=N` is dropped.
read_fasta <- function(path, label = "FASTA") {
  if (!file.exists(path) || file.access(path, 4L) != 0L) {
    fail(label, " is not readable: ", path)
  }
  seqs <- tryCatch(
    Biostrings::readDNAStringSet(path, format = "fasta"),
    error = function(e) fail("Cannot parse ", label, ": ", conditionMessage(e))
  )
  if (length(seqs) == 0L) fail(label, " contains no sequences")
  headers <- names(seqs)
  if (is.null(headers) || anyNA(headers) || any(!nzchar(headers))) {
    fail(label, " has an empty FASTA header")
  }
  names(seqs) <- parse_seq_ids(headers, label)
  seqs
}

## Strip labels down to bare sequence IDs, without requiring uniqueness.
##
## Use this for hit tables and other places where a label legitimately repeats.
strip_seq_ids <- function(labels) {
  sub(";$", "", sub("[[:space:];|].*$", "", labels))
}

## Strip a FASTA header down to a bare, unique sequence ID.
parse_seq_ids <- function(headers, label = "FASTA") {
  ids <- strip_seq_ids(headers)
  if (any(!nzchar(ids))) fail(label, " has an empty sequence identifier")
  duplicates <- unique(ids[duplicated(ids)])
  if (length(duplicates)) {
    fail(label, " has duplicate sequence IDs: ",
         paste(utils::head(duplicates, 10L), collapse = ", "))
  }
  ids
}

## Parse `size=N` abundance annotations out of FASTA headers, defaulting to 1.
parse_sizes <- function(headers) {
  captures <- regmatches(
    headers,
    regexec("(^|[;[:space:]])size=([0-9]+)($|[;[:space:]])", headers, perl = TRUE)
  )
  size <- vapply(
    captures,
    function(x) if (length(x)) as.numeric(x[[3L]]) else 1,
    numeric(1)
  )
  if (any(!is.finite(size)) || any(size < 1)) {
    fail("Every size=N abundance must be a positive integer")
  }
  size
}

## ---------------------------------------------------------------------------
## miscellaneous
## ---------------------------------------------------------------------------

file_md5 <- function(path) unname(tools::md5sum(path))

## Describe an input or output file for a JSON sidecar
##
## Only the base name is recorded: absolute paths leak user home directories
## and point into Nextflow work directories that do not outlive the run
file_record <- function(path, ...) {
  c(list(file = basename(path), md5 = file_md5(path)), list(...))
}

sanitize_id <- function(x, label = "identifier") {
  x <- gsub("[^A-Za-z0-9_.-]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  if (!nzchar(x)) fail(label, " contains no usable characters")
  x
}

## Write a set of files atomically: stage into a temporary sibling directory,
## then rename each into place.
##
## @param outputs final paths, in the same order as `writer` expects
## @param writer function of one argument (the staged paths) that writes them
## @param overwrite replace existing outputs
with_atomic_outputs <- function(outputs, writer, overwrite = FALSE) {
  parent <- unique(dirname(outputs))
  for (dir in parent) {
    if (!dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(dir)) fail("Cannot create output directory: ", dir)
  }
  conflicts <- outputs[file.exists(outputs)]
  if (length(conflicts) && !isTRUE(overwrite)) {
    fail("Output file(s) already exist; use --overwrite: ",
         paste(conflicts, collapse = ", "))
  }
  if (any(dir.exists(outputs))) fail("An output path is an existing directory")

  stage <- tempfile(pattern = ".optimotu-stage-", tmpdir = parent[[1L]])
  if (!dir.create(stage)) fail("Cannot create staging directory: ", stage)
  committed <- FALSE
  on.exit(
    if (!committed && dir.exists(stage)) {
      unlink(stage, recursive = TRUE, force = TRUE)
    },
    add = TRUE
  )

  staged <- file.path(stage, basename(outputs))
  writer(staged)

  unlink(conflicts, force = TRUE)
  for (i in seq_along(outputs)) {
    if (!file.rename(staged[[i]], outputs[[i]])) {
      fail("Could not move staged output into place: ", outputs[[i]])
    }
  }
  unlink(stage, recursive = TRUE, force = TRUE)
  committed <- TRUE
  invisible(outputs)
}

write_json_file <- function(x, path) {
  jsonlite::write_json(x, path, auto_unbox = TRUE, pretty = TRUE,
                       null = "null", digits = NA)
}

## ---------------------------------------------------------------------------
## reference taxonomy parsing
## ---------------------------------------------------------------------------
##
## Replaces `optimotu.pipeline::parse_reference_taxonomy()`, which is the only
## reason that package (and its ape / fastqindexr / targets dependency stack)
## would be needed here.

## Guess the taxonomy annotation format of a vector of FASTA headers.
detect_taxonomy_format <- function(headers, ranks = DEFAULT_RANKS) {
  probe <- utils::head(headers, 1000L)
  if (all(grepl(";tax=", probe, fixed = TRUE))) return("sintax")
  if (all(grepl("[a-z]__", probe))) {
    # UNITE puts the taxonomy in the second pipe-delimited field; Eukaryome
    # appends it to the ID with ';'.  Both are handled by the same token scan,
    # so the distinction is cosmetic.
    return(if (all(grepl("|", probe, fixed = TRUE))) "unite" else "eukaryome")
  }
  if (all(lengths(strsplit(probe, "|", fixed = TRUE)) >= 4L)) return("bold")
  # `ID;Kingdom;Phylum;...;Species` -- nothing labels the ranks, so they can only
  # be read positionally, and that is only safe when every header carries exactly
  # as many fields as there are ranks to fill.  A ragged file is a mistake worth
  # reporting rather than silently mis-assigning.
  n_fields <- lengths(strsplit(probe, ";", fixed = TRUE))
  if (length(unique(n_fields)) == 1L &&
      n_fields[[1L]] == length(setdiff(ranks, "root")) + 1L) {
    return("lineage")
  }
  fail("Could not recognize the reference taxonomy format. Supported: ",
       "eukaryome, unite, sintax, bold, lineage, tsv (use --taxonomy-format). ",
       "For a bare `ID;Kingdom;...;Species` header, every record must carry ",
       length(setdiff(ranks, "root")) + 1L, " ';'-separated fields to match ",
       "--ranks; this input has ",
       paste(sort(unique(n_fields)), collapse = "/"), ".")
}

## Parse `k__Fungi;p__Basidiomycota;...` style headers (Eukaryome, UNITE).
##
## Tokens are located by prefix rather than by position, so the surrounding
## delimiters and any extra header fields do not matter.
parse_prefixed_headers <- function(headers, ranks, rank_letters) {
  matches <- regmatches(headers, gregexpr("[A-Za-z]__[^;|\t]*", headers))
  tokens <- unlist(matches, use.names = FALSE)
  assemble_taxonomy(
    ids = parse_seq_ids(headers),
    ranks = ranks,
    record = rep.int(seq_along(matches), lengths(matches)),
    letter = tolower(substr(tokens, 1L, 1L)),
    taxon = trimws(sub("^[A-Za-z]__", "", tokens)),
    rank_letters = rank_letters
  )
}

## Parse `ID;tax=d:Fungi,p:Basidiomycota,...;` style headers (SINTAX).
parse_sintax_headers <- function(headers, ranks, rank_letters) {
  lineage <- sub(";.*$", "", sub("^.*;tax=", "", headers))
  matches <- regmatches(lineage, gregexpr("[A-Za-z]:[^,]*", lineage))
  tokens <- unlist(matches, use.names = FALSE)
  assemble_taxonomy(
    ids = parse_seq_ids(headers),
    ranks = ranks,
    record = rep.int(seq_along(matches), lengths(matches)),
    letter = tolower(substr(tokens, 1L, 1L)),
    taxon = trimws(sub("^[A-Za-z]:", "", tokens)),
    rank_letters = rank_letters
  )
}

## Parse `ID|marker|country|kingdom,phylum,...` style headers (BOLD).
parse_bold_headers <- function(headers, ranks) {
  bold_ranks <- c("kingdom", "phylum", "class", "order", "family", "subfamily",
                  "tribe", "genus", "species", "subspecies")
  fields <- strsplit(headers, "|", fixed = TRUE)
  lineage <- vapply(fields, function(x) if (length(x) >= 4L) x[[4L]] else "",
                    character(1))
  matches <- strsplit(lineage, ",", fixed = TRUE)
  position <- unlist(lapply(lengths(matches), seq_len), use.names = FALSE)
  assemble_taxonomy(
    ids = parse_seq_ids(headers),
    ranks = ranks,
    record = rep.int(seq_along(matches), lengths(matches)),
    letter = ifelse(position <= length(bold_ranks), bold_ranks[position],
                    NA_character_),
    taxon = trimws(unlist(matches, use.names = FALSE)),
    rank_letters = stats::setNames(bold_ranks, bold_ranks)
  )
}

## Parse `ID;Kingdom;Phylum;...;Species` style headers (positional lineage).
##
## Unlike the prefixed and SINTAX formats there is no per-token rank label, so
## field *position* carries the meaning: field 1 is the sequence ID and fields
## 2..n map onto `ranks` in order, skipping any synthetic root rank (which by
## definition is not in the header).
parse_lineage_headers <- function(headers, ranks) {
  target <- setdiff(ranks, "root")
  parts <- strsplit(headers, ";", fixed = TRUE)
  n <- lengths(parts)
  # 0 for the ID field, 1..length(target) for the rank fields
  position <- unlist(lapply(n, seq_len), use.names = FALSE) - 1L
  keep <- position >= 1L & position <= length(target)
  assemble_taxonomy(
    ids = parse_seq_ids(headers),
    ranks = ranks,
    record = rep.int(seq_along(parts), n)[keep],
    letter = target[position[keep]],
    taxon = trimws(unlist(parts, use.names = FALSE))[keep],
    rank_letters = stats::setNames(target, target)
  )
}

empty_taxonomy <- function(ids, ranks) {
  out <- data.frame(seq_id = ids, stringsAsFactors = FALSE)
  out[ranks] <- NA_character_
  out
}

## Scatter long-format (record, rank-key, taxon) triples into a wide table.
##
## Filling via matrix indexing rather than per-token column assignment keeps this
## linear; a naive loop copies a whole column per token, which turns parsing a
## million-record database into hours.
assemble_taxonomy <- function(ids, ranks, record, letter, taxon, rank_letters) {
  wide <- matrix(NA_character_, nrow = length(ids), ncol = length(ranks),
                 dimnames = list(NULL, ranks))
  if (length(record)) {
    mapped <- unname(rank_letters[letter])
    keep <- !is.na(mapped) & mapped %in% ranks & !is.na(taxon) & nzchar(taxon)
    if (any(keep)) {
      wide[cbind(record[keep], match(mapped[keep], ranks))] <- taxon[keep]
    }
  }
  data.frame(seq_id = ids, wide, stringsAsFactors = FALSE, check.names = FALSE)
}

## Parse a tabular taxonomy file.
##
## Accepts a header row naming `seq_id` (or any first column) plus one column per
## rank, or a two-column `id + delimited lineage` layout.
parse_taxonomy_tsv <- function(path, ranks, rank_letters) {
  # na.strings = "NA" reproduces read.delim's default exactly: a literal "NA"
  # cell becomes NA, an empty one stays an empty string
  tab <- data.table::fread(path, sep = "\t", header = TRUE, quote = "",
                           colClasses = "character", check.names = FALSE,
                           na.strings = "NA", showProgress = FALSE,
                           data.table = FALSE)
  if (!ncol(tab)) fail("Taxonomy table has no columns: ", path)
  names(tab)[[1L]] <- "seq_id"
  present <- intersect(ranks, names(tab))
  if (length(present)) {
    out <- empty_taxonomy(tab$seq_id, ranks)
    out[present] <- tab[present]
    return(out)
  }
  lineage_col <- if ("taxonomy" %in% names(tab)) "taxonomy" else names(tab)[[ncol(tab)]]
  lineage <- tab[[lineage_col]]
  if (all(grepl("[A-Za-z]__", utils::head(lineage, 1000L)))) {
    return(parse_prefixed_headers(paste0(tab$seq_id, ";", lineage), ranks,
                                  rank_letters))
  }
  if (all(grepl("^[A-Za-z]:", utils::head(lineage, 1000L)))) {
    return(parse_sintax_headers(paste0(tab$seq_id, ";tax=", lineage, ";"), ranks,
                                rank_letters))
  }
  # bare delimited lineage: assign positionally, skipping any root rank
  #
  # `tstrsplit()` transposes the split in one pass, so each rank is filled with a
  # single vector assignment.  The nested per-cell loop this replaces copied a
  # whole column for every cell it wrote
  target <- setdiff(ranks, "root")
  out <- empty_taxonomy(tab$seq_id, ranks)
  fields <- data.table::tstrsplit(lineage, "[,;|]", fixed = FALSE)
  for (j in seq_len(min(length(fields), length(target)))) {
    value <- trimws(fields[[j]])
    value[!nzchar(value)] <- NA_character_
    out[[target[[j]]]] <- value
  }
  out
}

## Read a reference taxonomy from FASTA headers or a tabular file.
##
## @param path FASTA or TSV file
## @param ranks target rank vector
## @param format one of `TAXONOMY_FORMATS`
## @param rank_letters letter -> rank mapping for prefixed / SINTAX formats
## @return data.frame with `seq_id` plus one column per rank
read_reference_taxonomy <- function(path, ranks, format = "auto",
                                    rank_letters = DEFAULT_RANK_LETTERS,
                                    headers = NULL) {
  if (!file.exists(path) || file.access(path, 4L) != 0L) {
    fail("Reference taxonomy is not readable: ", path)
  }
  is_fasta <- grepl("\\.(fa|fas|fna|fasta)(\\.gz)?$", path, ignore.case = TRUE)
  if (identical(format, "tsv") || (!is_fasta && identical(format, "auto"))) {
    return(parse_taxonomy_tsv(path, ranks, rank_letters))
  }
  if (is.null(headers)) {
    headers <- names(Biostrings::fasta.seqlengths(path, seqtype = "DNA"))
  }
  if (!length(headers)) fail("No FASTA headers found in ", path)
  if (identical(format, "auto")) format <- detect_taxonomy_format(headers, ranks)
  switch(
    format,
    eukaryome = parse_prefixed_headers(headers, ranks, rank_letters),
    unite     = parse_prefixed_headers(headers, ranks, rank_letters),
    sintax    = parse_sintax_headers(headers, ranks, rank_letters),
    bold      = parse_bold_headers(headers, ranks),
    lineage   = parse_lineage_headers(headers, ranks),
    fail("Unsupported --taxonomy-format for FASTA input: ", format)
  )
}

## Normalize, clean and validate a wide taxonomy table.
##
## Blanks become NA, `optimotu::clean_taxonomy()` NAs out placeholder names
## ("unclassified", "sp.", "incertae sedis", ...), lineages are truncated at the
## first gap so rank paths stay contiguous, and the root rank is checked to hold
## exactly one value -- an assertion `optimotu::calc_taxon_thresholds()` makes
## internally and would otherwise fail on much later.
clean_reference_taxonomy <- function(tab, ranks, label = "Reference taxonomy",
                                     require_single_root = TRUE,
                                     qualify = TRUE,
                                     keep_orphan_ranks = FALSE) {
  missing_columns <- setdiff(c("seq_id", ranks), names(tab))
  if (length(missing_columns)) {
    fail(label, " lacks required columns: ", paste(missing_columns, collapse = ", "))
  }
  # data.table::set() updates columns in place: the data.frame form copied a
  # whole column for every one of the ~2 * length(ranks) rewrites below
  tab <- data.table::as.data.table(tab)[, c("seq_id", ranks), with = FALSE]
  for (col in c("seq_id", ranks)) {
    data.table::set(tab, j = col, value = as.character(tab[[col]]))
  }

  if (anyNA(tab$seq_id) || any(!nzchar(trimws(tab$seq_id)))) {
    fail(label, " has an empty sequence ID")
  }
  duplicates <- unique(tab$seq_id[duplicated(tab$seq_id)])
  if (length(duplicates)) {
    fail(label, " has duplicate sequence IDs: ",
         paste(utils::head(duplicates, 10L), collapse = ", "))
  }

  for (rank in ranks) {
    value <- trimws(tab[[rank]])
    value[!nzchar(value)] <- NA_character_
    # ',' and ';' are structural in SINTAX headers and in our own TSV outputs
    data.table::set(tab, j = rank, value = gsub("[,;]", "_", value))
  }
  # optimotu::clean_taxonomy() asserts a data.frame and writes columns with
  # `[[<-`, which on a data.table warns about the shallow-copy selfref, so hand
  # it a plain data.frame and take the table back afterwards
  data.table::setDF(tab)
  tab <- optimotu::clean_taxonomy(tab, ranks = ranks)
  data.table::setDT(tab)

  if (keep_orphan_ranks) {
    # A placeholder in the middle of a lineage ("Helotiales_fam_Incertae_sedis")
    # would otherwise truncate away everything below it, discarding real species
    # names -- 4% of them in UNITE v10 Fungi.  Bridge such a gap with a synthetic
    # node named for the deepest genuinely assigned ancestor, but only where an
    # assigned rank actually exists further down; a lineage that simply stops
    # must stay unassigned below that point.
    #
    # This is a trade-off, not a free win: every orphan under one ancestor is
    # lumped into a single synthetic taxon, so unrelated lineages become
    # "cotaxonomic" at the bridged rank.  Hence opt-in.
    assigned <- !is.na(as.matrix(tab[, ranks, with = FALSE]))
    # has_below[, i] : does this record carry an assignment strictly below rank i?
    has_below <- matrix(FALSE, nrow = nrow(tab), ncol = length(ranks))
    acc <- rep(FALSE, nrow(tab))
    for (i in rev(seq_along(ranks))) {
      has_below[, i] <- acc
      acc <- acc | assigned[, i]
    }
    ancestor <- rep(NA_character_, nrow(tab))
    bridged <- stats::setNames(integer(length(ranks)), ranks)
    for (i in seq_along(ranks)) {
      rank <- ranks[[i]]
      here <- assigned[, i]
      ancestor[here] <- tab[[rank]][here]
      gap <- which(!here & !is.na(ancestor) & has_below[, i])
      if (length(gap)) {
        data.table::set(tab, i = gap, j = rank,
                        value = paste0(ancestor[gap], "_incertae_", rank))
        bridged[[rank]] <- length(gap)
      }
    }
    bridged <- bridged[bridged > 0L]
    if (length(bridged)) {
      note(label, ": bridged Incertae sedis gaps with synthetic taxa (",
           paste(sprintf("%s=%d", names(bridged), bridged), collapse = " "), ")")
    }
  }

  # truncate each lineage at its first gap so no rank sits below a missing parent
  for (i in seq_along(ranks)[-1L]) {
    parent <- ranks[[i - 1L]]
    child <- ranks[[i]]
    orphaned <- which(is.na(tab[[parent]]))
    if (length(orphaned)) {
      data.table::set(tab, i = orphaned, j = child, value = NA_character_)
    }
  }

  root_rank <- ranks[[1L]]
  root_taxa <- unique(tab[[root_rank]])
  if (anyNA(root_taxa)) {
    fail(label, " has records with no ", root_rank, " assignment")
  }
  if (require_single_root && length(root_taxa) != 1L) {
    fail(label, " must have exactly one distinct '", root_rank,
         "' value (found ", length(root_taxa),
         "); add a synthetic root rank via --ranks/--root-taxon")
  }

  if (qualify) {
    tab <- qualify_homonyms(tab, ranks)
    counts <- attr(tab, "qualified")
    counts <- counts[counts > 0L]
    if (length(counts)) {
      note(label, ": qualified ambiguous taxon names with their parent (",
           paste(sprintf("%s=%d", names(counts), counts), collapse = " "), ")")
    }
    attr(tab, "qualified") <- NULL
  }
  # every caller indexes the result with data.frame semantics
  # (`tab[i, , drop = FALSE]`, `tab[, ranks, drop = FALSE]`), so hand back a
  # plain data.frame and let each script opt into data.table where it needs to
  if (data.table::is.data.table(tab)) data.table::setDF(tab)
  tab
}

## Disambiguate taxon names that occur under more than one parent.
##
## `optimotu::summarize_by_rank()` groups reference sequences by the raw taxon
## name at each rank, and it enumerates every (superrank, subrank) pair -- so it
## scores "species within family" as well as "species within genus".  Databases
## that store bare epithets rather than binomials (Eukaryome is one) therefore
## merge homonyms from different genera into a single "true" taxon, which
## silently biases the optimized thresholds.  Prefixing only the ambiguous names
## with their parent restores a one-name-one-taxon mapping while leaving already
## unambiguous names untouched.
##
## @return the table, with an attribute `qualified` counting the rewrites
qualify_homonyms <- function(tab, ranks) {
  qualified <- stats::setNames(integer(length(ranks)), ranks)
  dt <- data.table::as.data.table(tab)
  # top-down, so that a parent is already unambiguous when its children are fixed
  for (i in seq_along(ranks)[-1L]) {
    parent <- ranks[[i - 1L]]
    child <- ranks[[i]]
    known <- !is.na(dt[[child]]) & !is.na(dt[[parent]])
    if (!any(known)) next
    # column names are deliberately not the rank names: `by = child` inside `[`
    # would resolve to the loop variable's value, not to the column
    pairs <- unique(data.table::data.table(ch = dt[[child]][known],
                                           pa = dt[[parent]][known]))
    homonyms <- pairs[, list(n_parents = .N), by = "ch"][n_parents > 1L, ch]
    if (!length(homonyms)) next
    hit <- which(known & dt[[child]] %in% homonyms)
    data.table::set(dt, i = hit, j = child,
                    value = paste(dt[[parent]][hit], dt[[child]][hit]))
    qualified[[child]] <- length(homonyms)
  }
  data.table::setDF(dt)
  attr(dt, "qualified") <- qualified
  dt
}

## Build SINTAX-formatted FASTA headers from a wide taxonomy table
##
## One melt plus one grouped paste, rather than a `vapply` over every record:
## the per-row form indexed the whole table once per sequence, which on a
## reference database of a million records was the slowest thing in this file
sintax_headers <- function(taxonomy, ranks, rank_letters = DEFAULT_RANK_LETTERS,
                           skip_ranks = "root") {
  used_ranks <- setdiff(ranks, skip_ranks)
  codes <- vapply(
    used_ranks,
    function(rank) {
      hit <- names(rank_letters)[rank_letters == rank]
      if (!length(hit)) fail("No SINTAX letter code known for rank '", rank, "'")
      hit[[1L]]
    },
    character(1)
  )
  long <- data.table::melt(
    data.table::as.data.table(taxonomy)[, c("seq_id", used_ranks), with = FALSE],
    id.vars = "seq_id", measure.vars = used_ranks,
    variable.name = "rank", value.name = "taxon",
    variable.factor = FALSE, na.rm = TRUE
  )
  # melt stacks the measure columns in order, so within a group the terms come
  # out most- to least-inclusive, exactly as the per-row version emitted them
  lineage <- long[, list(tax = paste0(codes[rank], ":", taxon, collapse = ",")),
                  by = "seq_id"]
  # records with nothing assigned at any used rank drop out of the melt
  tax <- lineage$tax[match(taxonomy$seq_id, lineage$seq_id)]
  tax[is.na(tax)] <- ""
  paste0(taxonomy$seq_id, ";tax=", tax, ";")
}
