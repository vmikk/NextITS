#!/usr/bin/env Rscript

## Usage example:
# extract_itsx_regions.R \
#   --fasta     Dereplicated.fa.gz \
#   --positions ITSx.positions.txt.gz \
#   --region    ITS \
#   --output    ITS.fasta.gz \
#   --report    ITS.extraction.tsv.gz

## Input data:
# - FASTA file with sequences
# - ITSx positions file: TSV, no header, eight columns
#   (SeqID, SeqLen, SSU, ITS1, 5.8S, ITS2, LSU, Info)

## Notes:
# - ITSx coordinates and Biostrings::subseq() are 1-based and end-inclusive
# - Coordinates must refer to the orientation of the supplied FASTA
#   (normally ITSx should be run with --complement F)


## Check time
start_time <- Sys.time()

cat("\nParsing input options and arguments...\n")

suppressPackageStartupMessages(require(optparse))

## Parse arguments
option_list <- list(
  make_option(c("-f", "--fasta"),     action = "store", default = NA, type = "character", help = "Input sequences in FASTA format"),
  make_option(c("-p", "--positions"), action = "store", default = NA, type = "character", help = "ITSx positions file"),
  make_option(c("-r", "--region"),    action = "store", default = NA, type = "character", help = "Region to extract: SSU, ITS, or LSU"),
  make_option(c("-o", "--output"),    action = "store", default = NA, type = "character", help = "Output FASTA file"),
  make_option(c("-a", "--report"),    action = "store", default = NA, type = "character", help = "Output audit table (TSV)")
)
opt <- parse_args(OptionParser(option_list = option_list))


############################################## Validate arguments

is_missing_argument <- function(x){
  is.null(x) || length(x) != 1L || is.na(x) || !nzchar(trimws(x))
}

required_options <- c("fasta", "positions", "region", "output", "report")
missing_options <- required_options[
  vapply(required_options, function(x) is_missing_argument(opt[[x]]), logical(1L))
]

if(length(missing_options) > 0L){
  stop(
    "Missing required option(s): ",
    paste0("--", missing_options, collapse = ", "),
    "\n",
    call. = FALSE)
}

FASTA     <- opt$fasta
POSITIONS <- opt$positions
REGION    <- toupper(trimws(opt$region))
OUTPUT    <- opt$output
REPORT    <- opt$report

if(!REGION %in% c("SSU", "ITS", "LSU")){
  stop("Region must be one of: SSU, ITS, LSU\n", call. = FALSE)
}

for(input_file in c(FASTA, POSITIONS)){
  if(!file.exists(input_file)){
    stop("Input file does not exist: ", input_file, "\n", call. = FALSE)
  }
  if(dir.exists(input_file)){
    stop("Expected an input file but found a directory: ", input_file, "\n", call. = FALSE)
  }
}

canonical_output_path <- function(filepath){
  filepath <- path.expand(filepath)
  output_dir <- dirname(filepath)

  if(dir.exists(filepath)){
    stop("Expected an output file but found a directory: ", filepath, "\n", call. = FALSE)
  }

  ## Resolve an existing output itself so symlinks to input files are detected
  if(file.exists(filepath)){
    return(normalizePath(filepath, winslash = "/", mustWork = TRUE))
  }

  if(!dir.exists(output_dir)){
    stop("Output directory does not exist: ", output_dir, "\n", call. = FALSE)
  }

  file.path(
    normalizePath(output_dir, winslash = "/", mustWork = TRUE),
    basename(filepath))
}

input_paths <- normalizePath(
  c(FASTA, POSITIONS),
  winslash = "/",
  mustWork = TRUE)
output_paths <- c(
  canonical_output_path(OUTPUT),
  canonical_output_path(REPORT))

if(anyDuplicated(output_paths)){
  stop("FASTA output and audit report must use different paths\n", call. = FALSE)
}
if(any(output_paths %in% input_paths)){
  stop("Output paths must not overwrite either input file\n", call. = FALSE)
}


## Log assigned variables
cat("\nParameters specified:\n")
cat(paste("Input sequences (FASTA): ", FASTA,     "\n", sep = ""))
cat(paste("ITSx positions: ",          POSITIONS, "\n", sep = ""))
cat(paste("Region to extract: ",       REGION,    "\n", sep = ""))
cat(paste("Output sequences: ",        OUTPUT,    "\n", sep = ""))
cat(paste("Output audit report: ",     REPORT,    "\n", sep = ""))
cat("\n")


############################################## Load packages

cat("Loading R packages...\n")

load_pckg <- function(pkg){
  suppressPackageStartupMessages(library(package = pkg, character.only = TRUE))
  cat(paste(pkg, packageVersion(pkg), "\n"))
}

load_pckg("data.table")
load_pckg("Biostrings")

cat("\n")


############################################## Helper functions

ITSX_COLUMNS <- c(
  "SeqID", "ITSxLengthRaw", "SSU_Raw", "ITS1_Raw",
  "S58_Raw", "ITS2_Raw", "LSU_Raw", "Info")

empty_positions_table <- function(){
  data.table(
    SeqID         = character(),
    ITSxLengthRaw = character(),
    SSU_Raw       = character(),
    ITS1_Raw      = character(),
    S58_Raw       = character(),
    ITS2_Raw      = character(),
    LSU_Raw       = character(),
    Info          = character())
}

read_itsx_positions <- function(filepath){
  positions <- tryCatch(
    suppressWarnings(
      fread(
        file         = filepath,
        header       = FALSE,
        sep          = "\t",
        quote        = "",
        fill         = FALSE,
        showProgress = FALSE)),
    error = function(e){
      error_message <- conditionMessage(e)
      if(grepl(
        "Input is either empty, fully whitespace",
        error_message,
        fixed = TRUE)){
        return(empty_positions_table())
      }

      stop(
        "Could not read ITSx positions file: ", filepath,
        "\n", error_message, "\n",
        call. = FALSE)
    })

  if(nrow(positions) == 0L && ncol(positions) == 0L){
    return(empty_positions_table())
  }

  if(ncol(positions) != length(ITSX_COLUMNS)){
    stop(
      "ITSx positions file must contain exactly eight tab-separated columns; found ",
      ncol(positions), "\n",
      call. = FALSE)
  }

  setnames(positions, ITSX_COLUMNS)
  positions[ , (ITSX_COLUMNS) := lapply(.SD, as.character), .SDcols = ITSX_COLUMNS ]
  positions
}

parse_itsx_length <- function(x){
  value <- trimws(x)
  value <- sub(
    pattern     = "[[:space:]]+bp\\.?$",
    replacement = "",
    x           = value,
    ignore.case = TRUE)

  valid <- !is.na(value) & grepl("^[0-9]+$", value)
  parsed <- rep(NA_real_, length(value))
  parsed[valid] <- suppressWarnings(as.numeric(value[valid]))
  valid <- valid & is.finite(parsed)

  res <- list(value = parsed, valid = valid)
  return(res)
}

MISSING_COORDINATE_STATUS <- c(
  "not found"   = "not_found",
  "no end"      = "no_end",
  "no start"    = "no_start",
  "overlap lsu" = "overlap_lsu",
  "overlap ssu" = "overlap_ssu")

parse_region_coordinates <- function(x, label){
  raw_value <- trimws(x)
  prefix <- paste0(label, ":")
  has_prefix <- !is.na(raw_value) & startsWith(raw_value, prefix)

  value <- rep(NA_character_, length(raw_value))
  value[has_prefix] <- trimws(substring(raw_value[has_prefix], nchar(prefix) + 1L))

  marker_status <- rep(NA_character_, length(raw_value))
  marker_status[has_prefix] <- unname(
    MISSING_COORDINATE_STATUS[tolower(value[has_prefix])])

  missing <- is.na(raw_value) | !is.na(marker_status)

  start <- rep(NA_real_, length(raw_value))
  end   <- rep(NA_real_, length(raw_value))
  malformed <- !missing & !has_prefix
  status <- rep("malformed", length(raw_value))
  status[is.na(raw_value) | !nzchar(raw_value)] <- "missing"
  status[!is.na(marker_status)] <- marker_status[!is.na(marker_status)]

  candidates <- which(!missing & has_prefix)
  if(length(candidates) > 0L){
    pattern <- "^([+-]?[0-9]+)[[:space:]]*-[[:space:]]*([+-]?[0-9]+)$"
    matches <- regmatches(
      value[candidates],
      regexec(pattern, value[candidates], perl = TRUE))
    matched <- lengths(matches) == 3L

    malformed[candidates[!matched]] <- TRUE

    if(any(matched)){
      matched_indices <- candidates[matched]
      matched_values <- matches[matched]
      start[matched_indices] <- suppressWarnings(as.numeric(
        vapply(matched_values, `[[`, character(1L), 2L)))
      end[matched_indices] <- suppressWarnings(as.numeric(
        vapply(matched_values, `[[`, character(1L), 3L)))
      status[matched_indices] <- "parsed"
    }
  }

  res <- list(
    Missing   = missing,
    Malformed = malformed,
    Status    = status,
    Start     = start,
    End       = end)
  return(res)
}

append_reason <- function(reasons, condition, reason){
  indices <- which(!is.na(condition) & condition)
  if(length(indices) == 0L){
    return(reasons)
  }

  already_set <- nzchar(reasons[indices])
  reasons[indices[!already_set]] <- reason
  reasons[indices[already_set]] <- paste0(
    reasons[indices[already_set]], ";", reason)
  return(reasons)
}

is_gzip_path <- function(filepath){
  grepl("\\.gz$", filepath, ignore.case = TRUE)
}


############################################## Load inputs

cat("Loading input sequences...\n")
seqs <- tryCatch(
  readDNAStringSet(filepath = FASTA, format = "fasta", use.names = TRUE),
  error = function(e){
    stop(
      "Could not read input FASTA file: ", FASTA,
      "\n", conditionMessage(e), "\n",
      call. = FALSE)
  })

seq_headers <- names(seqs)
if(length(seqs) > 0L &&
   (is.null(seq_headers) || anyNA(seq_headers) || any(!nzchar(seq_headers)))){
  stop("Every FASTA sequence must have a non-empty header\n", call. = FALSE)
}

## ITSx uses the first whitespace-delimited token as its sequence identifier
seq_ids <- sub(pattern = "[[:space:]].*$", replacement = "", x = seq_headers)

if(anyDuplicated(seq_ids)){
  duplicated_ids <- unique(seq_ids[duplicated(seq_ids)])
  stop(
    "FASTA contains duplicate ITSx identifier tokens. Examples: ",
    paste(head(duplicated_ids, 5L), collapse = ", "),
    "\n",
    call. = FALSE)
}

fasta_data <- data.table(
  SeqIndex    = seq_along(seqs),
  SeqID       = seq_ids,
  Header      = seq_headers,
  FastaLength = as.numeric(width(seqs)))

cat("Loading ITSx positions...\n")
positions <- read_itsx_positions(POSITIONS)

if(nrow(positions) > 0L){
  if(anyNA(positions$SeqID) || any(!nzchar(positions$SeqID))){
    stop("ITSx positions file contains an empty sequence identifier\n", call. = FALSE)
  }

  if(anyDuplicated(positions$SeqID)){
    duplicated_ids <- unique(positions$SeqID[duplicated(positions$SeqID)])
    stop(
      "ITSx positions file contains duplicate sequence identifiers. Examples: ",
      paste(head(duplicated_ids, 5L), collapse = ", "),
      "\n",
      call. = FALSE)
  }

  missing_fasta_ids <- setdiff(positions$SeqID, fasta_data$SeqID)
  if(length(missing_fasta_ids) > 0L){
    stop(
      "Could not match ", length(missing_fasta_ids),
      " ITSx position identifier(s) to the FASTA file. Examples: ",
      paste(head(missing_fasta_ids, 5L), collapse = ", "),
      "\n",
      call. = FALSE)
  }
}

positions[ , PositionRow := TRUE ]

## Preserve FASTA input order and add missing ITSx rows as non-detections
audit <- merge(
  x     = fasta_data,
  y     = positions,
  by    = "SeqID",
  all.x = TRUE,
  sort  = FALSE)
setorder(audit, SeqIndex)
audit[ is.na(PositionRow), PositionRow := FALSE ]


############################################## Validate ITSx records

audit[ , Reason := "" ]
audit[ , Region := REGION ]

audit[ , Reason := append_reason(Reason, !PositionRow, "no_positions") ]

itsx_length <- parse_itsx_length(audit$ITSxLengthRaw)
audit[ , ITSxLength := itsx_length$value ]
audit[ , Reason := append_reason(
  Reason,
  PositionRow & (!itsx_length$valid | ITSxLength <= 0),
  "invalid_itsx_length") ]
audit[ , Reason := append_reason(
  Reason,
  PositionRow & itsx_length$valid & ITSxLength != FastaLength,
  "sequence_length_mismatch") ]

is_chimeric <- grepl("chimeric", audit$Info, ignore.case = TRUE)
is_chimeric[is.na(is_chimeric)] <- FALSE
audit[ , Reason := append_reason(Reason, PositionRow & is_chimeric, "chimeric") ]

coordinate_fields <- c(
  SSU  = "SSU_Raw",
  ITS1 = "ITS1_Raw",
  S58  = "S58_Raw",
  ITS2 = "ITS2_Raw",
  LSU  = "LSU_Raw")
coordinate_labels <- c(
  SSU  = "SSU",
  ITS1 = "ITS1",
  S58  = "5.8S",
  ITS2 = "ITS2",
  LSU  = "LSU")

coordinates <- vector("list", length(coordinate_fields))
names(coordinates) <- names(coordinate_fields)
coordinate_is_valid <- vector("list", length(coordinate_fields))
names(coordinate_is_valid) <- names(coordinate_fields)

for(region_name in names(coordinate_fields)){
  parsed <- parse_region_coordinates(
    audit[[coordinate_fields[[region_name]]]],
    coordinate_labels[[region_name]])
  coordinates[[region_name]] <- parsed

  detected <- audit$PositionRow & !parsed$Missing
  finite_coordinates <- is.finite(parsed$Start) & is.finite(parsed$End)

  audit[ , Reason := append_reason(
    Reason,
    detected & parsed$Malformed,
    paste0("malformed_", region_name, "_coordinates")) ]
  audit[ , Reason := append_reason(
    Reason,
    detected & !parsed$Malformed & !finite_coordinates,
    paste0("nonfinite_", region_name, "_coordinates")) ]
  audit[ , Reason := append_reason(
    Reason,
    detected & !parsed$Malformed & finite_coordinates &
      (parsed$Start <= 0 | parsed$End <= 0),
    paste0("nonpositive_", region_name, "_coordinates")) ]
  audit[ , Reason := append_reason(
    Reason,
    detected & !parsed$Malformed & finite_coordinates &
      parsed$End <= parsed$Start,
    paste0("invalid_", region_name, "_length")) ]
  audit[ , Reason := append_reason(
    Reason,
    detected & !parsed$Malformed & finite_coordinates &
      (parsed$Start > FastaLength | parsed$End > FastaLength),
    paste0("out_of_bounds_", region_name, "_coordinates")) ]

  ## Assign one interval-local status using a deterministic precedence.
  ## Cross-region overlaps and ordering errors are recorded in `Reason` below.
  region_status <- parsed$Status
  region_status[!audit$PositionRow] <- "no_positions"

  parsed_numeric <- audit$PositionRow & parsed$Status == "parsed"
  region_status[parsed_numeric & !finite_coordinates] <- "nonfinite"
  region_status[
    parsed_numeric & finite_coordinates &
      (parsed$Start <= 0 | parsed$End <= 0)
  ] <- "nonpositive"
  region_status[
    parsed_numeric & finite_coordinates &
      parsed$Start > 0 & parsed$End > 0 &
      parsed$End <= parsed$Start
  ] <- "invalid_length"
  region_status[
    parsed_numeric & finite_coordinates &
      parsed$Start > 0 & parsed$End > parsed$Start &
      (parsed$Start > audit$FastaLength | parsed$End > audit$FastaLength)
  ] <- "out_of_bounds"
  region_status[
    parsed_numeric & finite_coordinates &
      parsed$Start > 0 & parsed$End > parsed$Start &
      parsed$Start <= audit$FastaLength & parsed$End <= audit$FastaLength
  ] <- "valid"

  region_length <- rep(NA_real_, nrow(audit))
  region_length[finite_coordinates] <-
    parsed$End[finite_coordinates] - parsed$Start[finite_coordinates] + 1

  normalized_columns <- paste0(
    region_name,
    c("_Status", "_Start", "_End", "_Length"))
  audit[ , (normalized_columns) := list(
    region_status,
    parsed$Start,
    parsed$End,
    region_length) ]

  coordinate_is_valid[[region_name]] <- region_status == "valid"
}

## Detected regions must be non-overlapping and follow their biological order
region_order <- names(coordinate_fields)
for(i in seq_len(length(region_order) - 1L)){
  for(j in seq.int(i + 1L, length(region_order))){
    left_name  <- region_order[[i]]
    right_name <- region_order[[j]]
    left  <- coordinates[[left_name]]
    right <- coordinates[[right_name]]
    both_valid <- coordinate_is_valid[[left_name]] &
      coordinate_is_valid[[right_name]]

    out_of_order <- both_valid & left$Start >= right$Start
    overlap <- both_valid & !out_of_order & left$End >= right$Start

    audit[ , Reason := append_reason(
      Reason,
      out_of_order,
      paste0("out_of_order_", left_name, "_", right_name)) ]
    audit[ , Reason := append_reason(
      Reason,
      overlap,
      paste0("overlap_", left_name, "_", right_name)) ]
  }
}


############################################## Select extraction coordinates

audit[ , `:=`(Start = NA_real_, End = NA_real_) ]

if(REGION == "SSU"){
  audit[ , Start := coordinates$SSU$Start ]
  audit[ , End   := coordinates$SSU$End ]
  audit[ , Reason := append_reason(
    Reason,
    PositionRow & coordinates$SSU$Missing,
    "SSU_not_detected") ]
}

if(REGION == "ITS"){
  audit[ , Start := coordinates$ITS1$Start ]
  audit[ , End   := coordinates$ITS2$End ]
  audit[ , Reason := append_reason(
    Reason,
    PositionRow & coordinates$ITS1$Missing,
    "ITS1_not_detected") ]
  audit[ , Reason := append_reason(
    Reason,
    PositionRow & coordinates$ITS2$Missing,
    "ITS2_not_detected") ]
}

if(REGION == "LSU"){
  audit[ , Start := coordinates$LSU$Start ]
  audit[ , End   := coordinates$LSU$End ]
  audit[ , Reason := append_reason(
    Reason,
    PositionRow & coordinates$LSU$Missing,
    "LSU_not_detected") ]
}

audit[ , Status := fifelse(nzchar(Reason), "excluded", "extracted") ]
audit[ , ExtractedLength := fifelse(
  Status == "extracted",
  End - Start + 1,
  NA_real_) ]


############################################## Extract and export sequences

extract_rows <- audit[ Status == "extracted" ]

cat("Extracting ", nrow(extract_rows), " sequence(s)...\n", sep = "")

if(nrow(extract_rows) > 0L){
  selected_seqs <- seqs[extract_rows$SeqIndex]
  result <- Biostrings::subseq(
    selected_seqs,
    start = extract_rows$Start,
    end   = extract_rows$End)

  ## Defensive assignment: retain complete input headers, including annotations
  names(result) <- names(selected_seqs)
} else {
  warning("No sequences passed validation; writing an empty FASTA file", call. = FALSE)
  result <- DNAStringSet()
}

cat("Writing extracted sequences...\n")
writeXStringSet(
  x        = result,
  filepath = OUTPUT,
  compress = is_gzip_path(OUTPUT),
  format   = "fasta",
  width    = 9999)

report_columns <- c(
  "SeqID", "Header", "FastaLength", "ITSxLength",
  "SSU_Status",  "SSU_Start",  "SSU_End",  "SSU_Length",
  "ITS1_Status", "ITS1_Start", "ITS1_End", "ITS1_Length",
  "S58_Status",  "S58_Start",  "S58_End",  "S58_Length",
  "ITS2_Status", "ITS2_Start", "ITS2_End", "ITS2_Length",
  "LSU_Status",  "LSU_Start",  "LSU_End",  "LSU_Length",
  "Info", "Region", "Start", "End", "ExtractedLength", "Status", "Reason")

report_data <- audit[ , ..report_columns ]
report_data[ Status == "extracted", Reason := NA_character_ ]

cat("Writing audit report...\n")
fwrite(
  x         = report_data,
  file      = REPORT,
  sep       = "\t",
  col.names = TRUE,
  na        = "NA",
  compress  = if(is_gzip_path(REPORT)) "gzip" else "none")


############################################## Summary

cat("\nExtraction summary:\n")
cat("Input sequences: ", nrow(audit), "\n", sep = "")
cat("Extracted: ", sum(audit$Status == "extracted"), "\n", sep = "")
cat("Excluded: ", sum(audit$Status == "excluded"), "\n", sep = "")

excluded_reasons <- audit[ Status == "excluded", Reason ]
if(length(excluded_reasons) > 0L){
  reason_counts <- sort(
    table(unlist(strsplit(excluded_reasons, split = ";", fixed = TRUE))),
    decreasing = TRUE)
  cat("\nExclusion reasons (a sequence may have more than one):\n")
  for(reason_name in names(reason_counts)){
    cat("- ", reason_name, ": ", reason_counts[[reason_name]], "\n", sep = "")
  }
}

cat("\nAll done.\n")


##################### Session info

end_time <- Sys.time()
elapsed_minutes <- as.numeric(difftime(end_time, start_time, units = "min"))
cat("\nElapsed time: ", elapsed_minutes, " minutes\n", sep = "")

cat("\nSession info:\n")
sessionInfo()
cat("\n")
