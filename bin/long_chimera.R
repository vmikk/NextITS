#!/usr/bin/env Rscript

## Check time
start_time <- Sys.time()

cat("\nParsing input options and arguments...\n")

suppressPackageStartupMessages(require(optparse))

## Parse arguments
option_list <- list(
  make_option(c("-t", "--hits"),    action = "store", default = NA,  type = "character", help = "cmsearch tabular output (--tblout)"),
  make_option(c("-c", "--cm"),      action = "store", default = NA,  type = "character", help = "Covariance model file used for the search (to get model lengths)"),
  make_option(c("-f", "--fasta"),   action = "store", default = NA,  type = "character", help = "Input sequences in FASTA format (uncompressed)"),
  make_option(c("-s", "--sample"),  action = "store", default = NA,  type = "character", help = "Sample ID (used as the output prefix)"),
  make_option(c("-b", "--minbits"), action = "store", default = 50,  type = "double",    help = "Minimum bit score of a strong hit"),
  make_option(c("-v", "--mincov"),  action = "store", default = 0.4, type = "double",    help = "Minimum fraction of the model covered by a strong hit"),
  make_option(c("-x", "--threads"), action = "store", default = 1L,  type = "integer",   help = "Number of CPU threads")
)
opt <- parse_args(OptionParser(option_list = option_list))


## Validation of the required arguments
is_missing_argument <- function(x){
  is.null(x) || length(x) != 1L || is.na(x) || !nzchar(trimws(x))
}

required_options <- c("hits", "cm", "fasta", "sample")
missing_options <- required_options[
  vapply(required_options, function(x) is_missing_argument(opt[[x]]), logical(1L))
]
if(length(missing_options) > 0L){
  stop("Missing required option(s): ", paste0("--", missing_options, collapse = ", "), "\n", call. = FALSE)
}

for(input_file in c(opt$hits, opt$cm, opt$fasta)){
  if(!file.exists(input_file)){
    stop("Input file does not exist: ", input_file, "\n", call. = FALSE)
  }
}

if(opt$minbits < 0){
  stop("--minbits must be non-negative\n", call. = FALSE)
}
if(opt$mincov <= 0 || opt$mincov > 1){
  stop("--mincov must be in (0, 1]\n", call. = FALSE)
}

## Assign variables
HITS    <- opt$hits
CMFILE  <- opt$cm
FASTA   <- opt$fasta
SAMPLE  <- opt$sample
MINBITS <- opt$minbits
MINCOV  <- opt$mincov
CPUTHREADS <- opt$threads

## Fixed constants
OPPOSITE_OVERLAP <- 0.5   # max overlap with a higher-scoring opposite-strand hit (fraction of the hit length)
IR_LEN <- 200             # min inverted repeat length (nt) of an inverted duplicate
IR_PID <- 90              # min identity (%) of the inverted repeat of an inverted duplicate

## Output files
OUT_IDS    <- "flagged_ids.txt"
OUT_FASTA  <- paste0(SAMPLE, "_LongChimera.fa.gz")
OUT_REPORT <- paste0(SAMPLE, "_LongChimera.tsv.gz")
OUT_COUNTS <- paste0(SAMPLE, "_LongChimera_counts.tsv")

## Log assigned variables
cat("\nParameters specified:\n")
cat(paste("Hits table: ",     HITS,    "\n", sep = ""))
cat(paste("Models: ",         CMFILE,  "\n", sep = ""))
cat(paste("Sequences: ",      FASTA,   "\n", sep = ""))
cat(paste("Sample ID: ",      SAMPLE,  "\n", sep = ""))
cat(paste("Min bit score: ",  MINBITS, "\n", sep = ""))
cat(paste("Min model coverage: ", MINCOV, "\n", sep = ""))
cat(paste("Number of CPU threads to use: ", CPUTHREADS, "\n", sep = ""))

cat("\n")

############################################## Load packages

load_pckg <- function(pkg = "data.table"){
  suppressPackageStartupMessages( library(package = pkg, character.only = TRUE) )
  cat(paste(pkg, packageVersion(pkg), "\n"))
}

cat("Loading packages:\n")
load_pckg("data.table")

setDTthreads(threads = CPUTHREADS)


############################################## Model lengths

## CLEN of each model (from the CM header; HMM filter sections have no CLEN)
cm_lines <- readLines(CMFILE)
cm_lines <- cm_lines[ grepl(pattern = "^(INFERNAL|HMMER|NAME|CLEN)", x = cm_lines) ]
CLEN <- list()
in_cm <- FALSE
for(ln in cm_lines){
  if(startsWith(ln, "INFERNAL")){ in_cm <- TRUE;  next }
  if(startsWith(ln, "HMMER"))   { in_cm <- FALSE; next }
  if(in_cm && startsWith(ln, "NAME")){ cur_name <- sub("^NAME\\s+", "", ln) }
  if(in_cm && startsWith(ln, "CLEN")){ CLEN[[ cur_name ]] <- as.numeric(sub("^CLEN\\s+", "", ln)) }
}
CLEN <- unlist(CLEN)

if(length(CLEN) == 0){
  stop("No model lengths (CLEN) found in the CM file\n", call. = FALSE)
}

## Short region names
region_name <- function(x){
  fcase(
    grepl("5_8S|5\\.8S", x), "5.8S",
    grepl("SSU",         x), "SSU",
    grepl("LSU",         x), "LSU",
    default = x)
}

cat("\nModel lengths:\n")
print(data.table(Model = names(CLEN), Region = region_name(names(CLEN)), CLEN = CLEN))


############################################## Output helpers

## Abundance from the sequence ID
get_size <- function(x){
  sz <- as.numeric( sub(pattern = ".*;size=([0-9]+).*", replacement = "\\1", x = x) )
  sz[ !grepl(";size=[0-9]+", x) ] <- 1
  sz
}

write_counts <- function(res){
  cnt <- data.table(
    SampleID       = SAMPLE,
    Inverted_Reads = res[ Group == "inverted", sum(Size) ],
    Inverted_Uniq  = res[ Group == "inverted", .N ],
    Repeat_Reads   = res[ Group == "repeat",   sum(Size) ],
    Repeat_Uniq    = res[ Group == "repeat",   .N ])
  fwrite(x = cnt, file = OUT_COUNTS, sep = "\t")
  cnt
}

## Exit with empty outputs (nothing to remove)
no_chimeras <- function(msg){
  cat("\n", msg, "\n", sep = "")
  file.create(OUT_IDS)
  cnt <- write_counts(data.table(Group = character(0), Size = numeric(0)))
  print(cnt)
  quit(save = "no", status = 0)
}


############################################## Load hits

cat("\nLoading cmsearch hits\n")

## Tabular output columns (1-based):
## 1 target name, 3 query name, 6 mdl from, 7 mdl to, 8 seq from, 9 seq to, 10 strand, 15 score
hits <- fread(
  cmd = paste0("awk '!/^#/ { print $1\"\\t\"$3\"\\t\"$6\"\\t\"$7\"\\t\"$8\"\\t\"$9\"\\t\"$10\"\\t\"$15 }' ", shQuote(HITS)),
  header = FALSE, sep = "\t", quote = "",
  col.names  = c("SeqID", "Model", "MdlFrom", "MdlTo", "SeqFrom", "SeqTo", "Strand", "Bits"),
  colClasses = c("character", "character", "integer", "integer", "integer", "integer", "character", "numeric"))

cat("..Number of hits: ", nrow(hits), "\n")

if(nrow(hits) == 0){
  no_chimeras("No rRNA-model hits found")
}

unknown_models <- setdiff(unique(hits$Model), names(CLEN))
if(length(unknown_models) > 0){
  stop("Hits to models absent from the CM file: ", paste(unknown_models, collapse = ", "), "\n", call. = FALSE)
}

hits[ , Region := region_name(Model) ]
hits[ , `:=`(
  Lo  = pmin(SeqFrom, SeqTo),
  Hi  = pmax(SeqFrom, SeqTo),
  Len = abs(SeqTo - SeqFrom) + 1L,
  ModelLen = CLEN[ Model ]) ]
hits[ , Coverage := (MdlTo - MdlFrom + 1) / ModelLen ]
hits[ , Strong := Bits >= MINBITS & Coverage >= MINCOV ]
hits[ , HitID := .I ]

cat("..Strong hits: ", sum(hits$Strong), "\n")


############################################## Mixed orientation

cat("\nSearching for sequences with rRNA in both orientations\n")

strong <- hits[ Strong == TRUE ]

## Strong hits overlapping a higher-scoring hit on the opposite strand are discounted
opp <- merge(
  x = strong[ , .(SeqID, HitID, Strand, Lo, Hi, Len, Bits) ],
  y = hits[ ,   .(SeqID, oStrand = Strand, oLo = Lo, oHi = Hi, oBits = Bits) ],
  by = "SeqID", allow.cartesian = TRUE)
opp <- opp[ Strand != oStrand & oBits > Bits ]
opp[ , Overlap := pmax(0L, pmin(Hi, oHi) - pmax(Lo, oLo) + 1L) / Len ]
discounted <- unique(opp[ Overlap > OPPOSITE_OVERLAP, HitID ])

valid <- strong[ ! HitID %in% discounted ]

mixed_ids <- valid[ , .(Plus = any(Strand == "+"), Minus = any(Strand == "-")), by = "SeqID"
                  ][ Plus & Minus, SeqID ]

cat("..Number of sequences with mixed orientation: ", length(mixed_ids), "\n")


############################################## Repeated regions

cat("\nSearching for repeated rRNA regions\n")

## Pairs of strong hits of the same model on the same strand
rp <- merge(
  x = strong[ , .(SeqID, Region, Strand, ModelLen, HitID, Lo, Hi, MdlFrom, MdlTo) ],
  y = strong[ , .(SeqID, Region, Strand, HitID2 = HitID, Lo2 = Lo, Hi2 = Hi, MdlFrom2 = MdlFrom, MdlTo2 = MdlTo) ],
  by = c("SeqID", "Region", "Strand"), allow.cartesian = TRUE)
rp <- rp[ HitID < HitID2 ]

## Non-overlapping in the sequence, sharing a substantial part of the model
rp[ , SeqOverlap := pmin(Hi, Hi2) - pmax(Lo, Lo2) + 1L ]
rp[ , SharedModel := pmax(0L, pmin(MdlTo, MdlTo2) - pmax(MdlFrom, MdlFrom2) + 1L) / ModelLen ]
rp <- rp[ SeqOverlap <= 0 & SharedModel >= MINCOV ]

repeats <- rp[ , .(Repeated = paste(sort(unique(Region)), collapse = ",")), by = "SeqID" ]

cat("..Number of sequences with repeated regions: ", nrow(repeats), "\n")


############################################## Combine

flagged <- unique(c(mixed_ids, repeats$SeqID))

if(length(flagged) == 0){
  no_chimeras("No long chimeras found")
}

res <- data.table(SeqID = flagged)
res[ , Group := fifelse(SeqID %chin% mixed_ids, "inverted", "repeat") ]
res[ , Size := get_size(SeqID) ]

## Best scores on each strand
best <- hits[ , .(
    BestPlusBits  = max(c(0, Bits[ Strand == "+" ])),
    BestMinusBits = max(c(0, Bits[ Strand == "-" ])) ),
  by = "SeqID" ]
res <- merge(res, best, by = "SeqID", all.x = TRUE)

## Hit layout (all hits above the reporting threshold, in sequence order)
layout <- hits[ SeqID %chin% flagged ][ order(SeqID, Lo) ][ ,
  .(Layout = paste0(Region, Strand, "[", MdlFrom, "-", MdlTo, "]:", Lo, "-", Hi, "(", round(Bits), ")", collapse = " ")),
  by = "SeqID" ]
res <- merge(res, layout, by = "SeqID", all.x = TRUE)
res <- merge(res, repeats, by = "SeqID", all.x = TRUE)


############################################## Classify the inverted sequences

## SSU and LSU strong on one strand only, 5.8S strong only on the other strand, between them
its_block <- function(h){
  for(sf in c("+", "-")){
    f <- h[ Strand == sf ]
    o <- h[ Strand != sf ]
    if(!all(c("SSU", "LSU") %chin% f$Region) || "5.8S" %chin% f$Region){ next }
    if(!identical(unique(o$Region), "5.8S")){ next }
    lo <- min(f[ Region %chin% c("SSU", "LSU"), Hi ])
    hi <- max(f[ Region %chin% c("SSU", "LSU"), Lo ])
    if(any(o$Lo > lo & o$Hi < hi)){ return(TRUE) }
  }
  FALSE
}

res[ , `:=`(InvRepeatLen = NA_integer_, InvRepeatIdentity = NA_real_, Block = FALSE) ]

## Load only the flagged sequences (via FASTA index)
cat("\nLoading flagged sequences\n")
load_pckg("Biostrings")
load_pckg("pwalign")

fai <- fasta.index(FASTA, seqtype = "DNA")
fai <- fai[ sub(pattern = "\\s.*", replacement = "", x = fai$desc) %chin% flagged, ]
seqs <- readDNAStringSet(fai)
names(seqs) <- sub(pattern = "\\s.*", replacement = "", x = names(seqs))

if(!setequal(names(seqs), flagged)){
  stop("Some flagged sequences are absent from the FASTA file\n", call. = FALSE)
}

inv_ids <- res[ Group == "inverted", SeqID ]
if(length(inv_ids) > 0){

  cat("..Estimating inverted repeats for ", length(inv_ids), " sequences\n")

  ## Inverted repeat = best local alignment of the sequence with its reverse complement
  smat <- nucleotideSubstitutionMatrix(match = 2, mismatch = -3, baseOnly = FALSE)
  ir_fun <- function(id){
    s <- seqs[[ id ]]
    a <- pairwiseAlignment(reverseComplement(s), s, type = "local",
      substitutionMatrix = smat, gapOpening = 5, gapExtension = 2)
    data.table(SeqID = id,
      InvRepeatLen = nchar(pattern(a)),
      InvRepeatIdentity = round(pid(a, type = "PID1"), 1))
  }
  ## NB. a local alignment of a 4-5 kb sequence takes ~0.3 GB of RAM, so at most 2 alignments run in parallel
  ir <- rbindlist(parallel::mclapply(inv_ids, ir_fun, mc.cores = max(1L, min(CPUTHREADS, 2L))))

  blk <- valid[ SeqID %chin% inv_ids, .(Block = its_block(.SD)), by = "SeqID" ]

  res[ ir, on = "SeqID", `:=`(InvRepeatLen = i.InvRepeatLen, InvRepeatIdentity = i.InvRepeatIdentity) ]
  res[ blk, on = "SeqID", Block := i.Block ]
}

res[ , Class := fcase(
  Group == "inverted" & InvRepeatLen >= IR_LEN & InvRepeatIdentity >= IR_PID, "inverted_duplicate",
  Group == "inverted" & Block %in% TRUE,                                       "inverted_ITS_block",
  Group == "inverted",                                                         "inverted_chimera",
  default = paste0("repeated_", gsub(",", "_", Repeated)) ) ]

cat("\nSummary of detected long chimeras:\n")
print( res[ , .(NumUniq = .N, NumReads = sum(Size)), by = .(Group, Class) ][ order(Group, Class) ] )


############################################## Export

cat("\nExporting results\n")

## IDs to remove
writeLines(text = res$SeqID, con = OUT_IDS)

## Removed sequences, annotated
seqs <- seqs[ res$SeqID ]
names(seqs) <- paste0(res$SeqID, ";sample=", SAMPLE, ";class=", res$Class, ";")
writeXStringSet(x = seqs, filepath = OUT_FASTA, compress = TRUE, format = "fasta", width = 99999)

## Per-sequence evidence
res[ , SampleID := SAMPLE ]
setorder(res, Group, Class, -Size, SeqID)
fwrite(
  x = res[ , .(SampleID, SeqID, Size, Group, Class,
    BestPlusBits, BestMinusBits, InvRepeatLen, InvRepeatIdentity, Layout) ],
  file = OUT_REPORT, sep = "\t", na = "NA", compress = "gzip")

## Per-sample counts
cnt <- write_counts(res)
print(cnt)


cat("\nAll done.\n")
cat("Elapsed time: ", format(round(difftime(Sys.time(), start_time, units = "secs"), 1)), "\n")
