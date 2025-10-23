#!/usr/bin/env Rscript

## Script to perform tag-jump removal

# Input arguments:
#   1. Sequence table in long format, no header (`Seq_tab_not_filtered.txt.gz`),
#       with columns: `SeqID`, `Abundance`, `SampleID`
#   2. Pre-clustered membership table (`TJPreclust.uc.parquet`)
#   2. f-parameter of UNCROSS (e.g., 0.01)
#   3. p-parameter (e.g., 1.0)

# Outputs:
#  - Tag-jump-filtered sequence table (`Seq_tab_TagJumpFiltered.txt.gz`)
#  - Table with tag-jump scores (`TagJump_scores.qs`)
#  - Plot (`TagJump_plot.pdf`)

cat("\nParsing input options and arguments...\n")

suppressPackageStartupMessages(require(optparse))

## Parse arguments
option_list <- list(
  make_option(c("-s", "--seqtab"), action="store", default="seqtab.txt.gz", type='character', help="Sequence table in long format"),
  make_option(c("-c", "--precls"), action="store", default="precls.txt.gz", type='character', help="Table with pre-clustered sequence membership"),
  make_option(c("-f", "--uncross_f"), action="store", default=0.01, type='numeric', help="f-parameter of UNCROSS"),
  make_option(c("-p", "--uncross_p"), action="store", default=1,    type='numeric', help="Additional p-parameter for UNCROSS")
)
opt <- parse_args(OptionParser(option_list=option_list))

## Function to convert text "NA"s to NA
to_na <- function(x){
  if(x %in% c("NA", "null", "Null")){ x <- NA }
  return(x)
}

## Replaces "null"s from Nextflow with NA
opt <- lapply(X = opt, FUN = to_na)

## Validation of the required arguments
if(is.na(opt$seqtab)){
  stop("Input file with sequence table is not specified\n")
}
if(is.na(opt$precls)){
  stop("Input file with pre-clustered membership table is not specified\n")
}

## Set default params if not specified
if(is.na(opt$uncross_f) | is.null(opt$uncross_f) | is.nan(opt$uncross_f) | !is.numeric(opt$uncross_f)){
  cat("f-parameter is not specified, using default value: 0.01\n")
  opt$uncross_f <- 0.01
}
if(is.na(opt$uncross_p) | is.null(opt$uncross_p) | is.nan(opt$uncross_p) | !is.numeric(opt$uncross_p)){
  cat("p-parameter is not specified, using default value: 1\n")
  opt$uncross_p <- 1
}

## Assign variables
SEQTAB <- opt$seqtab
PRECLS <- opt$precls
F      <- opt$uncross_f
P      <- opt$uncross_p

## Log assigned variables
cat("\nParameters specified:\n")
cat(paste("Sequence table: " , SEQTAB, "\n", sep = ""))
cat(paste("Pre-clustered membership table: ", PRECLS, "\n", sep = ""))
cat(paste("f-parameter of UNCROSS: ", F, "\n", sep = ""))
cat(paste("p-parameter of UNCROSS: ", P, "\n", sep = ""))

cat("\n")


############################################## Load packages

cat("Loading R packages...\n")

load_pckg <- function(pkg = "data.table"){
  suppressPackageStartupMessages( library(package = pkg, character.only = TRUE) )
  cat(paste(pkg, packageVersion(pkg), "\n"))
}

load_pckg("data.table")
load_pckg("arrow")
load_pckg("ggplot2")
# load_pckg("qs")

theme_set(theme_classic(base_size = 14))

############################################## Workflow

## Load sequence table
cat("..Loading sequence table\n")
SEQTAB <- fread(file = SEQTAB,
  sep = "\t", header = FALSE,
  col.names = c("SeqID", "SampleID", "Abundance"))

## Load sequence membership table
cat("..Loading sequence membership table\n")
PRECLS <- read_parquet(file = PRECLS)
setDT(PRECLS)
setnames(PRECLS, new = c("SeqID", "OTU"))

## Add cluster membership to the sequence table
SEQTAB <- merge(x = SEQTAB, y = PRECLS, by = "SeqID", all.x = TRUE)

if(any(is.na(SEQTAB$OTU))){
  cat("WARNING: Sequences without cluster membership detected\n")
  cat("WARNING: Excluding these records from the analysis\n")
  SEQTAB <- SEQTAB[ !is.na(OTU) ]
}

cat("...Number of unique sequences: ", length(unique(SEQTAB$SeqID)),    "\n")
cat("...Number of clusters: ",         length(unique(SEQTAB$OTU)),      "\n")
cat("...Number of samples: ",          length(unique(SEQTAB$SampleID)), "\n")

## Summarize by sequence clusters
cat("..Summarizing by sequence clusters\n")
OTUTAB <- SEQTAB[ , .(Abundance = sum(Abundance, na.rm = TRUE)), by = c("OTU", "SampleID") ]

## Estimate total abundance of sequence per plate
cat("..Estimating total OTU abundance\n")
OTUTAB[ , Total := sum(Abundance, na.rm = TRUE), by = "OTU" ]

## UNCROSS score (with original parameter - take a root from the exp in denominator, to make curves more steep)
uncross_score <- function(x, N, n, f = 0.01, tmin = 0.1, p = 1){
  # x = OTU abundance in a sample
  # N = total OTU abundance
  # n = number of samples
  # f = expected cross-talk rate, e.g. 0.01
  # tmin = min score to be considered as cross-talk
  # p = power to rise the exponent (default, 1; use 1/2 or 1/3 to make cureves more stepp)

  z <- f * N / n               # Expected treshold
  sc <- 2 / (1 + exp(x/z)^p)   # t-score
  res <- data.table(Score = sc, TagJump = sc >= tmin)
  return(res)
}

## Esimate UNCROSS score
cat("..Estimating UNCROSS score\n")
OTUTAB <- cbind(
  OTUTAB,
  uncross_score(
    x = OTUTAB$Abundance,
    N = OTUTAB$Total,
    n = length(unique(OTUTAB$SampleID)),
    f = as.numeric(F),
    p = as.numeric(P)
    )
  )

## Truncate singletons with total OTU abundance > 99 reads
# OTUTAB[ Abundance == 1 & Total > 99  , TagJump := TRUE ]
# OTUTAB[ Abundance == 2 & Total > 999 , TagJump := TRUE ]

cat("...Number of tag-jumps: ", sum(OTUTAB$TagJump, na.rm = TRUE), "\n")

## Export tag-jump scores
setcolorder(OTUTAB,
  c("OTU", "SampleID", "TagJump", "Score", "Abundance", "Total"))

setorder(OTUTAB, OTU, -Abundance, SampleID)

cat("..Exporting tag-jump scores\n")
qs::qsave(OTUTAB,
  "TagJump_scores.qs",
  preset = "custom", algorithm = "zstd", compress_level = 5L, nthreads = 1L)


## Plot
cat("..Making a plot\n")
PP <- ggplot(data = OTUTAB, aes(x = Total, y = Abundance, color = TagJump)) +
    geom_point() + scale_x_log10() + scale_y_log10() +
    scale_color_manual(values = c("#0C7C59", "#D64933")) +
    labs(x = "Total abundance of OTU, reads", y = "Abundance of OTU in a sample, reads")

cat("..Exporting a plot\n")
pdf(file = "TagJump_plot.pdf", width = 12, height = 9.5, useDingbats = FALSE)
  PP
dev.off()



## Exporting tag-jump data
# cat("..Exporting tag-jump data\n")
# JMPS <- OTUTAB[ TagJump == TRUE, .(OTU, SampleID) ]
# 
# saveRDS(object = JMPS,
#   file = "TagJump_OTUs.RData",
#   compress = "xz")


## Prepare filtered sequence table, remove tag-jumps
cat("..Removing tag-jumps\n")

## Add tag-jump info to the sequence table
n1 <- nrow(SEQTAB)

RES <- merge(
  x = SEQTAB,
  y = OTUTAB[ , .(OTU, SampleID, TagJump) ],
  by = c("OTU", "SampleID"), all.x = TRUE)

n2 <- nrow(RES)
if(n1 != n2){
  cat("WARNING: merging went wrong likely\n")
  cat("WARNING: There might be duplicated sequences\n")
}


## TJ stats
cat("..Calculating tag-jump summary\n")
TJ <- data.table(
    Total_reads              = sum(RES$Abundance),
    Number_of_TagJump_Events = sum(RES$TagJump),
    TagJump_reads = sum(RES[ TagJump == TRUE ]$Abundance, na.rm = T)
    )

TJ$ReadPercent_removed <- with(TJ, (TagJump_reads / Total_reads * 100))

fwrite(x = TJ, file = "TagJump_stats.txt", sep = "\t")


## Keep only non-tag-jump reads
RES <- RES[ TagJump == FALSE , .(SampleID, SeqID, Abundance) ]
setorder(RES, SampleID, -Abundance)

cat("..Exporting tag-jump filtered table\n")

fwrite(x = RES,
  file = "Seq_tab_TagJumpFiltered.txt.gz",
  sep = "\t", compress = "gzip")
