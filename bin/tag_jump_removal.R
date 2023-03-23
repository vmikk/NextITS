#!/usr/bin/env Rscript

## Script to perform tag-jump removal

# Input is given as positional arguments:
#   1. OTU table (`OTU_tab_not_filtered.txt.gz`)
#   2. f-parameter of UNCROSS (e.g., 0.01)
#   3. p-parameter (e.g., 1.0)

# Outputs:
#  - Tag-jumpfiltered OTU table (`OTU_tab_TagJumpFiltered.txt.gz`)
#  - Table with tag-jumps (`TagJump_OTUs.RData`)
#  - Plot (`TagJump_plot.pdf`)

args <- commandArgs(trailingOnly = TRUE)

suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
# library(openxlsx)

theme_set(theme_classic(base_size = 14))

## Load OTU table
cat("..Loading OTU table\n")
OTUTABW <- fread(
  file = args[1],
  sep = "\t", header = TRUE)

colnames(OTUTABW)[1] <- "OTU"

cat("...Number of OTUs: ", nrow(OTUTABW), "\n")
cat("...Number of samples: ", ncol(OTUTABW) - 1, "\n")

## Convert to long format
cat("..Converting OTU table to long format\n")
OTUTAB <- melt(data = OTUTABW, id.vars = "OTU",
  variable.name = "SampleID", value.name = "Abundance")

## Remove zero-OTUs
OTUTAB <- OTUTAB[ Abundance > 0 ]
cat("...Number of non-zero records: ", nrow(OTUTAB), "\n")


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
    f = as.numeric(args[2]),
    p = as.numeric(args[3])
    )
  )

## Truncate singletons with total OTU abundance > 99 reads
# OTUTAB[ Abundance == 1 & Total > 99  , TagJump := TRUE ]
# OTUTAB[ Abundance == 2 & Total > 999 , TagJump := TRUE ]

cat("...Number of tag-jumps: ", sum(OTUTAB$TagJump, na.rm = TRUE), "\n")


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


## TJ stats
cat("..Calculating tag-jump summary\n")
TJ <- data.table(
    Total_reads = sum(OTUTAB$Abundance),
    Number_of_TagJump_Events = sum(OTUTAB$TagJump),
    TagJump_reads = sum(OTUTAB[ TagJump == TRUE ]$Abundance, na.rm = T)
    )

TJ$ReadPercent_removed <- with(TJ, (TagJump_reads / Total_reads * 100))

fwrite(x = TJ, file = "TagJump_stats.txt", sep = "\t")


## Exporting tag-jump data
cat("..Exporting tag-jump data\n")
JMPS <- OTUTAB[ TagJump == TRUE, .(OTU, SampleID) ]

saveRDS(object = JMPS,
  file = "TagJump_OTUs.RData",
  compress = "xz")


## Prepare OTU tables, remove tag-jumps
cat("..Removing tag-jumps\n")

OTUTAB <- OTUTAB[ TagJump == FALSE ]

## Convert to wide format
RES <- dcast(data = OTUTAB,
  formula = OTU ~ SampleID,
  value.var = "Abundance", fill = 0)

## Sort rows (by total abundance)
clz <- colnames(RES)[-1]
otu_sums <- rowSums(RES[, ..clz], na.rm = TRUE)
RES <- RES[ order(otu_sums, decreasing = TRUE) ]


cat("..Exporting tag-jump filtered table\n")

fwrite(x = RES,
  file = "OTU_tab_TagJumpFiltered.txt.gz",
  sep = "\t", compress = "gzip")
