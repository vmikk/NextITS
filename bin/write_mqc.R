#!/usr/bin/env Rscript

## Prepare MultiQC-ready custom-content files from the NextITS summary tables
##
## These are a side output, not the primary report: they let groups that already run MultiQC across a project fold NextITS results in.
## MultiQC itself is not a pipeline dependency - this script only writes the files, which any stock MultiQC picks up from the `*_mqc.*` filename alone
##
## See https://docs.seqera.io/multiqc/custom_content

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

option_list <- list(
  make_option("--per-sample", type = "character", default = NULL, help = "per_sample.tsv (Step-1)"),
  make_option("--per-run",    type = "character", default = NULL, help = "per_run.tsv (Step-1)"),
  make_option("--otu-summary", type = "character", default = NULL, help = "Step2_summary.tsv"),
  make_option("--versions",   type = "character", default = NULL, help = "software_versions.yml"),
  make_option("--outdir",     type = "character", default = ".", help = "Output directory [%default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1L]
script_dir <- if (!is.na(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg), mustWork = FALSE)) else getwd()
lib <- c("report_lib.R", file.path(script_dir, "report_lib.R"))
lib <- lib[file.exists(lib)][1L]
if (is.na(lib)) stop("Could not find report_lib.R", call. = FALSE)
source(lib)

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
out <- function(f) file.path(opt$outdir, f)

## MultiQC reads a leading block of `#`-prefixed lines as one YAML document
write_mqc <- function(file, header, dt) {
  con <- file(out(file), open = "wt")
  on.exit(close(con))
  writeLines(paste0("# ", header), con)
  write.table(dt, con, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  cat("Wrote ", out(file), "\n", sep = "")
}

written <- 0L

if (is_usable_file(opt$`per-sample`)) {
  ps <- fread(opt$`per-sample`, sep = "\t", na.strings = c("NA", ""), showProgress = FALSE)
  id <- pick_col(ps, "Sample", "file", "SampleID")
  if (!is.null(id)) setnames(ps, id, "Sample")

  ## Long chimeras are present only if the `long_chimera` filter was enabled
  has_longchim <- all(c("LongChimera_Inverted_Reads", "LongChimera_Repeat_Reads") %in% names(ps))
  lc_inv <- num0(col_or(ps, 0, "LongChimera_Inverted_Reads"))
  lc_rep <- num0(col_or(ps, 0, "LongChimera_Repeat_Reads"))
  primer_ok <- num0(col_or(ps, 0, "PrimerChecked_Reads"))

  ## General Statistics: one row per sample, headers as an ordered list
  gs <- data.table(
    Sample     = ps$Sample,
    Demux      = num0(col_or(ps, 0, "Demultiplexed_Reads")),
    Artefacts  = round(num0(col_or(ps, 0, "PrimerArtefacts_Percent")), 2),
    LongChim   = round(ifelse(primer_ok > 0, (lc_inv + lc_rep) / primer_ok * 100, 0), 2),
    ITSx_yield = round(num0(col_or(ps, 0, "ITSx_Yield_Percent")), 2),
    Retained   = round(num0(col_or(ps, 0, "Percentage_Reads_Retained")), 2),
    Final      = num0(col_or(ps, 0, "SeqTable_NumReads")))
  if (!has_longchim) gs[, LongChim := NULL]
  write_mqc("nextits_step1_general_stats_mqc.tsv", c(
    "id: 'nextits_step1_genstats'",
    "plot_type: 'generalstats'",
    "pconfig:",
    "    - Sample: { title: 'Sample' }",
    "    - Demux: { title: 'Demultiplexed', description: 'Demultiplexed reads', format: '{:,.0f}', min: 0, scale: 'Blues' }",
    "    - Artefacts: { title: 'Artefacts', description: 'Multiprimer artefacts, % of primer-screened reads', suffix: '%', min: 0, max: 100, scale: 'OrRd' }",
    if (has_longchim) "    - LongChim: { title: 'Long chimeras', description: 'Long chimeras (inverted and repeat), % of primer-checked reads', suffix: '%', min: 0, max: 100, scale: 'OrRd' }",
    "    - ITSx_yield: { title: 'ITSx yield', description: 'ITSx-extracted reads, % of primer-checked reads', suffix: '%', min: 0, max: 100, scale: 'RdYlGn' }",
    "    - Retained: { title: 'Retained', description: 'Reads in the sequence table, % of demultiplexed reads', suffix: '%', min: 0, max: 100, scale: 'RdYlGn' }",
    "    - Final: { title: 'Final reads', description: 'Reads in the sequence table', format: '{:,.0f}', min: 0, scale: 'Greens' }"
  ), gs)
  written <- written + 1L

  ## Read fate as a stacked bar, mirroring the main report's per-sample panel
  demux     <- num0(col_or(ps, 0, "Demultiplexed_Reads"))
  artefacts <- num0(col_or(ps, 0, "PrimerArtefacts_Reads"))
  itsx_ok   <- if ("ITSx_Extracted_Reads" %in% names(ps)) num0(ps$ITSx_Extracted_Reads) else pmax(primer_ok - lc_inv - lc_rep, 0)
  retained  <- num0(col_or(ps, 0, "SeqTable_NumReads"))
  fate <- data.table(
    Sample              = ps$Sample,
    Retained            = retained,
    `Primer artefacts`  = artefacts,
    `Long chimera (inverted)` = lc_inv,
    `Long chimera (repeat)`   = lc_rep,
    `No ITSx detection` = pmax(primer_ok - lc_inv - lc_rep - itsx_ok, 0),
    `Other filtering`   = pmax(itsx_ok - retained, 0))
  if (!has_longchim) fate[, c("Long chimera (inverted)", "Long chimera (repeat)") := NULL]
  write_mqc("nextits_step1_read_fate_mqc.tsv", c(
    "id: 'nextits_step1_read_fate'",
    "section_name: 'NextITS: read fate'",
    "description: 'What happened to each demultiplexed read during Step-1 processing.'",
    "plot_type: 'bargraph'",
    "pconfig:",
    "    id: 'nextits_step1_read_fate_plot'",
    "    title: 'NextITS: read fate by sample'",
    "    ylab: 'Reads'",
    "    cpswitch_counts_label: 'Reads'"
  ), fate)
  written <- written + 1L
}

if (is_usable_file(opt$`per-run`)) {
  pr <- fread(opt$`per-run`, sep = "\t", na.strings = c("NA", ""), showProgress = FALSE)
  if (nrow(pr) > 0L) {
    keep <- intersect(c("Total_Number_Of_Reads", "Reads_Passed_QC", "Reads_Demultiplexed",
                        "Reads_PrimerChecked", "Reads_ITSx_Extracted", "SeqTable_NumReads"),
                      names(pr))
    if (length(keep) > 0L) {
      write_mqc("nextits_step1_run_funnel_mqc.tsv", c(
        "id: 'nextits_step1_run_funnel'",
        "section_name: 'NextITS: run funnel'",
        "description: 'Reads surviving each Step-1 stage, summed over the whole run.'",
        "plot_type: 'bargraph'",
        "pconfig:",
        "    id: 'nextits_step1_run_funnel_plot'",
        "    title: 'NextITS: reads per stage'",
        "    ylab: 'Reads'",
        "    cpswitch: false"
      ), cbind(Stage = "Whole run", pr[1L, ..keep]))
      written <- written + 1L
    }
  }
}

if (is_usable_file(opt$`otu-summary`)) {
  s2 <- fread(opt$`otu-summary`, sep = "\t", showProgress = FALSE)
  id <- pick_col(s2, "SampleID", "Sample")
  if (!is.null(id)) setnames(s2, id, "Sample")
  write_mqc("nextits_step2_general_stats_mqc.tsv", c(
    "id: 'nextits_step2_genstats'",
    "plot_type: 'generalstats'",
    "pconfig:",
    "    - Sample: { title: 'Sample' }",
    "    - Reads: { title: 'OTU reads', description: 'Reads in the pooled OTU table', format: '{:,.0f}', min: 0, scale: 'Blues' }",
    "    - OTUs: { title: 'OTUs', description: 'OTUs observed in the sample', format: '{:,.0f}', min: 0, scale: 'Greens' }",
    "    - Singleton_OTUs: { title: 'Singletons', description: 'OTUs represented by a single read across the whole run', format: '{:,.0f}', min: 0, scale: 'OrRd' }"
  ), s2[, .(Sample, Reads, OTUs, Singleton_OTUs)])
  written <- written + 1L
}

## Versions matching a MultiQC module name are shown beside that module
if (is_usable_file(opt$versions)) {
  v <- read_versions(opt$versions)
  if (!is.null(v) && requireNamespace("yaml", quietly = TRUE)) {
    yaml::write_yaml(v, out("nextits_software_mqc_versions.yml"))
    cat("Wrote ", out("nextits_software_mqc_versions.yml"), "\n", sep = "")
    written <- written + 1L
  }
}

if (written == 0L) {
  cat("No usable inputs; nothing written.\n")
}
