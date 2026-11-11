#!/usr/bin/env Rscript

## Self-contained HTML report for the NextITS `seqstats` sub-workflow (quick QC, demultiplexing and primer check)
## Reads per_sample.tsv / per_run.tsv from quick_stats.R plus the optional demultiplexing and seqkit files
## Panels drop out when their input is missing or was replaced by a `no_*` placeholder

start_time <- Sys.time()

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

option_list <- list(
  make_option("--per-sample",   type = "character", default = NULL, help = "per_sample.tsv (from quick_stats.R, required)"),
  make_option("--per-run",      type = "character", default = NULL, help = "per_run.tsv (from quick_stats.R, required)"),
  make_option("--params",       type = "character", default = NULL, help = "pipeline_params.tsv"),
  make_option("--versions",     type = "character", default = NULL, help = "software_versions.yml"),
  make_option("--command",      type = "character", default = NULL, help = "execution_command.txt"),
  make_option("--logo",         type = "character", default = NULL, help = "NextITS logo (SVG or PNG)"),
  make_option("--trace-stamp",  type = "character", default = NULL, help = "Timestamp Nextflow stamped into pipeline_info filenames"),
  make_option("--lima-summary", type = "character", default = NULL, help = "lima.lima.summary (PacBio)"),
  make_option("--lima-counts",  type = "character", default = NULL, help = "lima.lima.counts (PacBio)"),
  make_option("--illumina-demux",  type = "character", default = NULL, help = "Demux_totals.tsv (Illumina, cutadapt)"),
  make_option("--illumina-qcheck", type = "character", default = NULL, help = "Quality_check.tsv (Illumina)"),
  make_option("--counts-dir",   type = "character", default = NULL, help = "Directory with Counts_*.txt"),
  make_option("--assets",       type = "character", default = "assets", help = "Report asset directory [%default]"),
  make_option("--out",          type = "character", default = "Seqstats_report.html", help = "Output HTML [%default]"),
  make_option("--min-reads",        type = "double", default = 100, help = "Warn below this many primer-checked reads per sample [%default]"),
  make_option("--max-artefact-pct", type = "double", default = 20,  help = "Warn above this percentage of multiprimer artefacts per sample [%default]"),
  make_option("--min-demux-pct",    type = "double", default = 50,  help = "Warn if less than this percentage of raw reads (Illumina: read pairs) was demultiplexed [%default]"),
  make_option("--min-primer-pct",   type = "double", default = 40,  help = "Warn if less than this percentage of raw reads (Illumina: read pairs) passed the primer check [%default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

## Nextflow stages report_lib.R next to this script; a local run finds it beside the file.
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1L]
script_dir <- if (!is.na(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg), mustWork = FALSE)) else getwd()
lib <- c("report_lib.R", file.path(script_dir, "report_lib.R"))
lib <- lib[file.exists(lib)][1L]
if (is.na(lib)) stop("Could not find report_lib.R", call. = FALSE)
source(lib)

cat("Loading seqstats report inputs\n")

if (!is_usable_file(opt$`per-sample`) || !is_usable_file(opt$`per-run`)) {
  stop("--per-sample and --per-run are required", call. = FALSE)
}

per_sample <- fread(opt$`per-sample`, sep = "\t", na.strings = c("NA", ""), showProgress = FALSE)
per_run    <- fread(opt$`per-run`,    sep = "\t", na.strings = c("NA", ""), showProgress = FALSE)

if (!"Sample" %in% names(per_sample)) {
  id_col <- pick_col(per_sample, "Sample", "file", "SampleID")
  if (is.null(id_col)) stop("per_sample.tsv has no sample identifier column", call. = FALSE)
  setnames(per_sample, id_col, "Sample")
}
per_sample[, Sample := as.character(Sample)]
per_sample[, RunID := run_of(Sample)]

params   <- read_params_tsv(opt$params)
versions <- read_versions(opt$versions)

lima        <- parse_lima_summary(opt$`lima-summary`)
lima_counts <- read_lima_counts(opt$`lima-counts`)

## Illumina paired-end reads are demultiplexed, reoriented and merged before QC;
## PacBio reads pass QC before demultiplexing.
## quick_stats.R writes the reorientation columns only for Illumina data
is_illumina <- if ("Reoriented_Pairs" %in% names(per_sample) || "Pairs_Reoriented" %in% names(per_run)) {
  TRUE
} else if (!is.null(lima) || !is.null(lima_counts)) {
  FALSE
} else {
  identical(getp(params, "seqplatform", "PacBio"), "Illumina")
}

n_samples <- nrow(per_sample)
runs      <- sort(unique(stats::na.omit(per_sample$RunID)))
run       <- if (nrow(per_run) > 0L) per_run[1L] else data.table()

rv <- function(nm, default = NA_real_) {
  if (nm %in% names(run)) suppressWarnings(as.numeric(run[[nm]][[1L]])) else default
}

## ------------------------------------------------------- per-run figures

reads_raw   <- rv("Total_Number_Of_Reads")
reads_qc    <- rv("Reads_Passed_QC")
reads_demux <- rv("Reads_Demultiplexed", sum(num0(per_sample$Demultiplexed_Reads)))
reads_prim  <- rv("Reads_PrimerChecked", sum(num0(per_sample$PrimerChecked_Reads)))
reads_artef <- rv("PrimerArtefacts_Reads", sum(num0(per_sample$PrimerArtefacts_Reads)))

## Illumina: read pairs after reorientation, merged and joined reads
reads_reo    <- rv("Pairs_Reoriented", sum(num0(col_or(per_sample, 0, "Reoriented_Pairs"))))
reads_merged <- rv("Reads_Merged",     sum(num0(col_or(per_sample, 0, "Merged_Reads"))))
reads_joined <- rv("Reads_Joined",     sum(num0(col_or(per_sample, 0, "Joined_Reads"))))
pct_merged   <- rv("Percentage_Merged", safe_pct(reads_merged, reads_reo))

pct_demux <- rv("Percentage_Demultiplexed", safe_pct(reads_demux, reads_raw))
pct_prim  <- rv("Percentage_PrimerChecked", safe_pct(reads_prim, reads_raw))
pct_artef <- rv("PrimerArtefacts_Percent",  safe_pct(reads_artef, reads_artef + reads_prim))

## PacBio QC acts on raw reads; Illumina QC acts on merged reads
pct_qc <- if (is_illumina) safe_pct(reads_qc, reads_merged) else rv("Percentage_QC", safe_pct(reads_qc, reads_raw))

prim_per_sample <- num0(per_sample$PrimerChecked_Reads)
med_reads <- if (n_samples > 0L) stats::median(prim_per_sample) else NA_real_

## ------------------------------------------------------------ KPI cards

unit <- if (is_illumina) "pairs" else "reads"

## Card grades: below the threshold is "Poor", below the midpoint between the threshold and 100% is "Fair"
grade_min <- function(x, thr) grade(x, warn = thr + (100 - thr) / 2, bad = thr)

kpis <- kpi_row(
  kpi_card(paste("Raw", if (is_illumina) "read pairs" else "reads"), fmt_int(reads_raw)),
  if (is_illumina) kpi_card("Merged pairs", fmt_pct(pct_merged),
           paste0(fmt_int(reads_merged), " of ", fmt_int(reads_reo), " reoriented pairs"),
           grade(pct_merged, warn = 80, bad = 50)),
  kpi_card("Passed QC", fmt_pct(pct_qc),
           if (is_illumina) paste0(fmt_int(reads_qc), " of ", fmt_int(reads_merged), " merged reads")
           else paste0(fmt_int(reads_qc), " of ", fmt_int(reads_raw), " raw reads")),
  kpi_card("Demultiplexed", fmt_pct(pct_demux),
           paste0(fmt_int(reads_demux), " of ", fmt_int(reads_raw), " raw ", unit),
           grade_min(pct_demux, opt$`min-demux-pct`)),
  kpi_card("Primer-checked", fmt_pct(pct_prim),
           paste0(fmt_int(reads_prim), " reads · ", fmt_pct(safe_pct(reads_prim, reads_demux)), " of demultiplexed"),
           grade_min(pct_prim, opt$`min-primer-pct`)),
  kpi_card("Samples", fmt_int(n_samples),
           if (length(runs) > 0L) paste0(length(runs), " sequencing run(s)") else NULL),
  kpi_card("Median reads per sample", fmt_int(med_reads),
           if (n_samples > 0L) paste0(fmt_int(min(prim_per_sample)), " – ", fmt_int(max(prim_per_sample)),
                                      " primer-checked reads") else NULL,
           grade(med_reads, warn = 5 * opt$`min-reads`, bad = opt$`min-reads`)),
  kpi_card("Multiprimer artefacts", fmt_pct(pct_artef),
           paste0(fmt_int(reads_artef), " reads of ", fmt_int(reads_artef + reads_prim), " primer-screened"),
           grade(pct_artef, warn = opt$`max-artefact-pct` / 2, bad = opt$`max-artefact-pct`, higher_better = FALSE))
)

## ------------------------------------------------------------ QC verdicts

flag <- function(dt, cond, name, detail_fmt, threshold, status_bad = "warn") {
  hit <- dt[cond]
  if (nrow(hit) == 0L) {
    return(qc_rule(name, "ok", "no samples affected", threshold))
  }
  qc_rule(name, status_bad, sprintf(detail_fmt, nrow(hit)), threshold, sort(hit$Sample))
}

run_flag <- function(x, thr, name, what) {
  thr_txt <- paste0("< ", fmt_pct(thr, 0), " of raw ", unit)
  if (!is.finite(x)) return(qc_rule(name, "warn", paste0("the ", what, " rate could not be estimated"), thr_txt))
  if (x < thr) qc_rule(name, "warn", paste0(fmt_pct(x), " of raw ", unit, " ", what), thr_txt)
  else qc_rule(name, "ok", paste0(fmt_pct(x), " of raw ", unit, " ", what), thr_txt)
}

rules <- list(
  flag(per_sample, num0(per_sample$PrimerChecked_Reads) < opt$`min-reads`,
       "Library depth", "%d sample(s) below the depth floor",
       paste0("< ", fmt_int(opt$`min-reads`), " primer-checked reads")),
  flag(per_sample, num0(per_sample$PrimerArtefacts_Percent) >= opt$`max-artefact-pct`,
       "Multiprimer artefacts", "%d sample(s) above the artefact ceiling",
       paste0("≥ ", fmt_pct(opt$`max-artefact-pct`, 0))),
  run_flag(pct_demux, opt$`min-demux-pct`, "Demultiplexing rate", "demultiplexed"),
  run_flag(pct_prim,  opt$`min-primer-pct`, "Primer-check rate", "primer-checked")
)

## ------------------------------------------------------------- run sankey

## Trunk nodes share one colour so the eye follows the surviving reads;
## each loss branch keeps its own hue from FATE_COL
TRUNK <- "var(--fate-trunk)"
NODE_COL <- c(FATE_COL,
  "Raw reads" = TRUNK, "Raw pairs" = TRUNK, "Passed QC" = TRUNK,
  "Demultiplexed" = TRUNK, "Reoriented" = TRUNK, "Merged" = TRUNK, "Joined" = TRUNK,
  "Primer-checked" = TRUNK,
  "Not reoriented" = "var(--fate-noprimers)")

lk <- data.table(source = character(), target = character(), value = numeric())
add_link <- function(from, to, value) {
  if (is.finite(value) && value > 0) lk <<- rbind(lk, data.table(source = from, target = to, value = value))
}

## A stage cannot pass on more than it received; clamping keeps every loss branch non-negative
## even if the counts are inconsistent. NA means the stage count is unknown and the stage is skipped
take <- function(x, avail) {
  if (!is.finite(x)) return(NA_real_)
  if (!is.finite(avail)) return(max(x, 0))
  max(min(x, avail), 0)
}

## Primer screening: artefacts, primer-checked, and the rest (no primers)
primer_links <- function(from, avail) {
  prim  <- take(reads_prim, avail)
  artef <- take(reads_artef, if (is.finite(avail)) avail - num0(prim) else NA_real_)
  add_link(from, "Primer-checked", prim)
  add_link(from, "Primer artefacts", artef)
  if (is.finite(avail)) add_link(from, "No primers", avail - num0(prim) - num0(artef))
}

if (is_illumina) {
  ## Raw pairs -> demultiplexed -> reoriented -> merged (+ joined) -> passed QC -> primer-checked
  node  <- "Raw pairs"
  avail <- reads_raw
  demux <- take(reads_demux, avail)
  if (is.finite(demux)) {
    if (is.finite(avail)) add_link(node, "Undemultiplexed", avail - demux)
    add_link(node, "Demultiplexed", demux)
    node <- "Demultiplexed"; avail <- demux
  }
  reo <- take(reads_reo, avail)
  if (is.finite(reo)) {
    if (is.finite(avail)) add_link(node, "Not reoriented", avail - reo)
    add_link(node, "Reoriented", reo)
    node <- "Reoriented"; avail <- reo
  }
  merged <- take(reads_merged, avail)
  joined <- take(reads_joined, if (is.finite(avail)) avail - num0(merged) else NA_real_)
  add_link(node, "Merged", merged)
  add_link(node, "Joined", joined)
  if (is.finite(avail)) add_link(node, "Not merged", avail - num0(merged) - num0(joined))
  qc <- take(reads_qc, merged)
  add_link("Merged", "Passed QC", qc)
  if (is.finite(merged)) add_link("Merged", "Failed QC", merged - num0(qc))
  ## Joined pairs skip the merged-read QC
  add_link("Joined", "Passed QC", joined)
  primer_links("Passed QC", num0(qc) + num0(joined))
} else {
  ## Raw reads -> passed QC -> demultiplexed -> primer-checked
  node  <- "Raw reads"
  avail <- reads_raw
  qc <- take(reads_qc, avail)
  if (is.finite(qc)) {
    if (is.finite(avail)) add_link(node, "Failed QC", avail - qc)
    add_link(node, "Passed QC", qc)
    node <- "Passed QC"; avail <- qc
  }
  demux <- take(reads_demux, avail)
  if (is.finite(avail)) add_link(node, "Undemultiplexed", avail - num0(demux))
  add_link(node, "Demultiplexed", demux)
  primer_links("Demultiplexed", demux)
}

sankey_fig <- if (nrow(lk) > 0L) {
  nodes <- unique(c(lk$source, lk$target))
  echart(sankey_option(nodes, lk, colours = NODE_COL),
         height = 360, title = "Read fate across the run",
         caption = paste0(
           "Hover a flow or node for exact read counts. Widths are proportional to reads. ",
           if (is_illumina) "Illumina read pairs are demultiplexed, reoriented and merged before quality filtering."
           else "PacBio reads are quality-filtered before demultiplexing."))
}

## ------------------------------------------------------------ stage table

stage_tbl <- if (is_illumina) {
  rbindlist(list(
    data.table(Stage = "Demultiplexing (pairs)", In = reads_raw,   Out = reads_demux),
    data.table(Stage = "Reorientation (pairs)",  In = reads_demux, Out = reads_reo),
    data.table(Stage = "Read merging",           In = reads_reo,   Out = reads_merged + reads_joined),
    data.table(Stage = "Quality filtering",      In = reads_merged + reads_joined, Out = reads_qc + reads_joined),
    data.table(Stage = "Primer screening",       In = reads_qc + reads_joined,     Out = reads_prim)))
} else {
  rbindlist(list(
    data.table(Stage = "Quality filtering", In = reads_raw, Out = reads_qc),
    data.table(Stage = "Demultiplexing",    In = if (is.finite(reads_qc)) reads_qc else reads_raw, Out = reads_demux),
    data.table(Stage = "Primer screening",  In = reads_demux, Out = reads_prim)))
}
stage_tbl <- stage_tbl[is.finite(In) & is.finite(Out)]
stage_tbl[, Lost := pmax(In - Out, 0)]
stage_tbl[, `Lost %` := safe_pct(Lost, In)]
stage_tbl[, `Cumulative %` := safe_pct(Out, reads_raw)]

## ---------------------------------------------------------- per-sample fate

## Every demultiplexed read (pair) ends up in exactly one category
fate <- data.table(
  Sample           = per_sample$Sample,
  demux            = num0(per_sample$Demultiplexed_Reads),
  `Primer-checked` = num0(per_sample$PrimerChecked_Reads),
  `Primer artefacts` = num0(per_sample$PrimerArtefacts_Reads))

if (is_illumina) {
  reo    <- num0(col_or(per_sample, 0, "Reoriented_Pairs"))
  merged <- num0(col_or(per_sample, 0, "Merged_Reads"))
  joined <- num0(col_or(per_sample, 0, "Joined_Reads"))
  qc_ok  <- num0(col_or(per_sample, 0, "QC_Reads"))
  fate[, `Not merged` := pmax(reo - merged - joined, 0)]
  fate[, `Failed QC`  := pmax(merged - qc_ok, 0)]
  fate[, `No primers` := pmax(demux - reo, 0) + pmax(qc_ok + joined - `Primer artefacts` - `Primer-checked`, 0)]
} else {
  fate[, `No primers` := pmax(demux - `Primer artefacts` - `Primer-checked`, 0)]
}

SAMPLE_COL <- c("Primer-checked" = FATE_COL[["Retained"]], FATE_COL)
FATE_CATS <- c("Primer-checked", "Primer artefacts", "No primers",
               if (is_illumina) c("Not merged", "Failed QC"))
FATE_CATS <- FATE_CATS[vapply(FATE_CATS, function(k) sum(fate[[k]]) > 0 || k == "Primer-checked", logical(1))]

setorder(fate, -demux, Sample)

fate_fig <- if (n_samples > 0L) {
  fate_id <- new_chart_id()
  fate_variants <- list(
    counts  = stacked_bar_option(fate, FATE_CATS, colours = SAMPLE_COL, percent = FALSE,
                                 y_name = if (is_illumina) "Reads or read pairs" else "Reads"),
    percent = stacked_bar_option(fate, FATE_CATS, colours = SAMPLE_COL, percent = TRUE))
  echart(fate_variants$counts, id = fate_id, height = 400,
    title = "Read fate by sample",
    controls = switch_buttons(fate_id, c("Counts", "Percent"), c("counts", "percent")),
    variants = fate_variants,
    caption = paste0(
      "Samples are ordered by demultiplexed yield. ",
      if (is_illumina) "Read pairs that failed reorientation are counted as having no primers. " else "",
      "Switch to Percent to compare rates; a shallow sample can look like 100% artefacts while contributing almost no data."))
}

depth_fig <- if (n_samples > 1L) {
  bins <- as.integer(max(5L, min(40L, 2L * ceiling(sqrt(n_samples)))))
  o <- hist_option(prim_per_sample, bins = bins, x_name = "Primer-checked reads per sample",
                   y_name = "Samples", colour = "var(--viz-dist)",
                   marklines = bin_markline(prim_per_sample, bins, opt$`min-reads`,
                                            paste0("min reads = ", fmt_int(opt$`min-reads`))))
  if (!is.null(o)) echart(o, height = 300, title = "Reads per sample",
    caption = paste0("Distribution of primer-checked reads across samples. The median is ",
                     fmt_int(med_reads), " reads; the depth floor is ", fmt_int(opt$`min-reads`), " reads."))
}

sample_cols <- intersect(c(
  "Sample", "RunID", "Demultiplexed_Reads",
  "Reoriented_Pairs", "Merged_Reads", "NotMerged_Pairs", "Merged_Percent", "Joined_Reads", "QC_Reads",
  "PrimerChecked_Reads", "PrimerArtefacts_Reads", "PrimerArtefacts_Percent"),
  names(per_sample))

sample_view <- copy(per_sample)[, ..sample_cols]
setorderv(sample_view, "Demultiplexed_Reads", order = -1L, na.last = TRUE)
if (length(runs) < 2L && "RunID" %in% names(sample_view)) sample_view[, RunID := NULL]

sample_labels <- c(
  Sample = "Sample", RunID = "Run",
  Demultiplexed_Reads = if (is_illumina) "Demultiplexed pairs" else "Demultiplexed",
  Reoriented_Pairs = "Reoriented pairs", Merged_Reads = "Merged", NotMerged_Pairs = "Not merged pairs",
  Merged_Percent = "Merged, %", Joined_Reads = "Joined", QC_Reads = "Passed QC",
  PrimerChecked_Reads = "Primer-checked",
  PrimerArtefacts_Reads = "Artefact reads", PrimerArtefacts_Percent = "Artefacts")

## ------------------------------------------------------- demultiplexing

demux_bits <- if (is_illumina) {
  demux_section_illumina(read_illumina_demux(opt$`illumina-demux`),
                         read_illumina_qcheck(opt$`illumina-qcheck`), per_sample)
} else {
  demux_section_lima(lima, lima_counts)
}

## ----------------------------------------------------------- read lengths

counts_dir <- opt$`counts-dir`
len <- stage_lengths(counts_dir, list(
  c("Counts_1.RawData.txt",         "Raw data"),
  c("Counts_2.QC.txt",              if (is_illumina) "Merged, passed QC" else "Passed QC"),
  if (!is_illumina) c("Counts_3.Demux.txt", "Demultiplexed"),
  if (is_illumina)  c("Counts_3.Joined.txt", "Joined"),
  c("Counts_4.PrimerCheck.txt",     "Primer-checked"),
  c("Counts_4.PrimerArtefacts.txt", "Primer artefacts")))
len_bits <- stage_length_table(len)

## ------------------------------------------------------------- assemble

platform <- if (is_illumina) "Illumina" else "PacBio"

meta <- c(
  "NextITS"  = nextits_version_label(versions),
  "Nextflow" = version_label(versions, "Nextflow"),
  "Platform" = platform,
  "Samples"  = fmt_int(n_samples))
if (length(runs) > 0L) meta["Runs"] <- if (length(runs) <= 3L) paste(runs, collapse = ", ") else fmt_int(length(runs))
meta["Generated"] <- format(Sys.time(), "%Y-%m-%d %H:%M %Z")

sections <- list(
  sec("overview", "Run overview",
      kpis,
      qc_panel(rules),
      sankey_fig,
      subhead("overview-stages", "Stage-by-stage losses"),
      dt_table(stage_tbl,
        fmt = list(In = "int", Out = "int", Lost = "int", `Lost %` = "pct", `Cumulative %` = "pct"),
        bar_cols = "Out", search = FALSE, cols_menu = FALSE,
        download = "seqstats_stage_losses.tsv"),
      pipeline_links(stamp = opt$`trace-stamp`),
      subs = list(c("overview-stages", "Stage losses"))),

  sec("samples", "Per-sample",
      fate_fig, depth_fig,
      subhead("samples-table", "Per-sample counts"),
      dt_table(sample_view,
        labels = unname(ifelse(names(sample_view) %in% names(sample_labels),
                               sample_labels[names(sample_view)],
                               gsub("_", " ", names(sample_view)))),
        fmt = list(Merged_Percent = "pct", PrimerArtefacts_Percent = "pct"),
        bar_cols = "PrimerChecked_Reads",
        download = "seqstats_per_sample.tsv"),
      subs = list(c("samples-table", "Counts table"))),

  if (!is.null(demux_bits)) sec("demux",
      if (is_illumina) "Demultiplexing and read merging" else "Demultiplexing",
      demux_bits,
      note = if (is_illumina)
               "Illumina paired-end reads: tag assignment with cutadapt, read merging with USEARCH or VSEARCH."
             else "PacBio lima barcode assignment."),

  if (!is.null(len_bits)) sec("lengths", "Read length", len_bits),

  sec("settings", "Run settings",
      command_panel(opt$command, extra = c(
        "Sequencing platform" = getp(params, "seqplatform", platform),
        "Forward primer"      = getp(params, "primer_forward", "—"),
        "Reverse primer"      = getp(params, "primer_reverse", "—"))),
      p(class = "hint", HTML(paste0(
        "Every resolved parameter, including defaults, is listed in ",
        '<a href="../pipeline_info/pipeline_params.tsv">pipeline_params.tsv</a>.')))),

  sec("versions", "Software versions",
      versions_panel(versions))
)

report_page(
  title = "NextITS seqstats report",
  meta = meta,
  sections = sections,
  assets_dir = opt$assets,
  out = opt$out,
  logo = opt$logo)

cat("Wrote ", normalizePath(opt$out, mustWork = FALSE), "\n", sep = "")
cat("Elapsed minutes: ", round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 3), "\n", sep = "")
