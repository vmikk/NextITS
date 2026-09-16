#!/usr/bin/env Rscript

## Self-contained HTML run report for NextITS Step-1.
## Reads the summary tables from read_count_summary.R plus whatever optional
## diagnostic files the run produced; every panel drops out cleanly when its
## input is missing or was replaced by a `no_*` sentinel.

start_time <- Sys.time()

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

option_list <- list(
  make_option("--per-sample",   type = "character", default = NULL, help = "per_sample.tsv"),
  make_option("--per-run",      type = "character", default = NULL, help = "per_run.tsv"),
  make_option("--params",       type = "character", default = NULL, help = "pipeline_params.tsv"),
  make_option("--versions",     type = "character", default = NULL, help = "software_versions.yml"),
  make_option("--methods",      type = "character", default = NULL, help = "README_Step1_Methods.txt"),
  make_option("--command",      type = "character", default = NULL, help = "execution_command.txt"),
  make_option("--logo",         type = "character", default = NULL, help = "NextITS logo (SVG or PNG)"),
  make_option("--trace-stamp",  type = "character", default = NULL, help = "Timestamp Nextflow stamped into pipeline_info filenames"),
  make_option("--lima-summary", type = "character", default = NULL, help = "lima.lima.summary"),
  make_option("--lima-counts",  type = "character", default = NULL, help = "lima.lima.counts"),
  make_option("--counts-dir",   type = "character", default = NULL, help = "Directory with Counts_*.txt"),
  make_option("--itsx-dir",     type = "character", default = NULL, help = "Directory with ITSx *.summary.txt / *.positions.txt / *.problematic.txt"),
  make_option("--tagjump-scores", type = "character", default = NULL, help = "TagJump_scores.qs"),
  make_option("--tagjump-stats",  type = "character", default = NULL, help = "TagJump_stats.txt"),
  make_option("--seqtab",       type = "character", default = NULL, help = "Seqs.parquet"),
  make_option("--seq-quals",    type = "character", default = NULL, help = "SeqQualities.parquet"),
  make_option("--assets",       type = "character", default = "assets", help = "Report asset directory [%default]"),
  make_option("--out",          type = "character", default = "Step1_report.html", help = "Output HTML [%default]"),
  make_option("--min-reads",        type = "double", default = 100,  help = "Warn below this many demultiplexed reads [%default]"),
  make_option("--max-artefact-pct", type = "double", default = 20,   help = "Warn above this multiprimer artefact %% [%default]"),
  make_option("--min-itsx-pct",     type = "double", default = 50,   help = "Warn below this ITSx yield %% [%default]"),
  make_option("--min-retained-pct", type = "double", default = 30,   help = "Warn below this retained-read %% [%default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

## Nextflow stages report_lib.R next to this script; a local run finds it beside the file.
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1L]
script_dir <- if (!is.na(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg), mustWork = FALSE)) else getwd()
lib <- c("report_lib.R", file.path(script_dir, "report_lib.R"))
lib <- lib[file.exists(lib)][1L]
if (is.na(lib)) stop("Could not find report_lib.R", call. = FALSE)
source(lib)

cat("Loading Step-1 report inputs\n")

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
methods  <- read_methods_text(opt$methods)

n_samples <- nrow(per_sample)
runs      <- sort(unique(stats::na.omit(per_sample$RunID)))
run       <- if (nrow(per_run) > 0L) per_run[1L] else data.table()

rv <- function(nm, default = NA_real_) {
  if (nm %in% names(run)) suppressWarnings(as.numeric(run[[nm]][[1L]])) else default
}

has_itsx <- "ITSx_Extracted_Reads" %in% names(per_sample)

## ------------------------------------------------------- per-run figures

reads_raw   <- rv("Total_Number_Of_Reads")
reads_qc    <- rv("Reads_Passed_QC")
reads_demux <- rv("Reads_Demultiplexed", sum(num0(per_sample$Demultiplexed_Reads)))
reads_prim  <- rv("Reads_PrimerChecked", sum(num0(per_sample$PrimerChecked_Reads)))
reads_artef <- rv("PrimerArtefacts_Reads", sum(num0(per_sample$PrimerArtefacts_Reads)))
reads_itsx  <- rv("Reads_ITSx_Extracted", if (has_itsx) sum(num0(per_sample$ITSx_Extracted_Reads)) else NA_real_)
reads_final <- rv("SeqTable_NumReads", sum(num0(per_sample$SeqTable_NumReads)))

pct_artef <- rv("PrimerArtefacts_Percent", safe_pct(reads_artef, reads_artef + reads_prim))
pct_itsx  <- rv("ITSx_Yield_Percent", safe_pct(reads_itsx, reads_prim))

## read_count_summary.R computes this against raw reads per run but against
## demultiplexed reads per sample. Report both, each labelled.
pct_keep_raw   <- safe_pct(reads_final, reads_raw)
pct_keep_demux <- safe_pct(reads_final, reads_demux)

## ------------------------------------------------------------ KPI cards

kpis <- kpi_row(
  kpi_card("Demultiplexed reads", fmt_int(reads_demux),
           paste0(fmt_pct(safe_pct(reads_demux, reads_raw)), " of ", fmt_int(reads_raw), " raw reads")),
  kpi_card("Multiprimer artefacts", fmt_pct(pct_artef),
           paste0(fmt_int(reads_artef), " reads of ", fmt_int(reads_artef + reads_prim), " primer-screened"),
           grade(pct_artef, warn = opt$`max-artefact-pct` / 2, bad = opt$`max-artefact-pct`, higher_better = FALSE)),
  if (has_itsx) kpi_card("ITSx yield", fmt_pct(pct_itsx),
           paste0(fmt_int(reads_itsx), " of ", fmt_int(reads_prim), " primer-checked reads"),
           grade(pct_itsx, warn = 80, bad = opt$`min-itsx-pct`)),
  kpi_card("Reads retained", fmt_pct(pct_keep_raw),
           paste0(fmt_int(reads_final), " reads · ", fmt_pct(pct_keep_demux), " of demultiplexed"),
           grade(pct_keep_raw, warn = 60, bad = opt$`min-retained-pct`))
)

## ------------------------------------------------------------ QC verdicts

flag <- function(dt, cond, name, detail_fmt, threshold, status_bad = "warn") {
  hit <- dt[cond]
  if (nrow(hit) == 0L) {
    return(qc_rule(name, "ok", "no samples affected", threshold))
  }
  qc_rule(name, status_bad, sprintf(detail_fmt, nrow(hit)), threshold, sort(hit$Sample))
}

rules <- list(
  flag(per_sample, num0(per_sample$Demultiplexed_Reads) < opt$`min-reads`,
       "Library depth", "%d sample(s) below the depth floor",
       paste0("< ", fmt_int(opt$`min-reads`), " reads")),
  flag(per_sample, num0(per_sample$PrimerArtefacts_Percent) >= opt$`max-artefact-pct`,
       "Multiprimer artefacts", "%d sample(s) above the artefact ceiling",
       paste0("≥ ", fmt_pct(opt$`max-artefact-pct`, 0))),
  if (has_itsx) flag(per_sample, num0(per_sample$ITSx_Yield_Percent) < opt$`min-itsx-pct`,
       "ITSx yield", "%d sample(s) below the yield floor",
       paste0("< ", fmt_pct(opt$`min-itsx-pct`, 0))),
  flag(per_sample, num0(per_sample$Percentage_Reads_Retained) < opt$`min-retained-pct`,
       "Reads retained", "%d sample(s) below the retention floor",
       paste0("< ", fmt_pct(opt$`min-retained-pct`, 0), " of demultiplexed"))
)

## --------------------------------------------------------- read-fate model

fate <- data.table(
  Sample     = per_sample$Sample,
  demux      = num0(per_sample$Demultiplexed_Reads),
  artefacts  = num0(per_sample$PrimerArtefacts_Reads),
  primer_ok  = num0(per_sample$PrimerChecked_Reads),
  retained   = num0(per_sample$SeqTable_NumReads)
)
fate[, itsx_ok := if (has_itsx) num0(per_sample$ITSx_Extracted_Reads) else primer_ok]
fate[, no_itsx := pmax(primer_ok - itsx_ok, 0)]

chim_denovo <- num0(col_or(per_sample, 0, "DeNovoChimeras_NumReads"))
chim_ref    <- num0(col_or(per_sample, 0, "ReferenceBasedChimera_Reads"))
recovered   <- num0(col_or(per_sample, 0, "Recovered_ReferenceBasedChimea_Reads",
                                          "Recovered_ReferenceBasedChimera_Reads"))
fate[, chimeric := pmax(chim_denovo + chim_ref - recovered, 0)]
fate[, tagjump  := num0(col_or(per_sample, 0, "TagJump_Reads"))]

## Losses must not exceed what entered the stage, and must not double-count.
fate[, budget := pmax(itsx_ok - retained, 0)]
fate[, chimeric := pmin(chimeric, budget)]
fate[, tagjump  := pmin(tagjump, pmax(budget - chimeric, 0))]
fate[, other    := pmax(budget - chimeric - tagjump, 0)]
fate[, budget := NULL]

setnames(fate,
  c("artefacts", "no_itsx", "chimeric", "tagjump", "other", "retained"),
  c("Primer artefacts", "No ITSx detection", "Chimeric", "Tag jumps", "Other filtering", "Retained"))

FATE_CATS <- c("Retained", "Primer artefacts", "No ITSx detection",
               "Chimeric", "Tag jumps", "Other filtering")
FATE_CATS <- FATE_CATS[vapply(FATE_CATS, function(k) sum(fate[[k]]) > 0 || k == "Retained", logical(1))]

setorder(fate, -demux, Sample)

## ------------------------------------------------------------- run sankey

## Trunk nodes share one colour so the eye follows the surviving reads
## each loss branch keeps its own hue from FATE_COL
NODE_COL <- c(FATE_COL,
  "Raw reads"      = "#2e7d4f", "Passed QC"      = "#2e7d4f",
  "Demultiplexed"  = "#2e7d4f", "Primer-checked" = "#2e7d4f",
  "ITSx extracted" = "#2e7d4f")

lk <- data.table(source = character(), target = character(), value = numeric())
add_link <- function(from, to, value) {
  if (is.finite(value) && value > 0) lk <<- rbind(lk, data.table(source = from, target = to, value = value))
}

entry <- "Raw reads"
if (is.finite(reads_raw) && is.finite(reads_qc)) {
  add_link("Raw reads", "Passed QC", reads_qc)
  add_link("Raw reads", "Failed QC", reads_raw - reads_qc)
  entry <- "Passed QC"
} else if (is.finite(reads_raw)) {
  entry <- "Raw reads"
}
add_link(entry, "Demultiplexed", reads_demux)
add_link(entry, "Undemultiplexed", (if (is.finite(reads_qc)) reads_qc else reads_raw) - reads_demux)
add_link("Demultiplexed", "Primer artefacts", reads_artef)
add_link("Demultiplexed", "Primer-checked", reads_prim)

last <- "Primer-checked"
if (has_itsx && is.finite(reads_itsx)) {
  add_link("Primer-checked", "ITSx extracted", reads_itsx)
  add_link("Primer-checked", "No ITSx detection", reads_prim - reads_itsx)
  last <- "ITSx extracted"
}
for (k in setdiff(FATE_CATS, c("Retained", "Primer artefacts", "No ITSx detection"))) {
  add_link(last, k, sum(fate[[k]]))
}
add_link(last, "Retained", reads_final)

sankey_fig <- if (nrow(lk) > 0L) {
  nodes <- unique(c(lk$source, lk$target))
  echart(sankey_option(nodes, lk, colours = NODE_COL),
         height = 380, title = "Read fate across the run",
         caption = "Hover a flow or node for exact read counts. Widths are proportional to reads.")
}

## ------------------------------------------------------------ stage table

stage_tbl <- rbindlist(list(
  data.table(Stage = "Quality filtering", In = reads_raw,   Out = reads_qc),
  data.table(Stage = "Demultiplexing",    In = reads_qc,    Out = reads_demux),
  data.table(Stage = "Primer screening",  In = reads_demux, Out = reads_prim),
  if (has_itsx) data.table(Stage = "ITSx extraction", In = reads_prim, Out = reads_itsx),
  data.table(Stage = "Chimera / tag-jump filtering",
             In = if (has_itsx) reads_itsx else reads_prim, Out = reads_final)
), fill = TRUE)
stage_tbl <- stage_tbl[is.finite(In) & is.finite(Out)]
stage_tbl[, Lost := pmax(In - Out, 0)]
stage_tbl[, `Lost %` := safe_pct(Lost, In)]
stage_tbl[, `Cumulative %` := safe_pct(Out, reads_raw)]

## -------------------------------------------------------- per-sample plots

fate_id <- new_chart_id()
fate_variants <- list(
  counts  = stacked_bar_option(fate, FATE_CATS, percent = FALSE),
  percent = stacked_bar_option(fate, FATE_CATS, percent = TRUE))

fate_fig <- echart(fate_variants$counts, id = fate_id, height = 400,
  title = "Read fate by sample",
  controls = switch_buttons(fate_id, c("Counts", "Percent"), c("counts", "percent")),
  variants = fate_variants,
  caption = paste0(
    "Samples are ordered by demultiplexed yield. Switch to Percent to compare rates; ",
    "a shallow sample can look like 100% artefacts while contributing almost no data."))

scatter_fig <- {
  keep_pct <- num0(per_sample$Percentage_Reads_Retained)[match(fate$Sample, per_sample$Sample)]
  pts <- lapply(seq_len(nrow(fate)), function(i) list(
    x = max(fate$demux[[i]], 1), y = keep_pct[[i]],
    name = fate$Sample[[i]], size = max(fate$Retained[[i]], 1)))
  sopt <- scatter_option(pts,
    "Demultiplexed reads", "Reads retained, % of demultiplexed",
    x_log = TRUE, size_name = "Final reads",
    marklines = markline_y(opt$`min-retained-pct`, "retention floor"))
  echart(sopt, height = 360, title = "Depth versus retention",
    caption = "The x-axis is log-scaled and point area scales with the final read count. Samples in the lower-left are both shallow and lossy.")
}

spread_fig <- {
  groups <- list()
  groups[["Artefacts, %"]] <- num0(per_sample$PrimerArtefacts_Percent)
  if (has_itsx) groups[["ITSx yield, %"]] <- num0(per_sample$ITSx_Yield_Percent)
  groups[["Retained, %"]] <- num0(per_sample$Percentage_Reads_Retained)
  o <- box_option(groups, y_name = "Percent")
  if (!is.null(o)) echart(o, height = 300, title = "Spread of per-sample rates",
    caption = "Box shows the quartiles across samples; whiskers reach the minimum and maximum.")
}

sample_cols <- intersect(c(
  "Sample", "RunID", "Demultiplexed_Reads", "PrimerChecked_Reads",
  "PrimerArtefacts_Reads", "PrimerArtefacts_Percent",
  "ITSx_Extracted_Reads", "ITSx_Yield_Percent",
  "DeNovoChimeras_NumReads", "TagJump_Reads", "TagJump_Events",
  "SeqTable_NumReads", "SeqTable_NumUniqSeqs", "Percentage_Reads_Retained"),
  names(per_sample))

sample_view <- copy(per_sample)[, ..sample_cols]
setorderv(sample_view, "Demultiplexed_Reads", order = -1L, na.last = TRUE)
if (length(runs) < 2L && "RunID" %in% names(sample_view)) sample_view[, RunID := NULL]

sample_fmt <- list(
  PrimerArtefacts_Percent   = "pct",
  ITSx_Yield_Percent        = "pct",
  Percentage_Reads_Retained = "pct")

sample_labels <- c(
  Sample = "Sample", RunID = "Run",
  Demultiplexed_Reads = "Demultiplexed", PrimerChecked_Reads = "Primer-checked",
  PrimerArtefacts_Reads = "Artefact reads", PrimerArtefacts_Percent = "Artefacts",
  ITSx_Extracted_Reads = "ITSx reads", ITSx_Yield_Percent = "ITSx yield",
  DeNovoChimeras_NumReads = "Chimeric reads", TagJump_Reads = "Tag-jump reads",
  TagJump_Events = "Tag-jump events", SeqTable_NumReads = "Final reads",
  SeqTable_NumUniqSeqs = "Unique seqs", Percentage_Reads_Retained = "Retained")

## --------------------------------------------------------- demultiplexing

parse_lima_summary <- function(path) {
  if (!is_usable_file(path)) return(NULL)
  lines <- readLines(path, warn = FALSE)
  kv <- function(pat) {
    hit <- grep(pat, lines, value = TRUE, ignore.case = TRUE)[1L]
    if (is.na(hit)) return(NA_real_)
    suppressWarnings(as.numeric(sub("^[^:]*:\\s*([0-9.]+).*$", "\\1", hit)))
  }
  ## The "ZMW marginals for (C)" block explains why reads were dropped.
  start <- grep("ZMW marginals", lines)[1L]
  marg <- NULL
  if (!is.na(start)) {
    blk <- lines[seq(start + 1L, length(lines))]
    blk <- blk[seq_len(max(0L, which(!nzchar(trimws(blk)))[1L] - 1L))]
    hit <- regmatches(blk, regexec("^\\s*(.+?)\\s*:\\s*([0-9]+)", blk))
    hit <- Filter(function(m) length(m) == 3L, hit)
    if (length(hit) > 0L) {
      marg <- data.table(
        Reason = vapply(hit, function(m) trimws(m[[2L]]), ""),
        ZMWs   = as.numeric(vapply(hit, function(m) m[[3L]], "")))
    }
  }
  list(input = kv("ZMWs input"), pass = kv("above all thresholds"),
       fail = kv("below any threshold"), marginals = marg)
}

lima <- parse_lima_summary(opt$`lima-summary`)

lima_counts <- if (is_usable_file(opt$`lima-counts`)) {
  x <- tryCatch(fread(opt$`lima-counts`, sep = "\t", showProgress = FALSE), error = function(e) NULL)
  if (!is.null(x) && "Counts" %in% names(x)) x else NULL
} else NULL

demux_section <- {
  bits <- list()
  if (!is.null(lima)) {
    bits <- c(bits, list(kpi_row(
      kpi_card("ZMWs input", fmt_int(lima$input)),
      kpi_card("Above all thresholds", fmt_int(lima$pass),
               fmt_pct(safe_pct(lima$pass, lima$input))),
      kpi_card("Below any threshold", fmt_int(lima$fail),
               fmt_pct(safe_pct(lima$fail, lima$input)),
               grade(safe_pct(lima$fail, lima$input), warn = 10, bad = 25, higher_better = FALSE)))))
    if (!is.null(lima$marginals) && nrow(lima$marginals) > 0L) {
      m <- lima$marginals[ZMWs > 0][order(-ZMWs)]
      if (nrow(m) > 0L) {
        o <- bar_option(m$Reason, m$ZMWs, y_name = "ZMWs", colour = "#d98324")
        o$xAxis$axisLabel$rotate <- 30
        o$grid$bottom <- 30
        bits <- c(bits, list(
          subhead("demux-reasons", "Why reads were discarded"),
          echart(o, height = 300,
            caption = "A ZMW can fail several thresholds at once, so these bars overlap and do not sum to the total.")))
      }
    }
  }
  if (!is.null(lima_counts)) {
    nm <- pick_col(lima_counts, "IdxFirstNamed", "IdxCombinedNamed")
    cnt <- lima_counts[order(-Counts)]
    o <- bar_option(cnt[[nm]], cnt$Counts, y_name = "Reads")
    mean_reads <- mean(cnt$Counts)
    o$series[[1L]]$markLine <- markline_y(round(mean_reads, 1), "mean")
    bits <- c(bits, list(
      subhead("demux-barcodes", "Reads per barcode"),
      echart(o, height = 320,
        caption = "Even barcode representation is expected; a strong skew points at unbalanced pooling or a failing barcode.")))
    if ("MeanScore" %in% names(lima_counts)) {
      bits <- c(bits, list(dt_table(
        cnt[, .(Barcode = get(nm), Reads = Counts, `Mean score` = MeanScore)],
        fmt = list(Reads = "int", `Mean score` = "num"),
        bar_cols = "Reads", download = "lima_barcode_counts.tsv", cols_menu = FALSE)))
    }
  }
  if (length(bits) == 0L) NULL else bits
}

## ------------------------------------------------------------------ ITSx

itsx_dir <- opt$`itsx-dir`
itsx_files <- if (!is.null(itsx_dir) && dir.exists(itsx_dir)) {
  list(summary     = list.files(itsx_dir, pattern = "\\.summary\\.txt$",     full.names = TRUE),
       positions   = list.files(itsx_dir, pattern = "\\.positions\\.txt$",   full.names = TRUE),
       problematic = list.files(itsx_dir, pattern = "\\.problematic\\.txt$", full.names = TRUE))
} else list(summary = character(), positions = character(), problematic = character())

parse_itsx_summary <- function(path) {
  lines <- readLines(path, warn = FALSE)
  get1 <- function(pat) {
    hit <- grep(pat, lines, value = TRUE)[1L]
    if (is.na(hit)) return(NA_real_)
    suppressWarnings(as.numeric(trimws(sub("^.*:\\s*", "", hit))))
  }
  start <- grep("by preliminary origin", lines)[1L]
  origins <- NULL
  if (!is.na(start)) {
    blk <- lines[seq(start + 1L, length(lines))]
    stop_at <- which(grepl("^-{5,}", blk))[1L]
    if (!is.na(stop_at)) blk <- blk[seq_len(stop_at - 1L)]
    hit <- regmatches(blk, regexec("^\\s+(.+?):\\s*([0-9]+)\\s*$", blk))
    hit <- Filter(function(m) length(m) == 3L, hit)
    if (length(hit) > 0L) {
      origins <- data.table(
        Origin = vapply(hit, function(m) trimws(m[[2L]]), ""),
        N      = as.numeric(vapply(hit, function(m) m[[3L]], "")))
    }
  }
  list(sample = sub("\\.summary\\.txt$", "", basename(path)),
       input = get1("Number of sequences in input file"),
       detected = get1("Sequences detected as ITS by ITSx"),
       chimeric = get1("Sequences detected as chimeric by ITSx"),
       origins = origins)
}

itsx <- if (length(itsx_files$summary) > 0L) lapply(itsx_files$summary, parse_itsx_summary) else list()

itsx_section <- if (length(itsx) > 0L) {
  det <- rbindlist(lapply(itsx, function(s) data.table(
    Sample = s$sample, Input = s$input, Detected = s$detected, Chimeric = s$chimeric)), fill = TRUE)
  det[, `Not detected` := pmax(num0(Input) - num0(Detected) - num0(Chimeric), 0)]
  det[, `Detection rate` := safe_pct(Detected, Input)]
  setorder(det, -Input)

  origins <- rbindlist(lapply(itsx, function(s) {
    if (is.null(s$origins)) return(NULL)
    cbind(Sample = s$sample, s$origins)
  }), fill = TRUE)

  det_mat <- det[, .(Sample, Detected = num0(Detected), `Not detected`, Chimeric = num0(Chimeric))]
  det_cats <- c("Detected", "Not detected", "Chimeric")
  det_cats <- det_cats[vapply(det_cats, function(k) sum(det_mat[[k]]) > 0, logical(1))]

  bits <- list(
    callout(HTML(paste0(
      "ITSx runs on dereplicated sequences, so the counts in this section are ",
      "<b>unique sequences</b>, not reads."))),
    subhead("itsx-detection", "Detection rate"),
    echart(stacked_bar_option(det_mat, det_cats,
        colours = c(Detected = "#2a9d6e", `Not detected` = "#9b5de5", Chimeric = "#c9184a"),
        y_name = "Unique sequences"),
      height = 340, title = "ITSx detection per sample"))

  if (!is.null(origins) && nrow(origins) > 0L) {
    keep <- origins[, .(N = sum(N)), by = Origin][N > 0][order(-N)]$Origin
    if (length(keep) > 0L) {
      wide <- dcast(origins[Origin %in% keep], Sample ~ Origin, value.var = "N", fill = 0)
      setorderv(wide, keep[[1L]], order = -1L)
      pal <- setNames(c("#2a9d6e", PAL[-3L], rep(PAL, 3L))[seq_along(keep)], keep)
      if ("Fungi" %in% keep) pal[["Fungi"]] <- "#2a9d6e"
      oid <- new_chart_id()
      vars <- list(counts  = stacked_bar_option(wide, keep, colours = pal, percent = FALSE, y_name = "Unique sequences"),
                   percent = stacked_bar_option(wide, keep, colours = pal, percent = TRUE,  y_name = "Unique sequences"))
      bits <- c(bits, list(
        subhead("itsx-origin", "Preliminary taxonomic origin"),
        echart(vars$counts, id = oid, height = 360,
          controls = switch_buttons(oid, c("Counts", "Percent"), c("counts", "percent")),
          variants = vars,
          caption = paste0("ITSx assigns each sequence a preliminary origin from its HMM profile. ",
                           "Anything other than the target group is a contamination signal, not a taxonomic result."))))
    }
  }

  ## Region lengths, parsed from the per-sequence position lines.
  if (length(itsx_files$positions) > 0L) {
    pos <- rbindlist(lapply(itsx_files$positions, function(f) {
      ln <- readLines(f, warn = FALSE)
      if (length(ln) == 0L) return(NULL)
      m <- regmatches(ln, gregexpr("(SSU|ITS1|5\\.8S|ITS2|LSU): ([0-9]+)-([0-9]+)", ln))
      rbindlist(lapply(m, function(hits) {
        if (length(hits) == 0L) return(NULL)
        parts <- regmatches(hits, regexec("(SSU|ITS1|5\\.8S|ITS2|LSU): ([0-9]+)-([0-9]+)", hits))
        data.table(Region = vapply(parts, `[`, "", 2L),
                   Length = as.numeric(vapply(parts, `[`, "", 4L)) - as.numeric(vapply(parts, `[`, "", 3L)) + 1)
      }), fill = TRUE)
    }), fill = TRUE)
    if (!is.null(pos) && nrow(pos) > 0L) {
      lev <- intersect(c("SSU", "ITS1", "5.8S", "ITS2", "LSU"), unique(pos$Region))
      groups <- lapply(lev, function(r) pos[Region == r, Length])
      names(groups) <- lev
      o <- box_option(groups, y_name = "Length, bp")
      if (!is.null(o)) bits <- c(bits, list(
        subhead("itsx-regions", "Detected region lengths"),
        echart(o, height = 320,
          caption = "Length distribution of each rRNA region detected by ITSx, pooled across samples.")))
    }
  }

  if (length(itsx_files$problematic) > 0L) {
    prob <- rbindlist(lapply(itsx_files$problematic, function(f) {
      x <- tryCatch(fread(f, sep = "\t", header = FALSE, showProgress = FALSE), error = function(e) NULL)
      if (is.null(x) || ncol(x) < 2L) return(NULL)
      data.table(Reason = as.character(x[[2L]]))
    }), fill = TRUE)
    if (!is.null(prob) && nrow(prob) > 0L) {
      agg <- prob[, .(Sequences = .N), by = Reason][order(-Sequences)]
      bits <- c(bits, list(
        subhead("itsx-problems", "Problematic sequences"),
        dt_table(agg, fmt = list(Sequences = "int"), bar_cols = "Sequences",
                 search = FALSE, cols_menu = FALSE, download = "itsx_problematic.tsv")))
    }
  }

  c(bits, list(
    subhead("itsx-table", "Per-sample ITSx counts"),
    dt_table(det[, .(Sample, Input, Detected, `Not detected`, Chimeric, `Detection rate`)],
             fmt = list(Input = "int", Detected = "int", `Not detected` = "int",
                        Chimeric = "int", `Detection rate` = "pct"),
             bar_cols = "Input", download = "itsx_per_sample.tsv", cols_menu = FALSE)))
} else NULL

## ------------------------------------------------- chimeras and tag jumps

counts_dir <- opt$`counts-dir`
counts_file <- function(name) {
  if (is.null(counts_dir) || !dir.exists(counts_dir)) return(NULL)
  f <- file.path(counts_dir, name)
  if (is_usable_file(f)) f else NULL
}

chim_section <- {
  bits <- list()
  cutoff <- suppressWarnings(as.numeric(getp(params, "max_ChimeraScore", NA)))

  ## Prefer the full sequence table: it carries a score for every sequence, so
  ## the cut-off can be drawn against the distribution it actually acts on.
  all_scores <- if (is_usable_file(opt$seqtab) && requireNamespace("arrow", quietly = TRUE)) {
    tryCatch({
      x <- as.data.table(arrow::read_parquet(opt$seqtab, col_select = c("DeNovo_Chimera_Score")))
      v <- suppressWarnings(as.numeric(x[[1L]]))
      v[is.finite(v)]
    }, error = function(e) NULL)
  } else NULL

  if (!is.null(all_scores) && length(all_scores) > 0L) {
    o <- hist_option(all_scores, bins = 40L, x_name = "VSEARCH de novo chimera score",
                     y_name = "Sequences", colour = "#1d6fa5",
                     marklines = bin_markline(all_scores, 40L, cutoff,
                                              paste0("max_ChimeraScore = ", cutoff)),
                     log_y = TRUE)
    if (!is.null(o)) bits <- c(bits, list(
      subhead("chim-scores", "De novo chimera scores"),
      echart(o, height = 320, title = "Score distribution across all sequences",
        caption = paste0("Sequences scoring at or above ", if (is.finite(cutoff)) cutoff else "the cut-off",
                         " are discarded. Counts are on a log scale. If the cut-off does not fall in a trough ",
                         "between two modes, it is separating little."))))
  }

  f <- counts_file("Counts_7.ChimDenov.txt")
  if (!is.null(f)) {
    sc <- tryCatch(fread(f, sep = "\t", header = FALSE, showProgress = FALSE,
                         col.names = c("SeqID", "Score", "Sample")), error = function(e) NULL)
    if (!is.null(sc) && nrow(sc) > 0L) {
      if (is.null(all_scores)) {
        o <- hist_option(sc$Score, bins = 40L, x_name = "VSEARCH de novo chimera score",
                         y_name = "Sequences", colour = "#c9184a", log_y = TRUE)
        if (!is.null(o)) bits <- c(bits, list(
          subhead("chim-scores", "De novo chimera scores"),
          echart(o, height = 320, title = "Score distribution of flagged sequences",
            caption = "Only sequences that were called chimeric are listed in this file, so every value is above the cut-off.")))
      }
      per <- sc[, .(`Chimeric sequences` = .N, `Median score` = round(stats::median(Score), 3)), by = Sample][order(-`Chimeric sequences`)]
      bits <- c(bits, list(
        if (is.null(all_scores)) NULL else subhead("chim-per-sample", "Chimeras per sample"),
        dt_table(per,
          fmt = list(`Chimeric sequences` = "int", `Median score` = "num"),
          bar_cols = "Chimeric sequences", download = "chimera_per_sample.tsv", cols_menu = FALSE)))
    }
  }

  ref_reads <- sum(chim_ref)
  if (ref_reads > 0) {
    ref <- data.table(Sample = per_sample$Sample, Chimeric = chim_ref, Recovered = recovered)
    ref <- ref[Chimeric > 0 | Recovered > 0][order(-Chimeric)]
    bits <- c(bits, list(
      subhead("chim-ref", "Reference-based chimeras"),
      dt_table(ref, fmt = list(Chimeric = "int", Recovered = "int"),
               bar_cols = "Chimeric", download = "chimera_reference.tsv", cols_menu = FALSE)))
  }

  ## Tag jumps
  ## A real run has millions of sequence-by-sample occurrences, so report summary statistics rather than the distribution itself
  if (is_usable_file(opt$`tagjump-scores`) && requireNamespace("qs2", quietly = TRUE)) {
    tj <- tryCatch(as.data.table(qs2::qs_read(opt$`tagjump-scores`, nthreads = 1L)),
                   error = function(e) NULL)
    if (!is.null(tj) && nrow(tj) > 0L && "Score" %in% names(tj)) {
      tj[, Score := suppressWarnings(as.numeric(Score))]
      flagged <- if ("TagJump" %in% names(tj)) sum(tj$TagJump %in% TRUE) else NA_integer_
      reads_out <- if (all(c("TagJump", "Abundance") %in% names(tj))) {
        sum(num0(tj$Abundance[tj$TagJump %in% TRUE]))
      } else NA_real_
      reads_all <- if ("Abundance" %in% names(tj)) sum(num0(tj$Abundance)) else NA_real_

      tj_summary <- data.table(
        Metric = c("Sequence-by-sample occurrences", "Flagged as tag jumps",
                   "Reads removed", "Reads removed, % of pre-filter total",
                   "Median UNCROSS2 score", "Maximum UNCROSS2 score"),
        Value  = c(fmt_int(nrow(tj)), fmt_int(flagged), fmt_int(reads_out),
                   fmt_pct(safe_pct(reads_out, reads_all), 3),
                   fmt_num(stats::median(tj$Score, na.rm = TRUE), 4),
                   fmt_num(max(tj$Score, na.rm = TRUE), 4)))

      bits <- c(bits, list(
        subhead("tagjump", "Tag-jump filtering"),
        dt_table(tj_summary, fmt = list(Value = "chr"), search = FALSE, cols_menu = FALSE,
                 download = "tagjump_summary.tsv")))
    }
  }

  if (is_usable_file(opt$`tagjump-stats`)) {
    ts <- tryCatch(fread(opt$`tagjump-stats`, sep = "\t", showProgress = FALSE), error = function(e) NULL)
    if (!is.null(ts) && nrow(ts) > 0L) {
      bits <- c(bits, list(dt_table(ts,
        fmt = as.list(setNames(rep("int", ncol(ts)), names(ts))),
        search = FALSE, cols_menu = FALSE)))
    }
  }

  tj_events <- sum(num0(col_or(per_sample, 0, "TagJump_Events")))
  if (tj_events > 0) {
    tjt <- data.table(Sample = per_sample$Sample,
                      Events = num0(per_sample$TagJump_Events),
                      Reads  = num0(per_sample$TagJump_Reads))[Events > 0][order(-Reads)]
    bits <- c(bits, list(dt_table(tjt, fmt = list(Events = "int", Reads = "int"),
      bar_cols = "Reads", download = "tagjump_per_sample.tsv", cols_menu = FALSE)))
  }

  if (length(bits) == 0L) NULL else bits
}

## ------------------------------------------------------- length & quality

len_section <- {
  bits <- list()

  ## Per-stage read length, straight out of the seqkit stat tables.
  stages <- list(
    c("Counts_3.Demux.txt",           "Demultiplexed"),
    c("Counts_4.PrimerCheck.txt",     "Primer-checked"),
    c("Counts_4.PrimerArtefacts.txt", "Primer artefacts"))
  len <- rbindlist(lapply(stages, function(s) {
    f <- counts_file(s[[1L]])
    if (is.null(f)) return(NULL)
    x <- tryCatch(fread(f, sep = "\t", showProgress = FALSE), error = function(e) NULL)
    if (is.null(x) || !"avg_len" %in% names(x)) return(NULL)
    data.table(Stage = s[[2L]], Sample = clean_sample_name(x$file),
               Min = x$min_len, Avg = x$avg_len, Max = x$max_len, N = x$num_seqs)
  }), fill = TRUE)

  if (!is.null(len) && nrow(len) > 0L) {
    bits <- c(bits, list(
      subhead("len-stage", "Read length by stage"),
      p(class = "hint", paste0(
        "Length of the reads entering and leaving each stage, per sample. ",
        "Artefact reads that differ sharply in length point at primer-dimer or ",
        "concatemer products.")),
      dt_table(len[order(Stage, -N)],
        fmt = list(Min = "bp", Avg = "bp", Max = "bp", N = "int"),
        labels = c("Stage", "Sample", "Min length", "Mean length", "Max length", "Sequences"),
        download = "read_length_by_stage.tsv")))
  }

  ## Read length against quality, binned in R so the page stays small.
  if (is_usable_file(opt$`seq-quals`) && requireNamespace("arrow", quietly = TRUE)) {
    q <- tryCatch({
      ds <- arrow::read_parquet(opt$`seq-quals`,
                                col_select = c("Length", "AvgPhredScore"))
      as.data.table(ds)
    }, error = function(e) NULL)
    if (!is.null(q) && nrow(q) > 0L) {
      setnames(q, c("Length", "AvgPhredScore"), c("len", "phred"), skip_absent = TRUE)
      q <- q[is.finite(len) & is.finite(phred)]
      if (nrow(q) > 0L) {
        nb <- 44L
        lb <- seq(min(q$len), max(q$len), length.out = nb + 1L)
        pb <- seq(min(q$phred), max(q$phred), length.out = nb + 1L)
        if (diff(range(lb)) > 0 && diff(range(pb)) > 0) {
          q[, xi := pmin(nb, pmax(1L, findInterval(len, lb, rightmost.closed = TRUE)))]
          q[, yi := pmin(nb, pmax(1L, findInterval(phred, pb, rightmost.closed = TRUE)))]
          agg <- q[, .(value = .N), by = .(xi, yi)]
          xl <- as.character(round((lb[-1L] + lb[-length(lb)]) / 2))
          yl <- as.character(round((pb[-1L] + pb[-length(pb)]) / 2, 1))
          agg[, x := xl[xi]][, y := yl[yi]]
          bits <- c(bits, list(
            subhead("len-quality", "Read length versus quality"),
            echart(heatmap_option(agg, xl, yl, "Read length, bp", "Mean Phred score"),
              height = 400,
              caption = paste0("Density of ", fmt_int(nrow(q)), " reads. A long low-quality tail usually means the ",
                               "length or expected-error filters need tightening."))))
        }
      }
    }
  }

  ## MEEP, which drives the max_MEEP filter.
  if (is_usable_file(opt$seqtab) && requireNamespace("arrow", quietly = TRUE)) {
    st <- tryCatch(as.data.table(arrow::read_parquet(opt$seqtab, col_select = c("MEEP"))),
                   error = function(e) NULL)
    if (!is.null(st) && nrow(st) > 0L && any(is.finite(st$MEEP))) {
      meep_cut <- suppressWarnings(as.numeric(getp(params, "max_MEEP", NA)))
      v <- st$MEEP[is.finite(st$MEEP)]
      ml <- bin_markline(v, 40L, meep_cut, paste0("max_MEEP = ", meep_cut))
      o <- hist_option(v, bins = 40L, x_name = "MEEP (max expected errors per 100 bp)",
                       y_name = "Sequences", colour = "#1d6fa5", marklines = ml, log_y = TRUE)
      if (!is.null(o)) bits <- c(bits, list(
        subhead("len-meep", "Expected error rate"),
        echart(o, height = 300,
          caption = "MEEP of the sequences that reached the sequence table.")))
    }
  }

  if (length(bits) == 0L) NULL else bits
}

## ------------------------------------------------------------- assemble

meta <- c(
  "NextITS"  = nextits_version_label(versions),
  "Nextflow" = version_label(versions, "Nextflow"),
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
        download = "step1_stage_losses.tsv"),
      pipeline_links(stamp = opt$`trace-stamp`),
      subs = list(c("overview-stages", "Stage losses"))),

  sec("samples", "Per-sample",
      fate_fig, scatter_fig, spread_fig,
      subhead("samples-table", "Per-sample counts"),
      dt_table(sample_view,
        labels = unname(ifelse(names(sample_view) %in% names(sample_labels),
                               sample_labels[names(sample_view)],
                               gsub("_", " ", names(sample_view)))),
        fmt = sample_fmt, bar_cols = "Demultiplexed_Reads",
        download = "step1_per_sample.tsv"),
      subs = list(c("samples-table", "Counts table"))),

  if (!is.null(demux_section)) sec("demux", "Demultiplexing", demux_section,
      note = "PacBio lima barcode assignment."),

  if (!is.null(itsx_section)) sec("itsx", "ITSx", itsx_section,
      note = paste0("rRNA region extraction with ITSx (target group: ",
                    getp(params, "ITSx_tax", "unspecified"), ").")),

  if (!is.null(chim_section)) sec("chimeras", "Chimeras and tag jumps", chim_section),

  if (!is.null(len_section)) sec("quality", "Read length and quality", len_section),

  sec("settings", "Run settings",
      command_panel(opt$command, extra = c(
        "Sequencing platform" = getp(params, "seqplatform", "\u2014"),
        "rRNA region"         = getp(params, "its_region", "\u2014"),
        "Forward primer"      = getp(params, "primer_forward", "\u2014"),
        "Reverse primer"      = getp(params, "primer_reverse", "\u2014"))),
      p(class = "hint", HTML(paste0(
        "Every resolved parameter, including defaults, is listed in ",
        '<a href="../pipeline_info/pipeline_params.tsv">pipeline_params.tsv</a>.')))),

  sec("methods", "Methods and software",
      methods_panel(methods),
      if (!is.null(versions_panel(versions))) subhead("methods-versions", "Software versions"),
      versions_panel(versions),
      subs = list(c("methods-versions", "Software versions")))
)

report_page(
  title = "NextITS Step-1 report",
  meta = meta,
  sections = sections,
  assets_dir = opt$assets,
  out = opt$out,
  logo = opt$logo)

cat("Wrote ", normalizePath(opt$out, mustWork = FALSE), "\n", sep = "")
cat("Elapsed minutes: ", round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 3), "\n", sep = "")
