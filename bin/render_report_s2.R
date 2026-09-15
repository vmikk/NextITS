#!/usr/bin/env Rscript

## Self-contained HTML run report for NextITS Step-2.
## Built from the pooled OTU table plus the clustering lineage in
## UC_Pooled.parquet; optional panels drop out when their input is missing.

start_time <- Sys.time()

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

option_list <- list(
  make_option("--otutab",      type = "character", default = NULL, help = "OTU_table_long.txt.gz"),
  make_option("--otus",        type = "character", default = NULL, help = "OTUs.fa.gz"),
  make_option("--lulu-otutab", type = "character", default = NULL, help = "OTU_table_LULU.txt.gz"),
  make_option("--lulu-stats",  type = "character", default = NULL, help = "LULU_merging_statistics.txt.gz"),
  make_option("--uc-pooled",   type = "character", default = NULL, help = "UC_Pooled.parquet"),
  make_option("--params",      type = "character", default = NULL, help = "pipeline_params.tsv"),
  make_option("--versions",    type = "character", default = NULL, help = "software_versions.yml"),
  make_option("--methods",     type = "character", default = NULL, help = "README_Step2_Methods.txt"),
  make_option("--schema",      type = "character", default = NULL, help = "nextflow_schema.json"),
  make_option("--assets",      type = "character", default = "assets", help = "Report asset directory [%default]"),
  make_option("--out",         type = "character", default = "Step2_report.html", help = "Output HTML [%default]"),
  make_option("--summary-tsv", type = "character", default = "Step2_summary.tsv", help = "Per-sample summary TSV [%default]"),
  make_option("--min-reads",     type = "double", default = 100, help = "Warn below this library size [%default]"),
  make_option("--max-singleton-pct", type = "double", default = 50, help = "Warn above this singleton-OTU %% [%default]"),
  make_option("--rarefaction-max", type = "integer", default = 300L, help = "Skip rarefaction above this sample count [%default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1L]
script_dir <- if (!is.na(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg), mustWork = FALSE)) else getwd()
lib <- c("report_lib.R", file.path(script_dir, "report_lib.R"))
lib <- lib[file.exists(lib)][1L]
if (is.na(lib)) stop("Could not find report_lib.R", call. = FALSE)
source(lib)

if (!is_usable_file(opt$otutab)) stop("--otutab is required", call. = FALSE)

cat("Loading Step-2 OTU table\n")
otu <- fread(opt$otutab, sep = "\t", showProgress = FALSE)

## In dereplication-only mode (preclustering = clustering = "none") the table
## is keyed by DerepID rather than OTU; everything downstream treats the first
## column as the sequence unit, whatever it is called.
unit_col <- pick_col(otu, "OTU", "DerepID", "SeqID")
if (is.null(unit_col)) unit_col <- names(otu)[[1L]]
if (!identical(unit_col, "OTU")) setnames(otu, unit_col, "OTU")

missing <- setdiff(c("OTU", "SampleID", "Abundance"), names(otu))
if (length(missing) > 0L) {
  stop("OTU table is missing column(s): ", paste(missing, collapse = ", "),
       ". Found: ", paste(names(otu), collapse = ", "), call. = FALSE)
}
otu[, Abundance := as.numeric(Abundance)]
otu <- otu[is.finite(Abundance) & Abundance > 0]

params   <- read_params_tsv(opt$params)
versions <- read_versions(opt$versions)
methods  <- read_methods_text(opt$methods)

## What the pipeline actually produced depends on how it was configured.
clustering <- tolower(getp(params, "clustering", "none"))
preclust   <- tolower(getp(params, "preclustering", "none"))
unit <- if (clustering == "none" && preclust == "none") {
  "unique sequences"
} else if (grepl("dada2|unoise", paste(preclust, clustering))) {
  "ASVs / OTUs"
} else {
  "OTUs"
}

otu_tot <- otu[, .(Total = sum(Abundance), Samples = .N), by = OTU]
setorder(otu_tot, -Total)
n_otu     <- nrow(otu_tot)
n_samples <- uniqueN(otu$SampleID)
n_reads   <- sum(otu$Abundance)
singleton <- otu_tot[Total < 2, OTU]

per_sample <- otu[, .(
  Reads = sum(Abundance),
  OTUs  = uniqueN(OTU),
  Singleton_OTUs  = sum(OTU %chin% singleton),
  Singleton_reads = sum(Abundance[OTU %chin% singleton])
), by = SampleID]
per_sample[, Singleton_Percent := safe_pct(Singleton_OTUs, OTUs)]
per_sample[, RunID := run_of(SampleID)]
setorder(per_sample, -Reads, SampleID)
runs <- sort(unique(stats::na.omit(per_sample$RunID)))

fwrite(per_sample[, .(SampleID, Reads, OTUs, Singleton_OTUs, Singleton_reads)],
       opt$`summary-tsv`, sep = "\t")

depth <- as.numeric(per_sample$Reads)
dq <- stats::quantile(depth, c(0, 0.25, 0.5, 0.75, 1), names = FALSE)

## ------------------------------------------------------------- KPI cards

singleton_pct <- safe_pct(length(singleton), n_otu)

kpis <- kpi_row(
  kpi_card("Samples", fmt_int(n_samples),
           if (length(runs) > 0L) paste0(length(runs), " sequencing run(s)") else NULL),
  kpi_card(unit, fmt_int(n_otu), paste0(fmt_int(n_reads), " reads assigned")),
  kpi_card("Median depth", fmt_int(dq[[3L]]),
           paste0(fmt_int(dq[[1L]]), " – ", fmt_int(dq[[5L]]), " reads per sample"),
           grade(dq[[3L]], warn = 5 * opt$`min-reads`, bad = opt$`min-reads`)),
  kpi_card("Singletons", fmt_pct(singleton_pct),
           paste0(fmt_int(length(singleton)), " of ", fmt_int(n_otu), " ", unit,
                  " seen as a single read"),
           grade(singleton_pct, warn = opt$`max-singleton-pct` / 2,
                 bad = opt$`max-singleton-pct`, higher_better = FALSE))
)

overall_tbl <- data.table(
  Metric = c("Samples", "Sequencing runs", paste0("Total ", unit), "Total reads",
             "Minimum depth", "Lower quartile depth", "Median depth",
             "Upper quartile depth", "Maximum depth",
             paste0("Median ", unit, " per sample"), paste0("Singleton ", unit)),
  Value  = c(n_samples, max(length(runs), 1L), n_otu, n_reads,
             dq[[1L]], dq[[2L]], dq[[3L]], dq[[4L]], dq[[5L]],
             stats::median(per_sample$OTUs), length(singleton)))

## ------------------------------------------------------------ QC verdicts

shallow <- per_sample[Reads < opt$`min-reads`]
rules <- list(
  if (nrow(shallow) > 0L)
    qc_rule("Library size", "warn", sprintf("%d sample(s) below the depth floor", nrow(shallow)),
            paste0("< ", fmt_int(opt$`min-reads`), " reads"), sort(shallow$SampleID))
  else qc_rule("Library size", "ok", "no samples affected",
               paste0("< ", fmt_int(opt$`min-reads`), " reads")),
  if (is.finite(singleton_pct) && singleton_pct >= opt$`max-singleton-pct`)
    qc_rule("Singletons", "warn",
            paste0(fmt_pct(singleton_pct), " of ", unit, " occur as a single read"),
            paste0("≥ ", fmt_pct(opt$`max-singleton-pct`, 0)))
  else qc_rule("Singletons", "ok", paste0(fmt_pct(singleton_pct), " of ", unit),
               paste0("≥ ", fmt_pct(opt$`max-singleton-pct`, 0)))
)

## --------------------------------------------------------- clustering funnel

lulu_wide <- if (is_usable_file(opt$`lulu-otutab`)) {
  x <- tryCatch(fread(opt$`lulu-otutab`, sep = "\t", showProgress = FALSE), error = function(e) NULL)
  if (!is.null(x) && ncol(x) > 1L) x else NULL
} else NULL

n_lulu <- if (!is.null(lulu_wide)) nrow(lulu_wide) else NA_integer_

funnel <- NULL
if (is_usable_file(opt$`uc-pooled`) && requireNamespace("arrow", quietly = TRUE)) {
  uc <- tryCatch(as.data.table(arrow::read_parquet(opt$`uc-pooled`)), error = function(e) NULL)
  if (!is.null(uc) && nrow(uc) > 0L) {
    steps <- list(
      c("SeqID",        "Input sequences"),
      c("DerepID",      "Dereplicated"),
      c("PreclusterID", if (grepl("unoise", preclust)) "Denoised (UNOISE)" else if (grepl("dada2", preclust)) "Denoised (DADA2)" else "Pre-clustered"),
      c("OTU",          if (clustering == "none") "Final sequences" else paste0("Clustered (", clustering, ")")))
    funnel <- rbindlist(lapply(steps, function(s) {
      if (!s[[1L]] %in% names(uc)) return(NULL)
      data.table(Stage = s[[2L]], N = uniqueN(uc[[s[[1L]]]]))
    }), fill = TRUE)
    ## Pre-clustering and clustering can both be off, leaving identical counts.
    if (!is.null(funnel)) funnel <- funnel[!duplicated(funnel$N) | seq_len(.N) == 1L]
  }
}
if (!is.null(funnel) && nrow(funnel) > 0L && funnel$N[[nrow(funnel)]] != n_otu) {
  funnel <- rbind(funnel, data.table(Stage = "In OTU table", N = n_otu))
}
if (!is.null(funnel) && is.finite(n_lulu)) {
  funnel <- rbind(funnel, data.table(Stage = "LULU-curated", N = n_lulu))
}
if (!is.null(funnel) && nrow(funnel) > 0L) {
  funnel[, Retained := safe_pct(N, N[[1L]])]
  funnel[, Merged := c(NA_real_, -diff(N))]
}

funnel_fig <- if (!is.null(funnel) && nrow(funnel) > 1L) {
  o <- bar_option(funnel$Stage, funnel$N, y_name = paste0("Distinct sequences / ", unit))
  o$xAxis$axisLabel$rotate <- 0
  echart(o, height = 300, title = "Clustering funnel",
    caption = paste0("Distinct sequences remaining after each collapsing step, traced through the ",
                     "SeqID \u2192 DerepID \u2192 PreclusterID \u2192 OTU lineage in UC_Pooled.parquet."))
}

## ------------------------------------------------------------ per-sample

lib_fig <- {
  o <- bar_option(per_sample$SampleID, per_sample$Reads, y_name = "Reads")
  o$series[[1L]]$markLine <- markline_y(round(stats::median(depth)), "median")
  echart(o, height = 340, title = "Library size",
    caption = "Reads per sample after all Step-1 filtering and Step-2 pooling, sorted descending.")
}

rich_fig <- {
  pts <- lapply(seq_len(nrow(per_sample)), function(i) list(
    x = max(per_sample$Reads[[i]], 1), y = per_sample$OTUs[[i]],
    name = per_sample$SampleID[[i]], size = max(per_sample$OTUs[[i]], 1)))
  echart(scatter_option(pts, "Reads", paste0(unit, " observed"),
                        x_log = TRUE, size_name = NULL),
    height = 340, title = paste0("Depth versus ", unit),
    caption = "A curve that is still climbing steeply at the deepest samples means richness is limited by sequencing effort.")
}

## Rarefaction answers the question the scatter above only hints at.
rare_fig <- if (n_samples >= 2L && n_samples <= opt$`rarefaction-max` &&
                requireNamespace("vegan", quietly = TRUE)) {
  mat <- tryCatch({
    w <- dcast(otu, SampleID ~ OTU, value.var = "Abundance", fill = 0, fun.aggregate = sum)
    ids <- w$SampleID
    m <- as.matrix(w[, -1L])
    rownames(m) <- ids
    m
  }, error = function(e) NULL)
  if (!is.null(mat) && nrow(mat) > 0L) {
    ## One curve per sample, thinned to a manageable number of points.
    series <- lapply(seq_len(nrow(mat)), function(i) {
      n <- sum(mat[i, ])
      if (n < 2) return(NULL)
      steps <- unique(round(seq(1, n, length.out = min(60L, n))))
      y <- as.numeric(vegan::rarefy(mat[i, , drop = FALSE], sample = steps))
      list(name = rownames(mat)[[i]], x = steps, y = round(y, 2))
    })
    series <- Filter(Negate(is.null), series)
    if (length(series) > 0L) {
      o <- line_option(series, "Reads sampled", paste0(unit, " expected"),
                       show_legend = length(series) <= 12L, smooth = FALSE)
      o$tooltip <- list(trigger = "item", formatter = js_fn("scatterPoint", "reads", unit, ""))
      echart(o, height = 380, title = "Rarefaction curves",
        caption = "Expected richness as a function of sampling effort. Curves that have flattened are sequenced deeply enough; curves still rising are not.")
    }
  }
}

sample_tbl <- copy(per_sample)
if (length(runs) < 2L) sample_tbl[, RunID := NULL]

## ----------------------------------------------------------------- OTUs

rank_fig <- {
  o <- list(
    grid = list(left = 8, right = 20, top = 22, bottom = 10, containLabel = TRUE),
    tooltip = list(trigger = "item", formatter = js_fn("scatterPoint", "rank", "reads", "")),
    xAxis = list(type = "value", name = "Rank", nameLocation = "middle", nameGap = 30, min = 1),
    yAxis = list(type = "log", name = "Reads", nameLocation = "middle", nameGap = 50),
    series = list(list(type = "line", showSymbol = FALSE, smooth = FALSE,
                       lineStyle = list(width = 1.8, color = PAL[[1L]]),
                       areaStyle = list(color = PAL[[1L]], opacity = 0.12),
                       data = Map(function(a, b) list(a, b), seq_len(n_otu), otu_tot$Total))))
  if (n_otu > 200L) {
    o$dataZoom <- list(list(type = "inside"), list(type = "slider", height = 18, bottom = 4))
    o$grid$bottom <- 40
  }
  echart(o, height = 320, title = paste0("Rank abundance of all ", fmt_int(n_otu), " ", unit),
    caption = "Total reads per OTU against its abundance rank, on a log scale. A long flat tail of ones is the singleton pool.")
}

occ_fig <- {
  pts <- lapply(seq_len(nrow(otu_tot)), function(i) list(
    x = otu_tot$Samples[[i]], y = otu_tot$Total[[i]], name = otu_tot$OTU[[i]], size = 16))
  o <- scatter_option(pts, "Samples occupied", "Total reads",
                      size_name = NULL)
  o$yAxis$type <- "log"
  o$yAxis$min <- NULL
  o$xAxis$min <- 0
  o$xAxis$minInterval <- 1
  o$series[[1L]]$symbolSize <- 6
  o$series[[1L]]$large <- TRUE
  o$series[[1L]]$largeThreshold <- 2000
  echart(o, height = 340, title = "Occupancy versus abundance",
    caption = "Abundant OTUs found in only one sample are the classic signature of an artefact or a contaminant.")
}

len_fig <- if (is_usable_file(opt$otus) && requireNamespace("Biostrings", quietly = TRUE)) {
  sq <- tryCatch(Biostrings::readDNAStringSet(opt$otus), error = function(e) NULL)
  if (!is.null(sq) && length(sq) > 0L) {
    w <- as.numeric(Biostrings::width(sq))
    lo <- suppressWarnings(as.numeric(getp(params, "ampliconlen_min", NA)))
    hi <- suppressWarnings(as.numeric(getp(params, "ampliconlen_max", NA)))
    ml <- bin_markline(w, 40L, lo, paste0("ampliconlen_min = ", lo))
    if (is.null(ml)) ml <- bin_markline(w, 40L, hi, paste0("ampliconlen_max = ", hi))
    o <- hist_option(w, bins = 40L, x_name = "Sequence length, bp",
                     y_name = paste0(unit, " count"), marklines = ml)
    if (!is.null(o)) echart(o, height = 300, title = paste0(unit, " length distribution"),
      caption = paste0("Length of the ", fmt_int(length(sq)), " representative sequences. ",
                       "A secondary mode well away from the target amplicon usually means off-target amplification."))
  }
}

## ------------------------------------------------------------------ LULU

lulu_section <- if (!is.null(lulu_wide)) {
  long <- melt(lulu_wide, id.vars = names(lulu_wide)[[1L]],
               variable.name = "SampleID", value.name = "Abundance")
  setnames(long, 1L, "OTU")   # whatever the first column was called
  long[, Abundance := as.numeric(Abundance)]
  long <- long[is.finite(Abundance) & Abundance > 0]
  lulu_tot <- long[, .(Total = sum(Abundance)), by = OTU]

  cmp <- data.table(
    Stage = c("Before LULU", "After LULU"),
    OTUs = c(n_otu, nrow(lulu_tot)),
    `Singleton OTUs` = c(length(singleton), sum(lulu_tot$Total < 2)),
    Reads = c(n_reads, sum(lulu_tot$Total)))
  merged <- n_otu - nrow(lulu_tot)

  bits <- list(
    callout(paste0("LULU merged ", fmt_int(merged), " ", unit, " (",
                   fmt_pct(safe_pct(merged, n_otu)), " of the total) into their parents.")),
    dt_table(cmp, fmt = list(OTUs = "int", `Singleton OTUs` = "int", Reads = "int"),
             search = FALSE, cols_menu = FALSE))

  if (is_usable_file(opt$`lulu-stats`)) {
    st <- tryCatch(fread(opt$`lulu-stats`, sep = "\t", showProgress = FALSE), error = function(e) NULL)
    if (!is.null(st) && nrow(st) > 0L) {
      if ("status" %in% names(st)) {
        agg <- st[, .(N = .N), by = status][order(-N)]
        bits <- c(bits, list(
          subhead("lulu-status", "Merge decisions"),
          dt_table(agg, labels = c("Status", unit), fmt = list(N = "int"),
                   bar_cols = "N", search = FALSE, cols_menu = FALSE)))
      }
      keep <- intersect(c("parent_id", "daughter_id", "match", "spread", "rel_cooccurence",
                          "curationlevel", "status", "total", "parent_total"), names(st))
      if (length(keep) > 0L) {
        show <- st[, ..keep]
        if ("status" %in% names(show)) show <- show[order(status)]
        bits <- c(bits, list(
          subhead("lulu-detail", "Merging statistics"),
          dt_table(head(show, 2000L),
                   download = "lulu_merging_statistics.tsv",
                   caption = if (nrow(show) > 2000L)
                     paste0("Showing the first 2,000 of ", fmt_int(nrow(show)),
                            " rows; the full table is in 05.LULU/.") else NULL)))
      }
    }
  }
  bits
} else NULL

## -------------------------------------------------------------- assemble

meta <- c(
  "NextITS"  = version_label(versions, "NextITS"),
  "Nextflow" = version_label(versions, "Nextflow"),
  "Samples"  = fmt_int(n_samples),
  setNames(fmt_int(n_otu), unit))
if (length(runs) > 0L) meta["Runs"] <- if (length(runs) <= 3L) paste(runs, collapse = ", ") else fmt_int(length(runs))
meta["Generated"] <- format(Sys.time(), "%Y-%m-%d %H:%M %Z")

sections <- list(
  sec("overview", "Run overview",
      kpis,
      qc_panel(rules),
      funnel_fig,
      if (!is.null(funnel) && nrow(funnel) > 1L) tagList(
        subhead("overview-funnel", "Sequence collapsing"),
        dt_table(funnel, labels = c("Stage", "Distinct", "Retained %", "Merged away"),
                 fmt = list(N = "int", Retained = "pct", Merged = "int"),
                 bar_cols = "N", search = FALSE, cols_menu = FALSE,
                 download = "step2_clustering_funnel.tsv")),
      subhead("overview-stats", "Overall statistics"),
      dt_table(overall_tbl, fmt = list(Value = "int"), search = FALSE, cols_menu = FALSE,
               download = "step2_overall_stats.tsv"),
      pipeline_links(),
      subs = list(c("overview-stats", "Overall statistics"))),

  sec("samples", "Per-sample",
      lib_fig, rich_fig, rare_fig,
      subhead("samples-table", "Per-sample totals"),
      dt_table(sample_tbl,
        labels = c(SampleID = "Sample", RunID = "Run", Reads = "Reads", OTUs = unit,
                   Singleton_OTUs = "Singletons", Singleton_reads = "Singleton reads",
                   Singleton_Percent = "Singleton %")[names(sample_tbl)],
        fmt = list(Singleton_Percent = "pct"),
        bar_cols = "Reads", download = "step2_per_sample.tsv"),
      subs = list(c("samples-table", "Totals table"))),

  sec("otus", unit,
      rank_fig, occ_fig, len_fig),

  if (!is.null(lulu_section)) sec("lulu", "LULU curation", lulu_section),

  sec("settings", "Run settings",
      params_panel(params, opt$schema,
                   only_groups = c("step_2_specific_parameters", "common_parameters"))),

  sec("methods", "Methods and software",
      methods_panel(methods),
      if (!is.null(versions_panel(versions))) subhead("methods-versions", "Software versions"),
      versions_panel(versions),
      subs = list(c("methods-versions", "Software versions")))
)

report_page(
  title = "NextITS Step-2 report",
  meta = meta,
  sections = sections,
  assets_dir = opt$assets,
  out = opt$out,
  footer = "Generated by NextITS render_report_s2.R. Charts use Apache ECharts (Apache-2.0).")

cat("Wrote ", normalizePath(opt$out, mustWork = FALSE), "\n", sep = "")
cat("Wrote ", normalizePath(opt$`summary-tsv`, mustWork = FALSE), "\n", sep = "")
cat("Elapsed minutes: ", round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 3), "\n", sep = "")
