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
  make_option("--command",     type = "character", default = NULL, help = "execution_command.txt"),
  make_option("--logo",        type = "character", default = NULL, help = "NextITS logo (SVG or PNG)"),
  make_option("--trace-stamp", type = "character", default = NULL, help = "Timestamp Nextflow stamped into pipeline_info filenames"),
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

## `Max` is the largest single-sample count, used to tell one-sample dominants from widespread OTUs
otu_tot <- otu[, .(Total = sum(Abundance), Samples = uniqueN(SampleID), Max = max(Abundance)), by = OTU]
setorder(otu_tot, -Total, OTU)
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
    ## Each step carries the label for the sequences it folds away, so the sankey branches can be named after what actually merged them
    steps <- list(
      list("SeqID",        "Input sequences", NA_character_),
      list("DerepID",      "Dereplicated",    "Identical sequences"),
      list("PreclusterID",
           if (grepl("unoise", preclust)) "Denoised (UNOISE)"
           else if (grepl("dada2", preclust)) "Denoised (DADA2)"
           else "Pre-clustered",
           if (grepl("unoise|dada2", preclust)) "Merged by denoising" else "Merged by pre-clustering"),
      list("OTU",
           if (clustering == "none") "Final sequences" else paste0("Clustered (", clustering, ")"),
           if (clustering == "none") "Collapsed" else "Merged by clustering"))
    funnel <- rbindlist(lapply(steps, function(s) {
      if (!s[[1L]] %in% names(uc)) return(NULL)
      data.table(Stage = s[[2L]], N = uniqueN(uc[[s[[1L]]]]), Loss = s[[3L]])
    }), fill = TRUE)
    ## Pre-clustering and clustering can both be off, leaving identical counts.
    if (!is.null(funnel)) funnel <- funnel[!duplicated(funnel$N) | seq_len(.N) == 1L]
  }
}
if (!is.null(funnel) && nrow(funnel) > 0L && funnel$N[[nrow(funnel)]] != n_otu) {
  funnel <- rbind(funnel,
    data.table(Stage = "In OTU table", N = n_otu, Loss = "No reads left"), fill = TRUE)
}
if (!is.null(funnel) && is.finite(n_lulu)) {
  funnel <- rbind(funnel,
    data.table(Stage = "LULU-curated", N = n_lulu, Loss = "Merged by LULU"), fill = TRUE)
}
if (!is.null(funnel) && nrow(funnel) > 0L) {
  funnel[, Retained := safe_pct(N, N[[1L]])]
  funnel[, Merged := c(NA_real_, -diff(N))]
}

funnel_fig <- if (!is.null(funnel) && nrow(funnel) > 1L) {
  ## Each step keeps some sequences and folds the rest into their parents
  ## the folded ones branch off so the widths add up at every stage
  lk <- data.table(source = character(), target = character(), value = numeric())
  cols <- c(setNames(rep("var(--fate-trunk)", nrow(funnel)), funnel$Stage))
  for (i in seq_len(nrow(funnel) - 1L)) {
    from <- funnel$Stage[[i]]
    to   <- funnel$Stage[[i + 1L]]
    kept <- funnel$N[[i + 1L]]
    gone <- funnel$N[[i]] - kept
    lk <- rbind(lk, data.table(source = from, target = to, value = kept))
    if (gone > 0) {
      lab <- funnel$Loss[[i + 1L]]
      if (is.na(lab) || !nzchar(lab)) lab <- paste0("Removed before ", to)
      ## Sankey node names must be unique.
      while (lab %in% names(cols)) lab <- paste0(lab, " ")
      lk <- rbind(lk, data.table(source = from, target = lab, value = gone))
      cols[[lab]] <- "var(--fate-other)"
    }
  }
  if (nrow(lk) > 0L) {
    echart(sankey_option(unique(c(lk$source, lk$target)), lk, colours = cols),
      height = 340, title = "Clustering funnel",
      caption = paste0("Where each distinct sequence ended up, traced through the ",
                       "SeqID \u2192 DerepID \u2192 PreclusterID \u2192 OTU lineage in ",
                       "UC_Pooled.parquet. Grey branches are sequences folded into a parent."))
  }
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
  ## Rarefied richness depends only on a sample's nonzero counts, so each sample is handled on its own
  ## a dense sample-by-OTU matrix would not fit in memory for large runs
  counts <- otu[, .(Abundance = sum(Abundance)), by = .(SampleID, OTU)]
  counts <- split(counts$Abundance, counts$SampleID)
  counts <- counts[sort(names(counts), method = "radix")]   # C-locale order, independent of the system locale
  if (length(counts) > 0L) {
    ## One curve per sample, thinned to a manageable number of points.
    series <- lapply(names(counts), function(s) {
      v <- counts[[s]]
      n <- sum(v)
      if (n < 2) return(NULL)
      steps <- unique(round(seq(1, n, length.out = min(60L, n))))
      y <- as.numeric(vegan::rarefy(v, sample = steps))
      list(name = s, x = steps, y = round(y, 2))
    })
    series <- Filter(Negate(is.null), series)
    if (length(series) > 0L) {
      o <- line_option(series, "Reads sampled", paste0(unit, " expected"),
                       show_legend = length(series) <= length(PAL), smooth = FALSE)
      o$tooltip <- list(trigger = "item", formatter = js_fn("scatterPoint", "reads", unit, ""))
      echart(o, height = 380, title = "Rarefaction curves",
        caption = "Expected richness as a function of sampling effort. Curves that have flattened are sequenced deeply enough; curves still rising are not.")
    }
  }
}

sample_tbl <- copy(per_sample)
if (length(runs) < 2L) sample_tbl[, RunID := NULL]

## ----------------------------------------------------------------- OTUs

## Ranks by which half and nine tenths of all reads are accounted for
cum_reads <- cumsum(otu_tot$Total)
rank_q <- vapply(c(0.5, 0.9), function(q) which(cum_reads >= q * n_reads)[1L], integer(1))

rank_fig <- {
  ## Tied totals form flat steps
  ## the first and last rank of each step reproduce the curve exactly with a few thousand points instead of one per OTU
  r <- rle(otu_tot$Total)
  last  <- cumsum(r$lengths)
  first <- last - r$lengths + 1L
  rk <- c(rbind(first, last))
  rv <- rep(r$values, each = 2L)
  keep <- !duplicated(rk)
  o <- list(
    grid = list(left = 8, right = 20, top = 22, bottom = 10),
    tooltip = list(trigger = "item", formatter = js_fn("scatterPoint", "rank", "reads", "")),
    xAxis = list(type = "log", name = "Rank (log scale)", nameLocation = "middle", nameGap = 30, min = 1),
    yAxis = list(type = "log", name = "Reads", nameLocation = "middle", nameGap = 50),
    series = list(list(type = "line", showSymbol = FALSE, smooth = FALSE,
                       lineStyle = list(width = 2, color = PAL[[1L]]),
                       areaStyle = list(color = "var(--viz-wash)"),
                       data = Map(function(a, b) list(a, b), rk[keep], rv[keep]))))
  ## Reference lines only where they fall inside the plotted range and do not coincide
  marks <- unique(rank_q[is.finite(rank_q) & rank_q > 1L])
  if (length(marks) > 0L) {
    labs <- c("50% of reads", "90% of reads")[match(marks, rank_q)]
    o$series[[1L]]$markLine <- list(silent = TRUE, symbol = "none",
      lineStyle = list(color = "var(--viz-ref)", type = "dashed", width = 1.2),
      label = list(position = "insideEndTop", fontSize = 10, formatter = "{b}"),
      data = lapply(seq_along(marks), function(i) list(xAxis = marks[[i]], name = labs[[i]])))
  }
  if (n_otu > 200L) o$dataZoom <- list(list(type = "inside"))
  echart(o, height = 320, title = paste0("Rank abundance of all ", fmt_int(n_otu), " ", unit),
    caption = paste0(
      "Total reads per ", unit, " against its abundance rank, both on log scales. ",
      "Half of all reads are in the top ", fmt_int(rank_q[[1L]]), " ", unit,
      " (", fmt_pct(safe_pct(rank_q[[1L]], n_otu), 2), "), nine tenths in the top ",
      fmt_int(rank_q[[2L]]), " (", fmt_pct(safe_pct(rank_q[[2L]], n_otu), 2), "). ",
      "The flat tail at one read is the singleton pool."))
}

## Occupancy: the summary sentence, the chart and the table share these numbers
occ_one     <- otu_tot$Samples == 1L
wide_cut    <- max(2L, ceiling(0.1 * n_samples))
occ_callout <- callout(paste0(
  fmt_pct(safe_pct(sum(occ_one), n_otu)), " of ", unit, " (holding ",
  fmt_pct(safe_pct(sum(otu_tot$Total[occ_one]), n_reads)), " of reads) occur in a single sample; ",
  fmt_int(sum(otu_tot$Samples >= wide_cut)), " occur in at least ", fmt_int(wide_cut),
  " samples (10% or more)."))

occ_fig <- {
  ## Many OTUs share the same occupancy and total (most are one read in one sample), 
  ## so one point per distinct pair keeps the picture at a fraction of the size
  agg <- otu_tot[, .(N = .N, OTU = OTU[1L]), by = .(Samples, Total)]
  binned <- nrow(agg) > 5000L
  if (binned) {
    ## Very large runs: snap totals to 1/8-decade steps so the point count stays bounded
    agg <- otu_tot[, .(Total = 10^(round(log10(Total) * 8) / 8), Samples)][
      , .(N = .N, OTU = NA_character_), by = .(Samples, Total)]
    agg[, Total := signif(Total, 4)]
  }
  setorder(agg, -N)   # single-OTU points are drawn last, on top of the dense core
  data <- lapply(seq_len(nrow(agg)), function(i) {
    d <- list(value = list(agg$Samples[[i]], agg$Total[[i]], agg$N[[i]]))
    ## A point that is a single OTU keeps its ID for the tooltip
    if (agg$N[[i]] == 1L && !is.na(agg$OTU[[i]])) d$name <- agg$OTU[[i]]
    d
  })
  ## The palest sequential steps vanish against the light surface, so the scale starts at --seq-2
  brk <- list(
    list(min = 1, max = 1, label = "1"),
    list(min = 2, max = 9, label = "2–9"),
    list(min = 10, max = 99, label = "10–99"),
    list(min = 100, label = "≥ 100"))
  for (k in seq_along(brk)) brk[[k]]$color <- sprintf("var(--seq-%d)", k + 1L)
  o <- list(
    grid = list(left = 8, right = 22, top = 34, bottom = 10),
    tooltip = list(trigger = "item",
                   formatter = js_fn("occPoint", "Samples occupied", "Total reads", unit, n_samples)),
    visualMap = list(type = "piecewise", dimension = 2L, pieces = brk,
                     orient = "horizontal", top = 0, right = 10, itemWidth = 12, itemHeight = 10,
                     text = list("", paste0(unit, " per point:")), textGap = 6),
    xAxis = list(type = "log", name = "Samples occupied (log scale)", min = 1,
                 nameLocation = "middle", nameGap = 30),
    yAxis = list(type = "log", name = "Total reads", min = 1, nameLocation = "end", nameGap = 12,
                 nameTextStyle = list(align = "left")),
    series = list(list(type = "scatter", data = data, symbolSize = 7,
                       itemStyle = list(opacity = 0.85),
                       emphasis = list(focus = "self", itemStyle = list(opacity = 1)))))
  if (length(data) > 2000L) {
    o$series[[1L]]$large <- TRUE
    o$series[[1L]]$largeThreshold <- 2000
  }
  echart(o, height = 360, title = "Occupancy versus abundance",
    caption = paste0(
      "Each point is one combination of occupancy (number of samples) and total read count; ",
      "colour shows how many ", unit, " share it",
      if (binned) " (totals rounded to 1/8-decade steps)" else "", ". ",
      "The mass at the bottom left is unreplicated low-count sequences: PCR or sequencing noise, or genuinely rare taxa. ",
      "Low read counts spread across many samples can point to cross-contamination or tag-jumps. ",
      "Abundant ", unit, " confined to one sample are often genuine local dominants, ",
      "but deserve a look if that sample is a control."))
}

top_tbl <- {
  t <- head(otu_tot, 50L)
  data.table(
    OTU         = t$OTU,
    Reads       = t$Total,
    Reads_pct   = safe_pct(t$Total, n_reads),
    Samples     = t$Samples,
    Samples_pct = safe_pct(t$Samples, n_samples),
    Top_share   = safe_pct(t$Max, t$Total))
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
    bounds <- paste(na.omit(c(
      if (is.finite(lo)) paste0("ampliconlen_min = ", lo),
      if (is.finite(hi)) paste0("ampliconlen_max = ", hi))), collapse = ", ")
    if (!is.null(o)) echart(o, height = 300, title = paste0(unit, " length distribution"),
      caption = paste0("Length of the ", fmt_int(length(sq)), " representative sequences. ",
                      if (nzchar(bounds)) paste0(" Dashed line: ", bounds, ".") else ""))
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
    }
  }
  bits
} else NULL

## -------------------------------------------------------------- assemble

meta <- c(
  "NextITS"  = nextits_version_label(versions),
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
        dt_table(funnel, cols = c("Stage", "N", "Retained", "Merged"),
                 labels = c("Stage", "Distinct", "Retained %", "Merged away"),
                 fmt = list(N = "int", Retained = "pct", Merged = "int"),
                 bar_cols = "N", search = FALSE, cols_menu = FALSE,
                 download = "step2_clustering_funnel.tsv")),
      subhead("overview-stats", "Overall statistics"),
      dt_table(overall_tbl, fmt = list(Value = "int"), search = FALSE, cols_menu = FALSE,
               download = "step2_overall_stats.tsv"),
      pipeline_links(stamp = opt$`trace-stamp`),
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
      rank_fig, occ_callout, occ_fig,
      subhead("otus-top", paste0("Most abundant ", unit)),
      dt_table(top_tbl,
        labels = c("ID", "Reads", "% of reads", "Samples", "% of samples", "Largest single-sample share"),
        fmt = list(Reads_pct = "pct", Samples_pct = "pct", Top_share = "pct"),
        bar_cols = "Reads", search = FALSE, download = "step2_top_otus.tsv",
        caption = paste0(
          "The ", min(50L, n_otu), " ", unit, " with the most reads. Present in many samples with a low ",
          "single-sample share: widespread taxa, and the usual place to look for contaminants or tag-jump sources. ",
          "A share near 100% means nearly all reads come from one sample.")),
      len_fig,
      subs = list(c("otus-top", "Most abundant"))),

  if (!is.null(lulu_section)) sec("lulu", "LULU curation", lulu_section),

  sec("settings", "Run settings",
      command_panel(opt$command, extra = c(
        "Pre-clustering"  = getp(params, "preclustering", "\u2014"),
        "Clustering"      = getp(params, "clustering", "\u2014"),
        if (clustering == "optimotu") c("Reference set" = getp(params, "optimotu_refset.id", "\u2014"))
        else c("OTU identity" = getp(params, "otu_id", "\u2014")),
        "LULU curation"   = getp(params, "lulu", "\u2014"))),
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
  title = "NextITS Step-2 report",
  meta = meta,
  sections = sections,
  assets_dir = opt$assets,
  out = opt$out,
  logo = opt$logo)

cat("Wrote ", normalizePath(opt$out, mustWork = FALSE), "\n", sep = "")
cat("Wrote ", normalizePath(opt$`summary-tsv`, mustWork = FALSE), "\n", sep = "")
cat("Elapsed minutes: ", round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 3), "\n", sep = "")
