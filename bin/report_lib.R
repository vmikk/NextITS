## Shared building blocks for the NextITS Step-1 / Step-2 HTML run reports
## (sourced by `render_report_s1.R` and `render_report_s2.R`)
##
## The R side emits structure and data only
## Colours that must follow the light/dark theme, and all interactivity comes from `assets/report/report.js`

suppressPackageStartupMessages({
  library(htmltools)
  library(jsonlite)
  library(data.table)
})

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

## ---------------------------------------------------------------- inputs

## Nextflow stages stubs like `no_lima`, `no_lulu`, ... in place of optional inputs that were never produced
## Treat those as absent, along with empty files and directories
is_usable_file <- function(path) {
  if (is.null(path) || length(path) == 0L || is.na(path[[1L]])) return(FALSE)
  path <- as.character(path[[1L]])
  if (!nzchar(path) || !file.exists(path) || dir.exists(path)) return(FALSE)
  if (file.info(path)$size <= 0) return(FALSE)
  !grepl("^no_", basename(path), ignore.case = TRUE)
}

num0 <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x[is.na(x)] <- 0
  x
}

## First matching column name, so reports tolerate the per-sample/per-run naming drift in `read_count_summary.R`
pick_col <- function(dt, ...) {
  cand <- c(...)
  hit <- intersect(cand, names(dt))
  if (length(hit) == 0L) return(NULL)
  hit[[1L]]
}

col_or <- function(dt, default = NA_real_, ...) {
  nm <- pick_col(dt, ...)
  if (is.null(nm)) return(rep(default, max(nrow(dt), 1L)))
  suppressWarnings(as.numeric(dt[[nm]]))
}

## seqkit reports the staged file name, so strip the extensions and per-stage suffixes the pipeline appends (same list as `read_count_summary.R`)
clean_sample_name <- function(x) {
  x <- as.character(x)
  ## Compound suffixes first, otherwise ".full.fasta" loses only ".fasta".
  x <- sub("\\.full\\.fasta(\\.gz)?$", "", x)
  x <- sub("\\.(ITS1|ITS2|5_8S|SSU|LSU)\\.fasta(\\.gz)?$", "", x)
  x <- sub("\\.(fastq|fq|fa|fasta)\\.gz$", "", x)
  x <- sub("\\.(fastq|fq|fa|fasta)$", "", x)
  x <- sub("_PrimerChecked$", "", x)
  x <- sub("_PrimerArtefacts$", "", x)
  x <- sub("_Chimera$", "", x)
  x <- sub("_RescuedChimera$", "", x)
  x
}

## Run IDs are only recoverable from the `Run__Sample` naming convention
run_of <- function(sample_id) {
  out <- sub("__.*$", "", as.character(sample_id))
  out[out == as.character(sample_id)] <- NA_character_
  out
}

## ------------------------------------------------------------ formatting

fmt_int <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  out <- rep("—", length(x))
  ok <- is.finite(x)
  out[ok] <- formatC(round(x[ok]), format = "d", big.mark = ",")
  out
}

fmt_num <- function(x, digits = 2L) {
  x <- suppressWarnings(as.numeric(x))
  out <- rep("—", length(x))
  ok <- is.finite(x)
  out[ok] <- formatC(x[ok], format = "f", digits = digits, big.mark = ",")
  out
}

fmt_pct <- function(x, digits = 1L) {
  x <- suppressWarnings(as.numeric(x))
  out <- rep("—", length(x))
  ok <- is.finite(x)
  out[ok] <- paste0(formatC(x[ok], format = "f", digits = digits), "%")
  out
}

fmt_bp <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  out <- rep("—", length(x))
  ok <- is.finite(x)
  out[ok] <- paste0(formatC(round(x[ok]), format = "d", big.mark = ","), " bp")
  out
}

safe_pct <- function(num, den) {
  num <- suppressWarnings(as.numeric(num))
  den <- suppressWarnings(as.numeric(den))
  n <- max(length(num), length(den))
  if (n == 0L) return(numeric())
  ## Recycle first: subsetting a length-1 denominator by a length-n mask would otherwise pad with NA
  num <- rep_len(num, n)
  den <- rep_len(den, n)
  out <- rep(NA_real_, n)
  ok <- is.finite(num) & is.finite(den) & den > 0
  out[ok] <- num[ok] / den[ok] * 100
  out
}

## ------------------------------------------------------------ parameters

read_params_tsv <- function(path) {
  if (!is_usable_file(path)) {
    return(data.table(name = character(), value = character()))
  }
  dt <- fread(path, sep = "\t", header = FALSE, col.names = c("name", "value"),
              fill = TRUE, showProgress = FALSE, colClasses = "character")
  dt[, name := as.character(name)]
  dt[, value := as.character(value)]
  dt[]
}

getp <- function(p, pname, default = NA_character_) {
  if (is.null(p) || nrow(p) == 0L) return(default)
  hit <- p[name == pname]
  if (nrow(hit) == 0L) return(default)
  v <- hit$value[[1L]]
  if (is.na(v) || identical(v, "null") || !nzchar(v)) return(default)
  v
}

read_methods_text <- function(path) {
  if (!is_usable_file(path)) return(NULL)
  paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
}

read_versions <- function(path) {
  if (!is_usable_file(path) || !requireNamespace("yaml", quietly = TRUE)) return(NULL)
  tryCatch(yaml::read_yaml(path), error = function(e) NULL)
}

version_label <- function(versions, key) {
  if (is.null(versions) || !is.list(versions)) return("unknown")
  v <- versions[[key]]
  if (is.null(v) || is.null(v$version)) return("unknown")
  as.character(v$version)
}

## Pipeline version plus the short commit it was run from, when known
## software_versions.yml carries `revision` under NextITS
nextits_version_label <- function(versions) {
  ver <- version_label(versions, "NextITS")
  rev <- if (is.list(versions) && !is.null(versions$NextITS$revision)) {
    as.character(versions$NextITS$revision)
  } else NULL
  if (is.null(rev) || !nzchar(rev) || identical(rev, "null")) ver else paste0(ver, " (", rev, ")")
}

## --------------------------------------------------------------- palette

PAL <- c("#2e7d4f", "#d98324", "#0f766e", "#9b5de5", "#c9184a",
         "#6f9c1f", "#b45309", "#5f7a67", "#3f6fa5", "#6d4c8f")

## Read-fate colours: greens are kept reads, everything else is a named loss
FATE_COL <- c(
  "Retained"          = "#2e7d4f",
  "Undemultiplexed"   = "#8a9d90",
  "Failed QC"         = "#5f7a67",
  "Primer artefacts"  = "#d98324",
  "No ITSx detection" = "#9b5de5",
  "Chimeric"          = "#c9184a",
  "Tag jumps"         = "#b45309",
  "Other filtering"   = "#adbcb1"
)

## ----------------------------------------------------------------- cards

## status: one of "ok", "warn", "bad", or NA for a neutral card.
kpi_card <- function(label, value, sub = NULL, status = NA_character_) {
  cls <- paste("kpi", if (!is.na(status)) status else "")
  div(class = trimws(cls),
      div(class = "k-label", label),
      div(class = "k-value", value),
      if (!is.null(sub)) div(class = "k-sub", sub))
}

kpi_row <- function(...) {
  items <- Filter(Negate(is.null), list(...))
  if (length(items) == 0L) return(NULL)
  div(class = "kpis", items)
}

## Grade a value against warn/bad cut-offs. `higher_better` flips the test
grade <- function(x, warn, bad, higher_better = TRUE) {
  if (!is.finite(x)) return(NA_character_)
  if (higher_better) {
    if (x < bad) "bad" else if (x < warn) "warn" else "ok"
  } else {
    if (x > bad) "bad" else if (x > warn) "warn" else "ok"
  }
}

callout <- function(..., type = "note") {
  div(class = paste("callout", if (type != "note") type else ""), ...)
}

## -------------------------------------------------------------- QC rules

qc_rule <- function(name, status, detail, threshold = NULL, samples = NULL) {
  list(name = name, status = status, detail = detail,
       threshold = threshold, samples = samples)
}

qc_panel <- function(rules) {
  rules <- Filter(Negate(is.null), rules)
  if (length(rules) == 0L) return(NULL)
  ## Worst first, so a failure is never buried under a list of passes.
  ord <- order(match(vapply(rules, function(r) r$status, ""), c("bad", "warn", "ok")))
  div(class = "qc", lapply(rules[ord], function(r) {
    head_bits <- tagList(
      span(class = "qc-badge", switch(r$status, ok = "pass", warn = "warn", bad = "fail", r$status)),
      span(class = "qc-name", r$name),
      span(class = "qc-detail", r$detail),
      if (!is.null(r$threshold)) span(class = "qc-thr", r$threshold)
    )
    if (length(r$samples) > 0L) {
      tags$details(class = paste("qc-item", r$status),
        tags$summary(head_bits),
        div(class = "qc-samples", lapply(r$samples, tags$code)))
    } else {
      div(class = paste("qc-item", r$status), div(class = "qc-head", head_bits))
    }
  }))
}

## ------------------------------------------------------------ data table

## Renders a real <table>; report.js adds sort, search, column toggle and TSV download on top
## Rows are pasted as raw HTML because building tens of thousands of htmltools tags is far too slow
dt_table <- function(dt, cols = NULL, labels = NULL, fmt = NULL,
                     bar_cols = NULL, download = NULL, search = TRUE,
                     cols_menu = TRUE, caption = NULL) {

  if (is.null(dt) || nrow(dt) == 0L) return(p(class = "empty", "No data available."))
  dt <- as.data.table(dt)
  if (!is.null(cols)) {
    cols <- intersect(cols, names(dt))
    if (length(cols) == 0L) return(p(class = "empty", "No data available."))
    dt <- dt[, ..cols]
  }
  cols   <- names(dt)
  labels <- labels %||% gsub("_", " ", cols)
  labels <- rep(labels, length.out = length(cols))
  fmt    <- fmt %||% list()

  esc <- function(x) htmlEscape(as.character(x), attribute = FALSE)

  ## Per-column: display strings, the raw value used for sorting/export, and an optional 0-100 bar width
  cells <- lapply(seq_along(cols), function(j) {
    nm  <- cols[[j]]
    val <- dt[[j]]
    kind <- if (nm %in% names(fmt)) fmt[[nm]] else if (is.numeric(val)) "int" else "chr"
    disp <- switch(kind,
      int  = fmt_int(val),
      num  = fmt_num(val),
      pct  = fmt_pct(val),
      bp   = fmt_bp(val),
      esc(ifelse(is.na(val) | !nzchar(as.character(val)), "—", as.character(val))))
    if (kind != "chr") disp <- esc(disp)
    numeric <- kind %in% c("int", "num", "pct", "bp")
    raw <- if (numeric) {
      v <- suppressWarnings(as.numeric(val)); ifelse(is.finite(v), format(v, scientific = FALSE, trim = TRUE), "")
    } else {
      esc(ifelse(is.na(val), "", as.character(val)))
    }
    width <- NULL
    if (nm %in% bar_cols && numeric) {
      v  <- suppressWarnings(as.numeric(val))
      mx <- suppressWarnings(max(v[is.finite(v)], na.rm = TRUE))
      if (is.finite(mx) && mx > 0) width <- pmax(0, pmin(100, v / mx * 100))
    }
    list(disp = disp, raw = raw, numeric = numeric, width = width, name = nm)
  })

  td <- lapply(cells, function(cc) {
    cls <- paste0("td", if (cc$numeric) " class=\"num\"" else if (cc$name %in% c("Sample", "SampleID", "file", "OTU")) " class=\"sample\"" else "")
    inner <- if (is.null(cc$width)) cc$disp else {
      w <- ifelse(is.finite(cc$width), formatC(cc$width, format = "f", digits = 1), "0")
      paste0('<span class="cell-bar" style="width:', w, '%"></span><span class="cell-val">', cc$disp, '</span>')
    }
    paste0("<", cls, " data-v=\"", cc$raw, "\">", inner, "</td>")
  })

  body <- paste0("<tr>", Reduce(function(a, b) paste0(a, b), td), "</tr>", collapse = "")
  head <- paste0("<tr>", paste0(
    "<th", ifelse(vapply(cells, function(cc) cc$numeric, logical(1)), " class=\"num\"", ""), ">",
    esc(labels), "</th>", collapse = ""), "</tr>")

  toolbar <- if (search || cols_menu || !is.null(download)) {
    div(class = "tbl-bar",
      if (search) tags$input(type = "search", placeholder = "Filter rows…", `aria-label` = "Filter table rows"),
      if (cols_menu) div(class = "cols-menu",
        tags$button(type = "button", class = "btn", "Columns"),
        div(class = "cols-list")),
      if (!is.null(download)) tags$button(type = "button", class = "btn", `data-download` = download, "Download TSV"),
      span(class = "count"))
  }

  div(class = "tbl-wrap",
    toolbar,
    div(class = "tbl-scroll",
      HTML(paste0('<table class="dt"><thead>', head, '</thead><tbody>', body, '</tbody></table>'))),
    if (!is.null(caption)) div(style = "padding:8px 12px 10px", class = "empty", caption))
}

## ---------------------------------------------------------------- charts

.chart_seq <- new.env(parent = emptyenv())
.chart_seq$n <- 0L

new_chart_id <- function() {
  .chart_seq$n <- .chart_seq$n + 1L
  paste0("ch", .chart_seq$n)
}

## `variants` is a named list of alternative options the toolbar can swap in (counts vs percent, one dataset per sample, ...)
echart <- function(option, height = 340, title = NULL, caption = NULL,
                   controls = NULL, variants = NULL, id = NULL,
                   renderer = "canvas") {
  id <- id %||% new_chart_id()
  spec <- list(option = option, height = height, renderer = renderer)
  if (!is.null(variants)) spec$variants <- variants
  json <- toJSON(spec, auto_unbox = TRUE, null = "null", na = "null", digits = 6)
  ## `</script>` inside a JSON string would close the block early
  json <- gsub("</", "<\\/", json, fixed = TRUE)

  tags$figure(class = "fig",
    if (!is.null(title) || !is.null(controls))
      div(class = "fig-bar",
        if (!is.null(title)) span(class = "fig-title", title),
        controls),
    div(class = "chart", id = id, `data-chart` = paste0(id, "_data")),
    tags$script(type = "application/json", id = paste0(id, "_data"), HTML(json)),
    if (!is.null(caption)) tags$figcaption(caption))
}

## Toolbar controls. `chart_id` must match the `id` passed to echart()
switch_buttons <- function(chart_id, labels, keys, active = 1L) {
  div(class = "btn-group",
    lapply(seq_along(labels), function(i)
      tags$button(type = "button",
        class = paste("btn", if (i == active) "on" else ""),
        `data-fig-switch` = chart_id, `data-variant` = keys[[i]], labels[[i]])))
}

## ------------------------------------------------- ECharts option builders

## Common scaffolding; colours that depend on the theme are added in report.js
.base_opt <- function(legend = FALSE, grid = list(left = 8, right = 18, top = 34, bottom = 8)) {
  o <- list(grid = grid, tooltip = list(trigger = "axis", axisPointer = list(type = "shadow")))
  if (legend) o$legend <- list(top = 0, itemWidth = 11, itemHeight = 11, icon = "roundRect")
  o
}

## Horizontal stacked bar, one bar per sample. `mat` is a data.table with a `Sample` column and one numeric column per category
stacked_bar_option <- function(
  mat, categories, colours = FATE_COL,
  percent = FALSE, y_name = "Reads",
  zoom_above = 45L) {

  samples <- as.character(mat$Sample)
  tot <- Reduce(`+`, lapply(categories, function(k) num0(mat[[k]])))
  series <- lapply(categories, function(k) {
    v <- num0(mat[[k]])
    list(name = k, type = "bar", stack = "fate", barMaxWidth = 26,
         itemStyle = list(color = unname(colours[[k]] %||% "#8296ad")),
         emphasis = list(focus = "series"),
         data = if (percent) round(safe_pct(v, tot), 3) else v)
  })
  o <- .base_opt(legend = TRUE, grid = list(left = 8, right = 18, top = 46, bottom = if (length(samples) > zoom_above) 46 else 8))
  o$xAxis <- list(type = "category", data = samples,
                  axisLabel = list(rotate = if (length(samples) > 8) 45 else 0, hideOverlap = TRUE, fontSize = 10),
                  axisTick = list(alignWithLabel = TRUE))
  o$yAxis <- list(type = "value", min = 0,
                  name = if (percent) "% of demultiplexed reads" else y_name,
                  nameLocation = "end", nameGap = 12,
                  nameTextStyle = list(align = "left"),
                  max = if (percent) 100 else NULL,
                  axisLabel = list(formatter = if (percent) "{value}%" else js_fn("intValue")))
  o$series <- series
  o$tooltip <- list(trigger = "axis", axisPointer = list(type = "shadow"),
                    valueFormatter = if (percent) js_fn("pctValue", 2L) else js_fn("intValue"),
                    order = "seriesDesc")
  if (length(samples) > zoom_above) {
    span <- max(10, round(zoom_above / length(samples) * 100))
    o$dataZoom <- list(
      list(type = "slider", start = 0, end = span, height = 18, bottom = 8),
      list(type = "inside", start = 0, end = span))
  }
  o
}

## Simple vertical bar, one value per sample
bar_option <- function(labels, values, y_name = "Reads", colour = PAL[[1L]],
                       zoom_above = 45L, tooltip_suffix = "") {
  o <- .base_opt(grid = list(left = 8, right = 18, top = 34, bottom = if (length(labels) > zoom_above) 46 else 8))
  o$xAxis <- list(type = "category", data = as.character(labels),
                  axisLabel = list(rotate = if (length(labels) > 8) 45 else 0, hideOverlap = TRUE, fontSize = 10),
                  axisTick = list(alignWithLabel = TRUE))
  o$yAxis <- list(type = "value", min = 0, name = y_name,
                  nameLocation = "end", nameGap = 12,
                  nameTextStyle = list(align = "left"),
                  axisLabel = list(formatter = js_fn("intValue")))
  o$tooltip$valueFormatter <- js_fn("intValue")
  o$series <- list(list(type = "bar", barMaxWidth = 26, data = values,
                        itemStyle = list(color = colour, borderRadius = c(2, 2, 0, 0)),
                        emphasis = list(focus = "series")))
  if (length(labels) > zoom_above) {
    span <- max(10, round(zoom_above / length(labels) * 100))
    o$dataZoom <- list(
      list(type = "slider", start = 0, end = span, height = 18, bottom = 8),
      list(type = "inside", start = 0, end = span))
  }
  o
}

## Scatter with a per-point label shown on hover. `points` is a list of list(x, y, name, size)
scatter_option <- function(points, x_name, y_name, x_log = FALSE,
                           colour = PAL[[1L]], marklines = NULL, size_name = NULL) {
  data <- lapply(points, function(p) list(value = list(p$x, p$y, p$size %||% 8), name = p$name))
  o <- .base_opt(grid = list(left = 8, right = 22, top = 34, bottom = 10))
  o$xAxis <- list(type = if (x_log) "log" else "value", name = x_name,
                  nameLocation = "middle", nameGap = 30,
                  min = if (x_log) NULL else 0, scale = TRUE)
  o$yAxis <- list(type = "value", name = y_name, nameLocation = "end", nameGap = 12,
                  nameTextStyle = list(align = "left"), scale = TRUE)
  o$tooltip <- list(trigger = "item",
                    formatter = js_fn("scatterPoint", x_name, y_name, size_name %||% ""))
  o$series <- list(list(
    type = "scatter", data = data,
    symbolSize = js_fn("sqrtSize", 7, 26),
    itemStyle = list(color = colour, opacity = 0.78),
    emphasis = list(focus = "self", itemStyle = list(opacity = 1)),
    markLine = marklines))
  o
}

## A few ECharts fields take callbacks, which cannot survive JSON.parse
## Emit a placeholder that report.js resolves against its own registry of named formatters; see FNS in assets/report/report.js
js_fn <- function(name, ...) {
  args <- list(...)
  if (length(args) == 0L) list(`__fn` = name) else list(`__fn` = name, args = args)
}

## Histogram from raw values; binning happens in R so the page stays small
hist_option <- function(values, bins = 40L, x_name = "Value", y_name = "Count",
                        colour = PAL[[1L]], marklines = NULL, log_y = FALSE) {
  values <- values[is.finite(values)]
  if (length(values) == 0L) return(NULL)
  rng <- range(values)
  if (diff(rng) <= 0) rng <- c(rng[1L] - 0.5, rng[2L] + 0.5)
  brk <- seq(rng[1L], rng[2L], length.out = bins + 1L)
  h   <- hist(values, breaks = brk, plot = FALSE)
  mid <- signif(h$mids, 4)   # not round(): the values may all be very small
  o <- .base_opt(grid = list(left = 8, right = 18, top = 34, bottom = 10))
  o$xAxis <- list(type = "category", data = mid, name = x_name,
                  nameLocation = "middle", nameGap = 30,
                  axisLabel = list(hideOverlap = TRUE, fontSize = 10, interval = "auto"))
  o$yAxis <- list(type = if (log_y) "log" else "value", min = if (log_y) NULL else 0,
                  name = y_name, nameLocation = "end", nameGap = 12,
                  nameTextStyle = list(align = "left"))
  o$series <- list(list(type = "bar", data = h$counts, barCategoryGap = "0%",
                        itemStyle = list(color = colour), markLine = marklines))
  o
}

## hist_option renders a category axis of bin midpoints, so a cut-off line has to be placed by bin index on the same grid the histogram used
bin_markline <- function(values, bins, cutoff, label) {
  values <- values[is.finite(values)]
  if (!is.finite(cutoff) || length(values) == 0L) return(NULL)
  rng <- range(values)
  ## Outside the plotted range the line would be pinned to an edge bin and suggest a boundary that is not there
  if (cutoff < rng[1L] || cutoff > rng[2L]) return(NULL)
  if (diff(rng) <= 0) return(NULL)
  mids <- seq(rng[1L], rng[2L], length.out = bins + 1L)
  mids <- (mids[-1L] + mids[-length(mids)]) / 2
  idx  <- which.min(abs(mids - cutoff))
  ## An ECharts markLine label is centred on the line, so near either end it spills past the plot
  ## The caption carries the value, so drop the label rather than render it clipped
  edge <- max(1L, round(bins * 0.12))
  markline_x(idx - 1L, if (idx <= edge || idx > bins - edge) "" else label)
}

## A vertical reference line, e.g. a filtering cut-off
markline_x <- function(value, label) {
  list(silent = TRUE, symbol = "none",
       lineStyle = list(color = "#c9184a", type = "dashed", width = 1.4),
       label = list(formatter = label, position = "insideMiddleTop",
                    rotate = 0, fontSize = 10, padding = c(0, 0, 0, 6),
                    align = "left"),
       data = list(list(xAxis = value)))
}

markline_y <- function(value, label) {
  list(silent = TRUE, symbol = "none",
       lineStyle = list(color = "#c9184a", type = "dashed", width = 1.4),
       label = list(formatter = label, position = "insideEndTop",
                    rotate = 0, fontSize = 10),
       data = list(list(yAxis = value)))
}

## Box plot from a named list of numeric vectors
box_option <- function(groups, y_name = "Value", colour = PAL[[1L]]) {
  groups <- groups[vapply(groups, function(v) sum(is.finite(v)) > 0L, logical(1))]
  if (length(groups) == 0L) return(NULL)
  stats <- lapply(groups, function(v) {
    v <- v[is.finite(v)]
    q <- as.numeric(stats::quantile(v, c(0, 0.25, 0.5, 0.75, 1), names = FALSE))
    round(q, 4)
  })
  o <- .base_opt(grid = list(left = 8, right = 18, top = 34, bottom = 10))
  o$xAxis <- list(type = "category", data = names(groups),
                  axisLabel = list(hideOverlap = TRUE, fontSize = 10, rotate = if (length(groups) > 8) 30 else 0))
  o$yAxis <- list(type = "value", name = y_name, nameLocation = "end", nameGap = 12,
                  nameTextStyle = list(align = "left"), scale = TRUE)
  o$tooltip <- list(trigger = "item")
  o$series <- list(list(type = "boxplot", data = unname(stats),
                        itemStyle = list(color = colour, borderColor = colour, borderWidth = 1.4)))
  o
}

## Multi-line chart. `series` is a list of list(name, x, y)
line_option <- function(series, x_name, y_name, x_log = FALSE, y_log = FALSE,
                        show_legend = FALSE, smooth = TRUE) {
  o <- .base_opt(legend = show_legend, grid = list(left = 8, right = 22, top = if (show_legend) 44 else 34, bottom = 10))
  o$xAxis <- list(type = if (x_log) "log" else "value", name = x_name,
                  nameLocation = "middle", nameGap = 30, scale = TRUE,
                  min = if (x_log) NULL else 0)
  o$yAxis <- list(type = if (y_log) "log" else "value", name = y_name,
                  nameLocation = "end", nameGap = 12,
                  nameTextStyle = list(align = "left"), scale = TRUE,
                  min = if (y_log) NULL else 0)
  o$tooltip <- list(trigger = "item")
  o$series <- lapply(seq_along(series), function(i) {
    s <- series[[i]]
    list(name = s$name, type = "line", smooth = smooth, showSymbol = FALSE,
         lineStyle = list(width = 1.6, color = PAL[[(i - 1L) %% length(PAL) + 1L]]),
         itemStyle = list(color = PAL[[(i - 1L) %% length(PAL) + 1L]]),
         emphasis = list(focus = "series", lineStyle = list(width = 3)),
         data = Map(function(a, b) list(a, b), s$x, s$y))
  })
  o
}

## Sankey. `nodes` is a character vector of names (optionally with a colours lookup); `links` a data.table of source/target/value
sankey_option <- function(nodes, links, colours = FATE_COL, value_fmt = NULL) {
  node_data <- lapply(nodes, function(n)
    list(name = n, itemStyle = list(color = unname(colours[[n]] %||% PAL[[1L]]))))
  link_data <- lapply(seq_len(nrow(links)), function(i)
    list(source = links$source[[i]], target = links$target[[i]], value = links$value[[i]]))
  list(
    tooltip = list(trigger = "item", triggerOn = "mousemove"),
    series = list(list(
      type = "sankey", left = 8, right = 130, top = 12, bottom = 12,
      nodeGap = 14, nodeWidth = 13, draggable = FALSE,
      emphasis = list(focus = "adjacency"),
      label = list(position = "right", formatter = "{b}"),
      data = node_data, links = link_data)))
}

## Heatmap from a long data.table of x / y / value
heatmap_option <- function(dt, x_levels, y_levels, x_name, y_name,
                           value_name = "Count", log_scale = TRUE) {
  xi <- match(dt$x, x_levels) - 1L
  yi <- match(dt$y, y_levels) - 1L
  v  <- dt$value
  data <- Map(function(a, b, c) list(a, b, c), xi, yi, v)
  vmax <- max(v[is.finite(v)], na.rm = TRUE)
  list(
    grid = list(left = 8, right = 84, top = 14, bottom = 10, containLabel = TRUE),
    tooltip = list(trigger = "item"),
    xAxis = list(type = "category", data = x_levels, name = x_name,
                 nameLocation = "middle", nameGap = 30,
                 axisLabel = list(hideOverlap = TRUE, fontSize = 10, interval = "auto"),
                 splitArea = list(show = FALSE)),
    yAxis = list(type = "category", data = y_levels, name = y_name,
                 nameLocation = "middle", nameGap = 44,
                 axisLabel = list(hideOverlap = TRUE, fontSize = 10, interval = "auto"),
                 splitArea = list(show = FALSE)),
    ## Vertical and to the right, so it cannot sit on top of the x-axis label
    visualMap = list(min = 0, max = vmax, calculable = TRUE, orient = "vertical",
                     right = 6, top = "middle", itemHeight = 140, itemWidth = 12,
                     text = list(value_name, ""),
                     inRange = list(color = c("#eef6f0", "#a9d8bb", "#3f8f63", "#0c3d26"))),
    series = list(list(type = "heatmap", data = data, progressive = 2000,
                       emphasis = list(itemStyle = list(borderColor = "#1c2938", borderWidth = 1)))))
}

## ----------------------------------------------------------- copy blocks

## A bordered block of text with a Copy button in its header. `id` must be unique within the page: report.js copies that element's textContent
copy_block <- function(id, title, content, pre_class = NULL) {
  div(class = "cmd",
    div(class = "cmd-bar",
      span(class = "cmd-title", title),
      span(class = "copied"),
      tags$button(type = "button", class = "btn", `data-copy` = id, "Copy")),
    tags$pre(id = id, class = pre_class, content))
}

## -------------------------------------------------------- command panel

## The exact `nextflow run` invocation, so a reader can reproduce the run
## Nextflow writes it to pipeline_info/execution_command.txt
command_panel <- function(path, extra = NULL) {
  cmd <- if (is_usable_file(path)) {
    trimws(paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n"))
  } else NULL
  if (is.null(cmd) || !nzchar(cmd)) {
    return(p(class = "empty", "The run command was not recorded."))
  }
  tagList(
    copy_block("run-command", "Run command", cmd),
    if (!is.null(extra) && length(extra) > 0L)
      tags$dl(class = "kv", unlist(lapply(names(extra), function(k)
        list(tags$dt(k), tags$dd(extra[[k]]))), recursive = FALSE)))
}

## --------------------------------------------------- methods / versions

## document_s1.R / document_s2.R write a plain-text file with `Methods:` and `References:` headers; split it so the citations become a real list
methods_panel <- function(txt) {
  if (is.null(txt)) return(NULL)
  lines <- strsplit(txt, "\n", fixed = TRUE)[[1L]]
  ref_at <- grep("^\\s*References:\\s*$", lines)
  if (length(ref_at) == 0L) {
    return(copy_block("methods-text", "Methods", txt, pre_class = "methods"))
  }
  body <- lines[seq_len(ref_at[[1L]] - 1L)]
  body <- body[!grepl("^\\s*Methods:\\s*$", body)]
  refs <- lines[seq(ref_at[[1L]] + 1L, length(lines))]
  refs <- trimws(sub("^\\s*-\\s*", "", refs))
  refs <- refs[nzchar(refs)]
  tagList(
    copy_block("methods-text", "Methods",
               trimws(paste(body, collapse = "\n")), pre_class = "methods"),
    if (length(refs) > 0L) tagList(
      h3(class = "sub", "References"),
      tags$ol(class = "refs", lapply(refs, tags$li))))
}

## software_versions.yml is keyed by process, so render one group per stage
versions_panel <- function(versions) {
  if (is.null(versions) || !is.list(versions)) return(NULL)
  keys <- setdiff(names(versions), c("NextITS", "Nextflow"))
  if (length(keys) == 0L) return(NULL)
  rows <- rbindlist(lapply(keys, function(k) {
    v <- versions[[k]]
    if (!is.list(v) || length(v) == 0L) return(NULL)
    data.table(Stage = k, Tool = names(v), Version = vapply(v, as.character, ""))
  }), fill = TRUE)
  if (is.null(rows) || nrow(rows) == 0L) return(NULL)
  dt_table(rows, labels = c("Pipeline stage", "Tool", "Version"),
           search = TRUE, cols_menu = FALSE, download = "software_versions.tsv")
}

## ------------------------------------------------------ pipeline_info links

## Nextflow writes its execution report, timeline, trace and DAG only after the last process finishes, so this process cannot see them
## Their names are known in advance though: nextflow.config stamps one `trace_timestamp` into all four, and that value is handed to the report
pipeline_links <- function(rel = "../pipeline_info", stamp = NULL, dir = "pipeline_info") {
  named <- list(
    c("execution_report_%s.html",   "Execution report"),
    c("execution_timeline_%s.html", "Timeline"),
    c("execution_trace_%s.txt",     "Trace"),
    c("pipeline_dag_%s.svg",        "Workflow DAG"))

  items <- list()
  if (!is.null(stamp) && nzchar(stamp) && !identical(stamp, "null")) {
    items <- lapply(named, function(w)
      tags$a(href = paste0(rel, "/", sprintf(w[[1L]], stamp)), w[[2L]]))
  } else if (dir.exists(dir)) {
    ## Standalone render next to a results tree: glob instead.
    items <- Filter(Negate(is.null), lapply(named, function(w) {
      hit <- sort(Sys.glob(file.path(dir, sprintf(w[[1L]], "*"))), decreasing = TRUE)
      if (length(hit) == 0L) return(NULL)
      tags$a(href = paste0(rel, "/", basename(hit[[1L]])), w[[2L]])
    }))
  }

  items <- c(
    list(tags$a(href = paste0(rel, "/"), "Open pipeline_info")),
    items,
    list(tags$a(href = paste0(rel, "/pipeline_params.tsv"), "Parameters"),
         tags$a(href = paste0(rel, "/software_versions.yml"), "Software versions")))

  tagList(
    p(class = "hint", paste0(
      "Nextflow execution artefacts, relative to this file inside the results ",
      "directory. The links resolve only where the report sits in its published location.")),
    div(class = "links", items))
}

## --------------------------------------------------------------- page

## A section registers itself in the sidebar. `subs` adds second-level nav entries; each must match an id used inside `body`
sec <- function(id, title, ..., note = NULL, subs = NULL) {
  body <- Filter(Negate(is.null), list(...))
  if (length(body) == 0L) return(NULL)
  list(id = id, title = title, subs = subs,
       tag = tags$section(class = "panel", id = id,
         h2(title),
         if (!is.null(note)) p(class = "panel-note", note),
         body))
}

subhead <- function(id, title) h3(class = "sub", id = id, title)

read_asset <- function(assets_dir, name) {
  path <- file.path(assets_dir, name)
  if (!file.exists(path)) stop("Missing report asset: ", path, call. = FALSE)
  paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
}

## The logo goes in as a data URI rather than inline markup: the SVG carries its own `.fil0`-style class names, which would leak into the page stylesheet
logo_tag <- function(path) {
  if (!is_usable_file(path)) return(NULL)
  mime <- if (grepl("\\.svg$", path, ignore.case = TRUE)) "image/svg+xml" else "image/png"
  raw  <- readBin(path, "raw", file.info(path)$size)
  tags$img(src = paste0("data:", mime, ";base64,", jsonlite::base64_enc(raw)),
           alt = "NextITS")
}

## Built as one HTML string: htmltools indents tag children onto separate lines, and the resulting newline renders as a space before the period
default_footer <- function() {
  HTML(paste0(
    'Generated by <a href="https://github.com/vmikk/NextITS">NextITS</a>. ',
    'Charts use <a href="https://echarts.apache.org/">Apache ECharts</a>'))
}

## meta: named character vector rendered as the header strapline.
report_page <- function(title, meta, sections, assets_dir, out,
                        footer = NULL, logo = NULL) {
  sections <- Filter(Negate(is.null), sections)

  nav <- tags$nav(unlist(lapply(sections, function(s) {
    c(list(tags$a(href = paste0("#", s$id), s$title)),
      lapply(s$subs %||% list(), function(x)
        tags$a(class = "sub", href = paste0("#", x[[1L]]), x[[2L]])))
  }), recursive = FALSE))

  ## htmltools pulls <head> tags out of a tree for its dependency machinery and drops them from as.character(), so the head is assembled as text instead
  head_html <- paste0(
    "<head>\n",
    '<meta charset="utf-8">\n',
    '<meta name="viewport" content="width=device-width, initial-scale=1">\n',
    "<title>", htmlEscape(title), "</title>\n",
    "<style>\n", read_asset(assets_dir, "report.css"), "\n</style>\n",
    "</head>")

  body <- tags$body(
    tags$button(type = "button", class = "btn theme-toggle", "☾ Dark"),
    div(class = "layout",
      tags$aside(class = "sidebar",
        div(class = "brand", logo_tag(logo) %||% div(class = "brand-text", "NextITS")),
        div(class = "brand-sub", title),
        nav),
      tags$main(class = "main",
        div(class = "page-head",
          h1(title),
          div(class = "meta", lapply(names(meta), function(k)
            span(k, " ", tags$b(meta[[k]]))))),
        lapply(sections, function(s) s$tag),
        tags$footer(class = "page-foot", footer %||% default_footer()))),
    tags$script(HTML(read_asset(assets_dir, "echarts.min.js"))),
    tags$script(HTML(read_asset(assets_dir, "report.js"))))

  writeLines(
    paste0('<!DOCTYPE html>\n<html lang="en" class="no-js">\n',
           head_html, "\n", as.character(body), "\n</html>"),
    out, useBytes = TRUE)
  invisible(out)
}
