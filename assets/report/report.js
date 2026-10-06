/* NextITS run report - nav, theme, interactive tables, ECharts lifecycle
   All state lives in the DOM; there is no build step and no framework */
(function () {
  "use strict";

  var root = document.documentElement;
  root.classList.remove("no-js");

  var $  = function (sel, ctx) { return (ctx || document).querySelector(sel); };
  var $$ = function (sel, ctx) { return Array.prototype.slice.call((ctx || document).querySelectorAll(sel)); };

  /* ------------------------------------------------------------- theme */

  /* The <head> script in report_lib.R sets data-theme before first paint,
     using the same key; this keeps it in sync and handles the toggle */
  var THEME_KEY = "nextits-report-theme";
  var darkQuery = window.matchMedia ? window.matchMedia("(prefers-color-scheme: dark)") : null;

  function storedTheme() {
    try { return localStorage.getItem(THEME_KEY); } catch (e) { return null; }
  }
  function isDark() {
    var t = storedTheme();
    return t ? t === "dark" : !!(darkQuery && darkQuery.matches);
  }
  function themeIsDark() { return root.getAttribute("data-theme") === "dark"; }
  function reducedMotion() {
    return !!(window.matchMedia && window.matchMedia("(prefers-reduced-motion: reduce)").matches);
  }

  function setThemeAttr(dark) {
    root.setAttribute("data-theme", dark ? "dark" : "light");
    var btn = $(".theme-toggle");
    if (btn) {
      var label = dark ? "Switch to light theme" : "Switch to dark theme";
      btn.setAttribute("aria-label", label);
      btn.setAttribute("title", label);
    }
  }
  function applyTheme(dark) {
    setThemeAttr(dark);
    redrawCharts();
  }

  function initTheme() {
    setThemeAttr(isDark());
    var btn = $(".theme-toggle");
    if (btn) {
      btn.addEventListener("click", function () {
        var next = !themeIsDark();
        try { localStorage.setItem(THEME_KEY, next ? "dark" : "light"); } catch (e) { /* private mode */ }
        if (reducedMotion()) { applyTheme(next); return; }
        /* The page colours cross-fade (CSS) while the charts fade out; they are
           repainted in the new theme once hidden, then fade back in.
           The icon spins out and back in on the same beat */
        root.classList.add("theme-xfade");
        btn.classList.add("spin");
        setThemeAttr(next);
        setTimeout(function () {
          redrawCharts();
          btn.classList.remove("spin");
          root.classList.remove("theme-xfade");
        }, 170);
      });
    }
    /* Follow the OS setting until the reader picks a theme */
    if (darkQuery) {
      var follow = function (e) { if (!storedTheme()) applyTheme(e.matches); };
      if (darkQuery.addEventListener) darkQuery.addEventListener("change", follow);
      else if (darkQuery.addListener) darkQuery.addListener(follow);
    }
  }

  /* --------------------------------------------------------------- nav */

  function initNav() {
    var links = $$(".sidebar nav a");
    if (!links.length) return;
    var targets = links
      .map(function (a) { return { a: a, el: document.getElementById(a.hash.slice(1)) }; })
      .filter(function (t) { return t.el; });
    if (!targets.length) return;

    function mark(el) {
      links.forEach(function (a) {
        var on = a.hash === "#" + el.id;
        a.classList.toggle("active", on);
        if (on) a.setAttribute("aria-current", "true"); else a.removeAttribute("aria-current");
      });
    }

    if ("IntersectionObserver" in window) {
      var visible = {};
      var io = new IntersectionObserver(function (entries) {
        entries.forEach(function (e) { visible[e.target.id] = e.isIntersecting; });
        for (var i = 0; i < targets.length; i++) {
          if (visible[targets[i].el.id]) { mark(targets[i].el); break; }
        }
      }, { rootMargin: "-8% 0px -70% 0px", threshold: 0 });
      targets.forEach(function (t) { io.observe(t.el); });
    }
    links.forEach(function (a) {
      a.addEventListener("click", function () { setTimeout(function () { mark(document.getElementById(a.hash.slice(1)) || {}); }, 0); });
    });
  }

  /* ------------------------------------------------------------ tables */

  function cellValue(td) {
    var v = td.getAttribute("data-v");
    if (v === null) return td.textContent.trim().toLowerCase();
    if (v === "") return null;
    var n = parseFloat(v);
    return isNaN(n) ? v.toLowerCase() : n;
  }

  function initTable(wrap) {
    var table = $("table.dt", wrap);
    if (!table) return;
    var tbody = $("tbody", table);
    var rows = $$("tr", tbody);
    var headers = $$("thead th", table);
    var counter = $(".count", wrap);
    var scroller = $(".tbl-scroll", wrap);
    var total = rows.length;

    /* Read the labels before the sort icons go in: they feed the column
       menu and the TSV header */
    var labels = headers.map(function (th) { return th.textContent.trim(); });

    function updateCount(n) {
      if (counter) counter.textContent = n === total ? total + " rows" : n + " of " + total + " rows";
    }
    updateCount(total);

    /* sort ------------------------------------------------------------ */
    headers.forEach(function (th, idx) {
      if (th.hasAttribute("data-nosort")) return;
      var arrow = document.createElement("span");
      arrow.className = "arrow";
      arrow.setAttribute("aria-hidden", "true");
      th.appendChild(arrow);
      th.tabIndex = 0;

      function sortBy() {
        var asc = th.getAttribute("aria-sort") !== "ascending";
        headers.forEach(function (h) { h.removeAttribute("aria-sort"); });
        th.setAttribute("aria-sort", asc ? "ascending" : "descending");

        rows.sort(function (ra, rb) {
          var a = cellValue(ra.cells[idx]), b = cellValue(rb.cells[idx]);
          if (a === null && b === null) return 0;
          if (a === null) return 1;          /* blanks always last */
          if (b === null) return -1;
          if (a < b) return asc ? -1 : 1;
          if (a > b) return asc ? 1 : -1;
          return 0;
        });
        var frag = document.createDocumentFragment();
        rows.forEach(function (r) { frag.appendChild(r); });
        tbody.appendChild(frag);
      }
      th.addEventListener("click", sortBy);
      th.addEventListener("keydown", function (e) {
        if (e.key === "Enter" || e.key === " ") { e.preventDefault(); sortBy(); }
      });
    });

    /* Lift the sticky header off the rows once they scroll under it */
    if (scroller) {
      scroller.addEventListener("scroll", function () {
        scroller.classList.toggle("scrolled", scroller.scrollTop > 0);
      }, { passive: true });
    }

    /* search ---------------------------------------------------------- */
    var search = $("input[type=search]", wrap);
    if (search) {
      var filter = function () {
        var q = search.value.trim().toLowerCase();
        var shown = 0;
        rows.forEach(function (r) {
          var hit = !q || r.textContent.toLowerCase().indexOf(q) !== -1;
          r.hidden = !hit;
          if (hit) shown++;
        });
        updateCount(shown);
      };
      search.addEventListener("input", filter);
      search.addEventListener("keydown", function (e) {
        if (e.key === "Escape" && search.value) { e.preventDefault(); search.value = ""; filter(); }
      });
    }

    /* column show / hide ---------------------------------------------- */
    var menu = $(".cols-menu", wrap);
    if (menu) {
      var list = $(".cols-list", menu);
      var toggle = $(".btn", menu);
      toggle.setAttribute("aria-haspopup", "true");
      toggle.setAttribute("aria-expanded", "false");
      headers.forEach(function (th, idx) {
        var label = document.createElement("label");
        var cb = document.createElement("input");
        cb.type = "checkbox";
        cb.checked = true;
        cb.addEventListener("change", function () {
          th.hidden = !cb.checked;
          rows.forEach(function (r) { if (r.cells[idx]) r.cells[idx].hidden = !cb.checked; });
        });
        label.appendChild(cb);
        label.appendChild(document.createTextNode(labels[idx]));
        list.appendChild(label);
      });
      var setOpen = function (open) {
        list.classList.toggle("open", open);
        toggle.setAttribute("aria-expanded", open ? "true" : "false");
      };
      toggle.addEventListener("click", function (e) { e.stopPropagation(); setOpen(!list.classList.contains("open")); });
      document.addEventListener("click", function (e) { if (!menu.contains(e.target)) setOpen(false); });
      menu.addEventListener("keydown", function (e) {
        if (e.key === "Escape" && list.classList.contains("open")) { setOpen(false); toggle.focus(); }
      });
    }

    /* TSV download ----------------------------------------------------- */
    var dl = $("[data-download]", wrap);
    if (dl) {
      dl.addEventListener("click", function () {
        var lines = [labels.join("\t")];
        rows.forEach(function (r) {
          if (r.hidden) return;
          lines.push(Array.prototype.map.call(r.cells, function (td) {
            var v = td.getAttribute("data-v");
            return (v === null ? td.textContent.trim() : v).replace(/[\t\r\n]/g, " ");
          }).join("\t"));
        });
        var blob = new Blob([lines.join("\n") + "\n"], { type: "text/tab-separated-values" });
        var url = URL.createObjectURL(blob);
        var a = document.createElement("a");
        a.href = url;
        a.download = dl.getAttribute("data-download");
        document.body.appendChild(a);
        a.click();
        document.body.removeChild(a);
        setTimeout(function () { URL.revokeObjectURL(url); }, 0);
      });
    }
  }

  /* ------------------------------------------------------------ charts */

  var charts = [];   /* {el, spec, inst} */

  function fmtNum(v) {
    if (v === null || v === undefined || v === "") return "—";
    var n = +v;
    if (!isFinite(n)) return String(v);
    return Math.abs(n - Math.round(n)) < 1e-9
      ? Math.round(n).toLocaleString()
      : n.toLocaleString(undefined, { maximumFractionDigits: 3 });
  }

  /* ECharts takes callbacks for a handful of fields, but a report embeds its
     options as JSON, and JSON has no functions. The R side emits
     {"__fn": "name", "args": [...]} placeholders which are swapped for real
     functions here */
  var FNS = {
    intValue: function () {
      return function (v) { return fmtNum(v); };
    },
    pctValue: function (digits) {
      return function (v) {
        return (v === null || v === undefined || v === "") ? "—" : (+v).toFixed(digits || 0) + "%";
      };
    },
    sqrtSize: function (lo, hi) {
      return function (v) {
        return Math.max(lo, Math.min(hi, Math.sqrt(v[2] || 1)));
      };
    },
    scatterPoint: function (xl, yl, zl) {
      return function (p) {
        var v = p.value || [];
        var s = "<b>" + (p.name || p.seriesName || "") + "</b>";
        s += "<br/>" + xl + ": " + fmtNum(v[0]);
        s += "<br/>" + yl + ": " + fmtNum(v[1]);
        if (zl && v.length > 2) s += "<br/>" + zl + ": " + fmtNum(v[2]);
        return s;
      };
    }
  };

  /* Colours that follow the theme arrive from R as "var(--token)" strings and
     are looked up in the active theme's CSS custom properties.
     One pass rebuilds the option tree, swapping both those strings and the
     function placeholders, so the spec stored for each chart is never mutated */
  var TOKEN = /^var\((--[\w-]+)\)$/;

  function tokenReader() {
    var cs = getComputedStyle(root);
    var cache = {};
    return function (name) {
      if (Object.prototype.hasOwnProperty.call(cache, name)) return cache[name];
      var v = cs.getPropertyValue(name).trim();
      if (!v) {
        if (window.console) console.warn("NextITS report: unknown colour token " + name);
        v = cs.getPropertyValue("--ink-3").trim() || "#888888";
      }
      cache[name] = v;
      return v;
    };
  }

  function resolve(o, tok) {
    if (typeof o === "string") {
      var m = TOKEN.exec(o);
      return m ? tok(m[1]) : o;
    }
    if (Array.isArray(o)) {
      var arr = new Array(o.length);
      for (var i = 0; i < o.length; i++) arr[i] = resolve(o[i], tok);
      return arr;
    }
    if (o && typeof o === "object") {
      if (typeof o.__fn === "string" && FNS[o.__fn]) {
        return FNS[o.__fn].apply(null, o.args || []);
      }
      var out = {};
      for (var k in o) {
        if (Object.prototype.hasOwnProperty.call(o, k)) out[k] = resolve(o[k], tok);
      }
      return out;
    }
    return o;
  }

  function themeColors(tok) {
    var dark = themeIsDark();
    var viz = [];
    for (var i = 1; i <= 8; i++) viz.push(tok("--viz-" + i));
    return {
      dark: dark,
      font: getComputedStyle(document.body).fontFamily,
      ink: tok("--ink"), ink2: tok("--ink-2"), ink3: tok("--ink-3"),
      line: tok("--line"), line2: tok("--line-2"),
      surface: tok("--surface"), surface2: tok("--surface-2"),
      accent: tok("--accent"), accentStrong: tok("--accent-strong"),
      grid: tok("--viz-grid"), axis: tok("--viz-axis"),
      zoomFill: tok("--zoom-fill"), tooltipBg: tok("--tooltip-bg"),
      viz: viz,
      onAccent: dark ? "#0c0f0d" : "#ffffff",
      pointer: dark ? "rgba(255,255,255,0.05)" : "rgba(25,32,27,0.05)",
      tipShadow: dark ? "0 12px 32px -8px rgba(0,0,0,.65)" : "0 12px 32px -10px rgba(25,32,27,.28)"
    };
  }

  /* Stamp the theme's chrome onto a resolved option tree.
     The R side emits data, layout and colour tokens; everything else that
     must follow the theme is set here */
  function themeOption(o, c) {
    o.backgroundColor = "transparent";
    o.color = c.viz;
    o.textStyle = Object.assign({ color: c.ink2, fontFamily: c.font }, o.textStyle || {});
    o.animation = false;

    /* Every chart gets the same toolbox, so readers can lift a figure out of
       the report or inspect the numbers behind it without leaving the page */
    o.toolbox = Object.assign({
      right: 10,
      top: 4,
      itemSize: 14,
      itemGap: 10,
      iconStyle: { borderColor: c.ink3, borderWidth: 1.2 },
      emphasis: { iconStyle: { borderColor: c.accent, textFill: c.ink2, textBackgroundColor: c.surface, textPadding: [2, 4] } },
      feature: {
        saveAsImage: { title: "Save as PNG", pixelRatio: 2, backgroundColor: c.surface },
        dataView: {
          title: "View data", readOnly: true, lang: ["Chart data", "Close", "Refresh"],
          backgroundColor: c.surface, textColor: c.ink,
          textareaColor: c.surface2, textareaBorderColor: c.line,
          buttonColor: c.accent, buttonTextColor: c.onAccent
        },
        restore: { title: "Reset" }
      }
    }, o.toolbox || {});

    if (o.title) {
      o.title = [].concat(o.title).map(function (t) {
        t.textStyle = Object.assign({ color: c.ink, fontWeight: 600, fontSize: 13 }, t.textStyle || {});
        t.subtextStyle = Object.assign({ color: c.ink3, fontSize: 11 }, t.subtextStyle || {});
        return t;
      });
    }
    if (o.legend) {
      o.legend = [].concat(o.legend).map(function (l) {
        l.textStyle = Object.assign({ color: c.ink2, fontSize: 12 }, l.textStyle || {});
        l.inactiveColor = c.line2;
        l.inactiveBorderColor = c.line2;
        if (l.itemGap === undefined) l.itemGap = 14;
        l.pageIconColor = c.ink2;
        l.pageIconInactiveColor = c.line2;
        l.pageTextStyle = { color: c.ink3 };
        return l;
      });
    }
    if (o.tooltip) {
      o.tooltip = [].concat(o.tooltip).map(function (t) {
        t.backgroundColor = c.tooltipBg;
        t.borderColor = c.line;
        t.borderWidth = 1;
        t.padding = [8, 11];
        t.confine = true;
        t.textStyle = Object.assign({ color: c.ink, fontSize: 12, fontFamily: c.font }, t.textStyle || {});
        t.extraCssText = "border-radius:10px;box-shadow:" + c.tipShadow +
          ";backdrop-filter:blur(8px);-webkit-backdrop-filter:blur(8px);";
        if (t.axisPointer && t.axisPointer.type === "shadow") {
          t.axisPointer.shadowStyle = Object.assign({ color: c.pointer }, t.axisPointer.shadowStyle || {});
        }
        return t;
      });
    }
    ["xAxis", "yAxis"].forEach(function (k) {
      if (!o[k]) return;
      o[k] = [].concat(o[k]).map(function (ax) {
        var cat = ax.type === "category";
        ax.axisLine  = Object.assign({ lineStyle: { color: c.axis } }, ax.axisLine || {});
        ax.axisTick  = Object.assign({ show: cat, length: 4, lineStyle: { color: c.axis } }, ax.axisTick || {});
        ax.axisLabel = Object.assign({ color: c.ink3, fontSize: 11 }, ax.axisLabel || {});
        ax.nameTextStyle = Object.assign({ color: c.ink3, fontSize: 11, fontWeight: 500 }, ax.nameTextStyle || {});
        ax.splitLine = Object.assign({ lineStyle: { color: c.grid, type: "solid", width: 1 } }, ax.splitLine || {});
        return ax;
      });
    });
    if (o.dataZoom) {
      o.dataZoom = [].concat(o.dataZoom).map(function (dz) {
        if (dz.type === "slider") {
          dz.borderColor = c.line;
          dz.backgroundColor = "transparent";
          dz.fillerColor = c.zoomFill;
          dz.textStyle = { color: c.ink3, fontSize: 10 };
          dz.handleStyle = { color: c.surface, borderColor: c.line2 };
          dz.moveHandleStyle = { color: c.line2, opacity: 1 };
          dz.dataBackground = { lineStyle: { color: c.line2, width: 1 }, areaStyle: { color: c.surface2, opacity: 1 } };
          dz.selectedDataBackground = { lineStyle: { color: c.accentStrong, width: 1 }, areaStyle: { color: c.zoomFill, opacity: 1 } };
          dz.emphasis = { handleStyle: { borderColor: c.accentStrong }, moveHandleStyle: { color: c.accentStrong } };
        }
        return dz;
      });
    }
    if (o.visualMap) {
      o.visualMap = [].concat(o.visualMap).map(function (vm) {
        vm.textStyle = Object.assign({ color: c.ink2, fontSize: 11 }, vm.textStyle || {});
        vm.handleStyle = Object.assign({ borderColor: c.ink3 }, vm.handleStyle || {});
        return vm;
      });
    }
    if (o.series) {
      o.series = [].concat(o.series).map(function (s) {
        if (s.type === "bar") {
          s.blur = Object.assign({ itemStyle: { opacity: 0.28 } }, s.blur || {});
        }
        if (s.type === "line") {
          /* Lines drawn without symbols only react to hover with this on */
          s.triggerLineEvent = true;
          s.blur = Object.assign({ lineStyle: { opacity: 0.15 }, areaStyle: { opacity: 0.05 } }, s.blur || {});
        }
        if (s.type === "scatter") {
          s.itemStyle = Object.assign({ borderColor: c.surface, borderWidth: 1.5 }, s.itemStyle || {});
        }
        if (s.type === "sankey") {
          s.label = Object.assign({ color: c.ink, fontSize: 12, fontWeight: 500 }, s.label || {});
          s.lineStyle = Object.assign({ color: "gradient", opacity: c.dark ? 0.26 : 0.22, curveness: 0.5 }, s.lineStyle || {});
          s.itemStyle = Object.assign({ borderColor: c.surface, borderWidth: 1 }, s.itemStyle || {});
          s.emphasis = Object.assign({}, s.emphasis || {}, { lineStyle: { opacity: 0.55 } });
        }
        if (s.markLine) {
          s.markLine.label = Object.assign({
            color: c.ink2, fontSize: 10.5,
            backgroundColor: c.surface, padding: [2, 5], borderRadius: 4
          }, s.markLine.label || {});
        }
        return s;
      });
    }
    /* Keep axis labels and axis names inside the grid's own rectangle.
       `outerBoundsContain: "all"` also reserves room for the axis names, which
       the deprecated `containLabel` never did - without it the x-axis name is
       drawn past the bottom of the canvas and never seen */
    if (o.grid) {
      o.grid = [].concat(o.grid).map(function (g) {
        return Object.assign({ outerBoundsMode: "same", outerBoundsContain: "all" }, g);
      });
    }
    return o;
  }

  /* opts.animate: morph from the previous render (variant switches);
     opts.print: draw everything in one pass, so nothing is left for a later frame */
  function renderChart(rec, opts) {
    if (!window.echarts) return;
    opts = opts || {};
    if (!rec.inst) {
      rec.inst = echarts.init(rec.el, null, { renderer: rec.spec.renderer || "canvas" });
    }
    var tok = tokenReader();
    var o = themeOption(resolve(rec.spec.option, tok), themeColors(tok));
    if (opts.animate && !reducedMotion()) {
      o.animation = true;
      o.animationDuration = 0;
      o.animationDurationUpdate = 350;
      o.animationEasingUpdate = "cubicOut";
    }
    if (opts.print && o.series) {
      o.series.forEach(function (s) { s.progressive = 0; });
    }
    rec.inst.setOption(o, true);
  }

  function redrawCharts(opts) {
    charts.forEach(function (rec) { if (rec.inst) renderChart(rec, opts); });
  }

  function initCharts() {
    if (!window.echarts) return;

    $$("[data-chart]").forEach(function (el) {
      var node = document.getElementById(el.getAttribute("data-chart"));
      if (!node) return;
      var spec;
      try { spec = JSON.parse(node.textContent); } catch (e) { return; }
      el.style.height = (spec.height || 340) + "px";
      charts.push({ el: el, spec: spec, inst: null });
    });

    /* Render everything up front
       A report carries on the order of ten small charts,
       and ECharts' canvas renderer handles those in a few hundred milliseconds;
       deferring them only makes printing, PDF export and screenshotting miss panels
       Heavy series carry `large`/`progressive` instead */
    charts.forEach(function (rec) { renderChart(rec); });

    var t = null;
    window.addEventListener("resize", function () {
      clearTimeout(t);
      t = setTimeout(function () {
        charts.forEach(function (rec) { if (rec.inst) rec.inst.resize(); });
      }, 120);
    });

    /* Print on the light theme, at the printed width.
       `html.printing` applies the print layout on screen for the moment the
       charts are re-measured; the canvases do not reflow on their own */
    var printDark = null;
    window.addEventListener("beforeprint", function () {
      printDark = themeIsDark();
      root.classList.add("printing");
      if (printDark) setThemeAttr(false);
      charts.forEach(function (rec) {
        if (!rec.inst) return;
        rec.inst.resize();
        renderChart(rec, { print: true });
      });
    });
    window.addEventListener("afterprint", function () {
      root.classList.remove("printing");
      if (printDark) setThemeAttr(true);
      printDark = null;
      charts.forEach(function (rec) {
        if (!rec.inst) return;
        rec.inst.resize();
        renderChart(rec);
      });
    });
  }

  /* ------------------------------------------------- per-figure switches */

  /* A figure may declare alternative datasets in its spec (`variants`),
     e.g., counts vs percent; the toolbar buttons just swap which one is active */
  function recFor(id) {
    return charts.filter(function (r) { return r.el.id === id; })[0];
  }

  function initFigureControls() {
    $$("[data-fig-switch]").forEach(function (btn) {
      btn.setAttribute("aria-pressed", btn.classList.contains("on") ? "true" : "false");
      btn.addEventListener("click", function () {
        var rec = recFor(btn.getAttribute("data-fig-switch"));
        if (!rec || !rec.spec.variants) return;
        var key = btn.getAttribute("data-variant");
        var variant = rec.spec.variants[key];
        if (!variant || rec.spec.option === variant) return;
        $$('[data-fig-switch="' + btn.getAttribute("data-fig-switch") + '"]').forEach(function (b) {
          b.classList.toggle("on", b === btn);
          b.setAttribute("aria-pressed", b === btn ? "true" : "false");
        });
        rec.spec.option = variant;
        renderChart(rec, { animate: true });
      });
    });
  }

  /* ---------------------------------------------------------- copy button */

  function initCopy() {
    $$("[data-copy]").forEach(function (btn) {
      var timer = null;
      btn.addEventListener("click", function () {
        var src = document.getElementById(btn.getAttribute("data-copy"));
        if (!src) return;
        var text = src.textContent;
        var done = function (ok) {
          var note = btn.parentNode.querySelector(".copied");
          if (note) note.textContent = ok ? "Copied" : "Press Ctrl+C";
          btn.classList.toggle("is-copied", ok);
          clearTimeout(timer);
          timer = setTimeout(function () {
            if (note) note.textContent = "";
            btn.classList.remove("is-copied");
          }, 1800);
        };
        /* navigator.clipboard is unavailable on file:// in some browsers */
        if (navigator.clipboard && navigator.clipboard.writeText) {
          navigator.clipboard.writeText(text).then(function () { done(true); }, function () { selectText(src); done(false); });
        } else {
          selectText(src);
          done(false);
        }
      });
    });
  }

  function selectText(el) {
    var range = document.createRange();
    range.selectNodeContents(el);
    var sel = window.getSelection();
    sel.removeAllRanges();
    sel.addRange(range);
  }

  /* --------------------------------------------------------------- boot */

  function boot() {
    initTheme();
    initNav();
    $$(".tbl-wrap").forEach(initTable);
    initCharts();
    initFigureControls();
    initCopy();
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", boot);
  } else {
    boot();
  }
})();
