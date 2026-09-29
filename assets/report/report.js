/* NextITS run report - nav, theme, interactive tables, ECharts lifecycle
   All state lives in the DOM; there is no build step and no framework */
(function () {
  "use strict";

  document.documentElement.classList.remove("no-js");

  var $  = function (sel, root) { return (root || document).querySelector(sel); };
  var $$ = function (sel, root) { return Array.prototype.slice.call((root || document).querySelectorAll(sel)); };

  /* ------------------------------------------------------------- theme */

  var THEME_KEY = "nextits-report-theme";

  function systemDark() {
    return window.matchMedia && window.matchMedia("(prefers-color-scheme: dark)").matches;
  }
  function storedTheme() {
    try { return localStorage.getItem(THEME_KEY); } catch (e) { return null; }
  }
  function isDark() {
    var t = storedTheme();
    return t ? t === "dark" : systemDark();
  }
  function applyTheme(dark) {
    document.documentElement.setAttribute("data-theme", dark ? "dark" : "light");
    var btn = $(".theme-toggle");
    if (btn) {
      btn.textContent = dark ? "☀ Light" : "☾ Dark";
      btn.setAttribute("aria-label", dark ? "Switch to light theme" : "Switch to dark theme");
    }
    redrawCharts();
  }

  function initTheme() {
    applyTheme(isDark());
    var btn = $(".theme-toggle");
    if (!btn) return;
    btn.addEventListener("click", function () {
      var next = !isDark();
      try { localStorage.setItem(THEME_KEY, next ? "dark" : "light"); } catch (e) { /* private mode */ }
      applyTheme(next);
    });
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
      links.forEach(function (a) { a.classList.toggle("active", a.hash === "#" + el.id); });
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
    var total = rows.length;

    function updateCount(n) {
      if (counter) counter.textContent = n === total ? total + " rows" : n + " of " + total + " rows";
    }
    updateCount(total);

    /* sort ------------------------------------------------------------ */
    headers.forEach(function (th, idx) {
      if (th.hasAttribute("data-nosort")) return;
      th.insertAdjacentHTML("beforeend", ' <span class="arrow">▴▾</span>');
      th.addEventListener("click", function () {
        var asc = th.getAttribute("aria-sort") !== "ascending";
        headers.forEach(function (h) { h.removeAttribute("aria-sort"); h.querySelector(".arrow") && (h.querySelector(".arrow").textContent = "▴▾"); });
        th.setAttribute("aria-sort", asc ? "ascending" : "descending");
        var arrow = th.querySelector(".arrow");
        if (arrow) arrow.textContent = asc ? "▴" : "▾";

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
      });
    });

    /* search ---------------------------------------------------------- */
    var search = $("input[type=search]", wrap);
    if (search) {
      search.addEventListener("input", function () {
        var q = search.value.trim().toLowerCase();
        var shown = 0;
        rows.forEach(function (r) {
          var hit = !q || r.textContent.toLowerCase().indexOf(q) !== -1;
          r.hidden = !hit;
          if (hit) shown++;
        });
        updateCount(shown);
      });
    }

    /* column show / hide ---------------------------------------------- */
    var menu = $(".cols-menu", wrap);
    if (menu) {
      var list = $(".cols-list", menu);
      var toggle = $(".btn", menu);
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
        label.appendChild(document.createTextNode(th.textContent.replace(/[▴▾]/g, "").trim()));
        list.appendChild(label);
      });
      toggle.addEventListener("click", function (e) { e.stopPropagation(); list.classList.toggle("open"); });
      document.addEventListener("click", function (e) { if (!menu.contains(e.target)) list.classList.remove("open"); });
    }

    /* TSV download ----------------------------------------------------- */
    var dl = $("[data-download]", wrap);
    if (dl) {
      dl.addEventListener("click", function () {
        var names = headers.map(function (h) { return h.textContent.replace(/[▴▾]/g, "").trim(); });
        var lines = [names.join("\t")];
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
    if (v === null || v === undefined || v === "") return "\u2014";
    var n = +v;
    if (!isFinite(n)) return String(v);
    return Math.abs(n - Math.round(n)) < 1e-9
      ? Math.round(n).toLocaleString()
      : n.toLocaleString(undefined, { maximumFractionDigits: 3 });
  }

  /* ECharts takes callbacks for a handful of fields, but a report embeds its
     options as JSON, and JSON has no functions. The R side emits
     {"__fn": "name", "args": [...]} placeholders which are swapped for real
     functions here, after theming (which round-trips through JSON) */
  var FNS = {
    intValue: function () {
      return function (v) { return fmtNum(v); };
    },
    pctValue: function (digits) {
      return function (v) {
        return (v === null || v === undefined || v === "") ? "\u2014" : (+v).toFixed(digits || 0) + "%";
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

  function resolveFns(o) {
    if (Array.isArray(o)) return o.map(resolveFns);
    if (o && typeof o === "object") {
      if (typeof o.__fn === "string" && FNS[o.__fn]) {
        return FNS[o.__fn].apply(null, o.args || []);
      }
      var out = {};
      for (var k in o) {
        if (Object.prototype.hasOwnProperty.call(o, k)) out[k] = resolveFns(o[k]);
      }
      return out;
    }
    return o;
  }

  function themeColors() {
    var cs = getComputedStyle(document.documentElement);
    var get = function (n) { return cs.getPropertyValue(n).trim(); };
    return {
      ink:  get("--ink"),
      ink2: get("--ink-2"),
      ink3: get("--ink-3"),
      line: get("--line"),
      surface: get("--surface"),
      dark: document.documentElement.getAttribute("data-theme") === "dark"
    };
  }

  /* Recursively stamp theme-dependent colours onto an option tree
     The R side emits only data and layout;
     every colour that must follow the theme is left out there and injected here */
  function themeOption(opt) {
    var c = themeColors();
    var o = JSON.parse(JSON.stringify(opt));

    o.backgroundColor = "transparent";
    o.textStyle = Object.assign({ color: c.ink2, fontFamily: getComputedStyle(document.body).fontFamily }, o.textStyle || {});
    o.animation = false;

    /* Every chart gets the same toolbox, so readers can lift a figure out of
       the report or inspect the numbers behind it without leaving the page */
    o.toolbox = Object.assign({
      right: 8,
      top: 2,
      itemSize: 13,
      itemGap: 8,
      iconStyle: { borderColor: c.ink3 },
      emphasis: { iconStyle: { borderColor: c.ink } },
      feature: {
        saveAsImage: { title: "Save as PNG", pixelRatio: 2, backgroundColor: c.surface },
        dataView: { title: "View data", readOnly: true, lang: ["Chart data", "Close", "Refresh"] },
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
        l.textStyle = Object.assign({ color: c.ink2 }, l.textStyle || {});
        l.inactiveColor = c.ink3;
        return l;
      });
    }
    if (o.tooltip) {
      o.tooltip = [].concat(o.tooltip).map(function (t) {
        t.backgroundColor = c.surface;
        t.borderColor = c.line;
        t.textStyle = Object.assign({ color: c.ink, fontSize: 12 }, t.textStyle || {});
        t.extraCssText = "box-shadow:0 4px 18px rgba(0,0,0,.18);border-radius:6px;";
        return t;
      });
    }
    ["xAxis", "yAxis"].forEach(function (k) {
      if (!o[k]) return;
      o[k] = [].concat(o[k]).map(function (ax) {
        ax.axisLine  = Object.assign({ lineStyle: { color: c.line } }, ax.axisLine || {});
        ax.axisTick  = Object.assign({ lineStyle: { color: c.line } }, ax.axisTick || {});
        ax.axisLabel = Object.assign({ color: c.ink2, fontSize: 11 }, ax.axisLabel || {});
        ax.nameTextStyle = Object.assign({ color: c.ink3, fontSize: 11 }, ax.nameTextStyle || {});
        ax.splitLine = Object.assign({ lineStyle: { color: c.line, type: "dashed" } }, ax.splitLine || {});
        return ax;
      });
    });
    if (o.dataZoom) {
      o.dataZoom = [].concat(o.dataZoom).map(function (dz) {
        if (dz.type === "slider") {
          dz.borderColor = c.line;
          dz.textStyle = { color: c.ink3, fontSize: 10 };
          dz.handleStyle = { color: c.surface, borderColor: c.ink3 };
          dz.fillerColor = c.dark ? "rgba(90,171,221,.14)" : "rgba(29,111,165,.10)";
        }
        return dz;
      });
    }
    if (o.visualMap) {
      o.visualMap = [].concat(o.visualMap).map(function (vm) {
        vm.textStyle = Object.assign({ color: c.ink2, fontSize: 11 }, vm.textStyle || {});
        return vm;
      });
    }
    if (o.series) {
      o.series = [].concat(o.series).map(function (s) {
        if (s.type === "sankey") {
          s.label = Object.assign({ color: c.ink, fontSize: 11 }, s.label || {});
          s.lineStyle = Object.assign({ color: "gradient", opacity: c.dark ? 0.32 : 0.22, curveness: 0.5 }, s.lineStyle || {});
          s.itemStyle = Object.assign({ borderColor: c.surface, borderWidth: 1 }, s.itemStyle || {});
        }
        if (s.markLine) {
          s.markLine.label = Object.assign({ color: c.ink2, fontSize: 10 }, s.markLine.label || {});
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

  function renderChart(rec) {
    if (!window.echarts) return;
    if (!rec.inst) {
      rec.inst = echarts.init(rec.el, null, { renderer: rec.spec.renderer || "canvas" });
    }
    rec.inst.setOption(resolveFns(themeOption(rec.spec.option)), true);
  }

  function redrawCharts() {
    charts.forEach(function (rec) { if (rec.inst) renderChart(rec); });
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
    charts.forEach(renderChart);

    var t = null;
    window.addEventListener("resize", function () {
      clearTimeout(t);
      t = setTimeout(function () {
        charts.forEach(function (rec) { if (rec.inst) rec.inst.resize(); });
      }, 120);
    });
    /* Canvas charts do not reflow for the print stylesheet on their own */
    window.addEventListener("beforeprint", function () {
      charts.forEach(function (rec) { if (rec.inst) rec.inst.resize(); });
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
      btn.addEventListener("click", function () {
        var rec = recFor(btn.getAttribute("data-fig-switch"));
        if (!rec || !rec.spec.variants) return;
        var key = btn.getAttribute("data-variant");
        var variant = rec.spec.variants[key];
        if (!variant) return;
        $$('[data-fig-switch="' + btn.getAttribute("data-fig-switch") + '"]').forEach(function (b) {
          b.classList.toggle("on", b === btn);
        });
        rec.spec.option = variant;
        renderChart(rec);
      });
    });
  }

  /* ---------------------------------------------------------- copy button */

  function initCopy() {
    $$("[data-copy]").forEach(function (btn) {
      btn.addEventListener("click", function () {
        var src = document.getElementById(btn.getAttribute("data-copy"));
        if (!src) return;
        var text = src.textContent;
        var done = function (ok) {
          var note = btn.parentNode.querySelector(".copied");
          if (note) {
            note.textContent = ok ? "Copied" : "Press Ctrl+C";
            setTimeout(function () { note.textContent = ""; }, 2000);
          }
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
