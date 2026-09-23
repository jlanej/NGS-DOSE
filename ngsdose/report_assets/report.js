/* NGS-DOSE cohort report: charts drawn from the embedded data. Plain SVG, no dependencies.
   Chart specs live in data.charts; sample values in data.samples. */
(function () {
  "use strict";
  var DATA = JSON.parse(document.getElementById("report-data").textContent);
  var S = DATA.samples || [];
  var COLORS = ["var(--s1)", "var(--s2)", "var(--s3)"];
  var NS = "http://www.w3.org/2000/svg";
  var tip = document.createElement("div"); tip.className = "tip"; document.body.appendChild(tip);

  function el(tag, attrs, parent) {
    var e = document.createElementNS(NS, tag);
    for (var k in attrs) if (attrs[k] !== undefined && attrs[k] !== null) e.setAttribute(k, attrs[k]);
    if (parent) parent.appendChild(e);
    return e;
  }
  function text(parent, x, y, s, attrs) {
    var t = el("text", Object.assign({ x: x, y: y, fill: "var(--ink-2)" }, attrs || {}), parent);
    t.textContent = s; return t;
  }
  function fmt(v, d) {
    if (v === null || v === undefined || !isFinite(v)) return "–";
    if (d === undefined) d = Math.abs(v) >= 100 ? 0 : Math.abs(v) >= 10 ? 1 : Math.abs(v) >= 1 ? 2 : 3;
    return v.toLocaleString(undefined, { minimumFractionDigits: d, maximumFractionDigits: d });
  }
  function nice(lo, hi, n) {                       // clean ticks
    if (!(hi > lo)) { hi = lo + 1; lo = lo - 1; }
    var span = hi - lo, step = Math.pow(10, Math.floor(Math.log10(span / n))), err = span / n / step;
    step *= err >= 7.5 ? 10 : err >= 3.5 ? 5 : err >= 1.5 ? 2 : 1;
    var ticks = [], t = Math.ceil(lo / step) * step;
    for (; t <= hi + 1e-9; t += step) ticks.push(+t.toFixed(10));
    return ticks;
  }
  function showTip(ev, lines) {
    tip.textContent = "";
    lines.forEach(function (l) {
      var d = document.createElement("div");
      if (l.b) { var b = document.createElement("b"); b.textContent = l.b; d.appendChild(b); if (l.k) { var k = document.createElement("span"); k.className = "k"; k.textContent = " " + l.k; d.appendChild(k); } }
      else d.textContent = l.k || l;
      tip.appendChild(d);
    });
    tip.style.display = "block";
    var x = ev.pageX + 14, y = ev.pageY + 14;
    if (x + tip.offsetWidth > document.documentElement.scrollWidth - 8) x = ev.pageX - tip.offsetWidth - 10;
    tip.style.left = x + "px"; tip.style.top = y + "px";
  }
  function hideTip() { tip.style.display = "none"; }

  function values(spec, col) {                    // rows matching spec.where, with a finite value in col
    var out = [];
    for (var i = 0; i < S.length; i++) {
      var r = S[i], ok = true;
      if (spec.where) for (var k in spec.where) if (r[k] !== spec.where[k]) { ok = false; break; }
      if (!ok) continue;
      var v = r[col];
      if (typeof v === "number" && isFinite(v)) out.push({ r: r, v: v });
    }
    return out;
  }
  function frame(container, w, h, m) {
    var svg = el("svg", { viewBox: "0 0 " + w + " " + h, role: "img" }, container);
    return { svg: svg, g: el("g", { transform: "translate(" + m.l + "," + m.t + ")" }, svg), W: w - m.l - m.r, H: h - m.t - m.b, m: m };
  }
  function stepFmt(ticks) {                          // one number of decimals per axis, from the tick step: 9.5, 10.0, 10.5
    if (ticks.length < 2) return fmt;
    var d = Math.max(0, Math.ceil(-Math.log10(Math.abs(ticks[1] - ticks[0])) - 1e-9));
    return function (v) { return fmt(v, d); };
  }
  function axes(f, xs, ys, xlabel, ylabel, xfmt, yfmt, integerY) {
    var ticks = nice(xs.lo, xs.hi, 6), yt = nice(ys.lo, ys.hi, 5);
    xfmt = xfmt || stepFmt(ticks); yfmt = yfmt || stepFmt(yt);
    if (integerY) { yt = yt.filter(function (t) { return Number.isInteger(t); }); if (yt.length < 2) yt = [0, Math.ceil(ys.hi)]; }
    yt.forEach(function (t) { if (t < ys.lo - 1e-9 || t > ys.hi + 1e-9) return; var y = ys.map(t);
      el("line", { x1: 0, x2: f.W, y1: y, y2: y, stroke: "var(--grid)", "stroke-width": 1 }, f.g);
      if (!ylabel || y >= 10) text(f.g, -8, y + 4, (yfmt || fmt)(t), { "text-anchor": "end", fill: "var(--muted)" }); });   // a tick at the very top would sit on the axis title
    el("line", { x1: 0, x2: f.W, y1: f.H, y2: f.H, stroke: "var(--axis)", "stroke-width": 1 }, f.g);
    ticks.forEach(function (t) { if (t < xs.lo - 1e-9 || t > xs.hi + 1e-9) return; var x = xs.map(t);
      text(f.g, x, f.H + 18, (xfmt || fmt)(t), { "text-anchor": "middle", fill: "var(--muted)" }); });     // clear of the y axis's lowest label
    if (xlabel) text(f.g, f.W / 2, f.H + 36, xlabel, { "text-anchor": "middle" });
    if (ylabel) text(f.g, -f.m.l + 12, -8, ylabel, { "text-anchor": "start" });
  }
  function scale(lo, hi, a, b) { var s = { lo: lo, hi: hi }; s.map = function (v) { return a + (v - lo) / (hi - lo || 1) * (b - a); }; return s; }
  function refLine(f, xs, ys, ref, vertical) {
    if (vertical) { var x = xs.map(ref.x); el("line", { x1: x, x2: x, y1: 0, y2: f.H, stroke: "var(--ink-2)", "stroke-width": 1 }, f.g);
      if (ref.label) text(f.g, x + 4, 12, ref.label, { fill: "var(--ink-2)" }); }
    else { var y = ys.map(ref.y); el("line", { x1: 0, x2: f.W, y1: y, y2: y, stroke: "var(--ink-2)", "stroke-width": 1 }, f.g);
      if (ref.label) text(f.g, f.W - 4, y - 4, ref.label, { "text-anchor": "end", fill: "var(--ink-2)" }); }
  }
  function legend(container, items, always) {
    if (items.length < (always ? 1 : 2)) return;
    var d = document.createElement("div"); d.className = "legend";
    items.forEach(function (it) { var s = document.createElement("span"); var i = document.createElement("i");
      i.style.background = it.dash ? "repeating-linear-gradient(90deg, " + it.color + " 0 4px, transparent 4px 7px)" : it.color;
      if (it.line) i.className = "line"; if (it.tick) i.className = "tick";
      if (it.ring) { i.className = "ring"; i.style.background = "transparent"; i.style.borderColor = it.color; }
      s.appendChild(i); s.appendChild(document.createTextNode(it.label)); d.appendChild(s); });
    container.insertBefore(d, container.firstChild);
  }
  function groupsOf(spec) {
    if (!spec.group) return [{ key: null, label: "", where: {} }];
    return spec.group.levels.map(function (lv) { var w = {}; w[spec.group.col] = lv[0]; return { key: lv[0], label: lv[1], where: w }; });
  }
  function empty(container, msg) { var d = document.createElement("div"); d.className = "empty"; d.textContent = msg || "no data yet"; container.appendChild(d); }

  // ---------------- histogram (one or a few series, grouped bars per bin)
  function hist(container, spec) {
    var groups = groupsOf(spec), series = groups.map(function (g) { var sp = Object.assign({}, spec, { where: Object.assign({}, spec.where || {}, g.where) }); return { g: g, vals: values(sp, spec.col).map(function (x) { return x.v; }) }; });
    var all = [].concat.apply([], series.map(function (s) { return s.vals; }));
    if (!all.length) return empty(container);
    var lo = spec.xlim ? spec.xlim[0] : Math.min.apply(null, all), hi = spec.xlim ? spec.xlim[1] : Math.max.apply(null, all);
    if (spec.ref) spec.ref.forEach(function (r) { lo = Math.min(lo, r.x); hi = Math.max(hi, r.x); });
    if (hi === lo) { lo -= 1; hi += 1; }
    var pad = (hi - lo) * 0.04; lo -= pad; hi += pad;
    var nb = spec.bins || Math.max(10, Math.min(40, Math.round(Math.sqrt(all.length) * 1.5)));
    var edges = nice(lo, hi, nb); if (edges.length < 3) edges = [lo, (lo + hi) / 2, hi];
    var step = edges[1] - edges[0];
    if (edges[0] > lo) edges.unshift(edges[0] - step); if (edges[edges.length - 1] < hi) edges.push(edges[edges.length - 1] + step);
    var counts = series.map(function (s) { var c = new Array(edges.length - 1).fill(0); s.vals.forEach(function (v) { var i = Math.min(edges.length - 2, Math.max(0, Math.floor((v - edges[0]) / step))); c[i]++; }); return c; });
    var ymax = Math.max.apply(null, counts.map(function (c) { return Math.max.apply(null, c); })) || 1;
    var f = frame(container, 640, 260, { l: 48, t: 22, r: 12, b: 44 });
    var xs = scale(edges[0], edges[edges.length - 1], 0, f.W), ys = scale(0, ymax * 1.05, f.H, 0);
    axes(f, xs, ys, spec.xlabel, "samples", spec.xfmt !== undefined ? function (v) { return fmt(v, spec.xfmt); } : undefined, function (v) { return fmt(v, 0); }, true);
    var slot = xs.map(edges[1]) - xs.map(edges[0]), inner = slot - 2, bw = Math.min(24, Math.max(2, (inner - 2 * (series.length - 1)) / series.length));
    counts.forEach(function (c, si) {
      c.forEach(function (n, i) {
        if (!n) return;
        var x0 = xs.map(edges[i]) + 1 + si * (bw + 2), y0 = ys.map(n), h = f.H - y0;
        var r = Math.min(4, bw / 2, h);
        var p = "M" + x0 + "," + f.H + " v" + (-(h - r)) + " a" + r + "," + r + " 0 0 1 " + r + ",-" + r + " h" + (bw - 2 * r) + " a" + r + "," + r + " 0 0 1 " + r + "," + r + " v" + (h - r) + " z";
        var bar = el("path", { d: p, fill: COLORS[si] }, f.g);
        var hit = el("rect", { x: x0 - 1, y: 0, width: bw + 2, height: f.H, fill: "transparent" }, f.g);
        var lines = [{ b: n + (n === 1 ? " sample" : " samples"), k: series[si].g.label }, { k: fmt(edges[i], spec.xfmt) + " – " + fmt(edges[i + 1], spec.xfmt) }];
        hit.addEventListener("pointermove", function (ev) { bar.setAttribute("opacity", 0.75); showTip(ev, lines); });
        hit.addEventListener("pointerleave", function () { bar.removeAttribute("opacity"); hideTip(); });
      });
    });
    if (spec.ref) spec.ref.forEach(function (r) { refLine(f, xs, ys, r, true); });
    legend(container, series.map(function (s, i) { return { color: COLORS[i], label: s.g.label + " (n = " + s.vals.length + ")" }; }));
  }

  // ---------------- scatter with nearest-point hover
  function scatter(container, spec) {
    var groups = groupsOf(spec), pts = [];
    if (spec.points) spec.points.forEach(function (p) { pts.push({ x: p.x, y: p.y, label: p.label, si: p.si || 0, extra: p.extra }); });
    else groups.forEach(function (g, si) {
      var sp = Object.assign({}, spec, { where: Object.assign({}, spec.where || {}, g.where) });
      values(sp, spec.x).forEach(function (o) { var y = o.r[spec.y]; if (typeof y === "number" && isFinite(y)) pts.push({ x: o.v, y: y, label: o.r.sample, si: si, r: o.r }); });
    });
    if (spec.log) pts = pts.filter(function (p) { return p.x > 0 && p.y > 0; });
    if (pts.length < 1) return empty(container);
    var tx = function (v) { return spec.log ? Math.log10(v) : v; };
    var xv = pts.map(function (p) { return tx(p.x); }), yv = pts.map(function (p) { return tx(p.y); });
    var xlo = Math.min.apply(null, xv), xhi = Math.max.apply(null, xv), ylo = Math.min.apply(null, yv), yhi = Math.max.apply(null, yv);
    if (spec.identity) { xlo = ylo = Math.min(xlo, ylo); xhi = yhi = Math.max(xhi, yhi); }
    if (spec.xref !== undefined) { xlo = Math.min(xlo, tx(spec.xref)); xhi = Math.max(xhi, tx(spec.xref)); }
    if (spec.yref !== undefined) { ylo = Math.min(ylo, tx(spec.yref)); yhi = Math.max(yhi, tx(spec.yref)); }
    var px = (xhi - xlo || 1) * 0.06, py = (yhi - ylo || 1) * 0.08;
    var f = frame(container, 640, 320, { l: 56, t: 22, r: 16, b: 44 });
    var xs = scale(xlo - px, xhi + px, 0, f.W), ys = scale(ylo - py, yhi + py, f.H, 0);
    var lf = spec.log ? function (v) { return fmt(Math.pow(10, v)); } : undefined;
    axes(f, xs, ys, spec.xlabel, spec.ylabel, lf, lf);
    if (spec.identity) el("line", { x1: xs.map(xlo - px), y1: ys.map(xlo - px), x2: xs.map(xhi + px), y2: ys.map(xhi + px), stroke: "var(--axis)", "stroke-width": 1 }, f.g);
    // y = k x, clipped to the plot: k = 1 is equality without forcing the two axes onto one range, dashed for any other slope
    var refKeys = [];
    function slopeLine(k, dashed) {
      if (spec.log || !(k > 0)) return false;
      var a = Math.max(xs.lo, ys.lo / k), b = Math.min(xs.hi, ys.hi / k); if (!(b > a)) return false;
      el("line", { x1: xs.map(a), y1: ys.map(k * a), x2: xs.map(b), y2: ys.map(k * b), stroke: dashed ? "var(--ink-2)" : "var(--axis)", "stroke-width": dashed ? 1.5 : 1, "stroke-dasharray": dashed ? "5 4" : null }, f.g);
      return true;
    }
    if (spec.diagonal && slopeLine(1, false)) refKeys.push({ color: "var(--axis)", label: "equality", line: true });
    if (spec.slope_ref && slopeLine(spec.slope_ref.k, true)) refKeys.push({ color: "var(--ink-2)", label: spec.slope_ref.label, line: true, dash: true });
    if (spec.xref !== undefined) refLine(f, xs, ys, { x: tx(spec.xref) }, true);
    if (spec.yref !== undefined) refLine(f, xs, ys, { y: tx(spec.yref), label: spec.yref_label }, false);
    var fitLabels = [];
    if (spec.fit) {                                         // least squares on the (transformed) values, one line per series; a line through a handful of points misleads
      var nser = Math.max.apply(null, pts.map(function (p) { return p.si; })) + 1;
      for (var s = 0; s < nser; s++) {
        var sx = [], sy = [];
        pts.forEach(function (p, i) { if (p.si === s) { sx.push(xv[i]); sy.push(yv[i]); } });
        if (sx.length < 10) continue;
        var n = sx.length, mx = sx.reduce(function (a, b) { return a + b; }, 0) / n, my = sy.reduce(function (a, b) { return a + b; }, 0) / n, sxy = 0, sxx = 0;
        for (var i = 0; i < n; i++) { sxy += (sx[i] - mx) * (sy[i] - my); sxx += (sx[i] - mx) * (sx[i] - mx); }
        var b = sxy / (sxx || 1), a0 = my - b * mx, x1 = Math.min.apply(null, sx), x2 = Math.max.apply(null, sx);   // over the data, never past it
        el("line", { x1: xs.map(x1), y1: ys.map(a0 + b * x1), x2: xs.map(x2), y2: ys.map(a0 + b * x2), stroke: nser > 1 ? COLORS[s] : "var(--ink-2)", "stroke-width": 2, "stroke-linecap": "round", opacity: nser > 1 ? 0.9 : 0.6 }, f.g);
        if (spec.fit_labels && spec.fit_labels[s]) fitLabels.push({ x: xs.map(x2), y: ys.map(a0 + b * x2), s: spec.fit_labels[s] });
      }
    }
    var dots = pts.map(function (p) { p.cx = xs.map(tx(p.x)); p.cy = ys.map(tx(p.y)); return el("circle", { cx: p.cx, cy: p.cy, r: 4, fill: COLORS[p.si], stroke: "var(--surface)", "stroke-width": 2 }, f.g); });
    fitLabels.sort(function (a, b) { return a.y - b.y; }).forEach(function (l, i, all) {   // at the line's end, above it, kept apart, on a halo
      var y = Math.max(12, l.y - 8); if (i > 0) y = Math.max(y, all[i - 1].ty + 15); l.ty = y;
      text(f.g, l.x, y, l.s, { "text-anchor": "end", fill: "var(--ink)", "font-weight": 600, stroke: "var(--surface)", "stroke-width": 4, "paint-order": "stroke" });
    });
    var hot = null;
    f.svg.addEventListener("pointermove", function (ev) {
      var rect = f.svg.getBoundingClientRect(), k = 640 / rect.width, mx2 = (ev.clientX - rect.left) * k - f.m.l, my2 = (ev.clientY - rect.top) * k - f.m.t, best = -1, bd = 24 * 24;
      for (var i = 0; i < pts.length; i++) { var d = (pts[i].cx - mx2) * (pts[i].cx - mx2) + (pts[i].cy - my2) * (pts[i].cy - my2); if (d < bd) { bd = d; best = i; } }
      if (hot !== null) dots[hot].setAttribute("r", 4);
      if (best < 0) { hot = null; hideTip(); return; }
      hot = best; dots[best].setAttribute("r", 6);
      var p = pts[best], lines = [{ b: p.label || "", k: spec.legend ? spec.legend[p.si] : (groups[p.si] && groups[p.si].label) }, { b: fmt(p.y), k: spec.ylabel }, { b: fmt(p.x), k: spec.xlabel }];
      if (p.extra) p.extra.forEach(function (e) { lines.push({ k: e }); });
      if (p.r && p.r.pop) lines.push({ k: p.r.pop + (p.r.sex_inferred ? ", " + p.r.sex_inferred : "") });
      showTip(ev, lines);
    });
    f.svg.addEventListener("pointerleave", function () { if (hot !== null) dots[hot].setAttribute("r", 4); hot = null; hideTip(); });
    legend(container, (spec.legend ? spec.legend.map(function (l, i) { return { color: COLORS[i], label: l }; })
      : groups.map(function (g, i) { return { color: COLORS[i], label: g.label }; }).filter(function (g) { return g.label; })).concat(refKeys), refKeys.length > 0);
  }

  // ---------------- strip: values by group, with the median
  function strip(container, spec) {
    var order = spec.order || [], groups = {};
    values(spec, spec.col).forEach(function (o) { var g = o.r[spec.by]; if (!g) return; if (!groups[g]) groups[g] = []; groups[g].push(o); });
    if (!order.length) order = Object.keys(groups).sort();
    order = order.filter(function (g) { return groups[g]; });
    if (!order.length) return empty(container);
    var all = [].concat.apply([], order.map(function (g) { return groups[g].map(function (o) { return o.v; }); }));
    var lo = Math.min.apply(null, all), hi = Math.max.apply(null, all);
    if (spec.ref !== undefined) { lo = Math.min(lo, spec.ref); hi = Math.max(hi, spec.ref); }
    var pad = (hi - lo || 1) * 0.08;
    var w = Math.max(640, 40 * order.length), f = frame(container, w, 300, { l: 56, t: 22, r: 12, b: order.length > 8 ? 60 : 44 });
    var ys = scale(lo - pad, hi + pad, f.H, 0), band = f.W / order.length;
    var yt = nice(ys.lo, ys.hi, 5);
    yt.forEach(function (t) { var y = ys.map(t); el("line", { x1: 0, x2: f.W, y1: y, y2: y, stroke: "var(--grid)" }, f.g); text(f.g, -8, y + 4, fmt(t), { "text-anchor": "end", fill: "var(--muted)" }); });
    el("line", { x1: 0, x2: f.W, y1: f.H, y2: f.H, stroke: "var(--axis)" }, f.g);
    if (spec.ylabel) text(f.g, -f.m.l + 12, -8, spec.ylabel);
    if (spec.ref !== undefined) refLine(f, null, ys, { y: spec.ref, label: spec.ref_label }, false);
    var dots = [], pts = [];
    order.forEach(function (g, gi) {
      var xs0 = gi * band + band / 2, vals = groups[g].map(function (o) { return o.v; }).sort(function (a, b) { return a - b; });
      var med = vals.length % 2 ? vals[(vals.length - 1) / 2] : (vals[vals.length / 2 - 1] + vals[vals.length / 2]) / 2;
      var lab = (spec.labels && spec.labels[g]) || g;
      var tx = text(f.g, xs0, f.H + 16, lab, { "text-anchor": order.length > 8 ? "end" : "middle", fill: "var(--muted)" });
      if (order.length > 8) tx.setAttribute("transform", "rotate(-40 " + xs0 + " " + (f.H + 16) + ")");
      groups[g].forEach(function (o, i) {
        var h = 0, s = o.r.sample || String(i); for (var c = 0; c < s.length; c++) h = (h * 31 + s.charCodeAt(c)) % 1000;   // deterministic jitter
        var cx = xs0 + (h / 1000 - 0.5) * Math.min(band * 0.7, 28), cy = ys.map(o.v);
        pts.push({ cx: cx, cy: cy, o: o, g: lab }); dots.push(el("circle", { cx: cx, cy: cy, r: 3.5, fill: "var(--s1)", stroke: "var(--surface)", "stroke-width": 1.5, opacity: 0.85 }, f.g));
      });
      el("line", { x1: xs0 - Math.min(band * 0.4, 16), x2: xs0 + Math.min(band * 0.4, 16), y1: ys.map(med), y2: ys.map(med), stroke: "var(--ink)", "stroke-width": 2, "stroke-linecap": "round" }, f.g);
    });
    var hot = null;
    f.svg.addEventListener("pointermove", function (ev) {
      var rect = f.svg.getBoundingClientRect(), k = w / rect.width, mx = (ev.clientX - rect.left) * k - f.m.l, my = (ev.clientY - rect.top) * k - f.m.t, best = -1, bd = 24 * 24;
      for (var i = 0; i < pts.length; i++) { var d = (pts[i].cx - mx) * (pts[i].cx - mx) + (pts[i].cy - my) * (pts[i].cy - my); if (d < bd) { bd = d; best = i; } }
      if (hot !== null) dots[hot].setAttribute("r", 3.5);
      if (best < 0) { hot = null; hideTip(); return; }
      hot = best; dots[best].setAttribute("r", 6);
      var p = pts[best]; showTip(ev, [{ b: p.o.r.sample || "", k: p.g }, { b: fmt(p.o.v), k: spec.ylabel }]);
    });
    f.svg.addEventListener("pointerleave", function () { if (hot !== null) dots[hot].setAttribute("r", 3.5); hot = null; hideTip(); });
  }

  // ---------------- lines with a crosshair (the PC sweep)
  function lines(container, spec) {
    // options: series[].ci (colour index), bands [{x0, x1, label}] shaded behind, ticks {x: [...], label} marked on the axis, xfmt
    var ser = spec.series.filter(function (s) { return s.y && s.y.length; });
    if (!ser.length) return empty(container);
    var col = function (s, si) { return COLORS[s.ci !== undefined ? s.ci : si]; };
    var xs0 = [].concat.apply([], ser.map(function (s) { return s.x; })), ys0 = [].concat.apply([], ser.map(function (s) { return s.y.concat(s.lo || [], s.hi || []); })).filter(isFinite);
    var lo = Math.min.apply(null, ys0), hi = Math.max.apply(null, ys0); if (spec.ref !== undefined) { lo = Math.min(lo, spec.ref); hi = Math.max(hi, spec.ref); }
    var pad = (hi - lo || 1) * 0.1;
    var f = frame(container, 640, 280, { l: 56, t: 22, r: 16, b: 44 });
    var xs = scale(Math.min.apply(null, xs0), Math.max.apply(null, xs0), 0, f.W), ys = scale(lo - pad, hi + pad, f.H, 0);
    var xd = spec.xfmt === undefined ? 0 : spec.xfmt;
    var bandLabels = [];
    (spec.bands || []).forEach(function (bd) {                // shaded x-ranges behind everything, labelled where the label fits
      var x0 = xs.map(Math.max(bd.x0, xs.lo)), x1 = xs.map(Math.min(bd.x1, xs.hi)); if (!(x1 > x0)) return;
      el("rect", { x: x0, y: 0, width: x1 - x0, height: f.H, fill: "var(--grid)", opacity: 0.55 }, f.g);
      if (bd.label && x1 - x0 >= bd.label.length * 6.5 + 6) bandLabels.push([(x0 + x1) / 2, bd.label]);
    });
    axes(f, xs, ys, spec.xlabel, spec.ylabel, function (v) { return fmt(v, xd); });
    bandLabels.forEach(function (l) { text(f.g, l[0], 12, l[1], { "text-anchor": "middle", fill: "var(--ink-2)", "font-size": 11, stroke: "var(--surface)", "stroke-width": 3, "paint-order": "stroke" }); });
    if (spec.ref !== undefined) refLine(f, xs, ys, { y: spec.ref, label: spec.ref_label }, false);
    ser.forEach(function (s, si) {
      if (s.lo && s.hi) { var d = "M" + s.x.map(function (x, i) { return xs.map(x) + "," + ys.map(s.lo[i]); }).join(" L") + " L" + s.x.slice().reverse().map(function (x, i) { var j = s.x.length - 1 - i; return xs.map(x) + "," + ys.map(s.hi[j]); }).join(" L") + " z";
        el("path", { d: d, fill: col(s, si), opacity: 0.14 }, f.g); }
      el("path", { d: "M" + s.x.map(function (x, i) { return xs.map(x) + "," + ys.map(s.y[i]); }).join(" L"), fill: "none", stroke: col(s, si), "stroke-width": 2, "stroke-linejoin": "round", "stroke-linecap": "round" }, f.g);
      var n = s.x.length - 1; el("circle", { cx: xs.map(s.x[n]), cy: ys.map(s.y[n]), r: 4, fill: col(s, si), stroke: "var(--surface)", "stroke-width": 2 }, f.g);
    });
    if (spec.ticks) spec.ticks.x.forEach(function (x) {       // e.g. the anchor windows: short marks standing on the x axis
      el("line", { x1: xs.map(x), x2: xs.map(x), y1: f.H, y2: f.H - 10, stroke: "var(--s1)", "stroke-width": 2 }, f.g);
    });
    var cross = el("line", { y1: 0, y2: f.H, stroke: "var(--ink-2)", "stroke-width": 1, opacity: 0 }, f.g);
    f.svg.addEventListener("pointermove", function (ev) {
      var rect = f.svg.getBoundingClientRect(), k = 640 / rect.width, mx = (ev.clientX - rect.left) * k - f.m.l;
      var xv = xs.lo + (mx / f.W) * (xs.hi - xs.lo); if (xv < xs.lo || xv > xs.hi) { cross.setAttribute("opacity", 0); hideTip(); return; }
      var out = [], shown = null;
      ser.forEach(function (s) {                             // each series' nearest point to the pointer
        var best = 0; for (var i = 1; i < s.x.length; i++) if (Math.abs(s.x[i] - xv) < Math.abs(s.x[best] - xv)) best = i;
        if (shown === null) { shown = s.x[best]; out.push({ b: spec.xlabel + " " + fmt(shown, xd + 1) }); }
        out.push({ b: fmt(s.y[best], 3) + (s.lo ? " (" + fmt(s.lo[best], 2) + " to " + fmt(s.hi[best], 2) + ")" : ""), k: s.name });
      });
      cross.setAttribute("x1", xs.map(shown)); cross.setAttribute("x2", xs.map(shown)); cross.setAttribute("opacity", 0.6);
      showTip(ev, out);
    });
    f.svg.addEventListener("pointerleave", function () { cross.setAttribute("opacity", 0); hideTip(); });
    legend(container, ser.map(function (s, i) { return { color: col(s, i), label: s.name, line: true }; })
      .concat(spec.ticks && spec.ticks.label ? [{ color: "var(--s1)", label: spec.ticks.label, tick: true }] : []), !!spec.ticks);
  }

  // ---------------- meters: done of total per category
  function meters(container, spec) {
    var cats = spec.categories; if (!cats.length) return empty(container);
    var longest = Math.max.apply(null, cats.map(function (c) { return c.length; })), rowH = 26;
    var f = frame(container, 640, cats.length * rowH + 10, { l: Math.min(200, 12 + longest * 6.8), t: 4, r: 76, b: 4 });
    var tot = Math.max.apply(null, spec.totals), xs = scale(0, tot, 0, f.W);
    cats.forEach(function (c, i) {
      var y = i * rowH + 6;
      text(f.g, -8, y + 13, c, { "text-anchor": "end" });
      el("rect", { x: 0, y: y, width: xs.map(spec.totals[i]), height: 16, rx: 4, fill: "var(--s1)", opacity: 0.18 }, f.g);
      el("rect", { x: 0, y: y, width: Math.max(0, xs.map(spec.values[i])), height: 16, rx: 4, fill: "var(--s1)" }, f.g);
      text(f.g, xs.map(spec.totals[i]) + 8, y + 13, fmt(spec.values[i], 0) + " / " + fmt(spec.totals[i], 0), { fill: "var(--ink-2)" });
    });
  }

  // ---------------- forest: one estimate and its interval per row, grouped rows, vertical reference lines
  function forest(container, spec) {
    var rows = spec.rows || []; if (!rows.length) return empty(container);
    var longest = Math.max.apply(null, rows.map(function (r) { return r.length; })), rowH = 20, groupH = 24;
    var order = [], last = null;
    rows.forEach(function (r, i) { var g = spec.groups ? spec.groups[i] : null; if (g !== last) { order.push({ header: g }); last = g; } order.push({ i: i }); });
    var H = order.reduce(function (a, o) { return a + (o.header !== undefined ? groupH : rowH); }, 0);
    var f = frame(container, 640, H + 50, { l: Math.min(230, 12 + longest * 6.6), t: 8, r: 16, b: 42 });
    var num = function (v) { return typeof v === "number" && isFinite(v); };
    var vals = [].concat(spec.lo, spec.hi, spec.est, spec.est2 || [], spec.refs || []).filter(num);
    var lo = spec.xmin !== undefined ? spec.xmin : Math.min.apply(null, vals), hi = spec.xmax !== undefined ? spec.xmax : Math.max.apply(null, vals);
    var xs = scale(lo, hi, 0, f.W), clamp = function (v) { return Math.max(lo, Math.min(hi, v)); };
    if (spec.legend) legend(container, [{ color: COLORS[spec.ci || 0], label: spec.legend[0] }].concat(spec.est2 ? [{ color: COLORS[spec.ci || 0], label: spec.legend[1], ring: true }] : []), true);
    var ticks = nice(lo, hi, 6), xf = stepFmt(ticks), colour = COLORS[spec.ci || 0];
    ticks.forEach(function (t) { if (t < lo - 1e-9 || t > hi + 1e-9) return; var x = xs.map(t);
      el("line", { x1: x, x2: x, y1: 0, y2: f.H, stroke: "var(--grid)", "stroke-width": 1 }, f.g);
      text(f.g, x, f.H + 18, xf(t), { "text-anchor": "middle", fill: "var(--muted)" }); });
    if (spec.xlabel) text(f.g, f.W / 2, f.H + 36, spec.xlabel, { "text-anchor": "middle" });
    (spec.refs || []).forEach(function (r) { var x = xs.map(r); el("line", { x1: x, x2: x, y1: 0, y2: f.H, stroke: "var(--ink-2)", "stroke-width": 1 }, f.g); });
    var y = 0;
    order.forEach(function (o) {
      if (o.header !== undefined) {                              // a group's title, on the surface so the gridlines stop behind it
        var ht = text(f.g, -f.m.l + 4, y + 17, o.header, { "font-weight": 600, "font-size": 11.5, fill: "var(--ink)" });
        try { var bb = ht.getBBox(); f.g.insertBefore(el("rect", { x: bb.x - 4, y: bb.y - 2, width: bb.width + 8, height: bb.height + 4, fill: "var(--surface)" }), ht); } catch (e) { /* not laid out */ }
        y += groupH; return;
      }
      var i = o.i, cy = y + rowH / 2, e = spec.est[i], l = spec.lo[i], h = spec.hi[i];
      text(f.g, -8, cy + 4, rows[i], { "text-anchor": "end" });
      if (num(l) && num(h)) el("line", { x1: xs.map(clamp(l)), x2: xs.map(clamp(h)), y1: cy, y2: cy, stroke: colour, "stroke-width": 2, "stroke-linecap": "round" }, f.g);
      if (num(e)) el("circle", { cx: xs.map(clamp(e)), cy: cy, r: 4.5, fill: colour, stroke: "var(--surface)", "stroke-width": 2 }, f.g);
      var e2 = spec.est2 ? spec.est2[i] : null;                  // a second estimate of the same row: an open ring over the first
      if (num(e2)) el("circle", { cx: xs.map(clamp(e2)), cy: cy, r: 4.5, fill: "none", stroke: colour, "stroke-width": 1.75 }, f.g);
      var hit = el("rect", { x: -f.m.l, y: y, width: f.W + f.m.l, height: rowH, fill: "transparent" }, f.g);
      var names = spec.legend || [];
      var lines = [{ b: rows[i] }, { b: fmt(e, 2) + (num(l) && num(h) ? " (" + fmt(l, 2) + " to " + fmt(h, 2) + ")" : ""), k: names[0] || "" }];
      if (num(e2)) lines.push({ b: fmt(e2, 2), k: names[1] || "" });
      if (spec.extra && spec.extra[i]) lines.push({ k: spec.extra[i] });
      hit.addEventListener("pointermove", function (ev) { showTip(ev, lines); });
      hit.addEventListener("pointerleave", hideTip);
      y += rowH;
    });
  }

  // ---------------- heatmap: metrics by statistics, grouped rows, values printed in the cells
  function heatmap(container, spec) {
    var rows = spec.rows || [], cols = spec.cols || []; if (!rows.length || !cols.length) return empty(container);
    var longest = Math.max.apply(null, rows.map(function (r) { return r.length; })), rowH = 22, groupH = 26, lo = spec.lo === undefined ? 0 : spec.lo, hi = spec.hi === undefined ? 1 : spec.hi;
    var order = [], last = null;                        // a header row wherever the group changes
    rows.forEach(function (r, i) { var g = spec.groups ? spec.groups[i] : null; if (g !== last) { order.push({ header: g }); last = g; } order.push({ i: i }); });
    var f = frame(container, 640, order.reduce(function (a, o) { return a + (o.header !== undefined ? groupH : rowH); }, 0) + 40, { l: Math.min(230, 12 + longest * 6.6), t: 34, r: 8, b: 6 });
    var cw = f.W / cols.length;
    cols.forEach(function (c, j) { text(f.g, j * cw + cw / 2, -10, c, { "text-anchor": "middle", "font-size": 11 }); });
    var y = 0;
    order.forEach(function (o) {
      if (o.header !== undefined) { text(f.g, -f.m.l + 4, y + 18, o.header, { "font-weight": 600, "font-size": 11.5, fill: "var(--ink)" }); y += groupH; return; }
      var i = o.i; text(f.g, -8, y + rowH / 2 + 4, rows[i], { "text-anchor": "end" });
      cols.forEach(function (c, j) {
        var v = spec.values[i][j], ok = typeof v === "number" && isFinite(v), a = ok ? Math.max(0, Math.min(1, (v - lo) / (hi - lo || 1))) : 0;
        var cell = el("rect", { x: j * cw + 1, y: y + 1, width: cw - 2, height: rowH - 2, rx: 3, fill: "var(--s1)", opacity: ok ? 0.06 + 0.94 * a : 0.03 }, f.g);
        if (ok) text(f.g, j * cw + cw / 2, y + rowH / 2 + 4, fmt(v, spec.nd === undefined ? 2 : spec.nd), { "text-anchor": "middle", fill: a > 0.55 ? "#ffffff" : "var(--ink)", "font-size": 11.5, "pointer-events": "none" });
        var lines = [{ b: rows[i] }, { b: ok ? fmt(v, 3) : "–", k: spec.col_titles ? spec.col_titles[j] : c }]; if (spec.extra && spec.extra[i] && spec.extra[i][j]) lines.push({ k: spec.extra[i][j] });
        cell.addEventListener("pointermove", function (ev) { cell.setAttribute("stroke", "var(--ink)"); showTip(ev, lines); });
        cell.addEventListener("pointerleave", function () { cell.removeAttribute("stroke"); hideTip(); });
      });
      y += rowH;
    });
  }

  var TYPES = { hist: hist, scatter: scatter, strip: strip, lines: lines, meters: meters, heatmap: heatmap, forest: forest };
  Object.keys(DATA.charts || {}).forEach(function (id) {
    var c = document.getElementById("chart-" + id); if (!c) return;
    try { TYPES[DATA.charts[id].type](c, DATA.charts[id]); } catch (e) { empty(c, "chart failed: " + e.message); }
  });

  // ---------------- tables: sort on header click, filter box
  document.querySelectorAll("table.data").forEach(function (t) {
    var ths = t.querySelectorAll("thead th");
    ths.forEach(function (th, ci) {
      th.addEventListener("click", function () {
        var rows = Array.prototype.slice.call(t.querySelectorAll("tbody tr")), asc = th.getAttribute("data-asc") !== "1";
        ths.forEach(function (o) { o.removeAttribute("data-asc"); }); th.setAttribute("data-asc", asc ? "1" : "0");
        rows.sort(function (a, b) {
          var x = a.children[ci].textContent, y = b.children[ci].textContent, nx = parseFloat(x.replace(/,/g, "")), ny = parseFloat(y.replace(/,/g, ""));
          var r = (!isNaN(nx) && !isNaN(ny)) ? nx - ny : (x === "–" ? 1 : y === "–" ? -1 : x.localeCompare(y));
          return asc ? r : -r;
        });
        rows.forEach(function (r) { t.querySelector("tbody").appendChild(r); });
      });
    });
    var box = t.parentElement.previousElementSibling;
    if (box && box.classList.contains("filter")) box.addEventListener("input", function () {
      var q = box.value.toLowerCase(); t.querySelectorAll("tbody tr").forEach(function (r) { r.style.display = r.textContent.toLowerCase().indexOf(q) >= 0 ? "" : "none"; });
    });
  });

  // ---------------- theme toggle
  var btn = document.getElementById("theme");
  function setTheme(v) { if (v) document.documentElement.setAttribute("data-theme", v); else document.documentElement.removeAttribute("data-theme"); try { localStorage.setItem("ngsdose-theme", v || ""); } catch (e) {} }
  try { var saved = localStorage.getItem("ngsdose-theme"); if (saved) setTheme(saved); } catch (e) {}
  if (btn) btn.addEventListener("click", function () {
    var dark = document.documentElement.getAttribute("data-theme") === "dark" || (!document.documentElement.getAttribute("data-theme") && window.matchMedia("(prefers-color-scheme: dark)").matches);
    setTheme(dark ? "light" : "dark");
  });
})();
