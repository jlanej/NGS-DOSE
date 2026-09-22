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
  function axes(f, xs, ys, xlabel, ylabel, xfmt, yfmt, integerY) {
    var ticks = nice(xs.lo, xs.hi, 6), yt = nice(ys.lo, ys.hi, 5);
    if (integerY) { yt = yt.filter(function (t) { return Number.isInteger(t); }); if (yt.length < 2) yt = [0, Math.ceil(ys.hi)]; }
    yt.forEach(function (t) { if (t < ys.lo - 1e-9 || t > ys.hi + 1e-9) return; var y = ys.map(t);
      el("line", { x1: 0, x2: f.W, y1: y, y2: y, stroke: "var(--grid)", "stroke-width": 1 }, f.g);
      text(f.g, -8, y + 4, (yfmt || fmt)(t), { "text-anchor": "end", fill: "var(--muted)" }); });
    el("line", { x1: 0, x2: f.W, y1: f.H, y2: f.H, stroke: "var(--axis)", "stroke-width": 1 }, f.g);
    ticks.forEach(function (t) { if (t < xs.lo - 1e-9 || t > xs.hi + 1e-9) return; var x = xs.map(t);
      text(f.g, x, f.H + 16, (xfmt || fmt)(t), { "text-anchor": "middle", fill: "var(--muted)" }); });
    if (xlabel) text(f.g, f.W / 2, f.H + 34, xlabel, { "text-anchor": "middle" });
    if (ylabel) text(f.g, -f.m.l + 12, -8, ylabel, { "text-anchor": "start" });
  }
  function scale(lo, hi, a, b) { var s = { lo: lo, hi: hi }; s.map = function (v) { return a + (v - lo) / (hi - lo || 1) * (b - a); }; return s; }
  function refLine(f, xs, ys, ref, vertical) {
    if (vertical) { var x = xs.map(ref.x); el("line", { x1: x, x2: x, y1: 0, y2: f.H, stroke: "var(--ink-2)", "stroke-width": 1 }, f.g);
      if (ref.label) text(f.g, x + 4, 12, ref.label, { fill: "var(--ink-2)" }); }
    else { var y = ys.map(ref.y); el("line", { x1: 0, x2: f.W, y1: y, y2: y, stroke: "var(--ink-2)", "stroke-width": 1 }, f.g);
      if (ref.label) text(f.g, f.W - 4, y - 4, ref.label, { "text-anchor": "end", fill: "var(--ink-2)" }); }
  }
  function legend(container, items) {
    if (items.length < 2) return;
    var d = document.createElement("div"); d.className = "legend";
    items.forEach(function (it) { var s = document.createElement("span"); var i = document.createElement("i"); i.style.background = it.color; if (it.line) i.className = "line"; s.appendChild(i); s.appendChild(document.createTextNode(it.label)); d.appendChild(s); });
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
    if (spec.points) spec.points.forEach(function (p) { pts.push({ x: p.x, y: p.y, label: p.label, si: 0, extra: p.extra }); });
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
    if (spec.xref !== undefined) refLine(f, xs, ys, { x: tx(spec.xref) }, true);
    if (spec.yref !== undefined) refLine(f, xs, ys, { y: tx(spec.yref), label: spec.yref_label }, false);
    if (spec.fit && pts.length >= 10) {                     // least squares on the (transformed) values; a line through a handful of points misleads
      var n = pts.length, mx = xv.reduce(function (a, b) { return a + b; }, 0) / n, my = yv.reduce(function (a, b) { return a + b; }, 0) / n, sxy = 0, sxx = 0;
      for (var i = 0; i < n; i++) { sxy += (xv[i] - mx) * (yv[i] - my); sxx += (xv[i] - mx) * (xv[i] - mx); }
      var b = sxy / (sxx || 1), a0 = my - b * mx;
      el("line", { x1: xs.map(xlo), y1: ys.map(a0 + b * xlo), x2: xs.map(xhi), y2: ys.map(a0 + b * xhi), stroke: "var(--ink-2)", "stroke-width": 2, "stroke-linecap": "round", opacity: 0.6 }, f.g);
    }
    var dots = pts.map(function (p) { p.cx = xs.map(tx(p.x)); p.cy = ys.map(tx(p.y)); return el("circle", { cx: p.cx, cy: p.cy, r: 4, fill: COLORS[p.si], stroke: "var(--surface)", "stroke-width": 2 }, f.g); });
    var hot = null;
    f.svg.addEventListener("pointermove", function (ev) {
      var rect = f.svg.getBoundingClientRect(), k = 640 / rect.width, mx2 = (ev.clientX - rect.left) * k - f.m.l, my2 = (ev.clientY - rect.top) * k - f.m.t, best = -1, bd = 24 * 24;
      for (var i = 0; i < pts.length; i++) { var d = (pts[i].cx - mx2) * (pts[i].cx - mx2) + (pts[i].cy - my2) * (pts[i].cy - my2); if (d < bd) { bd = d; best = i; } }
      if (hot !== null) dots[hot].setAttribute("r", 4);
      if (best < 0) { hot = null; hideTip(); return; }
      hot = best; dots[best].setAttribute("r", 6);
      var p = pts[best], lines = [{ b: p.label || "", k: groups[p.si] && groups[p.si].label }, { b: fmt(p.y), k: spec.ylabel }, { b: fmt(p.x), k: spec.xlabel }];
      if (p.extra) p.extra.forEach(function (e) { lines.push({ k: e }); });
      if (p.r && p.r.pop) lines.push({ k: p.r.pop + (p.r.sex_inferred ? ", " + p.r.sex_inferred : "") });
      showTip(ev, lines);
    });
    f.svg.addEventListener("pointerleave", function () { if (hot !== null) dots[hot].setAttribute("r", 4); hot = null; hideTip(); });
    legend(container, groups.map(function (g, i) { return { color: COLORS[i], label: g.label }; }).filter(function (g) { return g.label; }));
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
    var ser = spec.series.filter(function (s) { return s.y && s.y.length; });
    if (!ser.length) return empty(container);
    var xs0 = [].concat.apply([], ser.map(function (s) { return s.x; })), ys0 = [].concat.apply([], ser.map(function (s) { return s.y.concat(s.lo || [], s.hi || []); })).filter(isFinite);
    var lo = Math.min.apply(null, ys0), hi = Math.max.apply(null, ys0); if (spec.ref !== undefined) { lo = Math.min(lo, spec.ref); hi = Math.max(hi, spec.ref); }
    var pad = (hi - lo || 1) * 0.1;
    var f = frame(container, 640, 280, { l: 56, t: 22, r: 16, b: 44 });
    var xs = scale(Math.min.apply(null, xs0), Math.max.apply(null, xs0), 0, f.W), ys = scale(lo - pad, hi + pad, f.H, 0);
    axes(f, xs, ys, spec.xlabel, spec.ylabel, function (v) { return fmt(v, 0); });
    if (spec.ref !== undefined) refLine(f, xs, ys, { y: spec.ref, label: spec.ref_label }, false);
    ser.forEach(function (s, si) {
      if (s.lo && s.hi) { var d = "M" + s.x.map(function (x, i) { return xs.map(x) + "," + ys.map(s.lo[i]); }).join(" L") + " L" + s.x.slice().reverse().map(function (x, i) { var j = s.x.length - 1 - i; return xs.map(x) + "," + ys.map(s.hi[j]); }).join(" L") + " z";
        el("path", { d: d, fill: COLORS[si], opacity: 0.12 }, f.g); }
      el("path", { d: "M" + s.x.map(function (x, i) { return xs.map(x) + "," + ys.map(s.y[i]); }).join(" L"), fill: "none", stroke: COLORS[si], "stroke-width": 2, "stroke-linejoin": "round", "stroke-linecap": "round" }, f.g);
      var n = s.x.length - 1; el("circle", { cx: xs.map(s.x[n]), cy: ys.map(s.y[n]), r: 4, fill: COLORS[si], stroke: "var(--surface)", "stroke-width": 2 }, f.g);
    });
    var cross = el("line", { y1: 0, y2: f.H, stroke: "var(--ink-2)", "stroke-width": 1, opacity: 0 }, f.g);
    f.svg.addEventListener("pointermove", function (ev) {
      var rect = f.svg.getBoundingClientRect(), k = 640 / rect.width, mx = (ev.clientX - rect.left) * k - f.m.l;
      var xv = Math.round(xs.lo + (mx / f.W) * (xs.hi - xs.lo)); if (xv < xs.lo || xv > xs.hi) { cross.setAttribute("opacity", 0); hideTip(); return; }
      cross.setAttribute("x1", xs.map(xv)); cross.setAttribute("x2", xs.map(xv)); cross.setAttribute("opacity", 0.6);
      var out = [{ b: spec.xlabel + " " + xv }];
      ser.forEach(function (s) { var i = s.x.indexOf(xv); if (i >= 0) out.push({ b: fmt(s.y[i], 3) + (s.lo ? " (" + fmt(s.lo[i], 2) + " to " + fmt(s.hi[i], 2) + ")" : ""), k: s.name }); });
      showTip(ev, out);
    });
    f.svg.addEventListener("pointerleave", function () { cross.setAttribute("opacity", 0); hideTip(); });
    legend(container, ser.map(function (s, i) { return { color: COLORS[i], label: s.name, line: true }; }));
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

  var TYPES = { hist: hist, scatter: scatter, strip: strip, lines: lines, meters: meters };
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
