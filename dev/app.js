/* HydratePRSV2 app logic — requires engine_port.js (E) and anchors.json (ANCHORS) */
(function () {
'use strict';
const G = E.GASES;                       // ['CH4','C2H6','C3H8','N2','CO2']
const LABEL = { CH4: 'CH₄', C2H6: 'C₂H₆', C3H8: 'C₃H₈', N2: 'N₂', CO2: 'CO₂' };
const FULLNAME = { CH4: 'methane', C2H6: 'ethane', C3H8: 'propane', N2: 'nitrogen', CO2: 'carbon dioxide' };
const PRESETS = {
  CH4: { CH4: 1 }, C2H6: { C2H6: 1 }, C3H8: { C3H8: 1 }, N2: { N2: 1 }, CO2: { CO2: 1 },
  NG: { CH4: 0.932, C2H6: 0.0425, C3H8: 0.0161, N2: 0.0043, CO2: 0.0051 },
  BIO: { CH4: 0.55, CO2: 0.40, N2: 0.05 }
};
const Q2 = { C2H6: [287.75, 3.39], CO2: [283.0, 4.5] };
const $ = id => document.getElementById(id);
const tick = () => new Promise(r => setTimeout(r, 15));

/* ---------- formatting ---------- */
function fmt(x, sig) {
  if (x === null || x === undefined || !isFinite(x)) return '—';
  sig = sig || 4;
  const a = Math.abs(x);
  if (a !== 0 && (a < 1e-3 || a >= 1e6)) return x.toExponential(Math.max(1, sig - 2));
  if (a >= 1000) return x.toFixed(Math.max(0, sig - 4));
  if (a >= 10) return x.toFixed(Math.max(1, sig - 2));
  return x.toPrecision(sig);
}
function niceTicks(lo, hi, n) {
  if (!(hi > lo)) { hi = lo + 1e-9; }
  const span = hi - lo, step0 = span / Math.max(1, n - 1);
  const mag = Math.pow(10, Math.floor(Math.log10(step0)));
  let step = mag;
  for (const m of [1, 2, 2.5, 5, 10]) { if (step0 <= m * mag) { step = m * mag; break; } }
  const t0 = Math.ceil(lo / step) * step, out = [];
  for (let v = t0; v <= hi + step * 1e-9; v += step) out.push(+v.toPrecision(12));
  return out;
}
function logTicks(lo, hi) {
  const out = [];
  const e0 = Math.floor(Math.log10(lo)), e1 = Math.ceil(Math.log10(hi));
  for (let e = e0; e <= e1; e++) {
    for (const m of [1, 2, 5]) { const v = m * Math.pow(10, e); if (v >= lo * 0.999 && v <= hi * 1.001) out.push(v); }
  }
  if (out.length < 3) return niceTicks(lo, hi, 6);
  return out;
}

/* ---------- SVG helpers ---------- */
const NS = 'http://www.w3.org/2000/svg';
function sv(tag, attrs, parent) {
  const el = document.createElementNS(NS, tag);
  for (const k in attrs) el.setAttribute(k, attrs[k]);
  if (parent) parent.appendChild(el);
  return el;
}
const STCOL = { sI: '#2563eb', sII: '#ea580c' };

function renderChart(container, pts, opts) {
  // pts: [{x,y,struct,tag,T,P}] ; opts: {logY, xlab, ylab, anchors:[{x,y}], onPick, xFrac:false}
  container.innerHTML = '';
  const W = container.clientWidth || 760, H = Math.max(340, Math.min(520, W * 0.55));
  const M = { l: 64, r: 18, t: 14, b: 44 };
  const svg = sv('svg', { width: W, height: H, viewBox: `0 0 ${W} ${H}` }, container);
  const good = pts.filter(p => isFinite(p.x) && isFinite(p.y) && p.y > 0);
  if (!good.length) { container.innerHTML = '<div class="hint">No equilibrium found in this range.</div>'; return; }
  let xlo = Math.min(...good.map(p => p.x)), xhi = Math.max(...good.map(p => p.x));
  if (xhi === xlo) { xlo -= 1; xhi += 1; }
  const xpad = (xhi - xlo) * 0.03; xlo -= xpad; xhi += xpad;
  let ylo, yhi, ymap;
  if (opts.logY) {
    ylo = Math.min(...good.map(p => p.y)) / 1.35; yhi = Math.max(...good.map(p => p.y)) * 1.35;
    if (opts.anchors && opts.anchors.length) { ylo = Math.min(ylo, Math.min(...opts.anchors.map(a => a.y)) / 1.35); yhi = Math.max(yhi, Math.max(...opts.anchors.map(a => a.y)) * 1.35); }
    ymap = v => Math.log10(v);
  } else {
    ylo = 0; yhi = Math.max(...good.map(p => p.y)) * 1.12;
    if (opts.anchors && opts.anchors.length) yhi = Math.max(yhi, Math.max(...opts.anchors.map(a => a.y)) * 1.12);
    ymap = v => v;
  }
  const X = x => M.l + (x - xlo) / (xhi - xlo) * (W - M.l - M.r);
  const Y = y => H - M.b - (ymap(y) - ymap(ylo)) / (ymap(yhi) - ymap(ylo)) * (H - M.t - M.b);
  // grid + axes
  const xt = niceTicks(xlo, xhi, 8), yt = opts.logY ? logTicks(ylo, yhi) : niceTicks(ylo, yhi, 6);
  for (const v of yt) {
    sv('line', { x1: M.l, x2: W - M.r, y1: Y(v), y2: Y(v), stroke: '#e8edf5', 'stroke-width': 1 }, svg);
    sv('text', { x: M.l - 7, y: Y(v) + 4, 'text-anchor': 'end', 'font-size': 11, fill: '#5b6779' }, svg)
      .textContent = opts.logY ? fmt(v, 2) : fmt(v, 3);
  }
  for (const v of xt) {
    sv('line', { x1: X(v), x2: X(v), y1: M.t, y2: H - M.b, stroke: '#eef1f7', 'stroke-width': 1 }, svg);
    sv('text', { x: X(v), y: H - M.b + 17, 'text-anchor': 'middle', 'font-size': 11, fill: '#5b6779' }, svg)
      .textContent = fmt(v, 4);
  }
  sv('line', { x1: M.l, x2: W - M.r, y1: H - M.b, y2: H - M.b, stroke: '#98a4b8' }, svg);
  sv('line', { x1: M.l, x2: M.l, y1: M.t, y2: H - M.b, stroke: '#98a4b8' }, svg);
  sv('text', { x: (M.l + W - M.r) / 2, y: H - 8, 'text-anchor': 'middle', 'font-size': 12, fill: '#334' }, svg).textContent = opts.xlab;
  const ylab = sv('text', { x: 16, y: (M.t + H - M.b) / 2, 'text-anchor': 'middle', 'font-size': 12, fill: '#334', transform: `rotate(-90 16 ${(M.t + H - M.b) / 2})` }, svg);
  ylab.textContent = opts.ylab;
  // 273.15 K marker (ice | liquid water)
  if (xlo < 273.15 && xhi > 273.15 && !opts.xFrac) {
    sv('line', { x1: X(273.15), x2: X(273.15), y1: M.t, y2: H - M.b, stroke: '#7c3aed', 'stroke-width': 1, 'stroke-dasharray': '4 4', opacity: .55 }, svg);
    sv('text', { x: X(273.15) + 4, y: M.t + 12, 'font-size': 10.5, fill: '#7c3aed' }, svg).textContent = 'ice ← 273.15 K → liquid water';
  }
  // anchors
  if (opts.anchors) for (const a of opts.anchors) {
    if (a.x < xlo || a.x > xhi || a.y < ylo || a.y > yhi) continue;
    sv('circle', { cx: X(a.x), cy: Y(a.y), r: 4.2, fill: '#fff', stroke: '#64748b', 'stroke-width': 1.6 }, svg);
  }
  // curve segments by structure (split where struct changes or gap in x)
  const sorted = good.slice().sort((p, q) => p.x - q.x);
  let seg = [];
  const flush = () => {
    if (seg.length > 1) {
      const d = seg.map((p, i) => `${i ? 'L' : 'M'}${X(p.x).toFixed(2)},${Y(p.y).toFixed(2)}`).join('');
      sv('path', { d, fill: 'none', stroke: STCOL[seg[0].struct] || '#334', 'stroke-width': 2.2, 'stroke-linejoin': 'round' }, svg);
    }
    seg = [];
  };
  for (let i = 0; i < sorted.length; i++) {
    const p = sorted[i];
    if (seg.length) {
      const prev = sorted[i - 1];
      const dx = Math.abs(p.x - prev.x) > (xhi - xlo) * 0.25;    // gap
      if (p.struct !== prev.struct || dx) flush();
    }
    seg.push(p);
  }
  flush();
  // points
  const tip = $('tooltip');
  for (const p of sorted) {
    const ice = p.T !== undefined && p.T < 273.15;
    const c = sv(ice ? 'path' : 'circle', ice
      ? { d: `M${X(p.x)},${Y(p.y) - 4.4} L${X(p.x) + 4.4},${Y(p.y)} L${X(p.x)},${Y(p.y) + 4.4} L${X(p.x) - 4.4},${Y(p.y)} Z`, fill: STCOL[p.struct] || '#334', stroke: '#fff', 'stroke-width': 1 }
      : { cx: X(p.x), cy: Y(p.y), r: 3.6, fill: STCOL[p.struct] || '#334', stroke: '#fff', 'stroke-width': 1.2 }, svg);
    const hit = sv('circle', { cx: X(p.x), cy: Y(p.y), r: 11, fill: 'transparent', style: 'cursor:pointer' }, svg);
    hit.addEventListener('mouseenter', ev => {
      const r = container.getBoundingClientRect();
      tip.style.display = 'block';
      tip.style.left = Math.min(W - 170, X(p.x) + 12) + 'px';
      tip.style.top = Math.max(4, Y(p.y) - 30) + 'px';
      tip.textContent = opts.tipFmt ? opts.tipFmt(p) : `T=${fmt(p.T, 5)} K\nP=${fmt(p.P / 1e6, 4)} MPa\n${p.struct}${p.tag === 'lc' ? ' · H–LW–LC' : (p.T < 273.15 ? ' · H–I–V' : ' · H–LW–V')}`;
    });
    hit.addEventListener('mouseleave', () => { tip.style.display = 'none'; });
    hit.addEventListener('click', () => { if (opts.onPick) opts.onPick(p); });
  }
}

/* ---------- composition helpers ---------- */
function readComp() {
  const gf = {}; let s = 0;
  for (const g of G) { const v = Math.max(0, parseFloat($('c_' + g).value) || 0); gf[g] = v; s += v; }
  if (s <= 0) return null;
  for (const g of G) gf[g] /= s;
  return gf;
}
function setComp(gf) {
  let s = 0; for (const g of G) s += (gf[g] || 0);
  for (const g of G) $('c_' + g).value = +(100 * (gf[g] || 0) / (s || 1)).toFixed(3);
}
function isPure(gf) { for (const g of G) if (gf[g] === 1) return g; return null; }
function normGasFrac(gf) { const out = {}; for (const g of G) if (gf[g] > 1e-12) out[g] = gf[g]; return out; }

/* ---------- progress ---------- */
function prog(fr, msg) {
  const el = $('progress'); el.style.display = 'block';
  const bar = el.firstElementChild || el.appendChild(document.createElement('div'));
  bar.style.width = (100 * fr).toFixed(1) + '%';
  if (msg !== undefined) { $('status').textContent = msg; $('status').className = 'status'; }
}
function progDone(msg) { $('progress').style.display = 'none'; $('status').textContent = msg || ''; }

/* ---------- computation modes ---------- */
let lastResults = null, lastMode = null;

async function runCurve() {
  const gf = readComp();
  if (!gf) { $('status').textContent = 'Enter a non-zero gas composition.'; $('status').className = 'status err'; return; }
  const Tmin = parseFloat($('tmin').value), Tmax = parseFloat($('tmax').value);
  const n = Math.max(3, Math.min(200, parseInt($('npts').value) || 41));
  if (!(Tmax > Tmin) || Tmin < 150 || Tmax > 340) { $('status').textContent = 'Temperature range invalid (need 150 K < Tmin < Tmax < 340 K).'; $('status').className = 'status err'; return; }
  const mode = $('structmode').value;
  const logY = $('logy').checked;
  const gfn = normGasFrac(gf);
  $('calc').disabled = true;
  const Tlist = Array.from({ length: n }, (_, i) => Tmin + (Tmax - Tmin) * i / (n - 1));
  const pts = [];
  let warm = null;
  try {
    for (let i = 0; i < n; i++) {
      if (i % 3 === 0) { prog(i / n, `Computing equilibrium… ${i + 1}/${n} T-points`); await tick(); }
      const T = Tlist[i];
      let res;
      if (mode === 'sI' || mode === 'sII') {
        const P = E.peqStruct(T, gfn, mode, warm);
        res = { P, struct: P ? mode : null, sI: mode === 'sI' ? P : null, sII: mode === 'sII' ? P : null };
      } else res = E.peqSelect(T, gfn, warm);
      if (res.P) {
        warm = res.P;
        const st = E.wgf(T, res.P, gfn);
        pts.push({ x: T, y: res.P / 1e6, T, P: res.P, struct: res.struct, tag: st ? st.tag : '', sI: res.sI, sII: res.sII });
      } else pts.push({ x: T, y: NaN, T, P: null, struct: null, tag: '' });
    }
  } catch (err) { console.error(err); progDone('Computation failed: ' + err.message); $('status').className = 'status err'; $('calc').disabled = false; return; }
  lastResults = pts; lastMode = 'curve';
  // anchors overlay for pure gas
  const pure = isPure(gf);
  let anchors = [];
  if (pure && $('showanch').checked && window.ANCHORS && window.ANCHORS[pure])
    anchors = window.ANCHORS[pure].map(a => ({ x: a[0], y: a[1] / 1e6 }));
  const okpts = pts.filter(p => p.P);
  renderChart($('chart'), okpts, {
    logY, xlab: 'Temperature  T  [K]', ylab: 'Dissociation pressure  P  [MPa]',
    anchors, onPick: p => showState(p.T, gfn, p.struct, p.P)
  });
  renderLegend(anchors.length);
  renderCurveTable(okpts);
  const regs = regionFlags(okpts, pure);
  progDone(`Done — ${okpts.length}/${n} points. ${regs}`);
  $('calc').disabled = false;
  if (pure) $('cmpRow').style.display = 'block'; else $('cmpRow').style.display = 'none';
}

function regionFlags(pts, pure) {
  const lc = pts.some(p => p.tag === 'lc');
  const ice = pts.some(p => p.T < 273.15);
  let s = '';
  if (ice) s += 'Includes H–I–V (ice) branch. ';
  if (lc) s += '⚠ Includes H–LW–LC states — see Validation notes (documented deviations). ';
  const structs = [];
  for (const p of pts) if (p.struct && structs[structs.length - 1] !== p.struct) structs.push(p.struct);
  if (structs.length === 1) s += `Structure: ${structs[0]}. `;
  else if (structs.length > 1) s += `Structure changes along the curve: ${structs.join(' → ')}. `;
  return s;
}

function renderLegend(hasAnch) {
  $('legend').innerHTML =
    `<span><span class="sw" style="background:${STCOL.sI}"></span>structure I</span>` +
    `<span><span class="sw" style="background:${STCOL.sII}"></span>structure II</span>` +
    `<span>◆ below 273.15 K = H–I–V (ice)</span>` +
    (hasAnch ? `<span><span class="sw" style="background:#fff;border:1.6px solid #64748b"></span>reference data anchors</span>` : '') +
    `<span class="small">click a point for full state details</span>`;
}

function renderCurveTable(pts) {
  let h = '<table class="dt"><tr><th>T [K]</th><th>P [MPa]</th><th>Structure</th><th>Region</th><th>P(sI) [MPa]</th><th>P(sII) [MPa]</th></tr>';
  for (const p of pts) {
    if (!p.P) continue;
    const reg = p.T < 273.15 ? 'H–I–V' : (p.tag === 'lc' ? 'H–LW–LC' : 'H–LW–V');
    h += `<tr><td class="num">${fmt(p.T, 5)}</td><td class="num">${fmt(p.P / 1e6, 4)}</td>` +
      `<td><span class="pill ${p.struct}">${p.struct}</span></td><td>${reg}</td>` +
      `<td class="num">${p.sI ? fmt(p.sI / 1e6, 4) : '—'}</td><td class="num">${p.sII ? fmt(p.sII / 1e6, 4) : '—'}</td></tr>`;
  }
  $('resultTable').innerHTML = h + '</table>';
}

async function runPx() {
  const A = $('pxA').value, B = $('pxB').value;
  if (A === B) { $('status').textContent = 'Choose two different gases.'; $('status').className = 'status err'; return; }
  const T = parseFloat($('pxT').value);
  if (!(T > 200 && T < 340)) { $('status').textContent = 'Temperature out of range.'; $('status').className = 'status err'; return; }
  const n = Math.max(10, Math.min(120, parseInt($('pxn').value) || 49));
  $('calc').disabled = true;
  // y_B grid: dense near ends (transitions often at ppm level)
  const fr = [];
  const ne = Math.max(6, Math.floor(n / 4));
  for (let i = 0; i <= ne; i++) fr.push(i === 0 ? 0 : Math.pow(10, -6 + 6 * i / ne));  // 0 .. 1 (log 1e-6..1)
  for (let i = 1; i < n - ne; i++) fr.push(i / (n - ne));
  fr.push(1);
  const uniq = [...new Set(fr.map(v => +v.toPrecision(6)))].sort((a, b) => a - b);
  const pts = [];
  let warm = null;
  try {
    for (let i = 0; i < uniq.length; i++) {
      if (i % 3 === 0) { prog(i / uniq.length, `Computing P–x diagram… ${i + 1}/${uniq.length}`); await tick(); }
      const yB = uniq[i];
      const gf = {}; gf[A] = 1 - yB; gf[B] = yB;
      const gfn = normGasFrac(gf);
      const res = E.peqSelect(T, gfn, warm);
      if (res.P) {
        warm = res.P;
        pts.push({ x: yB * 100, y: res.P / 1e6, T, P: res.P, struct: res.struct, yB, gfn });
      }
    }
  } catch (err) { console.error(err); progDone('Failed: ' + err.message); $('status').className = 'status err'; $('calc').disabled = false; return; }
  lastResults = pts; lastMode = 'px';
  renderChart($('chart'), pts, {
    logY: $('logy').checked, xFrac: true,
    xlab: `Vapor mole fraction y(${LABEL[B]})  [%]`, ylab: 'Equilibrium pressure  P  [MPa]',
    tipFmt: p => `y(${LABEL[B]}) = ${p.yB < 1e-4 ? p.yB.toExponential(2) : fmt(p.yB * 100, 3) + '%'}\nP = ${fmt(p.P / 1e6, 4)} MPa\nstructure ${p.struct}`,
    onPick: p => showState(T, p.gfn, p.struct, p.P)
  });
  $('legend').innerHTML =
    `<span><span class="sw" style="background:${STCOL.sI}"></span>sI</span><span><span class="sw" style="background:${STCOL.sII}"></span>sII</span>` +
    `<span class="small">T = ${fmt(T, 5)} K · click a point for details</span>`;
  // find transitions
  const trans = [];
  for (let i = 1; i < pts.length; i++) if (pts[i].struct !== pts[i - 1].struct)
    trans.push(`s${pts[i - 1].struct === 'sI' ? 'I' : 'II'}→s${pts[i].struct === 'sI' ? 'I' : 'II'} near y(${LABEL[B]}) ≈ ${fmt(pts[i].yB, 3)}`);
  let h = '<table class="dt"><tr><th>y(' + LABEL[B] + ')</th><th>P [MPa]</th><th>Structure</th></tr>';
  for (const p of pts) h += `<tr><td class="num">${p.yB === 0 ? '0' : (p.yB < 1e-4 ? p.yB.toExponential(2) : fmt(p.yB, 4))}</td><td class="num">${fmt(p.P / 1e6, 4)}</td><td><span class="pill ${p.struct}">${p.struct}</span></td></tr>`;
  $('resultTable').innerHTML = h + '</table>';
  progDone(trans.length ? 'Structure transitions: ' + trans.join('; ') : 'Single structure across the whole composition range.');
  $('cmpRow').style.display = 'none';
  $('calc').disabled = false;
}

async function runState() {
  const gf = readComp();
  if (!gf) { $('status').textContent = 'Enter a non-zero gas composition.'; $('status').className = 'status err'; return; }
  const T = parseFloat($('stateT').value);
  if (!(T > 150 && T < 340)) { $('status').textContent = 'Temperature out of range.'; $('status').className = 'status err'; return; }
  const mode = $('structmode').value;
  $('calc').disabled = true;
  prog(0.3, 'Computing equilibrium…'); await tick();
  const gfn = normGasFrac(gf);
  let res;
  if (mode === 'sI' || mode === 'sII') { const P = E.peqStruct(T, gfn, mode, null); res = { P, struct: P ? mode : null, sI: mode === 'sI' ? P : null, sII: mode === 'sII' ? P : null }; }
  else res = E.peqSelect(T, gfn, null);
  if (!res.P) { progDone('No three-phase equilibrium found at this temperature.'); $('status').className = 'status err'; $('chart').innerHTML = ''; $('resultTable').innerHTML = ''; $('statePanel').innerHTML = ''; $('calc').disabled = false; return; }
  showState(T, gfn, res.struct, res.P, res);
  // mini chart: both-structure bars
  $('chart').innerHTML = '';
  const div = document.createElement('div'); div.className = 'kv'; div.style.maxWidth = '420px';
  div.innerHTML = `<span class="k">Selected structure (lowest P)</span><span class="v"><span class="pill ${res.struct}">${res.struct}</span></span>` +
    `<span class="k">P(sI)</span><span class="v">${res.sI ? fmt(res.sI / 1e6, 4) + ' MPa' : '—'}</span>` +
    `<span class="k">P(sII)</span><span class="v">${res.sII ? fmt(res.sII / 1e6, 4) + ' MPa' : '—'}</span>`;
  $('chart').appendChild(div);
  $('resultTable').innerHTML = '';
  progDone('Done.');
  const pure = isPure(gf);
  $('cmpRow').style.display = pure ? 'block' : 'none';
  $('calc').disabled = false;
}

/* ---------- state details ---------- */
function showState(T, gfn, struct, P, pres) {
  const st = E.stateAt(T, gfn, struct, P);
  if (!st) { $('statePanel').innerHTML = '<div class="warnbox">State details unavailable (flash did not converge).</div>'; return; }
  const cages = E.CAGES[struct];
  const region = T < 273.15 ? 'H–I–V (ice)' : (st.tag === 'lc' ? 'H–LW–LC (aqueous + liquid hydrocarbon)' : 'H–LW–V (aqueous + vapor)');
  let h = `<h2>State details — T = ${fmt(T, 6)} K, P = ${fmt(P / 1e6, 5)} MPa <span class="pill ${struct}">${struct}</span> <span class="pill ${T < 273.15 ? 'ice' : (st.tag === 'lc' ? 'lc' : 'sI')}">${region}</span></h2>`;
  if (st.tag === 'lc') h += `<div class="warnbox">⚠ H–LW–LC region: the calibrated model carries documented large deviations here (see <b>Validation notes</b> tab). Pressure shown is the model's H–LW–LC equilibrium.</div>`;
  h += `<div class="cols2" style="margin-top:10px"><div>`;
  h += `<h3>Water fugacities</h3><div class="kv">` +
    `<span class="k">f<sub>w</sub><sup>L</sup> (fluid phase)</span><span class="v">${fmt(st.fLw, 5)} Pa</span>` +
    `<span class="k">f<sub>w</sub><sup>H</sup> (hydrate)</span><span class="v">${fmt(st.fH, 5)} Pa</span>` +
    `<span class="k">f<sub>w</sub><sup>β</sup> (empty lattice, P)</span><span class="v">${fmt(st.fbeta, 5)} Pa</span>` +
    `<span class="k">Σν ln(1−θ)</span><span class="v">${fmt(st.S, 5)}</span>` +
    `<span class="k">water activity a<sub>w</sub> = f<sub>L</sub>/f°<sub>w</sub></span><span class="v">${fmt(st.fLw / st.psat_w, 5)}</span>` +
    `<span class="k">hydration number (H₂O/guest)</span><span class="v">${fmt(st.hydrationNumber, 4)}</span>` +
    `</div>`;
  h += `<h3>Cage occupancies θ &amp; Langmuir constants C [Pa⁻¹]</h3>`;
  h += `<table class="dt"><tr><th>Guest</th><th>θ small (${fmt(cages[0][0] * 1e10, 3)} Å, z=${cages[0][1]})</th><th>θ large (${fmt(cages[1][0] * 1e10, 3)} Å, z=${cages[1][1]})</th><th>C small</th><th>C large</th><th>f<sub>gas</sub> [Pa]</th></tr>`;
  for (const g of Object.keys(st.fg)) {
    h += `<tr><td>${LABEL[g]}</td><td class="num">${fmt(st.theta[g][0], 4)}</td><td class="num">${fmt(st.theta[g][1], 4)}</td>` +
      `<td class="num">${fmt(st.C[g][0], 3)}</td><td class="num">${fmt(st.C[g][1], 3)}</td><td class="num">${fmt(st.fg[g], 4)}</td></tr>`;
  }
  h += `</table><div class="hint">ν(sI) = 2/46 small + 6/46 large; ν(sII) = 16/136 small + 8/136 large. θ = Cf/(1+ΣCf).</div>`;
  h += `</div><div>`;
  h += `<h3>Fluid phases at equilibrium</h3>`;
  if (st.y) {
    h += `<table class="dt"><tr><th>Comp.</th><th>vapor/HC-rich y</th><th>aqueous x</th><th>K = y/x</th></tr>`;
    for (let i = 0; i < 6; i++) {
      const nm = E.NAMES[i];
      const yv = st.y[i], xv = st.xaq ? st.xaq[i] : NaN;
      h += `<tr><td>${nm === 'H2O' ? 'H₂O' : LABEL[nm]}</td><td class="num">${fmt(yv, 4)}</td><td class="num">${st.xaq ? fmt(xv, 4) : '—'}</td><td class="num">${st.xaq && xv > 0 ? fmt(yv / xv, 3) : '—'}</td></tr>`;
    }
    h += `</table>`;
    h += `<div class="kv" style="margin-top:8px"><span class="k">water content of gas y(H₂O)</span><span class="v">${fmt(st.yw, 3)}</span>` +
      `<span class="k">≈ mg H₂O / m³ gas (ideal)</span><span class="v">${fmt(st.yw * 18.015 * P / (E.Rg * T) / 1e3 * 1e6, 3)}</span></div>`;
    const diss = Object.keys(gfn).map(g => `${LABEL[g]}: ${fmt(st.xaq ? st.xaq[E.NAMES.indexOf(g)] : NaN, 3)}`).join(' · ');
    h += `<div class="hint">Gas solubility in the aqueous phase (indicative — the frozen paper k<sub>ij</sub> underpredicts solubilities, see Validation notes): ${diss}</div>`;
  } else h += `<div class="hint">Ice branch: vapor composition y(H₂O) = ${fmt(st.yw !== null ? st.yw : (st.y ? st.y[E.IW] : NaN), 3)}</div>`;
  if (pres) h += `<div class="kv" style="margin-top:8px"><span class="k">P(sI)</span><span class="v">${pres.sI ? fmt(pres.sI / 1e6, 4) + ' MPa' : '—'}</span><span class="k">P(sII)</span><span class="v">${pres.sII ? fmt(pres.sII / 1e6, 4) + ' MPa' : '—'}</span></div>`;
  h += `</div></div>`;
  $('statePanel').innerHTML = h;
}

/* ---------- anchor comparison (pure gas) ---------- */
async function runCompare() {
  const gf = readComp(); const pure = gf && isPure(gf);
  if (!pure) return;
  const an = (window.ANCHORS && window.ANCHORS[pure]) || [];
  if (!an.length) return;
  $('cmpBtn').disabled = true;
  const rows = [];
  try {
    for (let i = 0; i < an.length; i++) {
      prog(i / an.length, `Comparing with reference data… ${i + 1}/${an.length}`); await tick();
      const [T, Pref, st] = an[i];
      const P = E.peqStruct(T, { [pure]: 1 }, st, Pref);
      let reg = T < 273.15 ? 'ice' : 'V';
      if (reg === 'V') { const w = E.wgf(T, Pref, { [pure]: 1 }); if (w && w.tag === 'lc') reg = 'LC'; }
      rows.push({ T, Pref, P, st, dev: P ? P / Pref - 1 : NaN, reg });
    }
  } catch (e) { progDone('Compare failed: ' + e.message); $('cmpBtn').disabled = false; return; }
  const aad = r => { const v = rows.filter(x => x.reg === r && isFinite(x.dev)); return v.length ? 100 * v.reduce((s, x) => s + Math.abs(x.dev), 0) / v.length : null; };
  const all = 100 * rows.filter(x => isFinite(x.dev)).reduce((s, x) => s + Math.abs(x.dev), 0) / rows.length;
  const paper = PAPER_AADP[pure];
  let h = `<h3>Model vs reference data — ${LABEL[pure]} (${FULLNAME[pure]})</h3>`;
  h += `<table class="dt"><tr><th>T [K]</th><th>P ref [MPa]</th><th>P model [MPa]</th><th>dev</th><th>region</th><th>struct</th></tr>`;
  for (const x of rows) {
    const cl = Math.abs(x.dev) < 0.05 ? 'ok' : (Math.abs(x.dev) < 0.15 ? 'warn' : 'bad');
    h += `<tr><td class="num">${fmt(x.T, 5)}</td><td class="num">${fmt(x.Pref / 1e6, 4)}</td><td class="num">${x.P ? fmt(x.P / 1e6, 4) : '—'}</td>` +
      `<td class="num" style="color:var(--${cl})">${isFinite(x.dev) ? (100 * x.dev).toFixed(1) + '%' : '—'}</td><td>${x.reg === 'ice' ? 'H–I–V' : (x.reg === 'LC' ? 'H–LW–LC' : 'H–LW–V')}</td><td>${x.st}</td></tr>`;
  }
  h += `</table><div class="infobox" style="margin-top:8px"><b>AADP (this model, anchor set):</b> overall ${fmt(all, 3)}%` +
    (aad('V') !== null ? ` · H–LW–V ${fmt(aad('V'), 3)}%` : '') + (aad('ice') !== null ? ` · H–I–V ${fmt(aad('ice'), 3)}%` : '') + (aad('LC') !== null ? ` · H–LW–LC ${fmt(aad('LC'), 3)}%` : '') +
    `<br><b>Paper (Table 3, larger data set):</b> ${paper}</div>`;
  h += `<div class="cols2" style="margin-top:10px"><div id="devPlot"></div><div id="parPlot"></div></div>`;
  $('cmpOut').innerHTML = h;
  devParityPlots(rows, pure);
  $('cmpCard').style.display = '';
  progDone('');
  $('cmpBtn').disabled = false;
}

/* ---------- paper reproduction helpers ---------- */
const PAPER_AADP = {
  CH4: 'overall 1.37% (H–I–V 0.7%, H–LW–V 2.05%); to 827 MPa: 2.86% (Table 4)',
  C2H6: 'overall 2.47% (H–I–V 3.14%, H–LW–V 1.15%, H–LW–LC 11.15%)',
  C3H8: 'overall 3.83% (H–I–V 2.13%, H–LW–V 3.41%)',
  N2: 'overall 6.84% (272–299.7 K: 3.14%, 300.6–314.5 K: 8.36%)',
  CO2: 'overall 2.89% (H–I–V 2.59%, H–LW–V 2.87%, H–LW–LC 10.2%)'
};

function loadPreset(name) {
  const gf = PRESETS[name];
  setComp(gf);
  document.querySelectorAll('.presets button.mini').forEach(b => b.classList.toggle('on', b.dataset.p === name));
}
function switchTab(id) {
  document.querySelectorAll('nav#tabs button').forEach(b => b.classList.toggle('on', b.dataset.t === id));
  document.querySelectorAll('.tabpane').forEach(p => p.classList.toggle('on', p.id === 'tab-' + id));
}
function setMode(m) {
  document.querySelectorAll('#modeRow button.mini').forEach(b => b.classList.toggle('on', b.dataset.m === m));
  $('modeCurve').style.display = (m === 'curve' || m === 'overview') ? 'block' : 'none';
  $('modePx').style.display = m === 'px' ? 'block' : 'none';
  $('modeState').style.display = m === 'state' ? 'block' : 'none';
  $('curveOpts').style.display = m === 'state' ? 'none' : 'block';
  window.CURMODE = m;
}

async function paperFig(which) {
  switchTab('calc');
  if (which === 'fig3') { loadPreset('C2H6'); setMode('curve'); $('tmin').value = 263.15; $('tmax').value = 287.5; $('npts').value = 33; $('showanch').checked = true; await runCurve(); }
  else if (which === 'fig4') { loadPreset('CO2'); setMode('curve'); $('tmin').value = 273.15; $('tmax').value = 291.5; $('npts').value = 30; $('showanch').checked = true; await runCurve(); }
  else if (which === 'fig5') { setMode('px'); $('pxA').value = 'CH4'; $('pxB').value = 'C3H8'; $('pxT').value = 277.6; $('pxn').value = 49; await runPx(); }
  else if (which === 'fig6') { setMode('px'); $('pxA').value = 'CH4'; $('pxB').value = 'C2H6'; $('pxT').value = 277.6; $('pxn').value = 49; await runPx(); }
  else if (which === 'fig7') { loadPreset('NG'); setMode('curve'); $('tmin').value = 275.15; $('tmax').value = 297.15; $('npts').value = 30; $('showanch').checked = false; await runCurve(); }
  else if (which === 'fig1') { switchTab('calc'); await runSolub(); }
  else if (which === 'ch4') { loadPreset('CH4'); setMode('curve'); $('tmin').value = 258.15; $('tmax').value = 303.15; $('npts').value = 46; $('showanch').checked = true; await runCurve(); }
}

async function occTable() {
  const specs = [
    ['CH4', 273.2, 0.969, 0.935, '0.967 / 0.947 (Hsieh), 0.971 / 0.92 (Sum)'],
    ['CH4', 276.7, 0.981, 0.917, '0.972 / 0.955 (Hsieh), 0.973 / 0.866 (Sum)'],
    ['CO2', 273.2, 0.986, 0.581, '0.985 / 0.388 (Hsieh), Davidson θL/θS = 8.33']
  ];
  let h = '<table class="dt"><tr><th>Guest</th><th>T [K]</th><th>P<sub>eq</sub> [MPa]</th><th>θ<sub>L</sub> this app</th><th>θ<sub>S</sub> this app</th><th>θ<sub>L</sub>/θ<sub>S</sub> paper</th><th>other models (paper)</th></tr>';
  for (const [g, T, thl, ths, others] of specs) {
    const P = E.peqStruct(T, { [g]: 1 }, 'sI', null);
    const st = P ? E.stateAt(T, { [g]: 1 }, 'sI', P) : null;
    h += `<tr><td>${LABEL[g]}</td><td class="num">${T}</td><td class="num">${P ? fmt(P / 1e6, 4) : '—'}</td>` +
      `<td class="num">${st ? fmt(st.theta[g][1], 4) : '—'}</td><td class="num">${st ? fmt(st.theta[g][0], 4) : '—'}</td>` +
      `<td class="num">${thl} / ${ths}</td><td class="small">${others}</td></tr>`;
  }
  $('occTable').innerHTML = h + '</table><div class="hint">Computed live by this app at its own equilibrium pressures; paper values from Tables 6–8 ("This work" column). Agreement ≤ 1.3% absolute in θ.</div>';
}


/* ---------- multi-series chart + overview / solubility / dev-parity plots ---------- */
const OVCOL = { CH4: '#1f77b4', C2H6: '#d62728', C3H8: '#2ca02c', N2: '#7f7f7f', CO2: '#9467bd' };
const SOLPT = {
  CH4: [[273.15, 2.6, 1.05e-3], [298.15, 0.101325, 2.5e-5], [273.15, 0.101325, 4.15e-5], [283.15, 5.0, 1.55e-3]],
  C2H6: [[298.15, 0.101325, 2.0e-5], [273.15, 0.101325, 4.6e-5], [277.6, 0.9, 1.3e-3]],
  C3H8: [[298.15, 0.101325, 2.6e-6], [273.15, 0.101325, 3.5e-6]],
  N2: [[298.15, 0.101325, 1.15e-5], [273.15, 0.101325, 1.35e-5], [277.15, 10, 1.0e-3]],
  CO2: [[278.15, 3.5, 1.2e-2], [298.15, 0.101325, 6.0e-4], [273.15, 0.101325, 1.34e-3]]
};
function renderMulti(container, series, opts) {
  container.innerHTML = '';
  const W = container.clientWidth || 760, H = Math.max(340, Math.min(520, W * 0.55));
  const M = { l: 70, r: 16, t: 14, b: 46 };
  const svg = sv('svg', { width: W, height: H, viewBox: `0 0 ${W} ${H}` }, container);
  const allp = series.flatMap(s => s.pts.filter(p => isFinite(p.x) && isFinite(p.y) && p.y > 0));
  const allm = series.flatMap(s => (s.marks || []).filter(p => isFinite(p.x) && isFinite(p.y) && p.y > 0));
  const both = allp.concat(allm);
  if (!both.length) { container.innerHTML = '<div class="hint">No data.</div>'; return; }
  let xlo = Math.min(...both.map(p => p.x)), xhi = Math.max(...both.map(p => p.x));
  let ylo = Math.min(...both.map(p => p.y)), yhi = Math.max(...both.map(p => p.y));
  const xm = opts.logX ? v => Math.log10(v) : v => v, ym = opts.logY ? v => Math.log10(v) : v => v;
  if (opts.logX) { xlo /= 1.3; xhi *= 1.3; } else { const p = (xhi - xlo) * 0.04 || 1; xlo -= p; xhi += p; }
  if (opts.logY) { ylo /= 2.2; yhi *= 2.2; } else { const p = (yhi - ylo) * 0.1 || 1; ylo -= p; yhi += p; }
  const X = x => M.l + (xm(x) - xm(xlo)) / (xm(xhi) - xm(xlo)) * (W - M.l - M.r);
  const Y = y => H - M.b - (ym(y) - ym(ylo)) / (ym(yhi) - ym(ylo)) * (H - M.t - M.b);
  const xt = opts.logX ? logTicks(xlo, xhi) : niceTicks(xlo, xhi, 8);
  const yt = opts.logY ? logTicks(ylo, yhi) : niceTicks(ylo, yhi, 6);
  for (const v of yt) {
    sv('line', { x1: M.l, x2: W - M.r, y1: Y(v), y2: Y(v), stroke: '#e8edf5' }, svg);
    sv('text', { x: M.l - 7, y: Y(v) + 4, 'text-anchor': 'end', 'font-size': 11, fill: '#5b6779' }, svg).textContent = fmt(v, 3);
  }
  for (const v of xt) {
    sv('line', { x1: X(v), x2: X(v), y1: M.t, y2: H - M.b, stroke: '#eef1f7' }, svg);
    sv('text', { x: X(v), y: H - M.b + 17, 'text-anchor': 'middle', 'font-size': 11, fill: '#5b6779' }, svg).textContent = fmt(v, 4);
  }
  (opts.hlines || []).forEach(hl => {
    if (hl.y < ylo || hl.y > yhi) return;
    sv('line', { x1: M.l, x2: W - M.r, y1: Y(hl.y), y2: Y(hl.y), stroke: hl.color || '#bbb', 'stroke-dasharray': '5,4' }, svg);
    sv('text', { x: W - M.r - 4, y: Y(hl.y) - 3, 'text-anchor': 'end', 'font-size': 10, fill: hl.color || '#888' }, svg).textContent = hl.label || '';
  });
  (opts.lines || []).forEach(ln => {
    const p1 = { x: ln.x0, y: ln.f(ln.x0) }, p2 = { x: ln.x1, y: ln.f(ln.x1) };
    sv('line', { x1: X(p1.x), y1: Y(p1.y), x2: X(p2.x), y2: Y(p2.y), stroke: ln.color, 'stroke-dasharray': ln.dash || '', 'stroke-width': 1.2 }, svg);
  });
  sv('line', { x1: M.l, x2: W - M.r, y1: H - M.b, y2: H - M.b, stroke: '#98a4b8' }, svg);
  sv('line', { x1: M.l, x2: M.l, y1: M.t, y2: H - M.b, stroke: '#98a4b8' }, svg);
  sv('text', { x: (M.l + W - M.r) / 2, y: H - 8, 'text-anchor': 'middle', 'font-size': 12, fill: '#334' }, svg).textContent = opts.xlab;
  sv('text', { x: 16, y: (M.t + H - M.b) / 2, 'text-anchor': 'middle', 'font-size': 12, fill: '#334', transform: `rotate(-90 16 ${(M.t + H - M.b) / 2})` }, svg).textContent = opts.ylab;
  for (const s of series) {
    const good = s.pts.filter(p => isFinite(p.x) && isFinite(p.y) && p.y > 0);
    if (good.length > 1) {
      const d = good.map((p, i) => (i ? 'L' : 'M') + X(p.x).toFixed(1) + ' ' + Y(p.y).toFixed(1)).join(' ');
      sv('path', { d, fill: 'none', stroke: s.color, 'stroke-width': 2, 'stroke-dasharray': s.dash || '' }, svg);
    }
    (s.marks || []).forEach(p => {
      if (!(isFinite(p.x) && isFinite(p.y) && p.y > 0)) return;
      sv('circle', { cx: X(p.x), cy: Y(p.y), r: 3.4, fill: '#fff', stroke: s.color, 'stroke-width': 1.6 }, svg);
    });
  }
}
function multiLegend(container, series) {
  container.innerHTML = series.map(s =>
    `<span><span class="sw" style="background:${s.color}"></span>${s.name}</span>` +
    (s.hasMarks ? `<span><span class="sw" style="background:#fff;border:1.6px solid ${s.color}"></span>${s.markName || 'data'}</span>` : '')
  ).join('') + `<span class="small">${NOTE || ''}</span>`;
}
let NOTE = '';
async function runOverview() {
  $('calc').disabled = true;
  const Tmin = parseFloat($('tmin').value) || 263.15, Tmax = parseFloat($('tmax').value) || 303.15;
  const n = Math.max(9, Math.min(60, parseInt($('npts').value) || 21));
  const series = [];
  for (let k = 0; k < G.length; k++) {
    const g = G[k];
    prog(k / G.length, `Computing ${LABEL[g]} curve…`); await tick();
    const pts = [];
    for (let i = 0; i < n; i++) {
      const T = Tmin + (Tmax - Tmin) * i / (n - 1);
      const r = E.peqSelect(T, { [g]: 1 }, null);
      if (r.P) pts.push({ x: T, y: r.P / 1e6 });
    }
    const marks = (window.ANCHORS && window.ANCHORS[g] || []).map(a => ({ x: a[0], y: a[1] / 1e6 }));
    series.push({ name: LABEL[g], color: OVCOL[g], pts, marks, hasMarks: true, markName: 'exp' });
  }
  NOTE = 'all five pure-gas dissociation curves vs reference data (auto structure selection)';
  renderMulti($('chart'), series, { xlab: 'Temperature  T  [K]', ylab: 'Dissociation pressure  P  [MPa]', logY: $('logy').checked });
  multiLegend($('legend'), series);
  $('resultTable').innerHTML = '<div class="hint">Overview mode: five pure-gas curves with reference-data markers. Use Single-gas presets for tables and state details.</div>';
  $('statePanel').innerHTML = '';
  progDone('Overview done — 5 curves.');
  $('calc').disabled = false;
}
async function runSolub() {
  $('calc').disabled = true;
  const T = 278.15;
  const series = [];
  for (let k = 0; k < G.length; k++) {
    const g = G[k];
    prog(k / G.length, `Computing ${LABEL[g]} solubility…`); await tick();
    const pts = [];
    for (let i = 0; i < 16; i++) {
      const P = 0.05e6 * Math.pow(160, i / 15);   // 0.05 -> 8 MPa log
      const z = new Array(6).fill(0); z[E.IW] = 0.99; z[E.NAMES.indexOf(g)] = 0.01;
      const fl = E.flash_TP(z, T, P, 'water_first');
      if (!fl || !fl.x || !fl.y) continue;
      const xaq = fl.x[E.IW] >= fl.y[E.IW] ? fl.x : fl.y;
      const xg = xaq[E.NAMES.indexOf(g)];
      if (xg > 0) pts.push({ x: P / 1e6, y: xg });
    }
    const marks = (SOLPT[g] || []).map(a => ({ x: a[1], y: a[2] }));
    series.push({ name: LABEL[g], color: OVCOL[g], pts, marks, hasMarks: true, markName: 'exp anchors' });
  }
  NOTE = 'gas solubility in water at 278.15 K (paper Fig. 1 style); markers at their own T';
  renderMulti($('chart'), series, { xlab: 'Pressure  P  [MPa]', ylab: 'mole fraction gas in aqueous phase', logX: true, logY: true });
  multiLegend($('legend'), series);
  $('resultTable').innerHTML = '<div class="hint">Solubility mode: model curves at 278.15 K; open circles are experimental anchors (each at its own temperature, see Validation notes).</div>';
  progDone('Solubility diagram done.');
  $('calc').disabled = false;
}
function devParityPlots(rows, pure) {
  const dev = rows.filter(r => isFinite(r.dev));
  const colOf = r => r.reg === 'ice' ? '#1f77b4' : (r.reg === 'LC' ? '#e67e22' : '#2ca02c');
  const dc = $('devPlot'), pc = $('parPlot');
  if (!dc || !pc) return;
  dc.innerHTML = ''; pc.innerHTML = '';
  // deviation plot
  {
    const W = dc.clientWidth || 420, H = 240, M = { l: 52, r: 12, t: 12, b: 34 };
    const svg = sv('svg', { width: W, height: H, viewBox: `0 0 ${W} ${H}` }, dc);
    const Tlo = Math.min(...dev.map(r => r.T)) - 2, Thi = Math.max(...dev.map(r => r.T)) + 2;
    let dmax = Math.max(30, ...dev.map(r => Math.abs(100 * r.dev))) * 1.15;
    const X = T => M.l + (T - Tlo) / (Thi - Tlo) * (W - M.l - M.r);
    const Y = d => H - M.b - (d + dmax) / (2 * dmax) * (H - M.t - M.b);
    [-10, 0, 10].forEach(v => {
      sv('line', { x1: M.l, x2: W - M.r, y1: Y(v), y2: Y(v), stroke: v === 0 ? '#98a4b8' : '#d8b4b4', 'stroke-dasharray': v === 0 ? '' : '5,4' }, svg);
      sv('text', { x: M.l - 6, y: Y(v) + 4, 'text-anchor': 'end', 'font-size': 10, fill: '#5b6779' }, svg).textContent = v + '%';
    });
    for (const T of niceTicks(Tlo, Thi, 6)) sv('text', { x: X(T), y: H - M.b + 15, 'text-anchor': 'middle', 'font-size': 10, fill: '#5b6779' }, svg).textContent = fmt(T, 4);
    dev.forEach(r => sv('circle', { cx: X(r.T), cy: Y(100 * r.dev), r: 3.6, fill: colOf(r), stroke: '#fff', 'stroke-width': 0.8 }, svg));
    sv('text', { x: (M.l + W - M.r) / 2, y: H - 5, 'text-anchor': 'middle', 'font-size': 11, fill: '#334' }, svg).textContent = 'T [K] — deviation of P_eq from reference (±10% bands)';
  }
  // parity plot
  {
    const W = pc.clientWidth || 420, H = 240, M = { l: 56, r: 12, t: 12, b: 36 };
    const svg = sv('svg', { width: W, height: H, viewBox: `0 0 ${W} ${H}` }, pc);
    const good = dev.filter(r => r.P);
    const lo = Math.min(...good.map(r => Math.min(r.Pref, r.P))) / 1e6 / 1.5, hi = Math.max(...good.map(r => Math.max(r.Pref, r.P))) / 1e6 * 1.5;
    const X = v => M.l + (Math.log10(v) - Math.log10(lo)) / (Math.log10(hi) - Math.log10(lo)) * (W - M.l - M.r);
    const Y = v => H - M.b - (Math.log10(v) - Math.log10(lo)) / (Math.log10(hi) - Math.log10(lo)) * (H - M.t - M.b);
    for (const v of logTicks(lo, hi)) {
      sv('line', { x1: X(v), x2: X(v), y1: M.t, y2: H - M.b, stroke: '#eef1f7' }, svg);
      sv('line', { x1: M.l, x2: W - M.r, y1: Y(v), y2: Y(v), stroke: '#eef1f7' }, svg);
      sv('text', { x: X(v), y: H - M.b + 14, 'text-anchor': 'middle', 'font-size': 10, fill: '#5b6779' }, svg).textContent = fmt(v, 3);
      sv('text', { x: M.l - 6, y: Y(v) + 4, 'text-anchor': 'end', 'font-size': 10, fill: '#5b6779' }, svg).textContent = fmt(v, 3);
    }
    [[1, '#98a4b8', ''], [1.1, '#d8b4b4', '4,3'], [0.9091, '#d8b4b4', '4,3'], [1.25, '#e8d4b4', '2,3'], [0.8, '#e8d4b4', '2,3']].forEach(L =>
      sv('line', { x1: X(lo), y1: Y(lo * L[0]), x2: X(hi), y2: Y(hi * L[0]), stroke: L[1], 'stroke-dasharray': L[2] }, svg));
    good.forEach(r => sv('circle', { cx: X(r.Pref / 1e6), cy: Y(r.P / 1e6), r: 3.6, fill: colOf(r), stroke: '#fff', 'stroke-width': 0.8 }, svg));
    sv('text', { x: (M.l + W - M.r) / 2, y: H - 5, 'text-anchor': 'middle', 'font-size': 11, fill: '#334' }, svg).textContent = 'P reference [MPa] (parity, ±10% / ±25% bands)';
  }
}

/* ---------- init ---------- */
function init() {
  document.querySelectorAll('nav#tabs button').forEach(b => b.addEventListener('click', () => {
    switchTab(b.dataset.t);
    if (b.dataset.t === 'paper' && !$('occTable').innerHTML) occTable();
  }));
  document.querySelectorAll('#modeRow button.mini').forEach(b => b.addEventListener('click', () => setMode(b.dataset.m)));
  document.querySelectorAll('.presets button.mini').forEach(b => b.addEventListener('click', () => { loadPreset(b.dataset.p); $('cmpRow').style.display = 'none'; }));
  $('calc').addEventListener('click', () => { if (window.CURMODE === 'px') runPx(); else if (window.CURMODE === 'state') runState(); else if (window.CURMODE === 'overview') runOverview(); else runCurve(); });
  $('cmpBtn').addEventListener('click', runCompare);
  document.querySelectorAll('[data-fig]').forEach(b => b.addEventListener('click', () => paperFig(b.dataset.fig)));
  setMode('curve');
  loadPreset('CH4');
  window.addEventListener('resize', () => { if (lastResults && lastMode === 'curve') { /* keep simple: no auto-rerender */ } });
  // initial run
  runCurve();
}
if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', init); else init();
})();
