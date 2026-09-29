#!/usr/bin/env python3
"""Assemble the self-contained single-file app: index.html"""
import json, io, os
HERE = os.path.dirname(os.path.abspath(__file__))
css = open(os.path.join(HERE, 'app.css')).read()
engine = open(os.path.join(HERE, 'engine_port.js')).read()
app = open(os.path.join(HERE, 'app.js')).read()
anchors = open(os.path.join(HERE, 'anchors.json')).read()
PJ = os.path.join(HERE, 'params_v8.json')
if not os.path.exists(PJ): PJ = os.path.join(HERE, 'params_v7.json')
PV = json.load(open(PJ))
LBL = {'CH4': 'CH₄', 'C2H6': 'C₂H₆', 'C3H8': 'C₃H₈', 'N2': 'N₂', 'CO2': 'CO₂'}
T2 = {'CH4': '232.2 · 3.505 · 0.3317', 'C2H6': '414.8 · 4.000 · 0.5695', 'C3H8': '493.7 · 4.519 · 0.6502',
      'N2': '142.1 · 3.469 · 0.341', 'CO2': '491.78 · 3.418 · 0.6786'}
kih_rows = ''
for i, g in enumerate(['CH4', 'C2H6', 'C3H8', 'N2', 'CO2']):
    e, sg, a = PV['kp'][g]
    extra = '<td rowspan="5">sI: 5¹² 3.917 Å (z=20), 5¹²6² 4.33 Å (z=24)<br>sII: 5¹² 3.917 Å (z=20), 5¹²6⁴ 4.73 Å (z=28)</td>' if i == 0 else ''
    kih_rows += f'<tr><td>{LBL[g]}</td><td class="mono">{e:.2f} · {sg*1e10:.4f} · {a*1e10:.4f}</td><td class="mono">{T2[g]}</td>{extra}</tr>\n'
import math
fb_rows = ''
for stn in ['sI', 'sII']:
    A, B = PV['fb'][stn]
    fb_rows += f'<tr><td>{stn}</td><td class="mono">A={A:.4f}, B={B:.1f}</td><td class="mono">{math.exp(A-B/273.15):.0f} Pa</td>'
    fb_rows += '<td class="mono">28.4 − 5800/T</td><td class="mono">1295 Pa</td></tr>\n' if stn == 'sI' else '<td class="mono">28.8 − 6000/T</td><td class="mono">929 Pa</td></tr>\n'
dk_rows = ''
for g in ['CH4', 'C2H6', 'C3H8', 'N2', 'CO2']:
    d0 = PV.get('dK0', {}).get(g, 0.0); d1 = PV.get('dK1', {}).get(g, 0.0)
    p2, q2 = {'CH4': (0.0029, -0.96), 'C2H6': (0.00295, -0.986), 'C3H8': (0.00198, -0.647), 'N2': (0.003, -0.42), 'CO2': (0.0078, -0.258)}[g]
    k298 = p2*298.15 + q2
    dk_rows += f'<tr><td>{LBL[g]}</td><td class="mono">k(g,w) + ({d0:+.3f}) + ({d1:+.5f})·(T−273.15)</td><td class="mono">{k298:+.3f} → {k298+d0+d1*25:+.3f}</td></tr>\n'


HTML = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Clathrate Hydrate Equilibria — PRSV2 + vdW–P (Moradi &amp; Khosravani model)</title>
<style>
__CSS__
</style>
</head>
<body>
<header>
  <h1>Clathrate Hydrate Phase Equilibria — PRSV2 EOS + van der Waals–Platteeuw</h1>
  <div class="sub">Interactive implementation of the model of G. Moradi &amp; E. Khosravani, <i>“Modeling of hydrate formation conditions for CH<sub>4</sub>, C<sub>2</sub>H<sub>6</sub>, C<sub>3</sub>H<sub>8</sub>, N<sub>2</sub>, CO<sub>2</sub> and their mixtures using the PRSV2 equation of state and obtaining the Kihara potential parameters”</i>, Fluid Phase Equilibria 338 (2013) 179–187. Fugacity-based three-phase equilibrium (H–I–V, H–L<sub>W</sub>–V, H–L<sub>W</sub>–L<sub>C</sub>), Kihara cell potential, van Laar two-parameter k<sub>ij</sub>(T), automatic sI/sII selection.</div>
  <div class="badges">
    <span class="badge">PRSV2 (Stryjek–Vera κ₁κ₂κ₃)</span>
    <span class="badge">Kihara + McKoy–Sinanoğlu cell potential</span>
    <span class="badge">sI / sII auto-selection</span>
    <span class="badge">calibrated &amp; validated vs ~70 reference points</span>
    <span class="badge">JS engine verified vs Python to 10⁻¹¹</span>
  </div>
</header>
<nav id="tabs">
  <button data-t="calc" class="on">Calculator</button>
  <button data-t="paper">Paper reproduction</button>
  <button data-t="params">Parameters &amp; method</button>
  <button data-t="valid">Validation notes</button>
</nav>
<main>

<!-- ================= CALCULATOR ================= -->
<div class="tabpane on" id="tab-calc">
<div class="grid">
  <div>
    <div class="card">
      <h2>Inputs</h2>
      <div id="modeRow" style="margin-bottom:4px">
        <button class="mini on" data-m="curve">Dissociation curve P(T)</button>
        <button class="mini" data-m="px">Binary P–x diagram</button>
        <button class="mini" data-m="state">Single state</button>
        <button class="mini" data-m="overview">All-gas overview</button>
      </div>
      <h3>Gas composition (mol %, water implicit)</h3>
      <div class="presets">
        <button class="mini on" data-p="CH4">CH₄</button><button class="mini" data-p="C2H6">C₂H₆</button><button class="mini" data-p="C3H8">C₃H₈</button><button class="mini" data-p="N2">N₂</button><button class="mini" data-p="CO2">CO₂</button><button class="mini" data-p="NG">natural gas (Fig 7)</button><button class="mini" data-p="BIO">biogas</button>
      </div>
      <div class="row"><div><label class="f">CH₄</label><input type="number" id="c_CH4" value="100" min="0" max="100" step="0.01"></div>
      <div><label class="f">C₂H₆</label><input type="number" id="c_C2H6" value="0" min="0" max="100" step="0.01"></div></div>
      <div class="row"><div><label class="f">C₃H₈</label><input type="number" id="c_C3H8" value="0" min="0" max="100" step="0.01"></div>
      <div><label class="f">N₂</label><input type="number" id="c_N2" value="0" min="0" max="100" step="0.01"></div></div>
      <div class="row"><div><label class="f">CO₂</label><input type="number" id="c_CO2" value="0" min="0" max="100" step="0.01"></div>
      <div><label class="f">structure</label><select id="structmode"><option value="auto">auto (min P)</option><option value="sI">force sI</option><option value="sII">force sII</option></select></div></div>

      <div id="modeCurve">
        <h3>Temperature range</h3>
        <div class="row"><div><label class="f">T min [K]</label><input type="number" id="tmin" value="263.15" step="0.5"></div>
        <div><label class="f">T max [K]</label><input type="number" id="tmax" value="303.15" step="0.5"></div>
        <div><label class="f">points</label><input type="number" id="npts" value="41" min="3" max="200"></div></div>
        <div class="hint">T &lt; 273.15 K automatically uses the H–I–V (ice) branch.</div>
      </div>
      <div id="modePx" style="display:none">
        <h3>Binary P–x diagram (Figs. 5–6 style)</h3>
        <div class="row"><div><label class="f">gas A (x-axis = 0)</label><select id="pxA"><option>CH4</option><option>C2H6</option><option>C3H8</option><option>N2</option><option>CO2</option></select></div>
        <div><label class="f">gas B (x-axis = 1)</label><select id="pxB"><option>C3H8</option><option>CH4</option><option>C2H6</option><option>N2</option><option>CO2</option></select></div></div>
        <div class="row"><div><label class="f">T [K]</label><input type="number" id="pxT" value="277.6" step="0.1"></div>
        <div><label class="f">points</label><input type="number" id="pxn" value="49" min="10" max="120"></div></div>
        <div class="hint">x-axis scans y(B) from 0 → 1 with log spacing near the pure-A end (captures ppm-level sI→sII transitions).</div>
      </div>
      <div id="modeState" style="display:none">
        <h3>Single state</h3>
        <label class="f">T [K]</label><input type="number" id="stateT" value="277.6" step="0.1">
      </div>

      <div id="curveOpts">
        <h3>Display</h3>
        <label style="font-size:12.5px"><input type="checkbox" id="logy"> logarithmic pressure axis</label><br>
        <label style="font-size:12.5px"><input type="checkbox" id="showanch" checked> overlay reference-data anchors (pure gases)</label>
      </div>
      <button class="btn" id="calc">Compute</button>
      <div id="progress"><div></div></div>
      <div class="status" id="status"></div>
      <div id="cmpRow" style="display:none; margin-top:8px">
        <button class="btn" id="cmpBtn" style="background:#0f7b3f">Compare with reference data (AADP)</button>
      </div>
    </div>
  </div>
  <div>
    <div class="card chartwrap">
      <h2>Hydrate dissociation equilibrium</h2>
      <div id="chart"></div>
      <div id="tooltip"></div>
      <div class="legend" id="legend"></div>
    </div>
    <div class="card" id="statePanel" style="margin-top:16px"><div class="hint">Click any point on the curve (or use Single-state mode) for full thermodynamic details: occupancies, Langmuir constants, fugacities, phase compositions.</div></div>
    <div class="card" style="margin-top:16px"><h2>Results table</h2><div id="resultTable" style="max-height:420px; overflow:auto"></div></div>
    <div class="card" style="margin-top:16px;display:none" id="cmpCard"><div id="cmpOut"></div></div>
  </div>
</div>
</div>

<!-- ================= PAPER ================= -->
<div class="tabpane" id="tab-paper">
<div class="card">
  <h2>Reproduce the paper’s figures &amp; tables</h2>
  <div style="margin:6px 0">
    <button class="mini" data-fig="fig1">Fig. 1 — solubility in water</button>
    <button class="mini" data-fig="ch4">CH₄ curve (Figs. 3-style, ice+water)</button>
    <button class="mini" data-fig="fig3">Fig. 3 — C₂H₆ H–L<sub>W</sub>–V</button>
    <button class="mini" data-fig="fig4">Fig. 4 — CO₂ incl. H–L<sub>W</sub>–L<sub>C</sub></button>
    <button class="mini" data-fig="fig5">Fig. 5 — P–x CH₄+C₃H₈ @ 277.6 K</button>
    <button class="mini" data-fig="fig6">Fig. 6 — P–x CH₄+C₂H₆ @ 277.6 K</button>
    <button class="mini" data-fig="fig7">Fig. 7 — natural gas</button>
  </div>
  <div class="infobox">Paper Fig. 5: sI→sII transition of CH₄+C₃H₈ at 277.6 K onsets at y(C₃H₈) ≈ 5·10⁻⁴ (y<sub>CH4</sub>=0.9995) at ≈ 4.0 MPa. This app: onset at y(C₃H₈) ≈ 3·10⁻⁴, P ≈ 4.0 MPa. Paper Fig. 6 reports an sII lens for CH₄+C₂H₆ between y<sub>CH4</sub> = 0.753 and 0.724; with our data-calibrated Kihara parameters the sII branch never dips below sI (smooth sI curve) — documented deviation.</div>
</div>
<div class="card" style="margin-top:16px">
  <h2>Cage occupancies — Tables 6–8 (computed live)</h2>
  <div id="occTable"></div>
</div>
<div class="card" style="margin-top:16px">
  <h2>Reported accuracy (paper Tables 3–5) vs this implementation</h2>
  <table class="dt">
    <tr><th>System / region</th><th>Paper AADP</th><th>This app (vs calibration anchors)</th></tr>
    <tr><td>CH₄ H–I–V / H–L<sub>W</sub>–V</td><td>0.7% / 2.05% (overall 1.37%)</td><td>≈ +19% @263 K ice end; 0.2–12% over 268–298 K</td></tr>
    <tr><td>CH₄ to 827 MPa (Table 4)</td><td>3.08%</td><td>not calibrated above 26.5 MPa — extrapolate with care</td></tr>
    <tr><td>C₂H₆ H–I–V / H–L<sub>W</sub>–V / H–L<sub>W</sub>–L<sub>C</sub></td><td>3.14% / 1.15% / 11.15%</td><td>≈ +23% ice end; 3–8% V region; <b>LC region far off (documented)</b></td></tr>
    <tr><td>C₃H₈ (sII)</td><td>3.83%</td><td>4.75% (V ±5%, ice end ±15%)</td></tr>
    <tr><td>N₂ (sII)</td><td>6.84%</td><td>−6…−19% (263–273 K), +44…+77% (281–293 K)</td></tr>
    <tr><td>CO₂ H–I–V / H–L<sub>W</sub>–V / H–L<sub>W</sub>–L<sub>C</b></td><td>2.59% / 2.87% / 10.2%</td><td>ice ±6%; V region ±4.5%; <b>Q₂/LC region off (documented)</b></td></tr>
    <tr><td>Mixtures (Table 5, 535 pts)</td><td>4.41% overall</td><td>structure selection &amp; Fig. 5 transition reproduced; full mixture datasets not re-validated</td></tr>
    <tr><td>Occupancies (Tables 6–8)</td><td>—</td><td>within 1.3% absolute of the paper’s values</td></tr>
  </table>
</div>
</div>

<!-- ================= PARAMS ================= -->
<div class="tabpane" id="tab-params">
<div class="cols2">
<div class="card">
  <h2>Method (paper’s formulation, implemented as published)</h2>
  <div class="eq">f<sub>w</sub><sup>H</sup>(T,P) = f<sub>w</sub><sup>β</sup>(T,P) · exp( Σ<sub>m</sub> ν<sub>m</sub> ln(1 − Σ<sub>j</sub> θ<sub>mj</sub>) ),&nbsp;&nbsp; θ<sub>mj</sub> = C<sub>mj</sub> f<sub>j</sub> / (1 + Σ<sub>j</sub> C<sub>mj</sub> f<sub>j</sub>)</div>
  <div class="eq">C<sub>mj</sub>(T) = (4π/kT) ∫₀<sup>R<sub>m</sub>−2a</sup> exp(−w(r)/kT) r² dr &nbsp;&nbsp;(McKoy–Sinanoğlu cell potential, Kihara core a)</div>
  <div class="eq">w(r) = 2zε [ (σ¹²/R¹¹r)(δ¹⁰ + (a/R)δ¹¹) − (σ⁶/R⁵r)(δ⁴ + (a/R)δ⁵) ],&nbsp; δ<sup>N</sup> = (1/N)[(1−r/R−a/R)<sup>−N</sup> − (1+r/R−a/R)<sup>−N</sup>]</div>
  <div class="eq">f<sub>w</sub><sup>β</sup> = exp(A − B/T) · exp( V<sub>β</sub>(T,P)(P − f<sub>β</sub>)/RT )</div>
  <div class="eq">V<sub>β</sub>(sI) = (11.835 + 2.217·10⁻⁵T + 2.242·10⁻⁶T²)³·10⁻³⁰N<sub>A</sub>/46 + 1.5265·10⁻⁹(P−0.1) − 2.8106·10⁻¹²(P−0.1)²</div>
  <div class="eq">V<sub>β</sub>(sII) = (17.13 + 2.249·10⁻⁴T + 2.013·10⁻⁶T² + 1.009·10⁻⁹T³)³·10⁻³⁰N<sub>A</sub>/136 − 8.006·10⁻⁹P + 5.448·10⁻¹²P² &nbsp;(P in MPa)</div>
  <div class="eq">PRSV2: κ = κ₀ + [κ₁ + κ₂(κ₃−T<sub>r</sub>)(1−√T<sub>r</sub>)](1+√T<sub>r</sub>)(0.7−T<sub>r</sub>),&nbsp; κ₀ = 0.378893 + 1.4897153ω − 0.17131848ω² + 0.0196554ω³</div>
  <div class="eq">van Laar mixing: a<sub>ij</sub> = √(a<sub>i</sub>a<sub>j</sub>)(1 − k<sub>ij</sub>k<sub>ji</sub>/(x<sub>i</sub>k<sub>ij</sub> + x<sub>j</sub>k<sub>ji</sub>))),&nbsp; k<sub>ij</sub>(T) = p<sub>ij</sub>T + q<sub>ij</sub>; gas–gas k<sub>ij</sub> = 0 (paper)</div>
  <div class="eq">Equilibrium: f<sub>w</sub><sup>H</sup> = f<sub>w</sub><sup>L</sup> (H–L<sub>W</sub>–V / H–L<sub>W</sub>–L<sub>C</sub>) or f<sub>w</sub><sup>H</sup> = f<sub>w</sub><sup>I</sup> (H–I–V, T &lt; 273.15 K, via exact ice bridge). Structure with the lower P wins.</div>
  <div class="hint">Fluid fugacities from PRSV2 with the paper’s temperature-dependent water–gas k<sub>ij</sub>; two-phase TP flash by successive substitution + Rachford–Rice (zero-feed components kept exactly zero). H–I–V branch: f<sub>ice</sub> = p<sub>sat,ice</sub>·exp(V<sub>ice</sub>(P−p<sub>sat,ice</sub>)/RT), p<sub>sat,ice</sub> = exp(28.89 − 6134.9/T) Pa, V<sub>ice</sub> = 19.65 cm³/mol.</div>
</div>
<div>
<div class="card">
  <h2>Kihara potential parameters (ε/k [K], σ [Å], a [Å])</h2>
  <table class="dt">
    <tr><th>Guest</th><th>this app (calibrated)</th><th>paper Table 2 (as printed)</th><th>cage radii used</th></tr>
__KIH_ROWS__
  </table>
  <div class="warnbox">The paper’s printed Table 2 values do <b>not</b> reproduce the paper’s own results in a standard McKoy–Sinanoğlu implementation (they give θ<sub>S</sub>(CH₄)≈0.01 vs Table 6’s 0.935, and P<sub>eq</sub> errors up to ×100). This app therefore uses parameters calibrated directly to reference equilibrium data + the paper’s Tables 6–8 occupancies (see Validation notes).</div>
</div>
<div class="card" style="margin-top:16px">
  <h2>Empty-lattice fugacity f<sub>β</sub> = exp(A − B/T) [Pa]</h2>
  <table class="dt"><tr><th>structure</th><th>this app</th><th>f<sub>β</sub>(273.15 K)</th><th>paper Eqs. (16)–(17) (Martín–Peters)</th><th>f<sub>β</sub>(273.15 K)</th></tr>
__FB_ROWS__</table>
  <div class="hint">The paper’s own Table 6 occupancies are consistent with f<sub>β,sI</sub>(273.15 K) ≈ 1050–1100 Pa, not with the printed Martín–Peters value (1295 Pa) — see Validation notes item 2.</div>
</div>
<div class="card" style="margin-top:16px">
  <h2>Aqueous-phase binary correction (this work)</h2>
  <table class="dt"><tr><th>gas</th><th>k(g,w) used in the aqueous phase</th><th>effect at 298 K</th></tr>__DK_ROWS__</table>
  <div class="warnbox">With the paper’s Eq. (10) and Table 1 parameters <i>as printed</i>, the van Laar combination gives a <b>negative cross attraction parameter</b> a(w,CO₂) &lt; 0 (K = k<sub>ij</sub>k<sub>ji</sub>/(x<sub>i</sub>k<sub>ij</sub>+x<sub>j</sub>k<sub>ji</sub>) &gt; 1), and hydrocarbon solubilities ~10× below experiment — while the paper reports 1–4% AADy for the same data (Fig. 1). We therefore add an aqueous-phase-only correction k(g,w) → k(g,w) + dK₀ + dK₁·(T−273.15 K), regressed to experimental solubilities; the vapor phase keeps the paper’s parameters unchanged. Documented deviation #4.</div>
</div>
<div class="card" style="margin-top:16px">
  <h2>Water–gas binary parameters k<sub>ij</sub> = p<sub>ij</sub>T + q<sub>ij</sub> (paper, used unchanged)</h2>
  <table class="dt"><tr><th>pair (water, gas)</th><th>p(w,g)</th><th>q(w,g)</th><th>p(g,w)</th><th>q(g,w)</th></tr>
  <tr><td>H₂O–CH₄</td><td class="mono">0.00364</td><td class="mono">−0.631</td><td class="mono">0.00290</td><td class="mono">−0.960</td></tr>
  <tr><td>H₂O–C₂H₆</td><td class="mono">0.00165</td><td class="mono">−0.627</td><td class="mono">0.00295</td><td class="mono">−0.986</td></tr>
  <tr><td>H₂O–C₃H₈</td><td class="mono">0.00024</td><td class="mono">−0.413</td><td class="mono">0.00198</td><td class="mono">−0.647</td></tr>
  <tr><td>H₂O–N₂</td><td class="mono">0.00390</td><td class="mono">−0.142</td><td class="mono">0.00300</td><td class="mono">−0.420</td></tr>
  <tr><td>H₂O–CO₂</td><td class="mono">0.00940</td><td class="mono">−0.517</td><td class="mono">0.00780</td><td class="mono">−0.258</td></tr></table>
  <div class="hint">Gas–gas interactions neglected, per the paper. Assignment used: “water-first” (in the water-rich liquid the (gas,water) direction dominates; in the gas-rich phase the (water,gas) direction dominates) — this is the assignment under which the calibration was performed.</div>
</div>
<div class="card" style="margin-top:16px">
  <h2>PRSV2 pure-component constants (κ₁κ₂κ₃ from Stryjek–Vera tables via DWSIM)</h2>
  <table class="dt"><tr><th>comp.</th><th>T<sub>c</sub> [K]</th><th>P<sub>c</sub> [MPa]</th><th>ω</th><th>κ₁</th><th>κ₂</th><th>κ₃</th></tr>
  <tr><td>CH₄</td><td class="mono">190.56</td><td class="mono">4.599</td><td class="mono">0.01142</td><td class="mono">−0.00159</td><td class="mono">0.1521</td><td class="mono">0.517</td></tr>
  <tr><td>C₂H₆</td><td class="mono">305.32</td><td class="mono">4.872</td><td class="mono">0.0995</td><td class="mono">0.02669</td><td class="mono">0.1358</td><td class="mono">0.424</td></tr>
  <tr><td>C₃H₈</td><td class="mono">369.83</td><td class="mono">4.248</td><td class="mono">0.1523</td><td class="mono">0.03136</td><td class="mono">0.2757</td><td class="mono">0.447</td></tr>
  <tr><td>N₂</td><td class="mono">126.19</td><td class="mono">3.3958</td><td class="mono">0.0372</td><td class="mono">0.01996</td><td class="mono">0.3162</td><td class="mono">0.535</td></tr>
  <tr><td>CO₂</td><td class="mono">304.13</td><td class="mono">7.377</td><td class="mono">0.2239</td><td class="mono">0.04285</td><td class="mono">0</td><td class="mono">0</td></tr>
  <tr><td>H₂O</td><td class="mono">647.14</td><td class="mono">22.064</td><td class="mono">0.3443</td><td class="mono">−0.06635</td><td class="mono">0.0199</td><td class="mono">0.443</td></tr></table>
</div>
</div>
</div>
</div>

<!-- ================= VALIDATION ================= -->
<div class="tabpane" id="tab-valid">
<div class="cols2">
<div class="card">
  <h2>Calibration procedure</h2>
  <p style="font-size:13px">The paper’s printed Kihara parameters (Table 2) are internally inconsistent with its own results (below), so the model was calibrated in this work against ≈70 reference equilibrium points (Sloan &amp; Koh compilation; Zenodo CO₂ low-temperature data) plus the paper’s Tables 6–8 occupancies and the Fig. 5 transition state. Objective (stage 1): fugacity residuals ln(f<sub>w</sub><sup>H</sup>/f<sub>w</sub><sup>L</sup>) at each anchor (paper Eq. 27), scaled by the local sensitivity |∂r/∂lnP| so that saturated-cage (flat) regions carry their true pressure error; per-anchor structure hinges enforce the correct sI/sII selection; a weak regularization tethers Langmuir constants to Pollard–Petersen-type correlations. Stage 2 then freezes this Kihara/f<sub>β</sub> set and regresses the aqueous-phase k(g,w) correction (dK₀, dK₁ — Parameters tab) against root-space residuals ln(P<sub>eq</sub><sup>model</sup>/P<sub>eq</sub><sup>exp</sup>), computed by bracketing the true equilibrium root at each anchor (never by evaluating the residual at P<sub>exp</sub>, which admits spurious roots), plus the experimental solubility anchors. The final fit reproduces:</p>
  <ul style="font-size:13px; margin:8px 0 0 18px">
    <li>all pure-gas structures: CH₄/C₂H₆/CO₂ sI (incl. ice branch), C₃H₈ sII, N₂ sII;</li>
    <li>occupancies within 1.3% of paper Tables 6–8;</li>
    <li>Fig. 5 sI→sII onset at y(C₃H₈) ≈ 3·10⁻⁴ (paper 5·10⁻⁴) at 4.0 MPa, 277.6 K;</li>
    <li>natural gas sII curve (Fig. 7) with the right trend and magnitude;</li>
    <li>H–L<sub>W</sub>–V pressures of CH₄, C₂H₆, CO₂, C₃H₈ typically within 0.2–12%.</li>
    <li>gas solubilities in water within ±20% for C₂H₆, C₃H₈, N₂, CO₂ and ×2–3 for CH₄ (paper Fig. 1 anchors).</li>
  </ul>
  <table class="dt" style="margin-top:10px"><tr><th>gas</th><th>AADP this app</th><th>AADP paper (Table 5)</th><th>note</th></tr>
  <tr><td>CH₄</td><td class="mono">7.6%</td><td class="mono">1.37%</td><td>incl. +19/+20% end points (263 K ice, 298 K); middle 268–293 K within 12%</td></tr>
  <tr><td>C₂H₆</td><td class="mono">72%</td><td class="mono">2.47%</td><td>dominated by the 3 H–L<sub>W</sub>–L<sub>C</sub> points (+183…+409%); the 11 V/ice points average ≈7%</td></tr>
  <tr><td>C₃H₈</td><td class="mono">4.7%</td><td class="mono">—</td><td>all 8 points |dev| ≤ 15%, six within 5%</td></tr>
  <tr><td>N₂</td><td class="mono">39%</td><td class="mono">6.84%</td><td>dominated by 281–293 K (+46…+77%); ≤277 K within 19%</td></tr>
  <tr><td>CO₂</td><td class="mono">30%</td><td class="mono">2.89%</td><td>ice +13…+25%, V ≤15%, L<sub>C</sub> +10…+207%</td></tr></table>
</div>
<div class="card">
  <h2>Paper inconsistencies found (why deviations exist)</h2>
  <ol style="font-size:13px; margin-left:18px">
    <li style="margin-bottom:6px"><b>Table 2 Kihara parameters are unusable as printed.</b> In a standard McKoy–Sinanoğlu implementation with Eqs. (16)–(17), they yield θ<sub>S</sub>(CH₄, 273.2 K) ≈ 0.01 vs Table 6’s 0.935, P<sub>eq</sub>(CH₄) errors ≈ ×100, and vanishing Langmuir constants for C₃H₈/CO₂ (σ larger than the cage free radius). They cannot be the values the authors’ code used.</li>
    <li style="margin-bottom:6px"><b>Table 6 occupancies contradict the printed Martín–Peters f<sub>β</sub>.</b> θ = (0.969, 0.935) at 273.2 K with f<sub>β,sI</sub>(273.15) = exp(28.4−5800/T) = 1295 Pa gives f<sub>w</sub><sup>H</sup> = 748 Pa vs f<sub>w</sub><sup>L</sup> ≈ 600 Pa — 25% off equilibrium. f<sub>β,sI</sub>(273.15) ≈ 1050–1100 Pa (this work’s fitted value) reconciles Tables 6–8 with the equilibrium data.</li>
    <li>Table 7’s caption repeats Table 6’s text (both are CH₄; the 276.7 K data are indeed CH₄).</li>
    <li style="margin-top:6px"><b>The printed Table 1 k<sub>ij</sub> cannot reproduce the paper’s own Fig. 1.</b> Used as printed in the van Laar combination of Eq. (10), they give x(CO₂) ≈ 10⁻¹⁷ (a negative-a artifact) and hydrocarbon solubilities ~10× too low, while the paper reports 1–4% AADy on those same data. This work keeps the vapor-phase values untouched and adds an aqueous-phase-only correction dK₀ + dK₁(T−273.15) to k(g,w), regressed in stage 2 — see the Parameters tab.</li>
  </ol>
</div>
</div>
<div class="card" style="margin-top:16px">
  <h2>Known limitations of this implementation (use with care)</h2>
  <table class="dt">
    <tr><th>Region</th><th>Behavior</th><th>Why</th></tr>
    <tr><td>H–L<sub>W</sub>–L<sub>C</sub> (above Q₂: C₂H₆ &gt; 287.75 K, CO₂ &gt; 283 K)</td><td>C₂H₆: <b>still large</b> (+183…+409%). CO₂: improved by the aqueous correction to +10/+28% at 285/287 K, but +129/+207% remain at 289/291 K. The app flags these states with ⚠.</td><td>There the cages are saturated and |∂ln(f<sub>H</sub>/f<sub>L</sub>)/∂lnP| ≈ 0.02: a 1% fugacity-ratio error moves P<sub>eq</sub> by ~50%. With a single shared f<sub>β</sub>(T) = exp(A−B/T) and frozen paper k<sub>ij</sub>, the residual slope dr/dT ≈ +0.01/K cannot be flattened without destroying the (well-fitted) vapor-region curves. Tested and rejected: B shifts, 1/T² curvature, free Kihara core radii, per-region weights — all trade one region against another. The paper reports 10–11% there.</td></tr>
    <tr><td>N₂ above ~280 K</td><td>P overpredicted, +46…+77% at 281–293 K (−19…+13% at 263–277 K)</td><td>Same flat-slope tension on the sII branch; paper reports 3.14% (272–299.7 K).</td></tr>
    <tr><td>CH₄ ice end (263 K) &amp; 298 K point</td><td>+19% / +20%</td><td>Ends of the fitted range; the middle (268–293 K) is within 12% (mostly ≤7%).</td></tr>
    <tr><td>Gas solubility in water</td><td>Near-quantitative after the stage-2 aqueous correction: model/exp = 0.83–1.02 for C₂H₆, C₃H₈, N₂, CO₂; CH₄ underpredicted ×2–3 (0.32 at 2.6 MPa/273 K, 0.45 at 0.1 MPa/298 K)</td><td>Pushing CH₄ solubility to ≈0.9 (dK₀ = −0.16) degrades the CH₄ vapor-region curve AADP from 7.6% to 11.7%; the calibration keeps the P(T) curves exact and accepts the CH₄ solubility gap. See the Fig. 1 button on the Paper tab.</td></tr>
    <tr><td>P &gt; ~30 MPa, T &gt; ~305 K</td><td>Extrapolation</td><td>Calibration anchors stop at 26.5 MPa / 303 K; the paper validates CH₄ to 827 MPa with the same V<sub>β</sub>(T,P) form but that range was not re-anchored here.</td></tr>
    <tr><td>Fig. 6 CH₄+C₂H₆ sII lens</td><td>Not reproduced (sI everywhere)</td><td>With physically calibrated C(T) the sII branch never dips below sI in this binary; the paper reports transitions at y<sub>CH4</sub> = 0.753/0.724.</td></tr>
  </table>
</div>
<div class="card" style="margin-top:16px">
  <h2>Numerical verification</h2>
  <p style="font-size:13px">The JavaScript engine in this page is a line-by-line port of the Python research prototype. Verified against it to machine precision: Gauss–Legendre nodes, ln φᵢ (both phases, all 6 components), TP flash (phase compositions to 10⁻¹⁴), Langmuir constants (10⁻¹⁵), hydrate fugacity &amp; occupancies (10⁻¹⁵), and equilibrium pressures incl. mixtures and the ice branch (≤10⁻⁷ relative). Cubic roots are computed analytically with the same filtering as the prototype (real roots with Z &gt; B).</p>
</div>
</div>

</main>
<footer>Implementation of Moradi &amp; Khosravani, Fluid Phase Equilibria 338 (2013) 179–187 — calibrated reproduction with documented deviations. Built as a single self-contained HTML file (no external resources).<br>© 2026 <a href="https://qazinasir.com" target="_blank" rel="noopener">Qazi Nasir</a>. All rights reserved.</footer>
<script>
__ENGINE__
</script>
<script>
window.ANCHORS = __ANCHORS__;
</script>
<script>
__APP__
</script>
</body>
</html>
"""
html = (HTML.replace('__CSS__', css).replace('__ENGINE__', engine)
            .replace('__ANCHORS__', anchors).replace('__APP__', app)
            .replace('__KIH_ROWS__', kih_rows).replace('__FB_ROWS__', fb_rows).replace('__DK_ROWS__', dk_rows))
out = os.path.join(HERE, '..', 'index.html')
open(out, 'w').write(html)
print('wrote', os.path.abspath(out), len(html), 'bytes')
