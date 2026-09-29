/* HydratePRSV2 engine — faithful JS port of engine.py + calib4.wgf + calibrate.Peq
   Moradi & Khosravani (FPE 338 (2013) 179-187) methodology, calibrated parameters (calib10).
   Components order: [CH4, C2H6, C3H8, N2, CO2, H2O]  */
'use strict';
const kB = 1.380649e-23, NA = 6.02214076e23, Rg = 8.314462618;
const NAMES = ['CH4', 'C2H6', 'C3H8', 'N2', 'CO2', 'H2O'];
const IW = 5;
const TC = [190.56, 305.32, 369.83, 126.19, 304.13, 647.14];
const PC = [4.599e6, 4.872e6, 4.248e6, 3.3958e6, 7.377e6, 22.064e6];
const OMG = [0.01142, 0.0995, 0.1523, 0.0372, 0.2239, 0.3443];
const K1 = [-0.00159, 0.02669, 0.03136, 0.01996, 0.04285, -0.06635];
const K2 = [0.1521, 0.1358, 0.2757, 0.3162, 0.0, 0.0199];
const K3 = [0.517, 0.424, 0.447, 0.535, 0.0, 0.443];
// paper water-gas van Laar parameters: pij qij pji qji, kij = p*T+q
const KIJ_P = {
  CH4: [0.00364, -0.631, 0.0029, -0.96],
  C2H6: [0.00165, -0.627, 0.00295, -0.986],
  C3H8: [0.00024, -0.413, 0.00198, -0.647],
  N2: [0.0039, -0.142, 0.003, -0.42],
  CO2: [0.0094, -0.517, 0.0078, -0.258]
};
// FINAL calibrated Kihara (calib10): eps/k [K], sigma [m], a [m]
const KIH = {
  CH4: [161.05387837375346, 3.1595316757233196e-10, 3.317e-11],
  C2H6: [180.45998837815412, 3.1638101945332206e-10, 4.799e-11],
  C3H8: [210.4072040439076, 3.45538154525678e-10, 6.502e-11],
  N2: [151.37955442363904, 2.959609283049243e-10, 3.41e-11],
  CO2: [172.3698592701493, 2.997294068555943e-10, 6.786e-11],
};
// FINAL calibrated empty-lattice fugacity: ln f_beta = A - B/T  [Pa]
const FBETA = {
  sI: [28.79326627652638, 5953.707727312784],
  sII: [28.605594106116573, 5922.567549323371],
};
const CAGES = { sI: [[3.917e-10, 20], [4.33e-10, 24]], sII: [[3.917e-10, 20], [4.73e-10, 28]] };
const NU = { sI: [2 / 46, 6 / 46], sII: [16 / 136, 8 / 136] };
const NW = { sI: 46, sII: 136 };
const V_ICE = 1.965e-5;
const GASES = ['CH4', 'C2H6', 'C3H8', 'N2', 'CO2'];

/* ---------- Gauss-Legendre 200 nodes (matches numpy leggauss to ~1e-14) ---------- */
function leggauss(n) {
  const x = new Array(n), w = new Array(n);
  const m = (n + 1) >> 1;
  for (let i = 0; i < m; i++) {
    let z = Math.cos(Math.PI * (i + 0.75) / (n + 0.5));   // initial guess
    let pp = 0, z1 = 0;
    for (let it = 0; it < 100; it++) {
      let p0 = 1, p1 = z;
      for (let j = 2; j <= n; j++) { const p2 = ((2 * j - 1) * z * p1 - (j - 1) * p0) / j; p0 = p1; p1 = p2; }
      pp = n * (z * p1 - p0) / (z * z - 1);
      z1 = z; z = z1 - p1 / pp;
      if (Math.abs(z - z1) < 3e-16) break;
    }
    x[i] = -z; x[n - 1 - i] = z;
    w[i] = 2 / ((1 - z * z) * pp * pp); w[n - 1 - i] = w[i];
  }
  return [x, w];
}
const [GLX, GLW] = leggauss(200);

/* ---------- PRSV2 EOS ---------- */
function kappa(i, T) {
  const Tr = T / TC[i], s = Math.sqrt(Tr);
  const k0 = 0.378893 + 1.4897153 * OMG[i] - 0.17131848 * OMG[i] ** 2 + 0.0196554 * OMG[i] ** 3;
  return k0 + (K1[i] + K2[i] * (K3[i] - Tr) * (1 - s)) * (1 + s) * (0.7 - Tr);
}
function pure_ab(i, T) {
  const kap = kappa(i, T);
  const alpha = (1 + kap * (1 - Math.sqrt(T / TC[i]))) ** 2;
  return [0.457235 * Rg * Rg * TC[i] * TC[i] / PC[i] * alpha, 0.077796 * Rg * TC[i] / PC[i]];
}
const DK_AQ = { CH4: [-0.1, 0.0], C2H6: [-0.11277309551359646, 0.002865740232777846], C3H8: [-0.027736336740983317, -0.0010563388736334642], N2: [-0.4968215535879624, 0.007745720584847979], CO2: [-1.9906188961566265, -0.00808850070947816] };  // aqueous shift of k(gas,water): dK0+dK1*(T-273.15)
function aqShift(g, T) { const s = DK_AQ[g]; return s[0] + s[1] * (T - 273.15); }
function kij_pair(i, j, T, assign, x) {
  if (i === j) return [0, 0];
  const gi = i !== IW ? NAMES[i] : null;
  const gj = j !== IW ? NAMES[j] : null;
  const aq = (x !== undefined && x !== null) && x[IW] > 0.5;
  if (gi !== null && Object.prototype.hasOwnProperty.call(KIJ_P, gi) && gj === null) {
    const [p1, q1, p2, q2] = KIJ_P[gi];
    const k2 = p2 * T + q2 + (aq ? aqShift(gi, T) : 0);
    return assign === 'gas_first' ? [p1 * T + q1, k2] : [k2, p1 * T + q1];
  }
  if (gj !== null && Object.prototype.hasOwnProperty.call(KIJ_P, gj) && gi === null) {
    const [p1, q1, p2, q2] = KIJ_P[gj];
    const k2 = p2 * T + q2 + (aq ? aqShift(gj, T) : 0);
    return assign === 'gas_first' ? [k2, p1 * T + q1] : [p1 * T + q1, k2];
  }
  return [0, 0];
}
function mix_params(x, T, assign) {
  const n = x.length;
  const ai = new Array(n), bi = new Array(n);
  for (let i = 0; i < n; i++) [ai[i], bi[i]] = pure_ab(i, T);
  const K = Array.from({ length: n }, () => new Array(n).fill(0));
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) {
    if (i === j) continue;
    const [kij, kji] = kij_pair(i, j, T, assign, x);
    const D = x[i] * kij + x[j] * kji;
    K[i][j] = Math.abs(D) < 1e-14 ? 0 : kij * kji / D;
  }
  const aij = Array.from({ length: n }, () => new Array(n).fill(0));
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++)
    aij[i][j] = i === j ? ai[i] : Math.sqrt(ai[i] * ai[j]) * (1 - K[i][j]);
  let am = 0, bm = 0;
  for (let i = 0; i < n; i++) { bm += x[i] * bi[i]; for (let j = 0; j < n; j++) am += x[i] * aij[i][j] * x[j]; }
  return { am, bm, ai, bi, aij };
}
function cubicRealRoots(c2, c1, c0) {  // Z^3 + c2 Z^2 + c1 Z + c0
  const p = c1 - c2 * c2 / 3;
  const q = 2 * c2 ** 3 / 27 - c2 * c1 / 3 + c0;
  const disc = (q / 2) ** 2 + (p / 3) ** 3;
  const off = -c2 / 3;
  if (disc > 0) {
    const sq = Math.sqrt(disc);
    const u = Math.cbrt(-q / 2 + sq), v = Math.cbrt(-q / 2 - sq);
    return [u + v + off];
  }
  if (Math.abs(p) < 1e-300) return [off];
  const r = Math.sqrt(-(p ** 3) / 27);
  const th = Math.acos(Math.max(-1, Math.min(1, -q / (2 * r))));
  const m = 2 * Math.sqrt(-p / 3);
  return [m * Math.cos(th / 3) + off, m * Math.cos((th + 2 * Math.PI) / 3) + off, m * Math.cos((th + 4 * Math.PI) / 3) + off];
}
function cubic_Z(A, B, root) {
  const c2 = -(1 - B), c1 = A - 3 * B * B - 2 * B, c0 = -(A * B - B * B - B ** 3);
  const rts = cubicRealRoots(c2, c1, c0);
  const real = rts.filter(z => z > B * (1 + 1e-9)).sort((a, b) => a - b);
  if (!real.length) return null;
  return root === 'v' ? real[real.length - 1] : real[0];
}
function lnphi(x, T, P, root, assign) {
  const { am, bm, ai, bi, aij } = mix_params(x, T, assign);
  const A = am * P / (Rg * T) ** 2, B = bm * P / (Rg * T);
  const Z = cubic_Z(A, B, root);
  if (Z === null) return null;
  const s2 = Math.SQRT2;
  if (Z - B <= 0 || Z + (1 - s2) * B <= 0) return null;
  const L = Math.log((Z + (1 + s2) * B) / (Z + (1 - s2) * B));
  const f = Z - 1 - Math.log(Z - B) - A / (2 * s2 * B) * L;
  const fA = -L / (2 * s2 * B);
  const fB = 1 / (Z - B) + A * L / (2 * s2 * B * B) - A / (2 * s2 * B) * ((1 + s2) / (Z + (1 + s2) * B) - (1 - s2) / (Z + (1 - s2) * B));
  const n = x.length;
  const dam = new Array(n);
  for (let i = 0; i < n; i++) {
    let s = 2 * x[i] * ai[i];
    for (let j = 0; j < n; j++) {
      if (j === i) continue;
      s += 2 * x[j] * aij[i][j];
      const [kij, kji] = kij_pair(i, j, T, assign, x);
      const D = x[i] * kij + x[j] * kji;
      if (Math.abs(D) > 1e-14) s += 2 * x[i] * x[j] * Math.sqrt(ai[i] * ai[j]) * kij * kji * kij / (D * D);
    }
    dam[i] = s;
  }
  let sdam = 0, sbm = 0;
  for (let i = 0; i < n; i++) { sdam += x[i] * dam[i]; sbm += x[i] * bi[i]; }
  const Abar = sdam * P / (Rg * T) ** 2, Bbar = sbm * P / (Rg * T);
  const out = new Array(n);
  for (let i = 0; i < n; i++) {
    const Ai_ = dam[i] * P / (Rg * T) ** 2, Bi_ = bi[i] * P / (Rg * T);
    out[i] = f + fA * (Ai_ - Abar) + fB * (Bi_ - Bbar);
  }
  return out;
}

/* ---------- Langmuir / vdW-P ---------- */
function langmuir(T, Rc, z, eps_k, sig, a) {
  const rmax = (Rc - 2 * a) * 0.9999999;
  if (rmax <= 1e-12) return 0;
  const jac = 0.5 * rmax;
  const eps = eps_k * kB;
  let sum = 0;
  for (let k = 0; k < GLX.length; k++) {
    const r = 0.5 * rmax * (GLX[k] + 1);
    const u = 1 - r / Rc - a / Rc, v = 1 + r / Rc - a / Rc;
    const d10 = (u ** -10 - v ** -10) / 10, d11 = (u ** -11 - v ** -11) / 11;
    const d4 = (u ** -4 - v ** -4) / 4, d5 = (u ** -5 - v ** -5) / 5;
    const w = 2 * z * eps * (sig ** 12 / (Rc ** 11 * r) * (d10 + (a / Rc) * d11) - sig ** 6 / (Rc ** 5 * r) * (d4 + (a / Rc) * d5));
    const arg = Math.max(-700, Math.min(700, w / (kB * T)));
    sum += GLW[k] * Math.exp(-arg) * r * r;
  }
  return 4 * Math.PI / (kB * T) * sum * jac;
}
function langmuir_set(T, gasList, struct) {
  const out = {};
  for (const g of gasList) {
    const [e, s, a] = KIH[g];
    out[g] = CAGES[struct].map(([Rc, z]) => langmuir(T, Rc, z, e, s, a));
  }
  return out;
}
function psat_beta(T, struct) { const [A, B] = FBETA[struct]; return Math.exp(A - B / T); }
function Vbeta(T, P_MPa, struct) {
  if (struct === 'sI') {
    const base = (11.835 + 2.217e-5 * T + 2.242e-6 * T * T) ** 3 * 1e-30 * NA / 46;
    return base + 1.5265e-9 * (P_MPa - 0.1) - 2.8106e-12 * (P_MPa - 0.1) ** 2;
  }
  const base = (17.13 + 2.249e-4 * T + 2.013e-6 * T * T + 1.009e-9 * T ** 3) ** 3 * 1e-30 * NA / 136;
  return base - 8.006e-9 * P_MPa + 5.448e-12 * P_MPa * P_MPa;
}
function psat_ice(T) { return Math.exp(28.89 - 6134.9 / T); }
function fH_w(T, P, fgas, struct) {  // fgas: {gas: fugacity Pa}
  const gases = Object.keys(fgas);
  const Cset = langmuir_set(T, gases, struct);
  const nu = NU[struct];
  let S = 0; const thetas = {};
  for (let m = 0; m < 2; m++) {
    let tot = 0;
    for (const g of gases) tot += Cset[g][m] * fgas[g];
    for (const g of gases) thetas[g + '|' + m] = Cset[g][m] * fgas[g] / (1 + tot);
    S += nu[m] * Math.log(1 - Math.min(tot / (1 + tot), 1 - 1e-300));
  }
  const psb = psat_beta(T, struct);
  const Vb = Vbeta(T, P / 1e6, struct);
  const fbeta = psb * Math.exp(Vb * (P - psb) / (Rg * T));
  return { fH: fbeta * Math.exp(S), thetas, S, fbeta, Cset };
}

/* ---------- flash ---------- */
function wilson_K(T, P) {
  return PC.map((pc, i) => pc / P * Math.exp(5.373 * (1 + OMG[i]) * (1 - TC[i] / T)));
}
function flash_TP(zin, T, P, assign, maxiter = 500, tol = 1e-10) {
  let s = 0; const z = zin.map(v => Math.max(v, 0));
  for (const v of z) s += v;
  for (let i = 0; i < 6; i++) z[i] /= s;
  const act = z.map(v => v > 0);
  let K = wilson_K(T, P).map(v => Math.max(1e-8, Math.min(1e8, v)));
  let beta = 0.5;
  const split = (K, beta) => {
    const x = new Array(6); let sx = 0;
    for (let i = 0; i < 6; i++) { x[i] = act[i] ? z[i] / (1 + beta * (K[i] - 1)) : 0; sx += x[i]; }
    if (sx <= 0 || !isFinite(sx)) return null;
    for (let i = 0; i < 6; i++) x[i] /= sx;
    const y = new Array(6); let sy = 0;
    for (let i = 0; i < 6; i++) { y[i] = act[i] ? K[i] * x[i] : 0; sy += y[i]; }
    if (sy <= 0 || !isFinite(sy)) return null;
    for (let i = 0; i < 6; i++) y[i] /= sy;
    return [x, y];
  };
  const rr = b => { let v = 0; for (let i = 0; i < 6; i++) { const d = 1 + b * (K[i] - 1); v += z[i] * (K[i] - 1) / d; } return isFinite(v) ? v : (b < 0.5 ? 1e6 : -1e6); };
  for (let it = 0; it < maxiter; it++) {
    const r0 = rr(1e-12), r1 = rr(1 - 1e-12);
    if (r0 > 0 && r1 < 0) {
      let lo = 1e-12, hi = 1 - 1e-12;
      for (let k = 0; k < 100; k++) { const mid = 0.5 * (lo + hi); (rr(mid) > 0) ? lo = mid : hi = mid; if (hi - lo < 1e-15) break; }
      beta = 0.5 * (lo + hi);
    } else beta = r0 <= 0 ? 1e-10 : 1 - 1e-10;
    const sp = split(K, beta);
    if (!sp) return null;
    const [x, y] = sp;
    let lphi = lnphi(x, T, P, 'l', assign), vphi = lnphi(y, T, P, 'v', assign);
    if (lphi === null && vphi === null) return null;
    if (lphi === null) lphi = vphi.slice();
    if (vphi === null) vphi = lphi.slice();
    let bad = false;
    for (let i = 0; i < 6; i++) if (act[i] && (!isFinite(lphi[i]) || !isFinite(vphi[i]))) bad = true;
    if (bad) return null;
    const Knew = K.slice();
    for (let i = 0; i < 6; i++) if (act[i]) Knew[i] = Math.exp(Math.max(-40, Math.min(40, lphi[i] - vphi[i])));
    let err = 0;
    for (let i = 0; i < 6; i++) if (act[i]) err = Math.max(err, Math.abs(Knew[i] / K[i] - 1));
    if (err > 1.0) { for (let i = 0; i < 6; i++) if (act[i]) K[i] = K[i] * Math.sqrt(Knew[i] / K[i]); }
    else K = Knew;
    if (err < tol) {
      const sp2 = split(K, beta);
      if (!sp2) return null;
      return { x: sp2[0], y: sp2[1], beta, conv: true };
    }
  }
  const sp2 = split(K, beta);
  if (!sp2) return null;
  return { x: sp2[0], y: sp2[1], beta, conv: false };
}

/* ---------- fluid fugacities (calib4.wgf, assign=water_first per calibration) ---------- */
function wgf(T, P, gasFrac, assign = 'water_first') {
  const gl = Object.keys(gasFrac).sort();
  if (T < 273.15) {
    const fi0 = psat_ice(T);
    const fi = fi0 * Math.exp(V_ICE * (P - fi0) / (Rg * T));
    const y = new Array(6).fill(0);
    let yw = Math.min(fi / P, 0.5);
    for (let it = 0; it < 12; it++) {
      y.fill(0); y[IW] = yw;
      for (const g of gl) y[NAMES.indexOf(g)] = (1 - yw) * gasFrac[g];
      const vp = lnphi(y, T, P, 'v', assign);
      if (vp === null) break;
      const ywNew = Math.max(1e-14, Math.min(0.6, fi / (P * Math.exp(vp[IW]))));
      if (Math.abs(ywNew - yw) < 1e-13) { yw = ywNew; break; }
      yw = 0.5 * yw + 0.5 * ywNew;
    }
    y.fill(0); y[IW] = yw;
    for (const g of gl) y[NAMES.indexOf(g)] = (1 - yw) * gasFrac[g];
    const vp = lnphi(y, T, P, 'v', assign);
    if (vp === null) return null;
    const fg = {};
    for (const g of gl) { const i = NAMES.indexOf(g); fg[g] = y[i] * Math.exp(vp[i]) * P; }
    return { fLw: fi, fg, tag: 'ice', y, x: null };
  }
  const z = new Array(6).fill(0); z[IW] = 0.8;
  for (const g of gl) z[NAMES.indexOf(g)] = 0.2 * gasFrac[g];
  const fl = flash_TP(z, T, P, assign, 150, 1e-9);
  if (!fl) return null;
  let xaq = fl.x, yg = fl.y;
  if (xaq[IW] < yg[IW]) { xaq = fl.y; yg = fl.x; }
  const lp = lnphi(xaq, T, P, 'l', assign), vp = lnphi(yg, T, P, 'v', assign);
  if (lp === null || vp === null) return null;
  const fLw = xaq[IW] * Math.exp(lp[IW]) * P;
  const fg = {};
  for (const g of gl) { const i = NAMES.indexOf(g); fg[g] = yg[i] * Math.exp(vp[i]) * P; }
  // tag: 'lc' when the gas-rich phase is a liquid (low compressibility factor)
  const mp = mix_params(yg, T, assign);
  const Zg = cubic_Z(mp.am * P / (Rg * T) ** 2, mp.bm * P / (Rg * T), 'v');
  const tag = (Zg !== null && Zg < 0.5) ? 'lc' : 'lw';
  return { fLw, fg, tag, x: xaq, y: yg, beta: fl.beta, conv: fl.conv };
}
/* residual r(lnP) = ln(fH/fLw) for a structure */
function rOf(T, P, gasFrac, struct, assign = 'water_first') {
  const st = wgf(T, P, gasFrac, assign);
  if (!st) return NaN;
  const { fH } = fH_w(T, P, st.fg, struct);
  const v = Math.log(fH) - Math.log(st.fLw);
  return isFinite(v) ? v : NaN;
}
/* ---------- equilibrium pressure ---------- */
function peqStruct(T, gasFrac, struct, Pguess = null, assign = 'water_first') {
  const lo = Pguess ? Pguess / 6 : 2e2, hi = Pguess ? Pguess * 6 : 4e8;
  const n = Pguess ? 14 : 60;
  const lns = [], vals = [];
  for (let k = 0; k < n; k++) { const L = Math.log(lo) + (Math.log(hi) - Math.log(lo)) * k / (n - 1); lns.push(L); vals.push(rOf(T, Math.exp(L), gasFrac, struct, assign)); }
  for (let k = 0; k < n - 1; k++) {
    if (isFinite(vals[k]) && isFinite(vals[k + 1]) && vals[k] * vals[k + 1] < 0) {
      let a = lns[k], b = lns[k + 1], fa = vals[k];
      for (let i = 0; i < 80; i++) {
        const m = 0.5 * (a + b), fm = rOf(T, Math.exp(m), gasFrac, struct, assign);
        if (!isFinite(fm)) break;
        if (fa * fm < 0) { b = m; } else { a = m; fa = fm; }
        if (b - a < 1e-10) break;
      }
      return Math.exp(0.5 * (a + b));
    }
  }
  if (Pguess) return peqStruct(T, gasFrac, struct, null, assign);
  return null;
}
function peqSelect(T, gasFrac, Pguess = null, assign = 'water_first') {
  const P1 = peqStruct(T, gasFrac, 'sI', Pguess, assign);
  const P2 = peqStruct(T, gasFrac, 'sII', Pguess, assign);
  const cand = [];
  if (P1) cand.push([P1, 'sI']);
  if (P2) cand.push([P2, 'sII']);
  if (!cand.length) return { P: null, struct: null, sI: P1, sII: P2 };
  cand.sort((a, b) => a[0] - b[0]);
  return { P: cand[0][0], struct: cand[0][1], sI: P1, sII: P2 };
}
/* full state details at equilibrium */
function stateAt(T, gasFrac, struct, P, assign = 'water_first') {
  const st = wgf(T, P, gasFrac, assign);
  if (!st) return null;
  const H = fH_w(T, P, st.fg, struct);
  const out = {
    T, P, struct, tag: st.tag, fLw: st.fLw, fg: st.fg, fH: H.fH, S: H.S, fbeta: H.fbeta,
    theta: {}, C: H.Cset, yw: null, xaq: st.x, yg: st.y, K: {}, psat_w: null, y: st.y
  };
  for (const g of Object.keys(st.fg)) {
    out.theta[g] = [H.thetas[g + '|0'], H.thetas[g + '|1']];
    out.K[g] = H.Cset[g].map((C, m) => C * st.fg[g]);  // cage loading Cf
  }
  if (st.y) out.yw = st.y[IW];
  out.psat_w = psat_water(T);
  // hydration number: water molecules per guest molecule
  let occ = 0;
  for (let m = 0; m < 2; m++) { let tot = 0; for (const g of Object.keys(st.fg)) tot += H.Cset[g][m] * st.fg[g]; occ += NU[struct][m] * tot / (1 + tot); }
  out.hydrationNumber = occ > 0 ? 1 / occ : Infinity;
  return out;
}
/* water psat correlation (Wagner-Pruss IAPWS-95 auxiliary eq, 273.16-647 K) */
function psat_water(T) {
  if (T < 273.16) return psat_ice(T);
  const th = T / 647.096;
  const A = -7.85951783, B = 1.84408259, C = -11.7866497, D = 22.6807411, E = -15.9618719, F = 1.80122502;
  return 22.064e6 * Math.exp((647.096 / T) * (A * (1 - th) + B * (1 - th) ** 1.5 + C * (1 - th) ** 3 + D * (1 - th) ** 3.5 + E * (1 - th) ** 4 + F * (1 - th) ** 7.5));
}
/* dissociation curve over T grid with warm starts */
function peqCurve(Tlist, gasFrac, structMode = 'auto', assign = 'water_first') {
  const out = [];
  let warm = null;
  for (const T of Tlist) {
    let res;
    if (structMode === 'sI' || structMode === 'sII') {
      const P = peqStruct(T, gasFrac, structMode, warm, assign);
      res = { P, struct: P ? structMode : null, sI: structMode === 'sI' ? P : null, sII: structMode === 'sII' ? P : null };
    } else res = peqSelect(T, gasFrac, warm, assign);
    if (res.P) warm = res.P;
    out.push({ T, ...res });
  }
  return out;
}
const E = {
  kB, NA, Rg, NAMES, IW, TC, PC, OMG, K1, K2, K3, KIJ_P, KIH, FBETA, CAGES, NU, NW, V_ICE, GASES,
  GLX, GLW, kappa, pure_ab, kij_pair, mix_params, cubic_Z, lnphi, langmuir, langmuir_set,
  psat_beta, Vbeta, psat_ice, fH_w, wilson_K, flash_TP, wgf, rOf, peqStruct, peqSelect, stateAt, psat_water, peqCurve
};

// expose for both Node (tests) and browser (app.js uses global E)
if (typeof module !== 'undefined' && module.exports) module.exports = E;
if (typeof window !== 'undefined') window.E = E;

