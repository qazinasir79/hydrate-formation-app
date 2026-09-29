const E = require("./engine_port.js");
const ref = require("/tmp/jsref.json");
let fails = 0;
function chk(name, a, b, rtol) {
  const ok = isFinite(a) && Math.abs(a / b - 1) < rtol;
  if (!ok) { fails++; console.log(`FAIL ${name}: js=${a} py=${b}`); }
  else console.log(`ok   ${name} rel=${Math.abs(a / b - 1).toExponential(1)}`);
}
const C1 = E.langmuir(273.15, E.CAGES.sI[0][0], E.CAGES.sI[0][1], ...E.KIH.CH4);
const C2 = E.langmuir(273.15, E.CAGES.sI[1][0], E.CAGES.sI[1][1], ...E.KIH.CH4);
chk("C_CH4_sI_small", C1, ref.C_CH4_sI[0], 1e-10);
chk("C_CH4_sI_large", C2, ref.C_CH4_sI[1], 1e-10);
const w = E.wgf(273.15, 2.6e6, { CH4: 1.0 });
const H = E.fH_w(273.15, 2.6e6, w.fg, "sI");
chk("fH", H.fH, ref.wgf.fH, 1e-10);
chk("S", H.S, ref.wgf.S, 1e-10);
chk("theta_L", H.thetas["CH4|1"], ref.wgf.theta_L, 1e-10);
chk("theta_S", H.thetas["CH4|0"], ref.wgf.theta_S, 1e-10);
const t0 = Date.now();
chk("Peq CH4 sI 273.15", E.peqStruct(273.15, { CH4: 1 }, "sI", 2.6e6), ref.peq["CH4_sI_273.15"], 3e-6);
chk("Peq CO2 sI 273.15", E.peqStruct(273.15, { CO2: 1 }, "sI", 1.28e6), ref.peq["CO2_sI_273.15"], 3e-6);
chk("Peq C3H8 sII 277.59", E.peqStruct(277.59, { C3H8: 1 }, "sII", 0.42e6), ref.peq["C3H8_sII_277.59"], 3e-6);
chk("Peq N2 sII 273.15", E.peqStruct(273.15, { N2: 1 }, "sII", 16.4e6), ref.peq["N2_sII_273.15"], 3e-6);
chk("Peq CH4 sI 263.15 ice", E.peqStruct(263.15, { CH4: 1 }, "sI", 1.55e6), ref.peq["CH4_sI_263.15"], 3e-6);
chk("Peq CH4 sI 298.15", E.peqStruct(298.15, { CH4: 1 }, "sI", 26.5e6), ref.peq["CH4_sI_298.15"], 3e-6);
const ng = E.peqSelect(281.15, { CH4: 0.932, C2H6: 0.0425, C3H8: 0.0161, N2: 0.0043, CO2: 0.0051 }, 1.6e6);
chk("NG 281.15 P", ng.P, ref["ng_281.15"][0], 3e-6);
console.log("NG struct:", ng.struct, "(py:", ref["ng_281.15"][1], ")");
const f5 = E.peqSelect(277.6, { CH4: 0.9997, C3H8: 0.0003 }, 4e6);
chk("Fig5 y3e-4 P", f5.P, ref["fig5_y3e4"][0], 3e-6);
console.log("Fig5 struct:", f5.struct, "(py:", ref["fig5_y3e4"][1], ")");
if (ref.solub) {
  const SP = { CH4: [273.15, 2.6e6], CO2: [278.15, 3.5e6], N2: [298.15, 0.101325e6] };
  for (const g of ["CH4", "CO2", "N2"]) {
    const z = new Array(6).fill(0); z[E.IW] = 0.99; z[E.NAMES.indexOf(g)] = 0.01;
    const fl = E.flash_TP(z, SP[g][0], SP[g][1], "water_first");
    const xaq = fl.x[E.IW] >= fl.y[E.IW] ? fl.x : fl.y;
    chk("solub_" + g, xaq[E.NAMES.indexOf(g)], ref.solub[g], 1e-6);
  }
}
console.log("elapsed", ((Date.now() - t0) / 1000).toFixed(1), "s ;", fails ? `${fails} FAILURES` : "ALL PASS");
