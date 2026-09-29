#!/usr/bin/env python3
"""fitv12: v8 calibration.
params (28): [ln eps, sigma]*5  +  [A,B,c]  +  dK_aq*5  +  dp_wg*5 + dp_gw*5
residuals: anchor fugacity residuals (w*r) + occupancy (Tables 6-8) + Fig5 transition
           + structure hinges + weak P-P regularization + solubility (ln x) anchors.
"""
import numpy as np, json, time, sys
from scipy.optimize import least_squares
sys.path.insert(0, '.')
import engine as E
import calib4
from calib4 import wgf, _cache, cpp, PPAB, TPP, WREG
from calibrate import ANCH

GASES = ['CH4', 'C2H6', 'C3H8', 'N2', 'CO2']
d0 = np.load('calib10.npy', allow_pickle=True).item()
KP7, FB7 = d0['kp'], d0['fb']
P_WG = {g: E.KIJ_P[g][0] for g in GASES}
P_GW = {g: E.KIJ_P[g][2] for g in GASES}

SOL = {
 'CH4': [(273.15, 2.6e6, 1.05e-3), (298.15, 0.101325e6, 2.50e-5), (273.15, 0.101325e6, 4.15e-5),
         (283.15, 5.0e6, 1.55e-3)],
 'C2H6': [(298.15, 0.101325e6, 2.00e-5), (273.15, 0.101325e6, 4.60e-5), (277.6, 0.9e6, 1.3e-3)],
 'C3H8': [(298.15, 0.101325e6, 2.60e-6), (273.15, 0.101325e6, 3.50e-6)],
 'N2':  [(298.15, 0.101325e6, 1.15e-5), (273.15, 0.101325e6, 1.35e-5), (277.15, 10e6, 1.0e-3)],
 'CO2': [(278.15, 3.5e6, 1.20e-2), (298.15, 0.101325e6, 6.00e-4), (273.15, 0.101325e6, 1.34e-3)],
}
OCC = [('CH4', 273.2, 2.60e6, 0.969, 0.935), ('CH4', 276.7, 3.44e6, 0.981, 0.917), ('CO2', 273.2, 1.28e6, 0.986, 0.581)]

def unpack(p):
    kp = {}; i = 0
    for g in GASES:
        kp[g] = (np.exp(p[i]), p[i+1], E.KIH[g][2]); i += 2
    A, B, c = p[i], p[i+1], p[i+2]; i += 3
    fb = {'sI': (A, B), 'sII': (A+c, B)}
    dK0 = {g: p[i+k] for k, g in enumerate(GASES)}; i += 5
    dK1 = {g: p[i+k] for k, g in enumerate(GASES)}; i += 5
    return kp, fb, dK0, dK1

def apply(p):
    kp, fb, dK0, dK1 = unpack(p)
    for g in GASES: E.DK_AQ[g] = (dK0[g], dK1[g])
    return kp, fb

def solub_x(g, T, P):
    z = np.zeros(6); z[E.IW] = 0.99; z[E.NAMES.index(g)] = 0.01
    fl = E.flash_TP(z, T, P, 'water_first', maxiter=200, tol=1e-10)
    if fl is None: return None
    x, y, beta, conv = fl
    xaq, yg = (x, y) if x[E.IW] >= y[E.IW] else (y, x)
    return float(xaq[E.NAMES.index(g)])

from scipy.optimize import brentq
def _r_of(T, gf, st, P, kp, fb):
    fLw, fg, tag = wgf(T, P, gf, use_cache=False)
    if fLw is None: return None
    fH = E.fH_w(T, P, fg, st, fb, kp)[0]
    return np.log(fH/fLw)

def root_dev(T, gf, st, Pexp, kp, fb):
    """ln(P_root/P_exp) using local bracket; penalty if no local crossing."""
    lo, hi = np.log(Pexp/2.5), np.log(Pexp*2.5)
    f = lambda lp: _r_of(T, gf, st, np.exp(lp), kp, fb)
    rl = f(lo); rh = f(hi)
    if rl is None or rh is None: return 3.0
    if rl*rh > 0:
        grid = np.linspace(lo, hi, 13)
        vals = [f(v) for v in grid]
        seg = None
        for a, b, va, vb in zip(grid[:-1], grid[1:], vals[:-1], vals[1:]):
            if va is None or vb is None: continue
            if va*vb < 0:
                seg = (a, b)
                if b >= np.log(Pexp): break
        if seg is None:
            return 2.0 + min(abs(rl), abs(rh))
        lo, hi = seg
    try:
        Pr = brentq(f, lo, hi, xtol=1e-9, rtol=1e-12)
    except Exception:
        return 3.0
    return Pr - np.log(Pexp)

def resid_all(p):
    _cache.clear()
    kp, fb = apply(p)
    out = []
    for gas in GASES:
        gf = {gas: 1.0}
        for (T, P, w, st) in sorted(ANCH[gas]):
            rd = root_dev(T, gf, st, P, kp, fb)
            r0 = _r_of(T, gf, st, P, kp, fb) or 3.0
            out.append(w*(rd + 0.3*np.clip(r0, -3, 3)))
    # occupancies at model Peq
    for (g, T, P, thl, ths) in OCC:
        fLw, fg, tag = wgf(T, P, {g: 1.0}, use_cache=False)
        if fLw is None: out += [3.0, 3.0]; continue
        fH, thetas, S, fbeta = E.fH_w(T, P, fg, 'sI', fb, kp)
        out.append(4.0*(thetas[(g, 1)] - thl))
        out.append(4.0*(thetas[(g, 0)] - ths))
    # P-P regularization
    for gas in GASES:
        if gas in PPAB:
            st = 'sII' if gas in ('N2',) else 'sI'
            for T in TPP:
                for m in (0, 1):
                    ct = cpp(gas, T, m)
                    if ct:
                        C = E.langmuir(T, E.CAGES[st][m][0], E.CAGES[st][m][1], *kp[gas])
                        out.append(WREG*np.log(C/ct))
    # Fig5 transition
    gf5 = {'CH4': 0.9995, 'C3H8': 0.0005}
    fLw, fg, tag = wgf(277.6, 4.0e6, gf5, use_cache=False)
    if fLw is not None:
        for st in ('sI', 'sII'):
            fH, _, _, _ = E.fH_w(277.6, 4.0e6, fg, st, fb, kp)
            out.append(1.5*(np.log(fH)-np.log(fLw)))
    # hinges
    for (T, P) in [(273.15, 2.60e6), (283.15, 7.20e6)]:
        fLw, fg, tag = wgf(T, P, {'CH4': 1.0}, use_cache=False)
        if fLw is None: continue
        fH2, _, _, _ = E.fH_w(T, P, fg, 'sII', fb, kp)
        out.append(3.0*min(0.0, np.log(fH2)-np.log(fLw)-0.03))
    for (T, P) in [(277.59, 0.42e6), (273.15, 0.172e6)]:
        fLw, fg, tag = wgf(T, P, {'C3H8': 1.0}, use_cache=False)
        if fLw is None: continue
        fH1, _, _, _ = E.fH_w(T, P, fg, 'sI', fb, kp)
        out.append(3.0*min(0.0, np.log(fH1)-np.log(fLw)-0.03))
    fLw, fg, tag = wgf(273.15, 16.4e6, {'N2': 1.0}, use_cache=False)
    if fLw is not None:
        fH1, _, _, _ = E.fH_w(273.15, 16.4e6, fg, 'sI', fb, kp)
        out.append(3.0*min(0.0, np.log(fH1)-np.log(fLw)-0.03))
    # solubility
    for g in GASES:
        for (T, P, xe) in SOL[g]:
            xc = solub_x(g, T, P)
            if xc is None or xc <= 0: out.append(0.8*3.0); continue
            out.append(0.8*np.log(xc/xe))
    return np.clip(np.array(out), -6.0, 6.0)

if __name__ == '__main__':
    t0 = time.time()
    p0 = []
    for g in GASES: p0 += [np.log(KP7[g][0]), KP7[g][1]]
    p0 += [FB7['sI'][0], FB7['sI'][1], FB7['sII'][0]-FB7['sI'][0]]
    p0 += [0.0]*5 + [0.0]*5
    p0 = np.array(p0)
    r0 = resid_all(p0)
    print(f"start cost={0.5*np.sum(r0**2):.4f} n={len(r0)} [{time.time()-t0:.0f}s]")
    lo = []; hi = []
    for g in GASES: lo += [np.log(80), 2.5e-10]; hi += [np.log(400), 4.0e-10]
    lo += [26.0, 5200.0, -0.30]; hi += [31.0, 6600.0, 0.10]
    lo += [-2.50]*5 + [-0.020]*5
    hi += [+0.10]*5 + [+0.020]*5
    sol = least_squares(resid_all, p0, bounds=(lo, hi), diff_step=1e-4,
                        xtol=1e-10, ftol=1e-10, max_nfev=600, verbose=1)
    print(f"fit done cost={sol.cost:.5f} nfev={sol.nfev} [{time.time()-t0:.0f}s]")
    kp, fb, dK0, dK1 = unpack(sol.x)
    for g in GASES:
        print(f"  {g}: eps/k={kp[g][0]:8.2f} sigma={kp[g][1]*1e10:.4f} | dK0={dK0[g]:+.3f} dK1={dK1[g]:+.5f}")
    print(f"  fbeta sI A={fb['sI'][0]:.4f} B={fb['sI'][1]:.2f}; sII A={fb['sII'][0]:.4f}")
    np.save('calib17.npy', {'kp': kp, 'fb': fb, 'dK0': dK0, 'dK1': dK1, 'p': sol.x}, allow_pickle=True)
    json.dump({'kp': {g: list(map(float, kp[g])) for g in GASES},
               'fb': {k: list(map(float, v)) for k, v in fb.items()},
               'dK0': {g: float(dK0[g]) for g in GASES},
               'dK1': {g: float(dK1[g]) for g in GASES},
               'assign': 'water_first'}, open('params_v8.json', 'w'), indent=1)
    print("saved calib17.npy / params_v8.json")
