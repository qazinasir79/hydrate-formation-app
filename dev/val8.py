#!/usr/bin/env python3
"""val8: full validation of a parameter set (v8 style: kp,fb,dK,dpwg,dpgw)."""
import numpy as np, json, sys, time
sys.path.insert(0, '.')
import engine as E
from calib4 import wgf
from calibrate import ANCH, Peq

def load(fname):
    d = np.load(fname, allow_pickle=True).item()
    return d
def apply_params(d):
    Z = {g: (0.0,0.0) for g in ['CH4','C2H6','C3H8','N2','CO2']}
    if 'dK0' in d:
        for g in Z: E.DK_AQ[g] = (d['dK0'][g], d['dK1'][g])
    elif 'dK' in d:
        for g in Z: E.DK_AQ[g] = (d['dK'][g], 0.0)
    return d['kp'], d['fb']

def region_of(g, T, P):
    if T < 273.15: return 'I'
    Q2 = {'C2H6': 287.75, 'CO2': 283.0}
    if g in Q2 and T > Q2[g]: return 'LC'
    return 'V'

def main(fname):
    t0 = time.time()
    d = load(fname)
    kp, fb = apply_params(d)
    print(f"=== {fname} ===")
    tot = {}
    for g in ['CH4','C2H6','C3H8','N2','CO2']:
        rows = []
        for (T, P, w, st) in sorted(ANCH[g]):
            Pm = Peq(T, g, st, kp, fb)
            dev = (Pm-P)/P*100 if Pm else np.nan
            rows.append((T, P/1e6, (Pm/1e6 if Pm else np.nan), dev, region_of(g, T, P)))
        aadp = np.nanmean([abs(r[3]) for r in rows])
        print(f"{g}: AADP={aadp:6.2f}%")
        for r in rows:
            flag = '  <<<' if abs(r[3]) > 25 else ''
            pm = f"{r[2]:8.4f}" if np.isfinite(r[2]) else "     fail"
            print(f"   T={r[0]:7.2f} Pexp={r[1]:8.4f} Pmod={pm} dev={r[3]:+7.1f}% {r[4]}{flag}")
        tot[g] = aadp
    # structures
    print("structures:", end=' ')
    for g, st in [('CH4','sI'),('C2H6','sI'),('CO2','sI'),('C3H8','sII'),('N2','sII')]:
        T = 277.15 if g != 'CH4' else 273.15
        s1 = Peq(T, g, 'sI', kp, fb); s2 = Peq(T, g, 'sII', kp, fb)
        sel = 'sI' if (s1 and (not s2 or s1 <= s2)) else 'sII'
        print(f"{g}:{sel}({st})", end=' ')
    print()
    # occupancies
    for (g, T, P, thl, ths) in [('CH4',273.2,2.60e6,0.969,0.935),('CH4',276.7,3.44e6,0.981,0.917),('CO2',273.2,1.28e6,0.986,0.581)]:
        fLw, fg, tag = wgf(T, P, {g:1.0}, use_cache=False)
        fH, thetas, S, fbeta = E.fH_w(T, P, fg, 'sI', fb, kp)
        print(f"occ {g}@{T}: thL={thetas[(g,1)]:.4f} ({thl}) thS={thetas[(g,0)]:.4f} ({ths})")
    # fig5 + NG
    from calib4 import wgf as _w
    def peqmix(T, gf):
        s1 = None; s2 = None
        try: s1 = Peq(T, gf, 'sI', kp, fb) if isinstance(gf, dict) else None
        except Exception: pass
        return s1
    # mixture Peq via brentq on r for each structure
    def peq_mix(T, gf, st):
        from scipy.optimize import brentq
        def rr(lnP):
            P = np.exp(lnP)
            fLw, fg, tag = _w(T, P, gf, use_cache=False)
            if fLw is None: return -1.0
            fH = E.fH_w(T, P, fg, st, fb, kp)[0]
            return np.log(fH/fLw)
        lo, hi = np.log(2e4), np.log(2e8)
        try:
            if rr(lo)*rr(hi) > 0: return None
            return np.exp(brentq(rr, lo, hi, xtol=1e-10))
        except Exception:
            return None
    gf5 = {'CH4':0.9995,'C3H8':0.0005}
    p1 = peq_mix(277.6, gf5, 'sI'); p2 = peq_mix(277.6, gf5, 'sII')
    print(f"Fig5 @277.6 y=0.9995: sI={p1 and p1/1e6:.3f} sII={p2 and p2/1e6:.3f} -> {'sII' if p2 and (not p1 or p2<p1) else 'sI'}")
    NG = {'CH4':0.932,'C2H6':0.0425,'C3H8':0.0161,'N2':0.0043,'CO2':0.0051}
    for T in (277.65, 287.65, 296.65):
        p1 = peq_mix(T, NG, 'sI'); p2 = peq_mix(T, NG, 'sII')
        sel = ('sI', p1) if (p1 and (not p2 or p1 <= p2)) else ('sII', p2)
        print(f"NG @{T}: {sel[0]} {sel[1]/1e6 if sel[1] else None:.3f} MPa")
    # solubility
    def solub_x(g, T, P):
        z=np.zeros(6); z[E.IW]=0.99; z[E.NAMES.index(g)]=0.01
        fl=E.flash_TP(z,T,P,'water_first',maxiter=200,tol=1e-10)
        if fl is None: return None
        x,y,beta,conv=fl
        xaq,yg=(x,y) if x[E.IW]>=y[E.IW] else (y,x)
        return float(xaq[E.NAMES.index(g)])
    SOL = {'CH4':[(273.15,2.6e6,1.05e-3),(298.15,0.101325e6,2.50e-5)],
           'C2H6':[(298.15,0.101325e6,2.00e-5),(277.6,0.9e6,1.3e-3)],
           'C3H8':[(298.15,0.101325e6,2.60e-6)],
           'N2':[(298.15,0.101325e6,1.15e-5),(277.15,10e6,1.0e-3)],
           'CO2':[(278.15,3.5e6,1.20e-2),(298.15,0.101325e6,6.00e-4)]}
    for g in SOL:
        for (T,P,xe) in SOL[g]:
            xc = solub_x(g,T,P)
            print(f"solub {g} {T}/{P/1e6:.2f}MPa: x={xc:.3e} exp={xe:.2e} ratio={xc/xe if xc else float('nan'):.2f}")
    print(f"[{time.time()-t0:.0f}s]")

if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'calib18.npy')
