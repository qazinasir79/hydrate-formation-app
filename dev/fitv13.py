# -*- coding: utf-8 -*-
"""dK-only root-space fit. Kihara + fbeta FROZEN at v7 (calib10).
p = [dK0 x5, dK1 x5]. Ice anchors excluded (dK inactive there)."""
import numpy as np, sys, math, time
sys.path.insert(0, '.')
import engine as E
from scipy.optimize import brentq
from calibrate import ANCH
from calib4 import wgf

IW = E.IW; NAMES = E.NAMES
d7 = np.load('calib10.npy', allow_pickle=True).item()
KP7, FB7 = d7['kp'], d7['fb']
G5 = ['CH4', 'C2H6', 'C3H8', 'N2', 'CO2']

def unpack(p):
    return {g: (p[i], p[5+i]) for i, g in enumerate(G5)}

def _r_of(T, gf, st, P, dk):
    E.DK_AQ = dk
    w = wgf(T, P, gf, use_cache=False)
    if w[0] is None: return None
    return math.log(E.fH_w(T, P, w[1], st, FB7, KP7)[0] / w[0])

def root_dev(T, gf, st, Pexp, dk):
    lo, hi = math.log(Pexp/1.8), math.log(Pexp*4.0)
    f = lambda lp: _r_of(T, gf, st, math.exp(lp), dk)
    rl = f(lo); rh = f(hi)
    if rl is None or rh is None: return 3.0
    if rl*rh > 0:
        grid = np.linspace(lo, hi, 9)
        vals = [f(v) for v in grid]
        seg = None
        for a, b, va, vb in zip(grid[:-1], grid[1:], vals[:-1], vals[1:]):
            if va is None or vb is None: continue
            if va*vb < 0:
                seg = (a, b)
                if b >= math.log(Pexp): break
        if seg is None: return 2.0 + min(abs(rl), abs(rh))
        lo, hi = seg
    try:
        Pr = brentq(f, lo, hi, xtol=1e-8, rtol=1e-11)
    except Exception:
        return 3.0
    return Pr - math.log(Pexp)

SOL = [
    ('CH4', 273.15, 2.60e6, 1.05e-3, 0.6), ('CH4', 298.15, 0.101325e6, 2.5e-5, 0.6),
    ('CH4', 283.15, 5.0e6, 1.30e-3, 0.6),
    ('C2H6', 298.15, 0.101325e6, 2.0e-5, 0.6), ('C2H6', 277.6, 0.90e6, 1.3e-3, 0.6),
    ('C3H8', 298.15, 0.101325e6, 2.6e-6, 0.6),
    ('N2', 298.15, 0.101325e6, 1.15e-5, 0.6), ('N2', 277.15, 10.0e6, 1.0e-3, 0.6),
    ('CO2', 278.15, 3.5e6, 1.2e-2, 0.6), ('CO2', 298.15, 0.101325e6, 6.0e-4, 0.6),
]

# anchor groups: LC/high-T problem anchors w=1; V guards w=0.3; NG/Fig5 w=0.3
PROB = {('C2H6',288.15),('C2H6',293.15),('C2H6',298.15),('C2H6',303.15),
        ('N2',285.15),('N2',289.15),('N2',293.15),
        ('CO2',285.15),('CO2',287.15),('CO2',289.15),('CO2',291.15)}
ANCHS = []
for g in G5:
    for (T, P, w0, st) in sorted(ANCH[g]):
        if T < 273.2: continue          # ice: dK inactive
        w = 1.0 if (g, T) in PROB else 1.0
        ANCHS.append((g, T, P, w, st))
NG = {'CH4': 0.932, 'C2H6': 0.0425, 'C3H8': 0.0161, 'N2': 0.0043, 'CO2': 0.0051}
MIX = [(277.65, NG, 0.0636e6, 'sII'), (287.65, NG, 0.3243e6, 'sII'), (296.65, NG, 1.2112e6, 'sII'),
       (277.6, {'CH4': 0.9997, 'C3H8': 0.0003}, 4.0e6, 'sI')]

def resid(p):
    dk = unpack(p)
    out = []
    for (g, T, Pexp, w, st) in ANCHS:
        rd = root_dev(T, {g: 1.0}, st, Pexp, dk)
        r0 = _r_of(T, {g: 1.0}, st, Pexp, dk)
        out.append(w*(rd + (0.3*max(-3, min(3, r0)) if r0 is not None else 0.9)))
    for (g, T, Pp, xexp, w) in SOL:
        z = np.zeros(6); z[IW] = 0.99; z[NAMES.index(g)] = 0.01
        E.DK_AQ = dk
        fl = E.flash_TP(z, T, Pp, 'water_first', maxiter=300, tol=1e-11)
        if fl[0] is None: out.append(w*6.0); continue
        x, y = fl[0], fl[1]
        xaq = x if x[IW] >= y[IW] else y
        out.append(w*max(-3.0, min(3.0, math.log(xaq[NAMES.index(g)]/xexp))))
    E.DK_AQ = {g: (0.0, 0.0) for g in G5}
    return np.array(out)

if __name__ == '__main__':
    import json
    d = json.load(open('params_v8.json'))
    p0 = np.array([d['dK0'][g] for g in G5] + [d['dK1'][g] for g in G5])
    t0 = time.time()
    r = resid(p0)
    print('nres', len(r), 'cost', 0.5*float(np.sum(r**2)), f'[{time.time()-t0:.1f}s]')
    np.save('p0_v13.npy', p0)
