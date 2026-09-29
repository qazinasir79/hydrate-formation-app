#!/usr/bin/env python3
"""patch_v8: inject final v8 params into engine_port.js and regenerate JS test refs."""
import json, re, subprocess, sys, os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)
P = json.load(open('params_v8.json'))
G = ['CH4', 'C2H6', 'C3H8', 'N2', 'CO2']
s = open('engine_port.js').read()

def fmtarr(v): return '[' + ', '.join(repr(float(x)) for x in v) + ']'

# KIH block
kih = 'const KIH = {\n' + ''.join(
    f"  {g}: {fmtarr([P['kp'][g][0], P['kp'][g][1], P['kp'][g][2]])},\n" for g in G) + '};'
s = re.sub(r'const KIH = \{.*?\};', kih, s, flags=re.S)
# FBETA block
fb = ("const FBETA = {\n"
      f"  sI: {fmtarr(P['fb']['sI'])},\n"
      f"  sII: {fmtarr(P['fb']['sII'])},\n}};")
s = re.sub(r'const FBETA = \{.*?\};', fb, s, flags=re.S)
# DK_AQ block
dk = ("const DK_AQ = { " + ', '.join(
    f"{g}: {fmtarr([P['dK0'][g], P['dK1'][g]])}" for g in G) + " };")
s = re.sub(r'const DK_AQ = \{.*?\};', dk, s, flags=re.S)
open('engine_port.js', 'w').write(s)
print('engine_port.js patched with v8')

# regenerate refs with python engine + v8
import engine as E
E.DK_AQ.update({g: (P['dK0'][g], P['dK1'][g]) for g in G})
sys.path.insert(0, '.')
from calib4 import wgf
from calibrate import Peq
kp = {g: tuple(P['kp'][g]) for g in G}
fb = {k: tuple(v) for k, v in P['fb'].items()}
ref = {}
ref['C_CH4_sI'] = [float(E.langmuir(273.15, *E.CAGES['sI'][0], *kp['CH4'])),
                   float(E.langmuir(273.15, *E.CAGES['sI'][1], *kp['CH4']))]
w = wgf(273.15, 2.6e6, {'CH4': 1.0}, use_cache=False)
fh, th, S, fbeta = E.fH_w(273.15, 2.6e6, w[1], 'sI', fb, kp)
ref['wgf'] = {'fH': float(fh), 'S': float(S), 'theta_L': float(th[('CH4', 1)]), 'theta_S': float(th[('CH4', 0)])}
ref['peq'] = {}
for g, T, st in [('CH4', 273.15, 'sI'), ('CH4', 298.15, 'sI'), ('CO2', 273.15, 'sI'),
                 ('C3H8', 277.59, 'sII'), ('N2', 273.15, 'sII'), ('CH4', 263.15, 'sI')]:
    key = f'{g}_{st}_{T}'
    ref['peq'][key] = float(Peq(T, g, st, kp, fb))
NG = {'CH4': 0.932, 'C2H6': 0.0425, 'C3H8': 0.0161, 'N2': 0.0043, 'CO2': 0.0051}
from scipy.optimize import brentq
def peqmix(T, gf, st):
    def rr(lp):
        Pp = np.exp(lp)
        ww = wgf(T, Pp, gf, use_cache=False)
        if ww[0] is None: return -1.0
        return np.log(E.fH_w(T, Pp, ww[1], st, fb, kp)[0] / ww[0])
    lo, hi = np.log(2e4), np.log(3e8)
    if rr(lo) * rr(hi) > 0: return None
    return float(np.exp(brentq(rr, lo, hi, xtol=1e-10)))
def peqsel(T, gf):
    p1 = peqmix(T, gf, 'sI'); p2 = peqmix(T, gf, 'sII')
    if p1 and (not p2 or p1 <= p2): return [p1, 'sI']
    if p2: return [p2, 'sII']
    return [None, None]
ref['ng_281.15'] = peqsel(281.15, NG)
ref['fig5_y3e4'] = peqsel(277.6, {'CH4': 0.9997, 'C3H8': 0.0003})
ref['solub'] = {}
for g, T, Pp in [('CH4', 273.15, 2.6e6), ('CO2', 278.15, 3.5e6), ('N2', 298.15, 0.101325e6)]:
    z = np.zeros(6); z[E.IW] = 0.99; z[E.NAMES.index(g)] = 0.01
    fl = E.flash_TP(z, T, Pp, 'water_first', maxiter=300, tol=1e-11)
    x, y = fl[0], fl[1]
    xaq = x if x[E.IW] >= y[E.IW] else y
    ref['solub'][g] = float(xaq[E.NAMES.index(g)])
json.dump(ref, open('/tmp/jsref_v8.json', 'w'), indent=1)
print('refs written', json.dumps(ref)[:400])
