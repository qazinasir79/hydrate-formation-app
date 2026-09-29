import numpy as np, json, sys, math, itertools
sys.path.insert(0,'.')
import fitv13 as F
import engine as E
from scipy.optimize import brentq
d = json.load(open('params_v8.json'))
base0 = np.array([d['dK0'][g] for g in F.G5]); base1 = np.array([d['dK1'][g] for g in F.G5])

def aadp(p):
    dk = F.unpack(p)
    res = {}
    for g in F.G5:
        devs = []
        for (T, P, w0, st) in sorted(E.__dict__.get('_x', []) or []) : pass
        from calibrate import ANCH
        for (T, Pexp, w0, st) in sorted(ANCH[g]):
            lo, hi = math.log(Pexp/4), math.log(Pexp*8)
            f = lambda lp: F._r_of(T, {g:1.0}, st, math.exp(lp), dk)
            rl, rh = f(lo), f(hi)
            Pm = None
            if rl is not None and rh is not None and rl*rh < 0:
                Pm = math.exp(brentq(f, lo, hi, xtol=1e-9))
            if Pm: devs.append(abs(Pm/Pexp-1)*100)
        res[g] = float(np.mean(devs)) if devs else float('nan')
    return res

for ch4, c2h6 in itertools.product([-0.04,-0.06,-0.10,-0.164], [-0.04,-0.06,-0.113]):
    p = np.concatenate([base0, base1]).copy()
    p[0] = ch4; p[1] = c2h6; p[5] = 0.0; p[6] = max(base1[1], 0.001)
    c, r = F.resid(p)[0], None
    rr = F.resid(p)
    cost = 0.5*float(np.sum(rr**2))
    a = aadp(p)
    print(f'CH4 dK0={ch4:+.3f} C2H6 dK0={c2h6:+.3f} cost={cost:.3f} | AADP ' + ' '.join(f'{g}:{a[g]:.1f}' for g in F.G5))
