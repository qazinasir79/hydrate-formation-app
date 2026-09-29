#!/usr/bin/env python3
"""gnfit: damped Gauss-Newton with monotone line search for v8 (28 params)."""
import numpy as np, json, sys, time
sys.path.insert(0, '.')
import fitv12 as F

G = F.GASES
STEPS = np.array([1e-4]*10 + [1e-3, 1e-2, 1e-3] + [1e-4]*5 + [1e-5]*5)
STEPS[1::2][:5] = 3e-14   # sigma steps (indices 1,3,5,7,9)

def cost(p):
    r = F.resid_all(p)
    return 0.5*float(r @ r), r

def jac(p, r0):
    n = len(p)
    J = np.zeros((len(r0), n))
    for i in range(n):
        pp = p.copy(); pp[i] += STEPS[i]
        # respect bounds
        pp[i] = min(max(pp[i], LO[i]), HI[i])
        ri = F.resid_all(pp)
        J[:, i] = (ri - r0)/STEPS[i]
    return J

import json as _json
_V7 = _json.load(open('params_v7.json'))
LO = []; HI = []
for g in G:
    e7, s7 = _V7['kp'][g][0], _V7['kp'][g][1]
    LO += [np.log(e7)-0.25, 0.92*s7]; HI += [np.log(e7)+0.25, 1.08*s7]
LO += [28.4, 5700.0, -0.30]; HI += [29.2, 6200.0, 0.10]
LO += [-2.50]*5 + [-0.020]*5
HI += [+0.10]*5 + [+0.020]*5
LO = np.array(LO); HI = np.array(HI)

if __name__ == '__main__':
    t0 = time.time()
    d7 = np.load('calib10.npy', allow_pickle=True).item()
    p = []
    for gg in G: p += [np.log(d7['kp'][gg][0]), d7['kp'][gg][1]]
    p += [d7['fb']['sI'][0], d7['fb']['sI'][1], d7['fb']['sII'][0]-d7['fb']['sI'][0]]
    p += [-0.20, -0.10, 0.00, -0.40, -2.10] + [0.0]*5
    p = np.clip(np.array(p), LO, HI)
    c, r = cost(p)
    print(f"start cost={c:.4f}", flush=True)
    lam = 1e-3
    for it in range(16):
        J = jac(p, r)
        JtJ = J.T @ J; Jtr = J.T @ r
        dg = np.diag(JtJ).copy(); dg[dg < 1e-12] = 1e-12
        ok = False
        for _ in range(6):
            try:
                delta = np.linalg.solve(JtJ + lam*np.diag(dg), -Jtr)
            except np.linalg.LinAlgError:
                lam *= 4; continue
            for alpha in (1.0, 0.5, 0.25, 0.125, 0.0625, 0.03):
                pn = np.clip(p + alpha*delta, LO, HI)
                cn, rn = cost(pn)
                if cn < c - 1e-7:
                    print(f"it{it}: cost {c:.4f} -> {cn:.4f} (alpha={alpha}, lam={lam:.1e}, |d|={np.linalg.norm(delta):.2e}) [{time.time()-t0:.0f}s]", flush=True)
                    p, c, r = pn, cn, rn
                    lam = max(lam*0.5, 1e-8)
                    ok = True
                    break
            if ok: break
            lam *= 4
        if not ok:
            print(f"it{it}: no descent (lam={lam:.1e}); stop", flush=True)
            break
        if np.linalg.norm(delta) < 1e-7:
            print("step tiny; stop", flush=True)
            break
    kp, fb, dK0, dK1 = F.unpack(p)
    for g in G:
        print(f"  {g}: eps/k={kp[g][0]:8.2f} sigma={kp[g][1]*1e10:.4f} | dK0={dK0[g]:+.3f} dK1={dK1[g]:+.5f}")
    print(f"  fbeta sI A={fb['sI'][0]:.4f} B={fb['sI'][1]:.2f}; sII A={fb['sII'][0]:.4f} B={fb['sII'][1]:.2f}")
    np.save('calib18.npy', {'kp': kp, 'fb': fb, 'dK0': dK0, 'dK1': dK1, 'p': p}, allow_pickle=True)
    json.dump({'kp': {g: list(map(float, kp[g])) for g in G}, 'fb': {k: list(map(float, v)) for k, v in fb.items()},
               'dK0': {g: float(dK0[g]) for g in G}, 'dK1': {g: float(dK1[g]) for g in G},
               'assign': 'water_first'}, open('params_v8.json', 'w'), indent=1)
    print(f"saved calib18/params_v8 final cost={c:.4f} [{time.time()-t0:.0f}s]")
