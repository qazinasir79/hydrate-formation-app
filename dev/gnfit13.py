# -*- coding: utf-8 -*-
import numpy as np, json, time
import fitv13 as F

LO = np.array([-2.6]*5 + [-0.02]*5)
HI = np.array([+0.1]*5 + [+0.02]*5)

def clip(p): return np.clip(p, LO, HI)

def cost(p):
    r = F.resid(p)
    return 0.5*float(np.sum(r**2)), r

def jac(p, r0, eps=2e-4):
    n = len(p); J = np.zeros((len(r0), n))
    for k in range(n):
        pp = p.copy()
        h = eps*max(1.0, abs(p[k])) if k < 5 else eps*0.5
        pp[k] = min(p[k]+h, HI[k])
        hm = pp[k]-p[k]
        _, r1 = cost(pp)
        J[:, k] = (r1-r0)/hm
    return J

p = clip(F.resid and np.load('p0_v13.npy'))
c, r = cost(p)
print('start cost=%.4f' % c, flush=True)
lam = 1e-3
t00 = time.time()
for it in range(22):
    J = jac(p, r)
    JTJ = J.T @ J; JTr = J.T @ r
    improved = False
    for _ in range(8):
        try:
            dp = np.linalg.solve(JTJ + lam*(np.diag(np.diag(JTJ))+1e-12*np.eye(len(p))), -JTr)
        except np.linalg.LinAlgError:
            lam *= 10; continue
        for alpha in [1.0, 0.5, 0.25]:
            pn = clip(p + alpha*dp)
            cn, rn = cost(pn)
            if cn < c:
                print(f'it{it}: cost {c:.4f} -> {cn:.4f} (alpha={alpha}, lam={lam:.1e}, |d|={np.linalg.norm(dp):.2e}) [{time.time()-t00:.0f}s]', flush=True)
                c, r, p = cn, rn, pn
                lam = max(lam*0.5, 1e-6); improved = True; break
        if improved: break
        lam *= 6
        if lam > 1e6: break
    if not improved:
        print(f'it{it}: no improvement (lam={lam:.1e}); stop', flush=True); break

dk = F.unpack(p)
out = {'kp': {g: list(map(float, F.KP7[g])) for g in F.G5},
       'fb': {s: list(map(float, F.FB7[s])) for s in ['sI', 'sII']},
       'dK0': {g: float(dk[g][0]) for g in F.G5},
       'dK1': {g: float(dk[g][1]) for g in F.G5}}
np.save('calib18.npy', out)
json.dump(out, open('params_v8.json', 'w'), indent=1)
print('saved. final cost=%.4f' % c, flush=True)
for g in F.G5: print(f'  {g}: dK0={dk[g][0]:+.4f} dK1={dk[g][1]:+.6f}')
