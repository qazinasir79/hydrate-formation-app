import numpy as np, time, json
from scipy.optimize import least_squares
import engine as E, calib4 as C
from fitv3 import resid2, unpack2, LO, HI, fitstart, GASES, HINGE, WREG_GAS
np.seterr(all='ignore')

t0=time.time()
# ---------- stage 1: block descent (per-gas with published fbeta, then fbeta) ----------
import calib3 as C3
kp=dict(E.KIH); kp['C3H8']=(220.0,3.40e-10,0.6502e-10)
fb={'sI':E.FBETA['sI'],'sII':E.FBETA['sII']}
for rnd in range(2):
    for gas in ['CH4','CO2','N2','C2H6','C3H8']:
        e,s,sol=C3.fit_gas(gas,kp,fb)
        kp[gas]=(e,s,kp[gas][2])
        print(f"bd r{rnd} {gas}: eps/k={e:8.2f} sigma={s*1e10:.4f} cost={sol.cost:.5f} [{time.time()-t0:.0f}s]",flush=True)
    fb,solb=C3.fit_fbeta(kp,fb)
    print(f"bd r{rnd} fbeta: sI=({fb['sI'][0]:.4f},{fb['sI'][1]:.1f}) sII=({fb['sII'][0]:.4f},{fb['sII'][1]:.1f}) cost={solb.cost:.4f} [{time.time()-t0:.0f}s]",flush=True)
# ---------- stage 2: joint fit with hinges ----------
p0=[]
for g in GASES: p0+=[float(np.log(kp[g][0])),float(kp[g][1])]
c0=fb['sII'][0]-fb['sI'][0]
p0+=[float(fb['sI'][0]),float(fb['sI'][1]),float(min(max(c0,LO[-1]),HI[-1]))]
p0=np.array(p0)
C._cache.clear(); R0=resid2(np.clip(p0,LO,HI))
print(f"joint start cost={0.5*np.sum(R0**2):.5f} nan={np.sum(~np.isfinite(R0))}")
best=None
for cstart in [LO[-1]+1e-6, -0.02, 0.0, 0.04]:
    p=np.clip(p0.copy(),LO,HI); p[-1]=cstart
    sol=fitstart(p,maxn=400)
    if best is None or sol.cost<best.cost: best=sol
    print(f" joint c0={cstart:+.3f}: cost={sol.cost:.5f} nfev={sol.nfev} c={sol.x[-1]:.4f} [{time.time()-t0:.0f}s]",flush=True)
sol=best; p1=sol.x; kp,fb=unpack2(p1)
print(f"\nFINAL cost={sol.cost:.5f}")
for g in GASES: print(f"  {g}: eps/k={kp[g][0]:8.2f} sigma={kp[g][1]*1e10:.4f} a={kp[g][2]*1e10:.4f}")
print(f"  fbeta sI A={fb['sI'][0]:.4f} B={fb['sI'][1]:.2f} c={p1[-1]:+.4f}")
print(f"  fbeta_sI(273.15)={E.psat_beta(273.15,'sI',fb):.1f} Pa fbeta_sII(273.15)={E.psat_beta(273.15,'sII',fb):.1f} Pa")
R=resid2(p1); i=0
for gas in GASES:
    n=len(C.ANCH[gas]); print(f"   {gas} anchors sumsq={np.sum(R[i:i+n]**2):.4f} max|r|={np.max(np.abs(R[i:i+n])):.3f}"); i+=n
for gas in GASES:
    if gas in C.PPAB:
        n=sum(1 for T in C.TPP for m in (0,1) if C.cpp(gas,T,m)); print(f"   {gas} reg sumsq={np.sum(R[i:i+n]**2):.4f}"); i+=n
print(f"   fig5 r={np.round(R[i:i+2],4)}"); i+=2
print(f"   hinges r={np.round(R[i:],4)}")
json.dump({'kp':{g:[float(kp[g][0]),float(kp[g][1]),float(kp[g][2])] for g in GASES},
           'fb':{'sI':[float(fb['sI'][0]),float(fb['sI'][1])],'sII':[float(fb['sII'][0]),float(fb['sII'][1])]}},
          open('params_final.json','w'),indent=1)
np.save('calibfinal.npy',{'kp':kp,'fb':fb,'p':p1},allow_pickle=True)
print(f"saved params_final.json [{time.time()-t0:.0f}s]")
