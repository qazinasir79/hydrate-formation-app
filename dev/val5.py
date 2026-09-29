import numpy as np, sys, time
import engine as E, calib4 as C4, calibrate as CB
np.seterr(all='ignore')
f=sys.argv[1]
d=np.load(f,allow_pickle=True).item(); kp,fb=d['kp'],d['fb']
t0=time.time()
summary={}
details={}
for gas in ['CH4','C2H6','C3H8','CO2','N2']:
    devs=[]; det=[]
    for (T,Pref,w,st) in sorted(CB.ANCH[gas]):
        Ps={}
        for s in ('sI','sII'):
            try: Ps[s]=CB.Peq(T,gas,s,kp,fb,Pguess=Pref)
            except Exception: Ps[s]=None
        cand=[(p,s) for s,p in Ps.items() if p]
        if not cand: det.append((T,Pref,None,None)); devs.append(1.0); continue
        P,ss=min(cand)
        devs.append(abs(P/Pref-1)); det.append((T,Pref,P,ss,st))
    summary[gas]=100*np.mean(devs); details[gas]=det
for gas in summary:
    print(f"--- {gas}: AADP={summary[gas]:.2f}% ---")
    for (T,Pref,P,ss,st) in details[gas]:
        if P is None: print(f"  T={T:7.2f} Pref={Pref/1e6:8.3f} FAIL"); continue
        flag="" if ss==st else f" <<exp {st}"
        print(f"  T={T:7.2f} Pref={Pref/1e6:8.3f} P={P/1e6:8.3f} dev={100*(P/Pref-1):+7.2f}% [{ss}]{flag}")
print("AADP:",{g:round(v,2) for g,v in summary.items()})
print("paper: CH4 1.37, C2H6 2.47(LC 11.15), C3H8 3.83, CO2 2.89(LC 10.2), N2 6.84")
print(f"[{time.time()-t0:.0f}s]")
