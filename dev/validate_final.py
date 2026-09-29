import numpy as np, time
import engine as E, calib4 as C4, calibrate as CB
np.seterr(all='ignore')
d=np.load('calib7.npy',allow_pickle=True).item()
kp,fb=d['kp'],d['fb']
t0=time.time()
def fgas_at(gas,T,P,struct):
    fLw,fg,tag=C4.wgf(T,P,{gas:1.0},use_cache=False)
    return fLw,fg,tag
print("=== Peq vs anchors (paper Table 3/4/5-style AADP) ===")
alldev=[]
for gas in ['CH4','C2H6','C3H8','CO2','N2']:
    devs=[]
    print(f"--- {gas} ---")
    for (T,Pref,w,st_exp) in sorted(CB.ANCH[gas]):
        Ps={}
        for st in ('sI','sII'):
            try: Ps[st]=CB.Peq(T,gas,st,kp,fb,Pguess=Pref)
            except Exception: Ps[st]=None
        cand=[(p,s) for s,p in Ps.items() if p]
        if not cand:
            print(f"  T={T:7.2f} Pref={Pref/1e6:8.3f} -> FAIL"); continue
        P,ss=min(cand)
        dev=P/Pref-1; devs.append(abs(dev))
        flag="" if ss==st_exp else f"  << expected {st_exp}"
        print(f"  T={T:7.2f} Pref={Pref/1e6:8.3f} Pcalc={P/1e6:8.3f} dev={100*dev:+6.2f}% [{ss}]{flag}")
    aad=100*np.mean(devs); alldev.append((gas,aad))
    print(f"  {gas} AADP = {aad:.2f}%")
print("\nAADP summary:", {g:round(a,2) for g,a in alldev})
print("paper:        CH4 1.37, C2H6 2.47, C3H8 3.83, CO2 2.89, N2 6.84")
print(f"[{time.time()-t0:.0f}s]")
print("\n=== occupancies vs paper Tables 6-8 ===")
for gas,struct in [('CH4','sI'),('CO2','sI'),('N2','sII')]:
    for T in [273.2,278.2,283.2]:
        try: P=CB.Peq(T,gas,struct,kp,fb,Pguess=None)
        except Exception: P=None
        if not P: continue
        fLw,fg,tag=fgas_at(gas,T,P,struct)
        if fg is None: continue
        fH,th,S,fbeta=E.fH_w(T,P,fg,struct,fb,kp)
        if struct=='sI': tl,ts=th[(gas,1)],th[(gas,0)]
        else: tl,ts=th[(gas,1)],th[(gas,0)]
        print(f"  {gas} [{struct}] T={T}: P={P/1e6:6.3f} MPa theta_large={tl:.4f} theta_small={ts:.4f} (fH/fLw={fH/fLw:.4f})")
print("  paper CH4 sI @273.2: 0.969 / 0.935 ; CO2 sI @273.2: 0.986 / 0.581")
print(f"\n=== Fig5 sII onset (CH4+C3H8): transition near y_C3H8=0.0005 @277.6K ===")
for y3 in [1e-5,1e-4,3e-4,5e-4,1e-3,3e-3,1e-2]:
    gf={'CH4':1-y3,'C3H8':y3}
    row=[]
    for st in ('sI','sII'):
        # solve Peq for mixture via bracketing on resid
        import scipy.optimize as so
        def r(lnP):
            P=np.exp(lnP)
            fLw,fg,tag=C4.wgf(277.6,P,gf,use_cache=False)
            if fLw is None or fg is None: return np.nan
            fH,_,_,_=E.fH_w(277.6,P,fg,st,fb,kp)
            v=np.log(fH)-np.log(fLw)
            return v if np.isfinite(v) else np.nan
        lns=np.linspace(np.log(1e5),np.log(3e7),40); vals=np.array([r(L) for L in lns])
        P=None
        for k in range(len(lns)-1):
            if np.isfinite(vals[k]) and np.isfinite(vals[k+1]) and vals[k]*vals[k+1]<0:
                P=float(np.exp(so.brentq(r,lns[k],lns[k+1],xtol=1e-8))); break
        row.append(P)
    Ps=[p for p in row if p]
    ss=('sI' if row[0] and row[1] and row[0]<row[1] else ('sII' if row[1] else 'sI'))
    print(f"  y_C3={y3:8.1e}: P_sI={row[0]/1e6 if row[0] else float('nan'):8.4f} P_sII={row[1]/1e6 if row[1] else float('nan'):8.4f} -> {ss}")
print(f"[{time.time()-t0:.0f}s]")
