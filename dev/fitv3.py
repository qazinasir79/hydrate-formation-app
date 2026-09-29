import numpy as np, time, json
from scipy.optimize import least_squares, brentq
import engine as E, calib4 as C
np.seterr(all='ignore')
GASES=C.GASES
WREG_GAS={'CH4':0.30,'CO2':0.30,'C2H6':0.30,'N2':0.12}
HINGE=[(273.15,2.60e6,'CH4','sII'),(283.15,7.20e6,'CH4','sII'),(293.15,19.5e6,'CH4','sII'),
       (277.59,0.42e6,'C3H8','sI'),(273.15,0.172e6,'C3H8','sI'),(268.15,0.10e6,'C3H8','sI'),
       (273.15,16.4e6,'N2','sI'),(281.15,33.0e6,'N2','sI'),(263.15,9.3e6,'N2','sI'),
       (273.15,1.28e6,'CO2','sII'),(273.15,0.55e6,'C2H6','sII')]
HW={'N2':6.0}
def unpack2(p):
    kp={}; i=0
    for g in GASES: kp[g]=(np.exp(p[i]),p[i+1],E.KIH[g][2]); i+=2
    A,B,c=p[i],p[i+1],p[i+2]
    return kp,{'sI':(A,B),'sII':(A+c,B)}
def resid2(p):
    kp,fb=unpack2(p); out=[]
    for gas in GASES:
        for (T,P,w,st) in sorted(C.ANCH[gas]):
            fLw,fg,tag=C.wgf(T,P,{gas:1.0})
            if fLw is None: out.append(w*3.0); continue
            fH,_,_,_=E.fH_w(T,P,fg,st,fb,kp)
            r=np.log(fH)-np.log(fLw)
            out.append(w*(r if np.isfinite(r) else 3.0))
    for gas in GASES:
        if gas in C.PPAB:
            st='sII' if gas=='N2' else 'sI'; wr=WREG_GAS[gas]
            for T in C.TPP:
                for m in (0,1):
                    ct=C.cpp(gas,T,m)
                    if ct:
                        Cm=E.langmuir(T,E.CAGES[st][m][0],E.CAGES[st][m][1],*kp[gas])
                        out.append(wr*np.log(Cm/ct))
    fLw,fg,tag=C.wgf(277.6,4.0e6,{'CH4':0.9995,'C3H8':0.0005})
    if fLw is not None:
        for st in ('sI','sII'):
            fH,_,_,_=E.fH_w(277.6,4.0e6,fg,st,fb,kp)
            out.append(1.5*(np.log(fH)-np.log(fLw)))
    for (T,P,g,st) in HINGE:
        fLw,fg,tag=C.wgf(T,P,{g:1.0})
        if fLw is None: continue
        fH,_,_,_=E.fH_w(T,P,fg,st,fb,kp)
        out.append(HW.get(g,3.0)*min(0.0,np.log(fH)-np.log(fLw)-0.03))
    return np.array(out)
LO=np.array([]);HI=np.array([])
_lo=[];_hi=[]
for g in GASES: _lo+=[np.log(80),2.5e-10]; _hi+=[np.log(500),4.2e-10]
_lo+=[26.0,5200.0,-0.06]; _hi+=[31.0,6600.0,0.10]
LO=np.array(_lo,float); HI=np.array(_hi,float)
def fitstart(p0,maxn=250):
    p0=np.clip(np.asarray(p0,float),LO,HI)
    assert np.all(np.isfinite(p0))
    return least_squares(resid2,p0,bounds=(LO,HI),diff_step=3e-4,x_scale='jac',
                         xtol=1e-10,ftol=1e-10,max_nfev=maxn)
if __name__=='__main__':
    d=np.load('calib5.npy',allow_pickle=True).item(); p1=np.asarray(d['p'],float)
    print("p1:",p1)
    print("LO:",LO); print("HI:",HI)
    bad=[(i,p1[i],LO[i],HI[i]) for i in range(len(p1)) if not (LO[i]<=p1[i]<=HI[i])]
    print("violations:",bad)
    C._cache.clear(); R=resid2(np.clip(p1,LO,HI))
    print(f"cache {len(C._cache)}; start cost={0.5*np.sum(R**2):.5f} nan={np.sum(~np.isfinite(R))}")
    t0=time.time(); best=None
    for c0 in [-0.059,-0.03,0.0]:
        p=list(np.clip(p1,LO,HI)); p[-1]=c0
        sol=fitstart(p)
        if best is None or sol.cost<best.cost: best=sol
        print(f" c0={c0}: cost={sol.cost:.5f} nfev={sol.nfev} c={sol.x[-1]:.4f} [{time.time()-t0:.0f}s]",flush=True)
    sol=best; p2=sol.x; kp,fb=unpack2(p2)
    print(f"\nBEST cost={sol.cost:.5f}")
    for g in GASES: print(f"  {g}: eps/k={kp[g][0]:8.2f} sigma={kp[g][1]*1e10:.4f} a={kp[g][2]*1e10:.4f}")
    print(f"  fbeta sI A={fb['sI'][0]:.4f} B={fb['sI'][1]:.2f} c={p2[-1]:.4f}")
    print(f"  fbeta_sI(273.15)={E.psat_beta(273.15,'sI',fb):.1f} fbeta_sII(273.15)={E.psat_beta(273.15,'sII',fb):.1f} Pa")
    R=resid2(p2); i=0
    for gas in GASES:
        n=len(C.ANCH[gas]); print(f"   {gas} anchors sumsq={np.sum(R[i:i+n]**2):.4f} max|r|={np.max(np.abs(R[i:i+n])):.3f}"); i+=n
    for gas in GASES:
        if gas in C.PPAB:
            n=sum(1 for T in C.TPP for m in (0,1) if C.cpp(gas,T,m)); print(f"   {gas} reg sumsq={np.sum(R[i:i+n]**2):.4f}"); i+=n
    print(f"   fig5 r={np.round(R[i:i+2],4)}"); i+=2
    print(f"   hinges r={np.round(R[i:],4)}")
    json.dump({'kp':{g:list(map(float,kp[g])) for g in GASES},'fb':{'sI':list(map(float,fb['sI'])),'sII':list(map(float,fb['sII']))}},open('params_v3.json','w'),indent=1)
    np.save('calib6.npy',{'kp':kp,'fb':fb,'p':p2},allow_pickle=True)
    print("saved params_v3.json / calib6.npy")
