import numpy as np, time, json, sys
from scipy.optimize import least_squares
import engine as E, calib4 as C4, calibrate as CB
np.seterr(all='ignore')
GASES=['CH4','C2H6','C3H8','N2','CO2']
ANCH=CB.ANCH; PPAB=CB.PPAB; TPP=CB.TPP; cpp=CB.cpp
WREG={'CH4':0.25,'CO2':0.25,'C2H6':0.25,'N2':0.10,'C3H8':0.0}
FLOOR=0.015
# paper Eqs 16-17 (Martin & Peters)
MP={'sI':(28.4,5800.0),'sII':(28.8,6000.0)}

def unpack5(p,mode):
    kp={}; i=0
    for g in GASES:
        kp[g]=(np.exp(p[i]),p[i+1],E.KIH[g][2]); i+=2
    if mode=='mp': fb={'sI':MP['sI'],'sII':MP['sII']}
    else: fb={'sI':(p[i],p[i+1]),'sII':(p[i+2],p[i+3])}
    return kp,fb

def r_of(T,P,gf,st,kp,fb):
    fLw,fg,tag=C4.wgf(T,P,gf)
    if fLw is None or fg is None: return None
    fH,_,_,_=E.fH_w(T,P,fg,st,fb,kp)
    v=np.log(fH)-np.log(fLw)
    return v if np.isfinite(v) else None

def resid5(p,mode):
    kp,fb=unpack5(p,mode)
    out=[]
    for gas in GASES:
        gf={gas:1.0}
        for (T,P,w,st) in sorted(ANCH[gas]):
            r0=r_of(T,P,gf,st,kp,fb)
            if r0 is None: out.append(w*2.0); continue
            rp=r_of(T,P*1.03,gf,st,kp,fb); rm=r_of(T,P/1.03,gf,st,kp,fb)
            if rp is not None and rm is not None:
                slope=(rp-rm)/(2*np.log(1.03))
                sabs=max(abs(slope),FLOOR) if slope<0 else FLOOR*3
            else: sabs=1.0
            out.append(w*r0/sabs)
            # per-anchor structure hinge: selected structure must win at this state
            oth='sII' if st=='sI' else 'sI'
            ro=r_of(T,P,gf,oth,kp,fb)
            if ro is not None:
                out.append(10.0*min(0.0,(ro-r0)-0.02))
    for gas in GASES:
        if gas in PPAB:
            st='sII' if gas=='N2' else 'sI'
            for T in TPP:
                for m in (0,1):
                    ct=cpp(gas,T,m)
                    if ct:
                        Cc=E.langmuir(T,E.CAGES[st][m][0],E.CAGES[st][m][1],*kp[gas])
                        out.append(WREG[gas]*np.log(Cc/ct))
    gf5={'CH4':0.9995,'C3H8':0.0005}
    for st in ('sI','sII'):
        r=r_of(277.6,4.0e6,gf5,st,kp,fb)
        out.append(1.5*(r if r is not None else 3.0))
    H=[('CH4',273.15,2.60e6,'sI','sII'),('CH4',283.15,7.20e6,'sI','sII'),('CH4',263.15,1.55e6,'sI','sII'),
       ('C3H8',277.59,0.42e6,'sII','sI'),('C3H8',273.15,0.172e6,'sII','sI'),('C3H8',268.15,0.10e6,'sII','sI'),
       ('C2H6',273.15,0.55e6,'sI','sII'),('C2H6',293.15,4.65e6,'sI','sII'),
       ('CO2',273.15,1.28e6,'sI','sII'),('CO2',244.50,0.364e6,'sI','sII'),
       ('N2',273.15,16.4e6,'sII','sI'),('N2',281.15,33.0e6,'sII','sI'),('N2',263.15,9.3e6,'sII','sI'),('N2',293.15,79.0e6,'sII','sI')]
    for (gas,T,P,sel,oth) in H:
        rs=r_of(T,P,{gas:1.0},sel,kp,fb); ro=r_of(T,P,{gas:1.0},oth,kp,fb)
        if rs is None or ro is None: continue
        out.append(3.0*min(0.0,(ro-rs)-0.02))
    return np.array(out)

LO=np.array([np.log(90),2.6e-10, np.log(90),2.6e-10, np.log(90),2.8e-10, np.log(90),2.6e-10, np.log(90),2.5e-10])
HI=np.array([np.log(420),4.0e-10, np.log(420),4.0e-10, np.log(320),3.65e-10, np.log(420),4.0e-10, np.log(420),3.5e-10])
LOf=np.concatenate([LO,[27.0,5400.0,27.0,5400.0]]); HIf=np.concatenate([HI,[30.0,6200.0,30.0,6400.0]])

def fitrun(p0,mode,maxn=600):
    lo,hi=(LO,HI) if mode=='mp' else (LOf,HIf)
    p0=np.clip(np.asarray(p0,float),lo,hi)
    return least_squares(lambda p: resid5(p,mode),p0,bounds=(lo,hi),diff_step=3e-4,xtol=1e-12,ftol=1e-12,max_nfev=maxn)

if __name__=='__main__':
    mode=sys.argv[1] if len(sys.argv)>1 else 'mp'
    t0=time.time()
    d=np.load('calib7.npy',allow_pickle=True).item(); kp4=d['kp']
    d2=np.load('calibfinal.npy',allow_pickle=True).item(); kpf=d2['kp']
    def pk(kp): return np.array([v for g in GASES for v in (np.log(kp[g][0]),kp[g][1])])
    starts=[pk(kp4),pk(kpf)]
    kp5=dict(kpf); kp5['CH4']=(230.0,3.42e-10,kpf['CH4'][2]); starts.append(pk(kp5))  # M-P level shift guess
    kp6=dict(kp4)
    for g in GASES: kp6[g]=(kp4[g][0]*1.25,kp4[g][1]*1.03,kp4[g][2])
    starts.append(pk(kp6))
    best=None
    for j,p0 in enumerate(starts):
        if mode=='free': p0=np.concatenate([p0,[28.4,5800.,28.8,6000.]])
        sol=fitrun(p0,mode)
        print(f" start{j}: cost={sol.cost:.5f} nfev={sol.nfev} [{time.time()-t0:.0f}s]",flush=True)
        if best is None or sol.cost<best.cost: best=sol
    sol=best; kp,fb=unpack5(sol.x,mode)
    print(f"\nBEST[{mode}] cost={sol.cost:.5f}")
    for g in GASES: print(f"  {g}: eps/k={kp[g][0]:8.2f} sigma={kp[g][1]*1e10:.4f}")
    print(f"  fbeta sI=({fb['sI'][0]:.4f},{fb['sI'][1]:.1f}) sII=({fb['sII'][0]:.4f},{fb['sII'][1]:.1f})")
    print(f"  fbsI(273)={E.psat_beta(273.15,'sI',fb):.1f} fbsII(273)={E.psat_beta(273.15,'sII',fb):.1f}")
    json.dump({'kp':{g:[float(kp[g][0]),float(kp[g][1]),float(kp[g][2])] for g in GASES},
               'fb':{'sI':[float(fb['sI'][0]),float(fb['sI'][1])],'sII':[float(fb['sII'][0]),float(fb['sII'][1])]}},
              open(f'params_v5_{mode}.json','w'),indent=1)
    np.save(f'calib8_{mode}.npy',{'kp':kp,'fb':fb,'p':sol.x},allow_pickle=True)
    print(f"saved [{time.time()-t0:.0f}s]")
