import numpy as np, time, json
from scipy.optimize import least_squares
import engine as E, calib4 as C4, calibrate as CB
np.seterr(all='ignore')
GASES=['CH4','C2H6','C3H8','N2','CO2']
ANCH=CB.ANCH; PPAB=CB.PPAB; TPP=CB.TPP; cpp=CB.cpp
WREG=0.03
CAP=15.0; FLOOR=1.0/CAP
HINGE_W=30.0; MARG=0.02
OCC=[('CH4',273.2,2.62e6,'sI',0.969,0.935),('CH4',276.7,3.70e6,'sI',0.981,0.917),
     ('CO2',273.2,1.28e6,'sI',0.986,0.581)]

def unpack7(p):
    kp={}; i=0
    for g in GASES:
        kp[g]=(np.exp(p[i]),p[i+1],E.KIH[g][2]); i+=2
    return kp,{'sI':(p[i],p[i+1]),'sII':(p[i+2],p[i+3])}

def r_of(T,P,gf,st,kp,fb):
    fLw,fg,tag=C4.wgf(T,P,gf)
    if fLw is None or fg is None: return None
    fH,_,_,_=E.fH_w(T,P,fg,st,fb,kp)
    v=np.log(fH)-np.log(fLw)
    return v if np.isfinite(v) else None

def resid7(p,scaled=True):
    kp,fb=unpack7(p)
    out=[]
    for gas in GASES:
        gf={gas:1.0}
        for (T,P,w,st) in sorted(ANCH[gas]):
            r0=r_of(T,P,gf,st,kp,fb)
            if r0 is None: out.append(w*2.0); continue
            if scaled:
                rp=r_of(T,P*1.03,gf,st,kp,fb); rm=r_of(T,P/1.03,gf,st,kp,fb)
                if rp is not None and rm is not None:
                    slope=(rp-rm)/(2*np.log(1.03))
                    sabs=max(abs(slope),FLOOR) if slope<0 else 3*FLOOR
                else: sabs=1.0
                out.append(w*r0/sabs)
            else:
                out.append(w*r0)
            oth='sII' if st=='sI' else 'sI'
            ro=r_of(T,P,gf,oth,kp,fb)
            if ro is not None:
                out.append(HINGE_W*min(0.0,(ro-r0)-MARG))
    # occupancies (paper Tables 6-8)
    for (gas,T,P,st,thl,ths) in OCC:
        fLw,fg,tag=C4.wgf(T,P,{gas:1.0})
        if fLw is None or fg is None: out.extend([2.0,2.0]); continue
        fH,th,S,fbeta=E.fH_w(T,P,fg,st,fb,kp)
        out.append((th[(gas,1)]-thl)/0.015)
        out.append((th[(gas,0)]-ths)/0.015)
    # weak P-P regularization
    for gas in GASES:
        if gas in PPAB:
            st='sII' if gas=='N2' else 'sI'
            for T in TPP:
                for m in (0,1):
                    ct=cpp(gas,T,m)
                    if ct:
                        Cc=E.langmuir(T,E.CAGES[st][m][0],E.CAGES[st][m][1],*kp[gas])
                        out.append(WREG*np.log(Cc/ct))
    gf5={'CH4':0.9995,'C3H8':0.0005}
    for st in ('sI','sII'):
        r=r_of(277.6,4.0e6,gf5,st,kp,fb)
        out.append(1.5*(r if r is not None else 3.0))
    return np.array(out)

LO=np.array([np.log(90),2.8e-10, np.log(90),2.8e-10, np.log(90),3.0e-10, np.log(90),2.7e-10, np.log(90),2.6e-10, 27.0,5400.0,27.0,5000.0])
HI=np.array([np.log(320),3.8e-10, np.log(320),3.9e-10, np.log(300),3.65e-10, np.log(320),3.7e-10, np.log(320),3.4e-10, 30.0,6250.0,30.0,6450.0])

def pk(kp,A1,B1,A2,B2):
    return np.array([v for g in GASES for v in (np.log(kp[g][0]),kp[g][1])]+[A1,B1,A2,B2])

if __name__=='__main__':
    t0=time.time()
    df=np.load('calibfinal.npy',allow_pickle=True).item()
    d9=np.load('calib9.npy',allow_pickle=True).item()
    d4=np.load('calib7.npy',allow_pickle=True).item()
    starts=[pk(df['kp'],df['fb']['sI'][0],df['fb']['sI'][1],df['fb']['sII'][0],df['fb']['sII'][1]),
            pk(d9['kp'],d9['fb']['sI'][0],d9['fb']['sI'][1],d9['fb']['sII'][0],d9['fb']['sII'][1]),
            pk(d4['kp'],d4['fb']['sI'][0],d4['fb']['sI'][1],d4['fb']['sII'][0],d4['fb']['sII'][1])]
    kpB=dict(df['kp']); kpB['CH4']=(150.0,3.14e-10,kpB['CH4'][2]); kpB['CO2']=(165.0,2.94e-10,kpB['CO2'][2])
    starts.append(pk(kpB,28.69,5933.,28.63,5933.))
    best=None
    for j,p0 in enumerate(starts):
        s1=least_squares(lambda p:resid7(p,False),np.clip(p0,LO,HI),bounds=(LO,HI),diff_step=1e-3,xtol=1e-12,ftol=1e-12,max_nfev=250)
        s2=least_squares(lambda p:resid7(p,True),s1.x,bounds=(LO,HI),diff_step=1e-3,xtol=1e-12,ftol=1e-12,max_nfev=500)
        print(f" start{j}: r-space={s1.cost:.4f} scaled={s2.cost:.5f} [{time.time()-t0:.0f}s]",flush=True)
        if best is None or s2.cost<best.cost: best=s2
    sol=best; kp,fb=unpack7(sol.x)
    print(f"\nBEST cost={sol.cost:.5f}")
    for g in GASES: print(f"  {g}: eps/k={kp[g][0]:8.2f} sigma={kp[g][1]*1e10:.4f}")
    print(f"  fbeta sI=({fb['sI'][0]:.4f},{fb['sI'][1]:.1f}) sII=({fb['sII'][0]:.4f},{fb['sII'][1]:.1f})")
    print(f"  fbsI(273.15)={E.psat_beta(273.15,'sI',fb):.1f} fbsII(273.15)={E.psat_beta(273.15,'sII',fb):.1f}")
    json.dump({'kp':{g:[float(kp[g][0]),float(kp[g][1]),float(kp[g][2])] for g in GASES},
               'fb':{'sI':[float(fb['sI'][0]),float(fb['sI'][1])],'sII':[float(fb['sII'][0]),float(fb['sII'][1])]}},
              open('params_v7.json','w'),indent=1)
    np.save('calib10.npy',{'kp':kp,'fb':fb,'p':sol.x},allow_pickle=True)
    print(f"saved [{time.time()-t0:.0f}s]")
