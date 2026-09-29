import numpy as np, time, json
from scipy.optimize import least_squares
import engine as E, calib4 as C4, calibrate as CB
np.seterr(all='ignore')
GASES=['CH4','C2H6','C3H8','N2','CO2']
ANCH=CB.ANCH; PPAB=CB.PPAB; TPP=CB.TPP; cpp=CB.cpp
WREG={'CH4':0.30,'CO2':0.30,'C2H6':0.30,'N2':0.12,'C3H8':0.0}
SLOPE_FLOOR=0.008

def unpack2(p):
    kp={}; i=0
    for g in GASES:
        kp[g]=(np.exp(p[i]),p[i+1],E.KIH[g][2]); i+=2
    A,B,c=p[i],p[i+1],p[i+2]
    return kp,{'sI':(A,B),'sII':(A+c,B)}

def r_of(T,P,gf,st,kp,fb):
    fLw,fg,tag=C4.wgf(T,P,gf)
    if fLw is None or fg is None: return None
    fH,_,_,_=E.fH_w(T,P,fg,st,fb,kp)
    v=np.log(fH)-np.log(fLw)
    return v if np.isfinite(v) else None

def resid4(p):
    kp,fb=unpack2(p)
    out=[]
    for gas in GASES:
        gf={gas:1.0}
        for (T,P,w,st) in sorted(ANCH[gas]):
            r0=r_of(T,P,gf,st,kp,fb)
            if r0 is None: out.append(w*3.0/SLOPE_FLOOR*0.05); continue
            rp=r_of(T,P*1.02,gf,st,kp,fb); rm=r_of(T,P/1.02,gf,st,kp,fb)
            if rp is not None and rm is not None:
                slope=abs((rp-rm)/(2*np.log(1.02)))
            else: slope=1.0
            slope=max(slope,SLOPE_FLOOR)
            out.append(w*r0/slope)
    # P-P regularization (r space)
    for gas in GASES:
        if gas in PPAB:
            st='sII' if gas=='N2' else 'sI'
            for T in TPP:
                for m in (0,1):
                    ct=cpp(gas,T,m)
                    if ct:
                        Cc=E.langmuir(T,E.CAGES[st][m][0],E.CAGES[st][m][1],*kp[gas])
                        out.append(WREG[gas]*np.log(Cc/ct))
    # Fig5 transition: both structures at equilibrium at 277.6 K, 4.0 MPa, y_CH4=0.9995
    gf5={'CH4':0.9995,'C3H8':0.0005}
    for st in ('sI','sII'):
        r=r_of(277.6,4.0e6,gf5,st,kp,fb)
        out.append(1.5*(r if r is not None else 3.0))
    # pairwise structure hinges: r_other - r_sel >= margin at key states
    H=[('CH4',273.15,2.60e6,'sI','sII',0.02),('CH4',283.15,7.20e6,'sI','sII',0.02),
       ('CH4',263.15,1.55e6,'sI','sII',0.02),
       ('C3H8',277.59,0.42e6,'sII','sI',0.02),('C3H8',273.15,0.172e6,'sII','sI',0.02),
       ('C3H8',268.15,0.10e6,'sII','sI',0.02),
       ('C2H6',273.15,0.55e6,'sI','sII',0.02),('C2H6',293.15,4.65e6,'sI','sII',0.02),
       ('CO2',273.15,1.28e6,'sI','sII',0.02),('CO2',244.50,0.364e6,'sI','sII',0.02),
       ('N2',273.15,16.4e6,'sII','sI',0.02),('N2',281.15,33.0e6,'sII','sI',0.02),
       ('N2',263.15,9.3e6,'sII','sI',0.02),('N2',293.15,79.0e6,'sII','sI',0.02)]
    for (gas,T,P,sel,oth,marg) in H:
        rs=r_of(T,P,{gas:1.0},sel,kp,fb); ro=r_of(T,P,{gas:1.0},oth,kp,fb)
        if rs is None or ro is None: continue
        out.append(3.0*min(0.0,(ro-rs)-marg))
    return np.array(out)

LO=np.array([np.log(90)]*1+[2.6e-10]+[np.log(90),2.6e-10]+[np.log(90),2.8e-10]+[np.log(90),2.6e-10]+[np.log(90),2.5e-10]+[26.0,5200.0,-0.25])
HI=np.array([np.log(400),4.0e-10]*2+[np.log(400),4.2e-10]+[np.log(400),4.0e-10]*2+[31.5,6600.0,0.10])

def fitstart(p0,maxn=500):
    p0=np.clip(np.asarray(p0,float),LO,HI)
    return least_squares(resid4,p0,bounds=(LO,HI),diff_step=3e-4,xtol=1e-12,ftol=1e-12,max_nfev=maxn)

if __name__=='__main__':
    t0=time.time()
    d=np.load('calibfinal.npy',allow_pickle=True).item()
    kp0,fb0=d['kp'],d['fb']
    def pack(kp,fb,c):
        p=[]
        for g in GASES: p+=[float(np.log(kp[g][0])),float(kp[g][1])]
        p+=[float(fb['sI'][0]),float(fb['sI'][1]),float(c)]
        return np.array(p)
    p_base=pack(kp0,fb0,fb0['sII'][0]-fb0['sI'][0])
    R=resid4(np.clip(p_base,LO,HI)); print(f"start cost={0.5*np.sum(R**2):.4f} [{time.time()-t0:.0f}s]")
    starts=[]
    for c0 in [-0.06,-0.12,-0.18,-0.22]:
        p=p_base.copy(); p[-1]=c0; starts.append(p)
    for de in [1.03,1.06]:   # N2 eps bumps
        p=p_base.copy(); iN2=2*GASES.index('N2'); p[iN2]=np.log(np.exp(p[iN2])*de); p[-1]=-0.14; starts.append(p)
    best=None
    for j,p in enumerate(starts):
        sol=fitstart(p)
        c=0.5*sol.cost*2
        print(f" start{j}: cost={sol.cost:.5f} c={sol.x[-1]:+.4f} nfev={sol.nfev} [{time.time()-t0:.0f}s]",flush=True)
        if best is None or sol.cost<best.cost: best=sol
    sol=best; p1=sol.x; kp,fb=unpack2(p1)
    print(f"\nBEST cost={sol.cost:.5f}")
    for g in GASES: print(f"  {g}: eps/k={kp[g][0]:8.2f} sigma={kp[g][1]*1e10:.4f}")
    print(f"  fbeta sI A={fb['sI'][0]:.4f} B={fb['sI'][1]:.2f} c={p1[-1]:+.4f}  fbsI(273)={E.psat_beta(273.15,'sI',fb):.1f} fbsII(273)={E.psat_beta(273.15,'sII',fb):.1f}")
    json.dump({'kp':{g:[float(kp[g][0]),float(kp[g][1]),float(kp[g][2])] for g in GASES},
               'fb':{'sI':[float(fb['sI'][0]),float(fb['sI'][1])],'sII':[float(fb['sII'][0]),float(fb['sII'][1])]}},
              open('params_v4.json','w'),indent=1)
    np.save('calib7.npy',{'kp':kp,'fb':fb,'p':p1},allow_pickle=True)
    print(f"saved [{time.time()-t0:.0f}s]")
