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
WOVER={('CH4',263.15):3.,('CH4',298.15):3.,
      ('C2H6',287.75):3.,('C2H6',293.15):3.,('C2H6',298.15):3.,('C2H6',303.15):3.,
      ('CO2',282.65):2.5,('CO2',283.0):2.5,('CO2',285.15):2.5,('CO2',287.15):2.5,('CO2',289.15):2.5,('CO2',291.15):2.5,
      ('N2',281.15):3.,('N2',285.15):3.,('N2',289.15):3.,('N2',293.15):3.}

def unpack9(p):
    kp={}; i=0
    for g in GASES:
        kp[g]=(np.exp(p[i]),p[i+1],p[i+2]); i+=3
    return kp,{'sI':(p[i],p[i+1],p[i+2]),'sII':(p[i+3],p[i+4],p[i+5])}

def r_of(T,P,gf,st,kp,fb):
    fLw,fg,tag=C4.wgf(T,P,gf)
    if fLw is None or fg is None: return None
    fH,_,_,_=E.fH_w(T,P,fg,st,fb,kp)
    v=np.log(fH)-np.log(fLw)
    return v if np.isfinite(v) else None

def resid9(p,scaled=True):
    kp,fb=unpack9(p)
    out=[]
    for gas in GASES:
        gf={gas:1.0}
        for (T,P,w,st) in sorted(ANCH[gas]):
            w=WOVER.get((gas,T),w)
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
    for (gas,T,P,st,thl,ths) in OCC:
        fLw,fg,tag=C4.wgf(T,P,{gas:1.0})
        if fLw is None or fg is None: out.extend([2.0,2.0]); continue
        fH,th,S,fbeta=E.fH_w(T,P,fg,st,fb,kp)
        out.append((th[(gas,1)]-thl)/0.015)
        out.append((th[(gas,0)]-ths)/0.015)
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

AB={'CH4':(0.15e-10,0.55e-10),'C2H6':(0.30e-10,0.75e-10),'C3H8':(0.45e-10,0.85e-10),
    'N2':(0.15e-10,0.55e-10),'CO2':(0.45e-10,0.85e-10)}
lo=[];hi=[]
for g in GASES:
    lo+=[np.log(90),2.6e-10,AB[g][0]]; hi+=[np.log(340),4.0e-10,AB[g][1]]
lo+=[26.0,5400.0,-4e4, 26.0,5000.0,-4e4]; hi+=[31.0,6400.0,4e4, 31.0,6450.0,4e4]
LO=np.array(lo); HI=np.array(hi)

def pk(kp,f1,f2):
    a=[]
    for g in GASES: a+=[np.log(kp[g][0]),kp[g][1],kp[g][2]]
    a+=[f1[0],f1[1],(f1[2] if len(f1)>2 else 0.0),f2[0],f2[1],(f2[2] if len(f2)>2 else 0.0)]
    return np.array(a,float)

if __name__=='__main__':
    t0=time.time()
    d8=np.load('calib11.npy',allow_pickle=True).item()
    d7=np.load('calib10.npy',allow_pickle=True).item()
    kpa=dict(d8['kp']); kpa['C2H6']=(kpa['C2H6'][0],kpa['C2H6'][1],0.5695e-10)
    starts=[pk(d8['kp'],d8['fb']['sI'],d8['fb']['sII']),
            pk(d7['kp'],d7['fb']['sI'],d7['fb']['sII']),
            pk(kpa,d8['fb']['sI'],d8['fb']['sII'])]
    best=None
    for j,p0 in enumerate(starts):
        s1=least_squares(lambda p:resid9(p,False),np.clip(p0,LO,HI),bounds=(LO,HI),diff_step=1e-3,xtol=1e-13,ftol=1e-13,max_nfev=300)
        s2=least_squares(lambda p:resid9(p,True),s1.x,bounds=(LO,HI),diff_step=1e-3,xtol=1e-13,ftol=1e-13,max_nfev=600)
        print(f" start{j}: r-space={s1.cost:.4f} scaled={s2.cost:.5f} [{time.time()-t0:.0f}s]",flush=True)
        if best is None or s2.cost<best.cost: best=s2
    sol=best; kp,fb=unpack9(sol.x)
    print(f"\nBEST cost={sol.cost:.5f}")
    for g in GASES: print(f"  {g}: eps/k={kp[g][0]:8.2f} sigma={kp[g][1]*1e10:.4f} a={kp[g][2]*1e10:.4f}")
    print(f"  fbeta sI=({fb['sI'][0]:.4f},{fb['sI'][1]:.1f},{fb['sI'][2]:.1f}) sII=({fb['sII'][0]:.4f},{fb['sII'][1]:.1f},{fb['sII'][2]:.1f})")
    json.dump({'kp':{g:[float(kp[g][0]),float(kp[g][1]),float(kp[g][2])] for g in GASES},
               'fb':{'sI':[float(v) for v in fb['sI']],'sII':[float(v) for v in fb['sII']]}},
              open('params_v9.json','w'),indent=1)
    np.save('calib12.npy',{'kp':kp,'fb':fb,'p':sol.x},allow_pickle=True)
    print(f"saved [{time.time()-t0:.0f}s]")
