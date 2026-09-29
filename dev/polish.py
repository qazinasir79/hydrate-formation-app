import numpy as np, time, json
from scipy.optimize import least_squares
import engine as E, calib4 as C4, calibrate as CB
np.seterr(all='ignore')
GASES=['CH4','C2H6','C3H8','N2','CO2']
ANCH=CB.ANCH
CAP=15.0; FLOOR=1.0/CAP
HINGE_W=30.0; MARG=0.02
WOVER={('CH4',263.15):3.,('CH4',298.15):3.,
      ('C2H6',287.75):3.,('C2H6',293.15):3.,('C2H6',298.15):3.,('C2H6',303.15):3.,
      ('CO2',282.65):2.5,('CO2',283.0):2.5,('CO2',285.15):2.5,('CO2',287.15):2.5,('CO2',289.15):2.5,('CO2',291.15):2.5,
      ('N2',281.15):3.,('N2',285.15):3.,('N2',289.15):3.,('N2',293.15):3.}
OCC=[('CH4',273.2,2.62e6,'sI',0.969,0.935),('CH4',276.7,3.70e6,'sI',0.981,0.917),('CO2',273.2,1.28e6,'sI',0.986,0.581)]
AB={'CH4':(0.25e-10,0.45e-10),'C2H6':(0.40e-10,0.62e-10),'C3H8':(0.55e-10,0.75e-10),
    'N2':(0.25e-10,0.45e-10),'CO2':(0.55e-10,0.78e-10)}
SB={'CH4':(2.8e-10,3.8e-10),'C2H6':(2.8e-10,3.9e-10),'C3H8':(3.0e-10,3.65e-10),
    'N2':(2.7e-10,3.7e-10),'CO2':(2.6e-10,3.4e-10)}

def r_of(T,P,gf,st,kp,fb):
    fLw,fg,tag=C4.wgf(T,P,gf)
    if fLw is None or fg is None: return None
    fH,_,_,_=E.fH_w(T,P,fg,st,fb,kp)
    v=np.log(fH)-np.log(fLw)
    return v if np.isfinite(v) else None

def gas_resid(p3,gas,kp,fb):
    kp=dict(kp); kp[gas]=(np.exp(p3[0]),p3[1],p3[2])
    out=[]
    gf={gas:1.0}
    for (T,P,w,st) in sorted(ANCH[gas]):
        w=WOVER.get((gas,T),w)
        r0=r_of(T,P,gf,st,kp,fb)
        if r0 is None: out.append(w*2.0); continue
        rp=r_of(T,P*1.03,gf,st,kp,fb); rm=r_of(T,P/1.03,gf,st,kp,fb)
        if rp is not None and rm is not None:
            slope=(rp-rm)/(2*np.log(1.03))
            sabs=max(abs(slope),FLOOR) if slope<0 else 3*FLOOR
        else: sabs=1.0
        out.append(w*r0/sabs)
        oth='sII' if st=='sI' else 'sI'
        ro=r_of(T,P,gf,oth,kp,fb)
        if ro is not None: out.append(HINGE_W*min(0.0,(ro-r0)-MARG))
    for (g2,T,P,st,thl,ths) in OCC:
        if g2!=gas: continue
        fLw,fg,tag=C4.wgf(T,P,{gas:1.0})
        if fLw is None or fg is None: out.extend([2.0,2.0]); continue
        fH,th,S,fbeta=E.fH_w(T,P,fg,st,fb,kp)
        out.append((th[(gas,1)]-thl)/0.015); out.append((th[(gas,0)]-ths)/0.015)
    if gas in ('CH4','C3H8'):
        gf5={'CH4':0.9995,'C3H8':0.0005}
        for st in ('sI','sII'):
            r=r_of(277.6,4.0e6,gf5,st,kp,fb)
            out.append(1.5*(r if r is not None else 3.0))
    if gas=='CH4':  # keep CH4/C3H8 cross-transition sane: fig6 lens check states
        gf6={'CH4':0.753,'C2H6':0.247}
        r1=r_of(277.6,3.6e6,gf6,'sI',kp,fb); r2=r_of(277.6,3.6e6,gf6,'sII',kp,fb)
        if r1 is not None and r2 is not None: out.append(1.0*(r1-r2))
    return np.array(out)

if __name__=='__main__':
    t0=time.time()
    d=np.load('calib10.npy',allow_pickle=True).item()
    kp=dict(d['kp']); fb=d['fb']
    for gas in ['C2H6','N2','CO2','CH4','C3H8']:
        lo=[np.log(90),SB[gas][0],AB[gas][0]]; hi=[np.log(340),SB[gas][1],AB[gas][1]]
        p0=np.array([np.log(kp[gas][0]),kp[gas][1],kp[gas][2]])
        if gas=='C2H6':
            for pa in (0.5695e-10,):
                s=least_squares(lambda p:gas_resid(p,gas,kp,fb),np.clip([np.log(kp[gas][0]),kp[gas][1],pa],lo,hi),bounds=(lo,hi),diff_step=1e-3,xtol=1e-13,ftol=1e-13,max_nfev=200)
                print(f"  {gas}(a-paper) cost={s.cost:.4f} -> eps={np.exp(s.x[0]):.2f} sig={s.x[1]*1e10:.4f} a={s.x[2]*1e10:.4f}",flush=True)
        s=least_squares(lambda p:gas_resid(p,gas,kp,fb),np.clip(p0,lo,hi),bounds=(lo,hi),diff_step=1e-3,xtol=1e-13,ftol=1e-13,max_nfev=250)
        c0=0.5*np.sum(gas_resid(np.clip(p0,lo,hi),gas,kp,fb)**2)
        print(f"{gas}: cost {c0:.4f} -> {s.cost:.4f}: eps/k={np.exp(s.x[0]):8.2f} sigma={s.x[1]*1e10:.4f} a={s.x[2]*1e10:.4f} [{time.time()-t0:.0f}s]",flush=True)
        if s.cost<c0: kp[gas]=(float(np.exp(s.x[0])),float(s.x[1]),float(s.x[2]))
    json.dump({'kp':{g:[float(kp[g][0]),float(kp[g][1]),float(kp[g][2])] for g in GASES},
               'fb':{'sI':[float(v) for v in fb['sI']],'sII':[float(v) for v in fb['sII']]}},
              open('params_v12.json','w'),indent=1)
    np.save('calib15.npy',{'kp':kp,'fb':fb},allow_pickle=True)
    print("saved calib15")
