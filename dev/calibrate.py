import numpy as np
from scipy.optimize import least_squares, brentq
import engine as E
np.seterr(all='ignore')

ASSIGN='water_first'

# ---------- anchors: gas -> list of (T, P_Pa, weight, struct) ----------
ANCH = {
 'CH4': [(273.15,2.60e6,3,'sI'),(275.15,3.24e6,2,'sI'),(277.59,4.03e6,3,'sI'),
         (280.15,5.50e6,1.5,'sI'),(283.15,7.20e6,2,'sI'),(286.15,10.4e6,1,'sI'),
         (289.15,14.2e6,1,'sI'),(293.15,19.5e6,0.7,'sI'),(298.15,26.5e6,0.4,'sI'),
         (263.15,1.55e6,1,'sI'),(268.15,2.05e6,0.5,'sI')],
 'C2H6':[(273.15,0.55e6,3,'sI'),(275.15,0.70e6,1.5,'sI'),(277.15,0.88e6,1.5,'sI'),
         (279.15,1.10e6,1,'sI'),(281.15,1.37e6,1,'sI'),(283.15,1.72e6,1,'sI'),
         (285.15,2.15e6,1,'sI'),(287.75,3.39e6,2,'sI'),
         (293.15,4.65e6,0.5,'sI'),(298.15,6.60e6,0.5,'sI'),(303.15,9.10e6,0.4,'sI'),
         (263.15,0.29e6,0.7,'sI'),(268.15,0.41e6,0.7,'sI'),(271.65,0.50e6,0.7,'sI')],
 'C3H8':[(273.15,0.172e6,3,'sII'),(274.15,0.215e6,1.5,'sII'),(275.15,0.26e6,1.5,'sII'),
         (276.15,0.33e6,1,'sII'),(277.59,0.42e6,3,'sII'),(278.15,0.47e6,0.7,'sII'),
         (268.15,0.10e6,0.5,'sII'),(271.15,0.145e6,0.5,'sII')],
 'CO2': [(273.15,1.28e6,3,'sI'),(275.15,1.60e6,1.5,'sI'),(277.15,1.99e6,1.5,'sI'),
         (279.15,2.45e6,1,'sI'),(281.15,3.20e6,1,'sI'),(282.65,4.30e6,1,'sI'),
         (283.0,4.50e6,2,'sI'),
         (285.15,5.60e6,0.6,'sI'),(287.15,6.75e6,0.6,'sI'),(289.15,8.05e6,0.5,'sI'),
         (291.15,9.55e6,0.4,'sI'),
         (263.15,0.80e6,0.7,'sI'),(268.15,1.02e6,0.7,'sI'),(253.15,0.45e6,0.4,'sI'),
         (244.5,0.364e6,0.4,'sI'),(238.5,0.280e6,0.4,'sI'),(234.3,0.232e6,0.3,'sI'),
         (226.1,0.160e6,0.3,'sI'),(217.8,0.1013e6,0.3,'sI')],
 'N2':  [(273.15,16.4e6,2,'sII'),(277.15,23.5e6,1,'sII'),(281.15,33.0e6,1,'sII'),
         (285.15,46.0e6,0.7,'sII'),(289.15,62.0e6,0.5,'sII'),(293.15,79.0e6,0.4,'sII'),
         (263.15,9.3e6,0.6,'sII'),(268.15,12.6e6,0.6,'sII')],
}

# ---------- P-P 1972 Langmuir A,B regularization (C in Pa^-1) ----------
PPAB = {'CH4':((3.7237e-3,2708.8),(1.8372e-2,2737.9)),
        'C2H6':(None,(0.6906e-2,3631.6)),
        'CO2':((1.1978e-3,2860.5),(0.8507e-2,3277.9)),
        'N2':((3.8087e-3,2205.5),(1.8420e-2,2301.3))}
TPP = [260,265,270,275,280,285,290,295,300]
WREG = 0.30   # weight per reg point (ln C space)

def cpp(gas,T,m):
    ab=PPAB[gas][m]
    if ab is None: return None
    A,B=ab
    return (A/T)*np.exp(B/T)/101325.0

# ---------- fast equilibrium with bracket around guess ----------
def Peq(T,gas,struct,kparams,fb,Pguess=None,assign=ASSIGN):
    zc=np.zeros(5); zc[E.NAMES.index(gas)]=1.0
    if Pguess:
        lo,hi=Pguess/6,Pguess*6; n=14
    else:
        lo,hi=2e2,4e8; n=60
    z=np.zeros(6); z[E.IW]=0.8; z[:5]=0.2*zc
    z/=z.sum()
    gas_idx=[E.NAMES.index(gas)]
    def resid(lnP):
        P=np.exp(lnP)
        try:
            if T<273.15:
                fi=E.psat_ice(T)*np.exp(E.V_ICE*(P-E.psat_ice(T))/(E.Rg*T))
                y=np.zeros(6); y[E.IW]=min(fi/P,0.5); y[gas_idx[0]]=1-y[E.IW]
                for _ in range(8):
                    vp=E.lnphi(y,T,P,'v',assign)
                    if vp is None: return np.nan
                    yw=min(max(fi/(P*np.exp(vp[E.IW])),1e-14),0.5)
                    if abs(yw-y[E.IW])<1e-13: y[E.IW]=yw; break
                    y[E.IW]=0.5*y[E.IW]+0.5*yw; y[gas_idx[0]]=1-y[E.IW]
                vp=E.lnphi(y,T,P,'v',assign)
                if vp is None: return np.nan
                fgas={gas:y[gas_idx[0]]*np.exp(vp[gas_idx[0]])*P}
                fH,_,_,_=E.fH_w(T,P,fgas,struct,fb,kparams)
                return np.log(fH)-np.log(fi)
            else:
                fl=E.flash_TP(z,T,P,assign,maxiter=200,tol=1e-9)
                if fl is None: return np.nan
                x,y,beta,conv=fl
                xaq,yg=(x,y) if x[E.IW]>=y[E.IW] else (y,x)
                lp=E.lnphi(xaq,T,P,'l',assign); vp=E.lnphi(yg,T,P,'v',assign)
                if lp is None or vp is None: return np.nan
                fLw=xaq[E.IW]*np.exp(lp[E.IW])*P
                fgas={gas:yg[gas_idx[0]]*np.exp(vp[gas_idx[0]])*P}
                fH,_,_,_=E.fH_w(T,P,fgas,struct,fb,kparams)
                return np.log(fH)-np.log(fLw)
        except Exception:
            return np.nan
    lns=np.linspace(np.log(lo),np.log(hi),n)
    vals=np.array([resid(L) for L in lns])
    for k in range(n-1):
        if np.isfinite(vals[k]) and np.isfinite(vals[k+1]) and vals[k]*vals[k+1]<0:
            try: return float(np.exp(brentq(resid,lns[k],lns[k+1],xtol=1e-9)))
            except Exception: pass
    # fallback wide
    if Pguess:
        return Peq(T,gas,struct,kparams,fb,None,assign)
    return None

def make_kp(gasvals):
    kp=dict(E.KIH)
    for g,(e,s) in gasvals.items(): kp[g]=(e,s,E.KIH[g][2])
    return kp

def fit_gas(gas,kp0,fb):
    """fit (ln eps, sigma) for one gas; others fixed at kp0"""
    anchors=ANCH[gas]
    regs=[]
    if gas in PPAB:
        for T in TPP:
            for m in (0,1):
                c=cpp(gas,T,m)
                if c: regs.append((T,m,c))
    st_used=anchors[0][3]
    def resid(p):
        gv={gas:(np.exp(p[0]),p[1])}
        kp=make_kp_from(gv,kp0)
        out=[]
        Pg=None
        for (T,P,w,st) in sorted(anchors):
            Pm=Peq(T,gas,st,kp,fb,Pguess=P)
            out.append(w*(np.log(Pm/P) if Pm else 3.0))
        for (T,m,ct) in regs:
            C=E.langmuir(T,E.CAGES[st_used][m][0],E.CAGES[st_used][m][1],*kp[gas])
            out.append(WREG*np.log(C/ct))
        return np.array(out)
    def make_kp_from(gv,kp0):
        kp=dict(kp0)
        for g,(e,s) in gv.items(): kp[g]=(e,s,kp0[g][2])
        return kp
    p0=np.array([np.log(kp0[gas][0]),kp0[gas][1]])
    sol=least_squares(resid,p0,method='lm',xtol=1e-4,ftol=1e-4,max_nfev=60)
    return np.exp(sol.x[0]),sol.x[1],sol

def fit_fbeta(kp,fb0):
    """fit A,B for both structures using all anchors"""
    def resid(p):
        fb={'sI':(p[0],p[1]),'sII':(p[2],p[3])}
        out=[]
        for gas,anchors in ANCH.items():
            for (T,P,w,st) in sorted(anchors):
                Pm=Peq(T,gas,st,kp,fb,Pguess=P)
                out.append(w*(np.log(Pm/P) if Pm else 3.0))
        return np.array(out)
    sol=least_squares(resid,np.array(fb0['sI']+fb0['sII']),method='lm',xtol=1e-3,ftol=1e-3,max_nfev=25)
    return {'sI':(sol.x[0],sol.x[1]),'sII':(sol.x[2],sol.x[3])},sol

if __name__=='__main__':
    import sys,time
    kp0=dict(E.KIH)
    kp0['C3H8']=(220.0,3.40e-10,0.6502e-10)
    fb={'sI':E.FBETA['sI'],'sII':E.FBETA['sII']}
    t0=time.time()
    # stage A: per-gas fits with published fbeta
    for gas in ['CH4','CO2','N2','C2H6','C3H8']:
        e,s,sol=fit_gas(gas,kp0,fb)
        kp0[gas]=(e,s,kp0[gas][2])
        print(f"{gas}: eps/k={e:.2f} sigma={s*1e10:.4f} A   cost={sol.cost:.4f}  [{time.time()-t0:.0f}s]",flush=True)
    # stage B: fbeta
    fb,solb=fit_fbeta(kp0,fb)
    print(f"fbeta: sI A={fb['sI'][0]:.4f} B={fb['sI'][1]:.1f}  sII A={fb['sII'][0]:.4f} B={fb['sII'][1]:.1f} cost={solb.cost:.3f} [{time.time()-t0:.0f}s]",flush=True)
    # stage A2: per-gas again with new fbeta
    for gas in ['CH4','CO2','N2','C2H6','C3H8']:
        e,s,sol=fit_gas(gas,kp0,fb)
        kp0[gas]=(e,s,kp0[gas][2])
        print(f"{gas}(2): eps/k={e:.2f} sigma={s*1e10:.4f} A cost={sol.cost:.4f} [{time.time()-t0:.0f}s]",flush=True)
    fb,solb=fit_fbeta(kp0,fb)
    print(f"fbeta(2): sI A={fb['sI'][0]:.4f} B={fb['sI'][1]:.1f} sII A={fb['sII'][0]:.4f} B={fb['sII'][1]:.1f} cost={solb.cost:.3f} [{time.time()-t0:.0f}s]",flush=True)
    np.save('calib_result.npy',{'kp':kp0,'fb':fb},allow_pickle=True)
    print("saved calib_result.npy")
