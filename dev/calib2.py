import numpy as np, time
from scipy.optimize import least_squares, brentq
import engine as E
np.seterr(all='ignore')
ASSIGN='water_first'
from calibrate import ANCH, PPAB, TPP, WREG, cpp

def flash2(z,T,P,assign=ASSIGN):
    """robust two-phase flash: Wilson init, fallback to forced split init"""
    fl=E.flash_TP(z,T,P,assign,maxiter=250,tol=1e-10)
    if fl is not None:
        x,y,beta,conv=fl
        K=y/x
        if not (np.all(np.abs(K-1)<1e-6)):
            return fl
    # forced split init: aqueous x, organic y
    K=np.zeros(6)
    K[E.IW]=0.002
    for i in range(5): K[i]=800.0
    def rr(b): return np.sum(z*(K-1)/(1+b*(K-1)))
    try:
        if rr(0)*rr(1)<0: beta=brentq(rr,1e-10,1-1e-10)
        else: beta=0.3
    except Exception: beta=0.3
    x=z/(1+beta*(K-1)); x=np.maximum(x,1e-15); x/=x.sum()
    y=K*x; y/=y.sum()
    for it in range(400):
        lp=E.lnphi(x,T,P,'l',assign); vp=E.lnphi(y,T,P,'v',assign)
        if lp is None or vp is None: return None
        Knew=np.exp(np.clip(lp-vp,-50,50))
        def rr2(b): return np.sum(z*(Knew-1)/(1+b*(Knew-1)))
        lo,hi=1e-10,1-1e-10
        try:
            if rr2(lo)*rr2(hi)<0: beta=brentq(rr2,lo,hi,xtol=1e-14)
            else: beta=lo if rr2(lo)>0 else hi
        except Exception: beta=0.3
        x=z/(1+beta*(Knew-1)); x=np.maximum(x,1e-16); x/=x.sum()
        y=Knew*x; y=np.maximum(y,1e-16); y/=y.sum()
        err=np.max(np.abs(Knew/K-1)); K=Knew
        if err<1e-11: return x,y,beta,True
    return x,y,beta,False

def water_and_gas_fug(T,P,gas_list,assign=ASSIGN):
    """returns (f_w_liq_or_ice, {gas: fugacity}, phase_tag) at T,P for excess-water feed"""
    if T<273.15:
        fi=E.psat_ice(T)*np.exp(E.V_ICE*(P-E.psat_ice(T))/(E.Rg*T))
        y=np.zeros(6); y[E.IW]=min(fi/P,0.5)
        for g in gas_list: y[E.NAMES.index(g)]=(1-y[E.IW])/len(gas_list)
        # distribute vapor among gases proportional to feed (pure gas case: single)
        for _ in range(10):
            vp=E.lnphi(y,T,P,'v',assign)
            if vp is None: return None,None,'fail'
            yw=min(max(fi/(P*np.exp(vp[E.IW])),1e-14),0.6)
            if abs(yw-y[E.IW])<1e-13: y[E.IW]=yw; break
            y[E.IW]=0.5*y[E.IW]+0.5*yw
            for g in gas_list: y[E.NAMES.index(g)]=(1-y[E.IW])/len(gas_list)
        vp=E.lnphi(y,T,P,'v',assign)
        if vp is None: return None,None,'fail'
        fg={g:y[E.NAMES.index(g)]*np.exp(vp[E.NAMES.index(g)])*P for g in gas_list}
        return fi,fg,'ice'
    z=np.zeros(6); z[E.IW]=0.8
    for g in gas_list: z[E.NAMES.index(g)]=0.2/len(gas_list)
    fl=flash2(z,T,P,assign)
    if fl is None: return None,None,'fail'
    x,y,beta,conv=fl
    if x[E.IW]>=y[E.IW]: xaq,yg=x,y
    else: xaq,yg=y,x
    lp=E.lnphi(xaq,T,P,'l',assign); vp=E.lnphi(yg,T,P,'v',assign)
    if lp is None or vp is None: return None,None,'fail'
    fLw=xaq[E.IW]*np.exp(lp[E.IW])*P
    fg={g:yg[E.NAMES.index(g)]*np.exp(vp[E.NAMES.index(g)])*P for g in gas_list}
    tag='lw' if yg[E.IW]>0.5 else ('lc' if xaq[E.IW]<0.9 else 'lw')
    return fLw,fg,tag

def anchor_resid(T,P,gas,struct,kp,fb):
    fLw,fg,tag=water_and_gas_fug(T,P,[gas])
    if fLw is None: return np.nan,tag
    fH,_,_,_=E.fH_w(T,P,fg,struct,fb,kp)
    return np.log(fH)-np.log(fLw),tag

def fit_gas(gas,kp0,fb,verbose=False):
    anchors=sorted(ANCH[gas]); st=anchors[0][3]
    regs=[]
    if gas in PPAB:
        for T in TPP:
            for m in (0,1):
                c=cpp(gas,T,m)
                if c: regs.append((T,m,c))
    def resid(p):
        kp=dict(kp0); kp[gas]=(np.exp(p[0]),p[1],kp0[gas][2])
        out=[]
        for (T,P,w,s) in anchors:
            r,_=anchor_resid(T,P,gas,s,kp,fb)
            out.append(w*(r if np.isfinite(r) else 3.0))
        for (T,m,ct) in regs:
            C=E.langmuir(T,E.CAGES[st][m][0],E.CAGES[st][m][1],*kp[gas])
            out.append(WREG*np.log(C/ct))
        return np.array(out)
    p0=np.array([np.log(kp0[gas][0]),kp0[gas][1]])
    sol=least_squares(resid,p0,method='lm',xtol=1e-5,ftol=1e-5,max_nfev=200)
    return np.exp(sol.x[0]),sol.x[1],sol

def fit_fbeta(kp,fb0):
    allanch=[(g,)+a for g in ANCH for a in ANCH[g]]
    def resid(p):
        fb={'sI':(p[0],p[1]),'sII':(p[2],p[3])}
        out=[]
        for (gas,T,P,w,s) in allanch:
            r,_=anchor_resid(T,P,gas,s,kp,fb)
            out.append(w*(r if np.isfinite(r) else 3.0))
        return np.array(out)
    sol=least_squares(resid,np.array(fb0['sI']+fb0['sII']),method='lm',xtol=1e-4,ftol=1e-4,max_nfev=80)
    return {'sI':(sol.x[0],sol.x[1]),'sII':(sol.x[2],sol.x[3])},sol

def full_Peq(T,gas,struct,kp,fb,Pguess=None):
    """true equilibrium P via bracketing on anchor_resid"""
    def resid(lnP):
        r,_=anchor_resid(T,np.exp(lnP),gas,struct,kp,fb)
        return r if np.isfinite(r) else np.nan
    if Pguess: lo,hi=np.log(Pguess/8),np.log(Pguess*8); n=25
    else: lo,hi=np.log(1e2),np.log(5e8); n=90
    lns=np.linspace(lo,hi,n)
    vals=np.array([resid(L) for L in lns])
    for k in range(n-1):
        if np.isfinite(vals[k]) and np.isfinite(vals[k+1]) and vals[k]*vals[k+1]<0:
            try: return float(np.exp(brentq(resid,lns[k],lns[k+1],xtol=1e-9)))
            except Exception: pass
    return None

if __name__=='__main__':
    t0=time.time()
    kp=dict(E.KIH); kp['C3H8']=(220.0,3.40e-10,0.6502e-10)
    fb={'sI':E.FBETA['sI'],'sII':E.FBETA['sII']}
    for rnd in range(3):
        for gas in ['CH4','CO2','N2','C2H6','C3H8']:
            e,s,sol=fit_gas(gas,kp,fb)
            kp[gas]=(e,s,kp[gas][2])
            print(f"r{rnd} {gas}: eps/k={e:8.2f} sigma={s*1e10:.4f} cost={sol.cost:.5f} [{time.time()-t0:.0f}s]",flush=True)
        fb,solb=fit_fbeta(kp,fb)
        print(f"r{rnd} fbeta: sI=({fb['sI'][0]:.4f},{fb['sI'][1]:.1f}) sII=({fb['sII'][0]:.4f},{fb['sII'][1]:.1f}) cost={solb.cost:.4f} [{time.time()-t0:.0f}s]",flush=True)
    np.save('calib2.npy',{'kp':kp,'fb':fb},allow_pickle=True)
    print("=== validation: true Peq vs anchors ===")
    tot=[]
    for gas in ANCH:
        st=sorted(ANCH[gas])[0][3]
        Pg=None
        for (T,P,w,s) in sorted(ANCH[gas]):
            Pm=full_Peq(T,gas,s,kp,fb,Pguess=P)
            if Pm:
                d=abs(Pm/P-1)*100; tot.append((w,d))
                print(f"  {gas} {s} T={T}: model {Pm/1e6:8.3f} vs {P/1e6:8.3f} MPa  ({d:5.1f}%)")
                Pg=Pm
            else:
                print(f"  {gas} {s} T={T}: NO ROOT (anchor {P/1e6} MPa)")
    aad=sum(d for w,d in tot)/len(tot)
    waad=sum(w*d for w,d in tot)/sum(w for w,_ in tot)
    print(f"AAD={aad:.2f}%  weightedAAD={waad:.2f}%  [{time.time()-t0:.0f}s]")
