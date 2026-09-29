import numpy as np, time
from scipy.optimize import least_squares, brentq
import engine as E
np.seterr(all='ignore')
ASSIGN='water_first'
from calibrate import ANCH, PPAB, TPP, WREG, cpp

_cache={}
def water_and_gas_fug(T,P,gas_list,assign=ASSIGN,use_cache=True):
    key=(round(T,6),round(P,3),tuple(gas_list),assign)
    if use_cache and key in _cache: return _cache[key]
    if T<273.15:
        fi=E.psat_ice(T)*np.exp(E.V_ICE*(P-E.psat_ice(T))/(E.Rg*T))
        y=np.zeros(6); y[E.IW]=min(fi/P,0.5)
        for g in gas_list: y[E.NAMES.index(g)]=(1-y[E.IW])/len(gas_list)
        for _ in range(12):
            vp=E.lnphi(y,T,P,'v',assign)
            if vp is None: break
            yw=min(max(fi/(P*np.exp(vp[E.IW])),1e-14),0.6)
            if abs(yw-y[E.IW])<1e-13: y[E.IW]=yw; break
            y[E.IW]=0.5*y[E.IW]+0.5*yw
            for g in gas_list: y[E.NAMES.index(g)]=(1-y[E.IW])/len(gas_list)
        vp=E.lnphi(y,T,P,'v',assign)
        if vp is None: res=(None,None,'fail')
        else:
            fg={g:y[E.NAMES.index(g)]*np.exp(vp[E.NAMES.index(g)])*P for g in gas_list}
            res=(fi,fg,'ice')
    else:
        z=np.zeros(6); z[E.IW]=0.8
        for g in gas_list: z[E.NAMES.index(g)]=0.2/len(gas_list)
        fl=E.flash_TP(z,T,P,assign,maxiter=150,tol=1e-9)
        if fl is None: res=(None,None,'fail')
        else:
            x,y,beta,conv=fl
            xaq,yg=(x,y) if x[E.IW]>=y[E.IW] else (y,x)
            lp=E.lnphi(xaq,T,P,'l',assign); vp=E.lnphi(yg,T,P,'v',assign)
            if lp is None or vp is None: res=(None,None,'fail')
            else:
                fLw=xaq[E.IW]*np.exp(lp[E.IW])*P
                fg={g:yg[E.NAMES.index(g)]*np.exp(vp[E.NAMES.index(g)])*P for g in gas_list}
                tag='lc' if yg[E.IW]<0.5 else 'lw'
                res=(fLw,fg,tag)
    if use_cache: _cache[key]=res
    return res

def anchor_resid(T,P,gas,struct,kp,fb):
    fLw,fg,tag=water_and_gas_fug(T,P,[gas])
    if fLw is None: return np.nan,tag
    fH,_,_,_=E.fH_w(T,P,fg,struct,fb,kp)
    return np.log(fH)-np.log(fLw),tag

def fit_gas(gas,kp0,fb):
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
    sol=least_squares(resid,p0,method='lm',diff_step=1e-4,xtol=1e-8,ftol=1e-8,max_nfev=300)
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
    sol=least_squares(resid,np.array(fb0['sI']+fb0['sII']),method='lm',diff_step=1e-5,xtol=1e-8,ftol=1e-8,max_nfev=100)
    return {'sI':(sol.x[0],sol.x[1]),'sII':(sol.x[2],sol.x[3])},sol

def full_Peq(T,gas,struct,kp,fb,Pguess=None):
    def resid(lnP):
        r,_=anchor_resid(T,np.exp(lnP),gas,struct,kp,fb,use_cache=False) if False else anchor_resid(T,np.exp(lnP),gas,struct,kp,fb)
        return r if np.isfinite(r) else np.nan
    if Pguess: lo,hi=np.log(Pguess/10),np.log(Pguess*10); n=40
    else: lo,hi=np.log(1e2),np.log(5e8); n=120
    lns=np.linspace(lo,hi,n)
    vals=np.array([resid(L) for L in lns])
    for k in range(len(vals)-1):
        if np.isfinite(vals[k]) and np.isfinite(vals[k+1]) and vals[k]*vals[k+1]<0:
            try: return float(np.exp(brentq(resid,lns[k],lns[k+1],xtol=1e-9)))
            except Exception: pass
    return None

def occupancy_report(gas,T,P,struct,kp,fb):
    fLw,fg,tag=water_and_gas_fug(T,P,[gas],use_cache=False)
    fH,thetas,S,fbeta=E.fH_w(T,P,fg,struct,fb,kp)
    return thetas,S,fg,tag

if __name__=='__main__':
    t0=time.time()
    # ---- diagnostics: C values of candidate param sets vs P-P ----
    print("=== Langmuir C diagnostics @273.15K (Pa^-1); P-P in parens ===")
    cands={'engine-KIH':{'CH4':E.KIH['CH4'],'CO2':E.KIH['CO2'],'N2':E.KIH['N2'],'C2H6':E.KIH['C2H6']},
           'classic':{'CH4':(154.38,3.0455e-10,0.3317e-10),'CO2':(169.09,2.9681e-10,0.6805e-10)},
           'paperT2':{'CH4':(232.2,3.505e-10,0.3317e-10),'CO2':(491.78,3.418e-10,0.6786e-10),'N2':(142.1,3.469e-10,0.341e-10)}}
    for name,sets in cands.items():
        for g,(e,s,a) in sets.items():
            Cs=E.langmuir(273.15,3.917e-10,20,e,s,a); Cl=E.langmuir(273.15,4.33e-10,24,e,s,a)
            ps=cpp(g,273.15,0); pl=cpp(g,273.15,1)
            ps=f"{ps:.3e}" if ps else "  --  "; pl=f"{pl:.3e}" if pl else "  --  "
            print(f"  {name:10s} {g:5s}: C_s={Cs:.3e} ({ps})  C_l={Cl:.3e} ({pl})")
    # ---- calibration ----
    kp=dict(E.KIH); kp['C3H8']=(220.0,3.40e-10,0.6502e-10)
    fb={'sI':E.FBETA['sI'],'sII':E.FBETA['sII']}
    for rnd in range(3):
        for gas in ['CH4','CO2','N2','C2H6','C3H8']:
            e,s,sol=fit_gas(gas,kp,fb)
            kp[gas]=(e,s,kp[gas][2])
            print(f"r{rnd} {gas}: eps/k={e:8.2f} sigma={s*1e10:.4f} cost={sol.cost:.5f} [{time.time()-t0:.0f}s]",flush=True)
        fb,solb=fit_fbeta(kp,fb)
        print(f"r{rnd} fbeta: sI=({fb['sI'][0]:.4f},{fb['sI'][1]:.1f}) sII=({fb['sII'][0]:.4f},{fb['sII'][1]:.1f}) cost={solb.cost:.4f} [{time.time()-t0:.0f}s]",flush=True)
    np.save('calib3.npy',{'kp':kp,'fb':fb},allow_pickle=True)
    print(f"saved. total {time.time()-t0:.0f}s")
