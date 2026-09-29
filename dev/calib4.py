import numpy as np, time, json
from scipy.optimize import least_squares, brentq
import engine as E
np.seterr(all='ignore')
ASSIGN='water_first'
from calibrate import ANCH, PPAB, TPP, WREG, cpp

_cache={}
def wgf(T,P,gas_frac,assign=ASSIGN,use_cache=True):
    """gas_frac: dict gas->mole frac of gas phase (sums ~1). Returns (f_w, {gas:f}, tag)"""
    gl=tuple(sorted(gas_frac))
    key=(round(T,6),round(P,3),tuple(round(gas_frac[g],8) for g in gl),assign)
    if use_cache and key in _cache: return _cache[key]
    if T<273.15:
        fi=E.psat_ice(T)*np.exp(E.V_ICE*(P-E.psat_ice(T))/(E.Rg*T))
        y=np.zeros(6); y[E.IW]=min(fi/P,0.5)
        for g in gl: y[E.NAMES.index(g)]=(1-y[E.IW])*gas_frac[g]
        for _ in range(12):
            vp=E.lnphi(y,T,P,'v',assign)
            if vp is None: break
            yw=min(max(fi/(P*np.exp(vp[E.IW])),1e-14),0.6)
            if abs(yw-y[E.IW])<1e-13: y[E.IW]=yw; break
            y[E.IW]=0.5*y[E.IW]+0.5*yw
            for g in gl: y[E.NAMES.index(g)]=(1-y[E.IW])*gas_frac[g]
        vp=E.lnphi(y,T,P,'v',assign)
        if vp is None: res=(None,None,'fail')
        else:
            fg={g:y[E.NAMES.index(g)]*np.exp(vp[E.NAMES.index(g)])*P for g in gl}
            res=(fi,fg,'ice')
    else:
        z=np.zeros(6); z[E.IW]=0.8
        for g in gl: z[E.NAMES.index(g)]=0.2*gas_frac[g]
        fl=E.flash_TP(z,T,P,assign,maxiter=150,tol=1e-9)
        if fl is None: res=(None,None,'fail')
        else:
            x,y,beta,conv=fl
            xaq,yg=(x,y) if x[E.IW]>=y[E.IW] else (y,x)
            lp=E.lnphi(xaq,T,P,'l',assign); vp=E.lnphi(yg,T,P,'v',assign)
            if lp is None or vp is None: res=(None,None,'fail')
            else:
                fLw=xaq[E.IW]*np.exp(lp[E.IW])*P
                fg={g:yg[E.NAMES.index(g)]*np.exp(vp[E.NAMES.index(g)])*P for g in gl}
                tag='lc' if yg[E.IW]<0.5 else 'lw'
                res=(fLw,fg,tag)
    if use_cache: _cache[key]=res
    return res

GASES=['CH4','C2H6','C3H8','N2','CO2']
def unpack(p):
    kp={}; i=0
    for g in GASES:
        kp[g]=(np.exp(p[i]),p[i+1],E.KIH[g][2]); i+=2
    A,B,c=p[i],p[i+1],p[i+2]
    fb={'sI':(A,B),'sII':(A+c*273.15/B*B,B)}  # fbeta_sII(T)=exp(A-B/T+c) -> express as (A2,B2): A2=A+c, B2=B
    fb={'sI':(A,B),'sII':(A+c,B)}
    return kp,fb
def resid_all(p,hinge=True):
    kp,fb=unpack(p)
    out=[]
    for gas in GASES:
        gf={gas:1.0}
        for (T,P,w,st) in sorted(ANCH[gas]):
            fLw,fg,tag=wgf(T,P,gf)
            if fLw is None: out.append(w*3.0); continue
            fH,_,_,_=E.fH_w(T,P,fg,st,fb,kp)
            r=np.log(fH)-np.log(fLw)
            out.append(w*(r if np.isfinite(r) else 3.0))
    # regularization to P-P
    for gas in GASES:
        if gas in PPAB:
            st='sII' if gas in ('N2',) else 'sI'
            for T in TPP:
                for m in (0,1):
                    ct=cpp(gas,T,m)
                    if ct:
                        C=E.langmuir(T,E.CAGES[st][m][0],E.CAGES[st][m][1],*kp[gas])
                        out.append(WREG*np.log(C/ct))
    # Fig5 transition state: 277.6 K, y_CH4=0.9995, P*=4.0 MPa: both structures at equilibrium
    gf5={'CH4':0.9995,'C3H8':0.0005}
    fLw,fg,tag=wgf(277.6,4.0e6,gf5)
    if fLw is not None:
        for st in ('sI','sII'):
            fH,_,_,_=E.fH_w(277.6,4.0e6,fg,st,fb,kp)
            out.append(1.5*(np.log(fH)-np.log(fLw)))
    if hinge:
        # CH4 must be sI: at its sI-equilibrium states, sII residual > 0
        for (T,P) in [(273.15,2.60e6),(283.15,7.20e6)]:
            fLw,fg,tag=wgf(T,P,{'CH4':1.0})
            if fLw is None: continue
            fH2,_,_,_=E.fH_w(T,P,fg,'sII',fb,kp)
            r2=np.log(fH2)-np.log(fLw)
            out.append(3.0*min(0.0,r2-0.03))
        # C3H8 must be sII: at its equilibrium, sI residual > 0
        for (T,P) in [(277.59,0.42e6),(273.15,0.172e6)]:
            fLw,fg,tag=wgf(T,P,{'C3H8':1.0})
            if fLw is None: continue
            fH1,_,_,_=E.fH_w(T,P,fg,'sI',fb,kp)
            r1=np.log(fH1)-np.log(fLw)
            out.append(3.0*min(0.0,r1-0.03))
        # N2 must be sII
        fLw,fg,tag=wgf(273.15,16.4e6,{'N2':1.0})
        if fLw is not None:
            fH1,_,_,_=E.fH_w(273.15,16.4e6,fg,'sI',fb,kp)
            out.append(3.0*min(0.0,np.log(fH1)-np.log(fLw)-0.03))
    return np.array(out)

def fit(p0):
    lo=[];hi=[]
    for g in GASES:
        lo+= [np.log(80),2.5e-10]; hi+=[np.log(400),4.0e-10]
    lo+=[26.0,5200.0,0.0]; hi+=[31.0,6600.0,0.10]
    sol=least_squares(resid_all,p0,bounds=(lo,hi),diff_step=1e-4,xtol=1e-10,ftol=1e-10,max_nfev=400,verbose=0)
    return sol

if __name__=='__main__':
    t0=time.time()
    d=np.load('calib3.npy',allow_pickle=True).item()
    kp3,fb3=d['kp'],d['fb']
    p0=[]
    for g in GASES:
        e,s,a=kp3[g]
        if g=='C3H8': e,s=220.0,3.40e-10   # reset diverged
        p0+=[np.log(e),s]
    A,B=fb3['sI']
    c=np.log(np.exp(fb3['sII'][0]-fb3['sII'][1]/273.15)/np.exp(A-B/273.15))
    c=min(max(c,0.0),0.10)
    p0+=[A,B,c]
    p0=np.array(p0)
    print("p0:",[f"{v:.4g}" for v in p0])
    # warm cache
    resid_all(p0)
    print(f"cache warm: {len(_cache)} entries, {time.time()-t0:.1f}s")
    sol=fit(p0)
    p1=sol.x
    kp,fb=unpack(p1)
    print(f"fit done cost={sol.cost:.5f} nfev={sol.nfev} [{time.time()-t0:.0f}s]")
    for g in GASES: print(f"  {g}: eps/k={kp[g][0]:8.2f} sigma={kp[g][1]*1e10:.4f} a={kp[g][2]*1e10:.4f}")
    print(f"  fbeta sI: A={fb['sI'][0]:.4f} B={fb['sI'][1]:.2f} ; c={p1[-1]:.4f} (sII: A={fb['sII'][0]:.4f} B={fb['sII'][1]:.2f})")
    print(f"  fbeta_sI(273.15)={E.psat_beta(273.15,'sI',fb):.1f} Pa, fbeta_sII(273.15)={E.psat_beta(273.15,'sII',fb):.1f} Pa")
    json.dump({'kp':{g:list(map(float,kp[g])) for g in GASES},'fb':{'sI':list(map(float,fb['sI'])),'sII':list(map(float,fb['sII']))},'assign':ASSIGN},open('params_v1.json','w'),indent=1)
    np.save('calib4.npy',{'kp':kp,'fb':fb,'p':p1},allow_pickle=True)
    print("saved params_v1.json, calib4.npy")
