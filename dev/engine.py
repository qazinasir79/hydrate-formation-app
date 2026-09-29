import numpy as np
from scipy.optimize import brentq
kB = 1.380649e-23; NA = 6.02214076e23; Rg = 8.314462618
CAL = 4.184

# ---------------- components ----------------
# order: [CH4, C2H6, C3H8, N2, CO2, H2O]
NAMES = ['CH4','C2H6','C3H8','N2','CO2','H2O']
IW = 5
# Tc(K), Pc(Pa), omega, PRSV kappa1 (kappa2=kappa3=0)
TC  = np.array([190.56, 305.32, 369.83, 126.19, 304.13, 647.14])
PC  = np.array([4.599e6, 4.872e6, 4.248e6, 3.3958e6, 7.377e6, 22.064e6])
OMG = np.array([0.01142, 0.0995, 0.1523, 0.0372, 0.2239, 0.3443])
# PRSV2 kappa1,kappa2,kappa3 (Stryjek-Vera 1986b, via DWSIM prsv2.dat)
K1  = np.array([-0.00159, 0.02669, 0.03136, 0.01996, 0.04285, -0.06635])
K2  = np.array([0.1521, 0.1358, 0.2757, 0.3162, 0.0, 0.0199])
K3  = np.array([0.517, 0.424, 0.447, 0.535, 0.0, 0.443])
# paper Table 1: solvent water(i?) solute gas(j?): pij qij pji qji   kij = pij*T+qij
KIJ_P = {'CH4':(0.00364,-0.631,0.0029,-0.96),
         'C2H6':(0.00165,-0.627,0.00295,-0.986),
         'C3H8':(0.00024,-0.413,0.00198,-0.647),
         'N2':(0.0039,-0.142,0.003,-0.42),
         'CO2':(0.0094,-0.517,0.0078,-0.258)}
# Kihara: eps/k (K), sigma (m), a (m)  -- initial: fitted to P-P A,B correlations
KIH = {'CH4':(156.31,3.2046e-10,0.3317e-10),
       'C2H6':(175.93,3.2661e-10,0.4799e-10),
       'C3H8':(200.0,3.45e-10,0.6502e-10),   # guess, to be fitted
       'N2':(131.50,3.2230e-10,0.341e-10),
       'CO2':(172.17,2.9856e-10,0.6786e-10)}

# ---------------- PRSV2 EOS ----------------
def kappa(i,T):
    Tr = T/TC[i]; s = np.sqrt(Tr)
    k0 = 0.378893 + 1.4897153*OMG[i] - 0.17131848*OMG[i]**2 + 0.0196554*OMG[i]**3
    return k0 + (K1[i] + K2[i]*(K3[i]-Tr)*(1-s))*(1+s)*(0.7-Tr)

def pure_ab(i,T):
    kap = kappa(i,T)
    alpha = (1+kap*(1-np.sqrt(T/TC[i])))**2
    a = 0.457235*Rg*Rg*TC[i]**2/PC[i]*alpha
    b = 0.077796*Rg*TC[i]/PC[i]
    return a,b

DK_AQ = {  # aqueous-phase-only correction of k(g,w): k -> k + dK0 + dK1*(T-273.15)  [v8]
    'CH4': (-0.100000, +0.000000),
    'C2H6': (-0.112773, +0.002866),
    'C3H8': (-0.027736, -0.001056),
    'N2': (-0.496822, +0.007746),
    'CO2': (-1.990619, -0.008089),
}  # aqueous shift of k(gas,water): dK0 + dK1*(T-273.15)
def _aq_shift(g, T):
    s = DK_AQ[g]
    return (s[0] + s[1]*(T-273.15)) if isinstance(s, tuple) else s
def kij_pair(i,j,T,assign='gas_first',x=None):
    # returns k_ij, k_ji for components i,j; x: phase composition (activates aqueous shift)
    if i==j: return 0.0,0.0
    gi = NAMES[i] if i!=IW else None
    gj = NAMES[j] if j!=IW else None
    aq = (x is not None) and x[IW] > 0.5
    if gi in KIJ_P and gj is None:   # i=gas, j=water
        p1,q1,p2,q2 = KIJ_P[gi]
        k2 = p2*T+q2 + (_aq_shift(gi,T) if aq else 0.0)
        if assign=='gas_first': return p1*T+q1, k2
        else: return k2, p1*T+q1
    if gj in KIJ_P and gi is None:   # i=water, j=gas
        p1,q1,p2,q2 = KIJ_P[gj]
        k2 = p2*T+q2 + (_aq_shift(gj,T) if aq else 0.0)
        if assign=='gas_first': return k2, p1*T+q1
        else: return p1*T+q1, k2
    return 0.0,0.0  # gas-gas neglected per paper

def mix_params(x,T,assign='gas_first',_x=None):
    n = len(x)
    ai = np.zeros(n); bi = np.zeros(n)
    for i in range(n):
        ai[i],bi[i] = pure_ab(i,T)
    K = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            if i==j: continue
            kij,kji = kij_pair(i,j,T,assign,x=x)
            D = x[i]*kij + x[j]*kji
            K[i,j] = 0.0 if abs(D)<1e-14 else kij*kji/D
    aij = np.sqrt(np.outer(ai,ai))*(1-K)
    np.fill_diagonal(aij, ai)
    am = x @ aij @ x
    bm = x @ bi
    return am,bm,ai,bi,aij,K

def cubic_Z(A,B,root):
    # Z^3 -(1-B)Z^2 + (A-3B^2-2B)Z -(AB-B^2-B^3)=0
    c = [1.0, -(1-B), A-3*B*B-2*B, -(A*B-B*B-B**3)]
    rts = np.roots(c)
    real = sorted(z.real for z in rts if abs(z.imag)<1e-9 and z.real>B*(1+1e-9))
    if not real: return None
    return real[-1] if root=='v' else real[0]

def lnphi(x,T,P,root,assign='gas_first'):
    am,bm,ai,bi,aij,K = mix_params(x,T,assign)
    A = am*P/(Rg*T)**2; B = bm*P/(Rg*T)
    Z = cubic_Z(A,B,root)
    if Z is None: return None
    if Z-B<=0 or Z+(1-np.sqrt(2.0))*B<=0: return None
    s2 = np.sqrt(2.0)
    L = np.log((Z+(1+s2)*B)/(Z+(1-s2)*B))
    f = Z-1-np.log(Z-B)-A/(2*s2*B)*L
    fA = -L/(2*s2*B)
    fB = 1/(Z-B) + A*L/(2*s2*B*B) - A/(2*s2*B)*((1+s2)/(Z+(1+s2)*B)-(1-s2)/(Z+(1-s2)*B))
    n = len(x)
    # d(am)/dx_i
    dam = np.zeros(n)
    for i in range(n):
        s = 2*x[i]*ai[i]
        for j in range(n):
            if j==i: continue
            s += 2*x[j]*aij[i,j]
            kij,kji = kij_pair(i,j,T,assign,x=x)
            D = x[i]*kij+x[j]*kji
            if abs(D)>1e-14:
                # d a_ij/d x_i = sqrt(ai aj) * k_ij k_ji * k_ij / D^2
                s += 2*x[i]*x[j]*np.sqrt(ai[i]*ai[j])*kij*kji*kij/D**2
        dam[i] = s
    dbm = bi.copy()
    Abar = float(x @ dam)*P/(Rg*T)**2
    Bbar = float(x @ dbm)*P/(Rg*T)
    out = np.zeros(n)
    for i in range(n):
        Ai_ = dam[i]*P/(Rg*T)**2; Bi_ = dbm[i]*P/(Rg*T)
        out[i] = f + fA*(Ai_-Abar) + fB*(Bi_-Bbar)
    return out

# ---------------- pure Psat via EOS (for validation) ----------------
def psat_pure(i,T):
    def resid(lnP):
        P=np.exp(lnP)
        lv=lnphi(np.array([1.0]*1*0+[1.0]) if False else np.eye(6)[i],T,P,'l')
        vv=lnphi(np.eye(6)[i],T,P,'v')
        if lv is None or vv is None: return None
        return lv[i]-vv[i]
    lo,hi=np.log(1.0),np.log(PC[i]*0.9999)
    lns=np.linspace(lo,hi,300); vals=[]
    for L in lns:
        v=resid(L); vals.append(v if v is not None else np.nan)
    vals=np.array(vals)
    for k in range(len(vals)-1):
        if np.isfinite(vals[k]) and np.isfinite(vals[k+1]) and vals[k]*vals[k+1]<0:
            return float(np.exp(brentq(resid,lns[k],lns[k+1],xtol=1e-10)))
    return None

# ---------------- Langmuir / vdW-P ----------------
GLX, GLW = np.polynomial.legendre.leggauss(200)
CAGES = {  # struct: [(R,z), ...] small,large
 'sI': [(3.917e-10,20),(4.33e-10,24)],
 'sII':[(3.917e-10,20),(4.73e-10,28)]}
NU = {'sI':(2/46.,6/46.), 'sII':(16/136.,8/136.)}
NW = {'sI':46., 'sII':136.}

def langmuir(T,Rc,z,eps_k,sig,a):
    rmax=(Rc-2*a)*0.9999999
    if rmax<=1e-12: return 0.0
    r=0.5*rmax*(GLX+1); jac=0.5*rmax
    u=1-r/Rc-a/Rc; v=1+r/Rc-a/Rc
    d10=(u**-10-v**-10)/10; d11=(u**-11-v**-11)/11
    d4=(u**-4-v**-4)/4; d5=(u**-5-v**-5)/5
    eps=eps_k*kB
    w=2*z*eps*(sig**12/(Rc**11*r)*(d10+(a/Rc)*d11)-sig**6/(Rc**5*r)*(d4+(a/Rc)*d5))
    arg=np.clip(w/(kB*T),-700,700)
    return 4*np.pi/(kB*T)*np.sum(GLW*np.exp(-arg)*r**2)*jac

def langmuir_set(T,gas_list,struct,params=None):
    out={}
    for g in gas_list:
        e,s,a = (params[g] if params else KIH[g])
        Cs=[]
        for (Rc,z) in CAGES[struct]:
            Cs.append(langmuir(T,Rc,z,e,s,a))
        out[g]=Cs
    return out

# empty hydrate fugacity (Martin-Peters form; A,B refittable)
FBETA = {'sI':(28.4,5800.0), 'sII':(28.8,6000.0)}
def psat_beta(T,struct,fb=None):
    par = fb[struct] if fb else FBETA[struct]
    A,B = par[0],par[1]
    D = par[2] if len(par)>2 else 0.0
    return np.exp(A-B/T+D/(T*T))
def Vbeta(T,P_MPa,struct,revised_sI=True):
    if struct=='sI':
        base=(11.835+2.217e-5*T+2.242e-6*T*T)**3*1e-30*NA/46.0
        if revised_sI: return base+1.5265e-9*(P_MPa-0.1)-2.8106e-12*(P_MPa-0.1)**2
        return base-8.006e-9*P_MPa+5.448e-12*P_MPa**2
    base=(17.13+2.249e-4*T+2.013e-6*T*T+1.009e-9*T**3)**3*1e-30*NA/136.0
    return base-8.006e-9*P_MPa+5.448e-12*P_MPa**2

def psat_ice(T): return np.exp(28.89-6134.9/T)
V_ICE=1.965e-5

LAM = {'sI':0.0,'sII':0.0}   # lattice non-ideality (guest-guest/strain) [K]; ln fH -= sum(nu*LAM*theta)/T
def fH_w(T,P,fgas,struct,fb=None,kparams=None,revised_sI=True):
    # fgas: dict gas->fugacity (Pa)
    Cset=langmuir_set(T,list(fgas.keys()),struct,kparams)
    nu=NU[struct]
    S=0.0; thetas={}
    for m,(Rc,z) in enumerate(CAGES[struct]):
        tot=0.0
        for g,f in fgas.items():
            tot+=Cset[g][m]*f
        for g,f in fgas.items():
            thetas[(g,m)]=Cset[g][m]*f/(1+tot)
        S+=nu[m]*np.log(1.0-min(tot/(1+tot),1-1e-300))
        S-=nu[m]*LAM[struct]*min(tot/(1+tot),1.0)/T
    psb=psat_beta(T,struct,fb)
    Vb=Vbeta(T,P/1e6,struct,revised_sI)
    fbeta=psb*np.exp(Vb*(P-psb)/(Rg*T))
    return fbeta*np.exp(S), thetas, S, fbeta

# ---------------- flash ----------------
def wilson_K(T,P):
    return PC/P*np.exp(5.373*(1+OMG)*(1-TC/T))

def flash_TP(z,T,P,assign='gas_first',maxiter=500,tol=1e-10):
    # two-phase VL; returns (x,y,beta,converged). Zero-feed components stay exactly zero.
    z=np.asarray(z,float); z=np.maximum(z,0.0); z=z/z.sum()
    act=z>0
    K=np.clip(wilson_K(T,P),1e-8,1e8)
    beta=0.5
    def split(K,beta):
        x=z/(1+beta*(K-1))
        x=np.where(act,x,0.0); x=x/x.sum()
        y=K*x
        y=np.where(act,y,0.0); s=y.sum()
        if s<=0 or not np.isfinite(s): return None,None
        return x,y/s
    for it in range(maxiter):
        def rr(b):
            v=np.sum(z*(K-1)/(1+b*(K-1)))
            return v if np.isfinite(v) else (1e6 if b<0.5 else -1e6)
        r0,r1=rr(1e-12),rr(1-1e-12)
        if r0>0 and r1<0:
            beta=brentq(rr,1e-12,1-1e-12,xtol=1e-14)
        else:
            beta=1e-10 if r0<=0 else 1-1e-10
        x,y=split(K,beta)
        if x is None: return None
        lphi=lnphi(x,T,P,'l',assign); vphi=lnphi(y,T,P,'v',assign)
        if lphi is None and vphi is None: return None
        if lphi is None: lphi=vphi.copy()
        if vphi is None: vphi=lphi.copy()
        if not (np.all(np.isfinite(lphi[act])) and np.all(np.isfinite(vphi[act]))): return None
        Knew=np.where(act,np.exp(np.clip(lphi-vphi,-40,40)),K)
        err=np.max(np.abs(Knew[act]/K[act]-1))
        K=K*np.sqrt(Knew/K) if err>1.0 else Knew   # geometric damping when far
        if err<tol:
            x,y=split(K,beta)
            if x is None: return None
            return x,y,beta,True
    x,y=split(K,beta)
    if x is None: return None
    return x,y,beta,False

# ---------------- equilibrium solvers ----------------
def equilib_hydrate(T,zgas,struct,assign='gas_first',fb=None,kparams=None,revised_sI=True,
                    Plo=2e2,Phi=4e8,nscan=120):
    z=np.zeros(6); z[IW]=0.8; 
    for gi,gn in enumerate(NAMES[:-1]): z[gi]=0.2*zgas[gi]
    z=z/z.sum()
    gas_idx=[i for i in range(5) if zgas[i]>0]
    def resid(lnP):
        P=np.exp(lnP)
        if T<273.15:
            # H-I-V: ice fugacity, vapor phase gas + water in equil with ice
            fi=psat_ice(T)*np.exp(V_ICE*(P-psat_ice(T))/(Rg*T))
            yw=min(fi/P,0.5)
            for _ in range(6):
                y=np.zeros(6); y[IW]=yw
                for i in gas_idx: y[i]=(1-yw)*zgas[i]
                vp=lnphi(y,T,P,'v',assign)
                if vp is None: break
                yw_new=min(max(fi/(P*np.exp(vp[IW])),1e-12),0.5)
                if abs(yw_new-yw)<1e-12: yw=yw_new; break
                yw=0.5*yw+0.5*yw_new
            fgas={}
            vp=lnphi(y,T,P,'v',assign)
            if vp is None: return 1e6
            for i in gas_idx: fgas[NAMES[i]]=y[i]*np.exp(vp[i])*P
            fH,_,_,_=fH_w(T,P,fgas,struct,fb,kparams,revised_sI)
            return np.log(fH)-np.log(fi)
        else:
            fl=flash_TP(z,T,P,assign)
            if fl is None: return 1e6
            x,y,beta,conv=fl
            # identify aqueous phase
            if x[IW]>=y[IW]: xaq,yg=x,y
            else: xaq,yg=y,x
            lp=lnphi(xaq,T,P,'l',assign); vp=lnphi(yg,T,P,'v',assign)
            if lp is None or vp is None: return 1e6
            fLw=xaq[IW]*np.exp(lp[IW])*P
            fgas={}
            for i in gas_idx:
                fgas[NAMES[i]]=yg[i]*np.exp(vp[i])*P
            fH,_,_,_=fH_w(T,P,fgas,struct,fb,kparams,revised_sI)
            return np.log(fH)-np.log(fLw)
    lns=np.linspace(np.log(Plo),np.log(Phi),nscan)
    vals=[]
    for L in lns:
        try: vals.append(resid(L))
        except Exception: vals.append(np.nan)
    vals=np.array(vals)
    for k in range(len(vals)-1):
        if np.isfinite(vals[k]) and np.isfinite(vals[k+1]) and vals[k]*vals[k+1]<0:
            try:
                r=brentq(resid,lns[k],lns[k+1],xtol=1e-8)
                return float(np.exp(r)),resid
            except Exception: pass
    return None,resid

def structure_select(T,zgas,**kw):
    out={}
    for st in ['sI','sII']:
        P,res=equilib_hydrate(T,zgas,st,**kw)
        out[st]=P
    if out['sI'] is None and out['sII'] is None: return None,out
    if out['sI'] is None: return 'sII',out
    if out['sII'] is None: return 'sI',out
    return ('sI' if out['sI']<=out['sII'] else 'sII'),out

if __name__=='__main__':
    print("=== pure Psat validation (PRSV, kappa1 published) ===")
    tests={'H2O':[(273.16,611.657),(298.15,3169.9),(323.15,12349),(373.15,101420)],
           'CH4':[(111.66,101325),(150,1.042e6)],
           'CO2':[(216.58,517.95e3),(253.15,1.786e6),(273.15,3.485e6),(293.15,5.729e6)],
           'C2H6':[(184.55,101325)],
           'C3H8':[(231.05,101325),(273.15,0.4745e6),(293.15,0.836e6)],
           'N2':[(77.35,101325),(100,0.7783e6)]}
    for nm,pts in tests.items():
        i=NAMES.index(nm)
        for T,Pt in pts:
            P=psat_pure(i,T)
            print(f"  {nm} T={T}: Psat={P if P is None else round(P,1)} (ref {Pt})  err={'' if P is None else f'{(P/Pt-1)*100:.2f}%'}")
