import numpy as np
from scipy.optimize import least_squares
kB = 1.380649e-23
xs, wg = np.polynomial.legendre.leggauss(400)
ATM = 101325.0

def langmuir(T, Rc, z, eps_k, sig, a):
    rmax = (Rc - 2*a)*0.9999999
    if rmax <= 0: return 1e-300
    r = 0.5*rmax*(xs+1); jac = 0.5*rmax
    u = 1 - r/Rc - a/Rc; v = 1 + r/Rc - a/Rc
    d10=(u**-10 - v**-10)/10; d11=(u**-11 - v**-11)/11
    d4=(u**-4 - v**-4)/4; d5=(u**-5 - v**-5)/5
    eps = eps_k*kB
    w = 2*z*eps*(sig**12/(Rc**11*r)*(d10+(a/Rc)*d11) - sig**6/(Rc**5*r)*(d4+(a/Rc)*d5))
    arg=np.clip(w/(kB*T),-700,700)
    return 4*np.pi/(kB*T)*np.sum(wg*np.exp(-arg)*r**2)*jac

sI_s=(3.917e-10,20); sI_l=(4.33e-10,24)
sII_s=(3.917e-10,20); sII_l=(4.73e-10,28)

# P-P A,B: C[atm^-1] = (A/T) exp(B/T);  small A x1e-3, large A x1e-2 (K/atm)
PP = {
 'CH4':  dict(struct='sI',  a=0.3317, S=(3.7237e-3,2708.8), L=(1.8372e-2,2737.9)),
 'CO2':  dict(struct='sI',  a=0.6786, S=(1.1978e-3,2860.5), L=(0.8507e-2,3277.9)),
 'N2':   dict(struct='sII', a=0.341,  S=(3.8087e-3,2205.5), L=(1.8420e-2,2301.3)),
 'C2H6': dict(struct='sI',  a=0.4799, S=None,               L=(0.6906e-2,3631.6)),
}
Ts = np.array([260,265,270,273.15,278,283,288,293,298,300.])
def Ccorr(AB,T): A,B=AB; return (A/T)*np.exp(B/T)/ATM   # 1/Pa

for gas,d in PP.items():
    cs,cl = (sI_s,sI_l) if d['struct']=='sI' else (sII_s,sII_l)
    a = d['a']*1e-10
    def resid(p):
        e, s = p[0], p[1]*1e-10
        out=[]
        if d['S'] is not None:
            tgt = Ccorr(d['S'],Ts); got=np.array([langmuir(T,*cs,e,s,a) for T in Ts])
            out += list(np.log(got/tgt))
        tgt = Ccorr(d['L'],Ts); got=np.array([langmuir(T,*cl,e,s,a) for T in Ts])
        out += list(np.log(got/tgt))
        return out
    best=None
    for e0 in (120,150,160,170,200,230):
        for s0 in (2.8,3.0,3.15,3.3,3.5):
            try:
                sol=least_squares(resid,[e0,s0],method='lm',xtol=1e-14,ftol=1e-14)
                if best is None or sol.cost<best.cost: best=sol
            except Exception: pass
    r=np.array(best.fun)
    print(f"{gas}: eps/k={best.x[0]:.2f} K, sigma={best.x[1]:.4f} A, a={d['a']} -> max|ln dev|={np.max(np.abs(r)):.3f} rms={np.sqrt(np.mean(r**2)):.3f}")
