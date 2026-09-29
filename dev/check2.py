import numpy as np
from scipy.optimize import brentq
kB = 1.380649e-23; NA = 6.02214076e23; Rg = 8.314462618
xs, wg = np.polynomial.legendre.leggauss(600)

def wcell(r, R, z, eps, sig, a):
    u = 1 - r/R - a/R; v = 1 + r/R - a/R
    d = lambda N: (u**(-N) - v**(-N))/N
    return 2*z*eps*(sig**12/(R**11*r)*(d(10) + (a/R)*d(11)) - sig**6/(R**5*r)*(d(4) + (a/R)*d(5)))

def langmuir(T, R, z, eps_k, sig, a):
    rmax = R - 2*a
    r = 0.5*rmax*(xs+1); jac = 0.5*rmax
    wv = wcell(r, R, z, eps_k*kB, sig, a)
    arg = np.clip(wv/(kB*T), -700, 700)
    return 4*np.pi/(kB*T)*np.sum(wg*np.exp(-arg)*r**2)*jac

sI_s=(3.917e-10,20); sI_l=(4.33e-10,24); sII_s=(3.917e-10,20); sII_l=(4.73e-10,28)
# Parrish-Prausnitz / Sloan-Koh classic Kihara: eps/k (K), sigma (m), a (m)
PP = {
 'CH4':(159.95,3.464e-10,0.3317e-10),
 'C2H6':(240.97,4.038e-10,0.4799e-10),
 'C3H8':(212.57,4.519e-10,0.6502e-10),
 'CO2':(171.44,3.418e-10,0.6786e-10),
 'N2':(76.4,3.469e-10,0.341e-10),
}
def th(C,f): x=C*f; return x/(1+x)

print("=== Occupancy check vs paper Tables 6-8 (P-P Kihara) ===")
T=273.2; f=0.87*2.6e6
Cs=langmuir(T,*sI_s,*PP['CH4']); Cl=langmuir(T,*sI_l,*PP['CH4'])
print(f"CH4 273.2K f=2.26MPa: thetaS={th(Cs,f):.4f} (target 0.935)  thetaL={th(Cl,f):.4f} (target 0.969)")
T=276.7; f=0.88*3.75e6
Cs=langmuir(T,*sI_s,*PP['CH4']); Cl=langmuir(T,*sI_l,*PP['CH4'])
print(f"CH4 276.7K f=3.3MPa: thetaS={th(Cs,f):.4f} (target 0.917)  thetaL={th(Cl,f):.4f} (target 0.981)")
T=273.2; f=0.90*1.28e6
Cs=langmuir(T,*sI_s,*PP['CO2']); Cl=langmuir(T,*sI_l,*PP['CO2'])
print(f"CO2 273.2K f=1.15MPa: thetaS={th(Cs,f):.4f} (target 0.581)  thetaL={th(Cl,f):.4f} (target 0.986)")

print()
print("=== Equilibrium curves with Martin-Peters f_beta (Pa) ===")
def psat_w(T): return 610.94*np.exp(17.625*(T-273.15)/(T-30.04))
def psat_ice(T): return np.exp(28.89-6134.9/T)
def psat_beta(T,struct): return np.exp(28.4-5800/T) if struct=='sI' else np.exp(28.8-6000/T)
def Vbeta(T,P_MPa,struct,revised_sI=True):
    if struct=='sI':
        base=(11.835+2.217e-5*T+2.242e-6*T*T)**3*1e-30*NA/46.0
        if revised_sI: return base + 1.5265e-9*(P_MPa-0.1) - 2.8106e-12*(P_MPa-0.1)**2
        return base - 8.006e-9*P_MPa + 5.448e-12*P_MPa**2
    return (17.13+2.249e-4*T+2.013e-6*T*T+1.009e-9*T**3)**3*1e-30*NA/136.0 - 8.006e-9*P_MPa + 5.448e-12*P_MPa**2
def cages(struct):
    return (sI_s,sI_l,2/46.,6/46.) if struct=='sI' else (sII_s,sII_l,16/136.,8/136.)

def Peq(T,gas,struct,phi_gas=0.9,fLmode='water'):
    (cs,cl,vs,vl)=cages(struct)
    Cs=langmuir(T,*cs,*PP[gas]); Cl=langmuir(T,*cl,*PP[gas])
    def resid(logP):
        P=np.exp(logP); f=phi_gas*P
        S=vs*np.log(1-th(Cs,f))+vl*np.log(1-th(Cl,f))
        fb=psat_beta(T,struct)*np.exp(Vbeta(T,P/1e6,struct)*(P-psat_beta(T,struct))/(Rg*T))
        if fLmode=='water':
            psw=psat_w(T) if T>=273.15 else psat_ice(T)
            Vw=1.8e-5 if T>=273.15 else 1.965e-5
            fl=psw*np.exp(Vw*(P-psw)/(Rg*T))
        return np.log(fb)+S-np.log(fl)
    lo,hi=np.log(5e2),np.log(5e8)
    ps=np.linspace(lo,hi,400); vals=np.array([resid(p) for p in ps])
    idx=np.where(np.sign(vals[:-1])*np.sign(vals[1:])<0)[0]
    if len(idx)==0: return None
    return np.exp(brentq(resid,ps[idx[0]],ps[idx[0]+1]))/1e6

for gas,struct,Ts,exp in [
  ('CH4','sI',[273.15,277.59,283.15,288.15,293.15],[2.60,4.03,7.2,12.4,19.5]),
  ('CO2','sI',[273.15,276.15,279.15,282.15],[1.28,1.79,2.45,3.6]),
  ('C2H6','sI',[273.15,277.15,281.15,285.15],[0.55,0.83,1.32,2.24]),
  ('C3H8','sII',[273.15,275.15,277.59],[0.172,0.26,0.42]),
  ('N2','sII',[273.15,277.15,281.15],[16.4,22.6,31.0]),
  ('CH4','sI',[255.15,265.15],[0.70,1.55]),  # H-I-V region rough expectations
]:
    print(f"--- {gas} ({struct}) ---")
    for T,e in zip(Ts,exp):
        p1=Peq(T,gas,struct,phi_gas=0.90)
        print(f"  T={T}: P_eq={p1 if p1 is None else round(p1,3)} MPa  (exp~{e})")
