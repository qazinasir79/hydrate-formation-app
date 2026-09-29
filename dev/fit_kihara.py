import numpy as np
from scipy.optimize import brentq, least_squares, fsolve
kB = 1.380649e-23; NA = 6.02214076e23; Rg = 8.314462618
xs, wg = np.polynomial.legendre.leggauss(400)

def langmuir(T, Rc, z, eps_k, sig, a, cutoff='R-2a'):
    rmax = (Rc - 2*a) if cutoff=='R-2a' else (Rc - a)
    if rmax <= 0: return 0.0
    r = 0.5*rmax*(xs+1); jac = 0.5*rmax
    u = 1 - r/Rc - a/Rc; v = 1 + r/Rc - a/Rc
    d10=(u**-10 - v**-10)/10; d11=(u**-11 - v**-11)/11
    d4=(u**-4 - v**-4)/4; d5=(u**-5 - v**-5)/5
    eps = eps_k*kB
    w = 2*z*eps*(sig**12/(Rc**11*r)*(d10+(a/Rc)*d11) - sig**6/(Rc**5*r)*(d4+(a/Rc)*d5))
    arg = np.clip(w/(kB*T), -700, 700)
    return 4*np.pi/(kB*T)*np.sum(wg*np.exp(-arg)*r**2)*jac

sI_s=(3.917e-10,20); sI_l=(4.33e-10,24)
def th(C,f): x=C*f; return x/(1+x)

# Fit CH4 (eps_k, sigma) to match paper Table 6: T=273.2, f=phi*2.6e6, thetaL=0.969, thetaS=0.935
# and Table 7: T=276.7, equilibrium P~3.75 MPa -> thetaL=0.981, thetaS=0.917
a_CH4 = 0.3317e-10
def resid(p):
    eps_k, sig = p[0], p[1]*1e-10
    out=[]
    T=273.2; f=0.87*2.6e6
    Cs=langmuir(T,*sI_s,eps_k,sig,a_CH4); Cl=langmuir(T,*sI_l,eps_k,sig,a_CH4)
    out += [th(Cs,f)-0.935, th(Cl,f)-0.969]
    T=276.7; f=0.88*3.75e6
    Cs=langmuir(T,*sI_s,eps_k,sig,a_CH4); Cl=langmuir(T,*sI_l,eps_k,sig,a_CH4)
    out += [th(Cs,f)-0.917, th(Cl,f)-0.981]
    return out
sol = least_squares(resid, [160.0, 3.46], method='lm', xtol=1e-12)
print("CH4 fit:", sol.x, "residuals:", sol.fun)
# What do the "remembered" P-P params give?
for (e,s,label) in [(159.95,3.464,'P-P memory'),(166.04,3.464,'alt'),(232.2,3.505,'paper T2')]:
    T=273.2; f=0.87*2.6e6
    Cs=langmuir(T,*sI_s,e,s*1e-10,a_CH4); Cl=langmuir(T,*sI_l,e,s*1e-10,a_CH4)
    print(f"{label}: eps/k={e} sig={s} -> thetaS={th(Cs,f):.4f} thetaL={th(Cl,f):.4f}  (C_S={Cs:.3e}, C_L={Cl:.3e} 1/Pa)")

# check cutoff sensitivity with fitted params
e,s = sol.x[0], sol.x[1]*1e-10
for cut in ['R-2a','R-a']:
    T=273.2; f=0.87*2.6e6
    Cs=langmuir(T,*sI_s,e,s,a_CH4,cut); Cl=langmuir(T,*sI_l,e,s,a_CH4,cut)
    print(f"cutoff {cut}: thetaS={th(Cs,f):.4f} thetaL={th(Cl,f):.4f}")

# Fit CO2 to Table 8: T=273.2, P_eq=1.28 MPa, phi~0.9 -> thetaL=0.986, thetaS=0.581
a_CO2=0.6786e-10
def resid2(p):
    eps_k, sig = p[0], p[1]*1e-10
    T=273.2; f=0.90*1.28e6
    Cs=langmuir(T,*sI_s,eps_k,sig,a_CO2); Cl=langmuir(T,*sI_l,eps_k,sig,a_CO2)
    return [th(Cs,f)-0.581, th(Cl,f)-0.986]
sol2 = least_squares(resid2, [171.0, 3.42], method='lm', xtol=1e-12)
print("CO2 fit:", sol2.x, "residuals:", sol2.fun)
for (e,s,label) in [(171.44,3.418,'P-P memory'),(491.78,3.418,'paper T2')]:
    T=273.2; f=0.90*1.28e6
    Cs=langmuir(T,*sI_s,e,s*1e-10,a_CO2); Cl=langmuir(T,*sI_l,e,s*1e-10,a_CO2)
    print(f"{label}: thetaS={th(Cs,f):.4f} thetaL={th(Cl,f):.4f} (C_S={Cs:.3e}, C_L={Cl:.3e})")
