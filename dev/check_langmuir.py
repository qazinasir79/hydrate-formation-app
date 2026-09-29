import numpy as np
kB = 1.380649e-23
NA = 6.02214076e23
Rg = 8.314462618

def wcell(r, R, z, eps, sig, a):
    d10 = ((1 - r/R - a/R)**(-10) - (1 + r/R - a/R)**(-10))/10.0
    d11 = ((1 - r/R - a/R)**(-11) - (1 + r/R - a/R)**(-11))/11.0
    d4  = ((1 - r/R - a/R)**(-4)  - (1 + r/R - a/R)**(-4))/4.0
    d5  = ((1 - r/R - a/R)**(-5)  - (1 + r/R - a/R)**(-5))/5.0
    return 2*z*eps*(sig**12/(R**11*r)*(d10 + (a/R)*d11) - sig**6/(R**5*r)*(d4 + (a/R)*d5))

xs, wg = np.polynomial.legendre.leggauss(500)

def langmuir(T, R, z, eps_k, sig, a):
    rmax = R - 2*a
    r = 0.5*rmax*(xs+1); jac = 0.5*rmax
    eps = eps_k*kB
    wv = wcell(r, R, z, eps, sig, a)
    arg = np.clip(wv/(kB*T), -700, 700)
    integrand = np.exp(-arg)*r**2
    return 4*np.pi/(kB*T)*np.sum(wg*integrand)*jac   # 1/Pa

# cages: (R in m, z)
sI_s = (3.917e-10, 20); sI_l = (4.33e-10, 24)
sII_s = (3.917e-10, 20); sII_l = (4.73e-10, 28)
# Kihara from paper Table 2: eps/k (K), sigma (A), a (A)
CH4 = (232.2, 3.505e-10, 0.3317e-10)
CO2 = (491.78, 3.418e-10, 0.6786e-10)
C2H6 = (414.8, 4.0e-10, 0.5695e-10)
C3H8 = (493.7, 4.519e-10, 0.6502e-10)
N2 = (142.1, 3.469e-10, 0.341e-10)

def theta(C, f):
    Cf = C*f
    return Cf/(1+Cf)

# Test 1: CH4 at 273.2 K, equilibrium P = 2.6 MPa; f_gas approx 0.85*2.6e6
T = 273.2
for name,(Rc,z) in [('sI_small',sI_s),('sI_large',sI_l)]:
    C = langmuir(T, Rc, z, *CH4)
    for phi in [1.0, 0.85]:
        f = phi*2.6e6
        print(f"CH4 {name}: C={C:.4e} 1/Pa  theta(f={f/1e6:.2f} MPa)={theta(C,f):.4f}")
print("Table 6 targets: thetaL=0.969 thetaS=0.935")

# Test 2: CO2 at 273.2 K equilibrium P ~ 1.28 MPa
for name,(Rc,z) in [('sI_small',sI_s),('sI_large',sI_l)]:
    C = langmuir(T, Rc, z, *CO2)
    for phi in [1.0, 0.88]:
        f = phi*1.28e6
        print(f"CO2 {name}: C={C:.4e} 1/Pa theta(f={f/1e6:.2f} MPa)={theta(C,f):.4f}")
print("Table 8 targets: thetaL=0.986 thetaS=0.581")

# Test 3: full equilibrium for CH4 sI at 273.2 K: solve fH = fL
def psat_w(T):  # Magnus, Pa
    return 610.94*np.exp(17.625*(T-273.15)/(T-30.04))
def psat_beta_sI(T):
    return np.exp(28.4 - 5800.0/T)
def Vbeta_sI(T,P_MPa, revised=True):
    base = (11.835 + 2.217e-5*T + 2.242e-6*T*T)**3 * 1e-30 * NA/46.0
    if revised:
        return base + 1.5265e-9*(P_MPa-0.1) - 2.8106e-12*(P_MPa-0.1)**2
    return base - 8.006e-9*P_MPa + 5.448e-12*P_MPa**2

def fH(T, P, phi_gas, struct='sI'):
    if struct=='sI':
        Cs = langmuir(T, *sI_s, *CH4); Cl = langmuir(T, *sI_l, *CH4)
        vs, vl = 2/46., 6/46.
        Pb = psat_beta_sI(T); Vb = Vbeta_sI(T, P/1e6)
    S = vs*np.log(1-theta(Cs, phi_gas*P)) + vl*np.log(1-theta(Cl, phi_gas*P))
    fbeta = Pb*np.exp(Vb*(P-Pb)/(Rg*T))
    return fbeta*np.exp(S), np.exp(S), theta(Cs,phi_gas*P), theta(Cl,phi_gas*P)

def fL(T,P):
    Vw = 1.8e-5
    psw = psat_w(T)
    return psw*np.exp(Vw*(P-psw)/(Rg*T))  # x_w ~ 1, phi_sat ~ 1

T = 273.2
from scipy.optimize import brentq
def resid(logP):
    P = np.exp(logP)
    fh,_,_,_ = fH(T,P,0.87)
    return np.log(fh) - np.log(fL(T,P))
try:
    lo, hi = np.log(1e4), np.log(2e8)
    # scan
    ps = np.linspace(lo,hi,300)
    vals = [resid(p) for p in ps]
    sign_changes = [(ps[i],ps[i+1]) for i in range(len(ps)-1) if vals[i]*vals[i+1]<0]
    print("sign change brackets (MPa):", [(np.exp(a)/1e6, np.exp(b)/1e6) for a,b in sign_changes])
    for a,b in sign_changes:
        root = brentq(resid, a, b)
        print("P_eq(CH4, 273.2K) =", np.exp(root)/1e6, "MPa   (experimental ~2.6)")
except ImportError:
    print("no scipy; manual")
