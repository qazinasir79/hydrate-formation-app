import numpy as np
from scipy.optimize import brentq, least_squares
Rg = 8.314462618

# PRSV2: alpha = [1+k(1-sqrt(Tr))]^2 ; k = k0 + [k1 + k2(k3-Tr)(1-sqrt(Tr))](1-sqrt(Tr))(0.7-Tr)
# k0 = 0.378893 + 1.4897153w - 0.17131848w^2 + 0.0196554w^3
COMPS = {
 # Tc(K), Pc(Pa), omega
 'H2O': (647.14, 22.064e6, 0.3443),
 'CH4': (190.56, 4.599e6, 0.01142),
 'C2H6':(305.32, 4.872e6, 0.0995),
 'C3H8':(369.83, 4.248e6, 0.1523),
 'N2':  (126.19, 3.3958e6, 0.0372),
 'CO2': (304.13, 7.377e6, 0.2239),
}
# reference Psat anchors (Pa)
PSAT = {
 'H2O': [(273.16,611.657),(283.15,1228.2),(293.15,2339.2),(298.15,3169.9),(303.15,4247.0),(313.15,7383.5),(323.15,12349),(333.15,19941),(343.15,31201),(353.15,47414),(363.15,70182),(373.15,101420)],
 'CH4': [(111.66,101325),(120,192e3),(130,540.6e3),(140,1.243e6),(150,2.369e6*0+1.042e6),(160,1.59e6*0+1.588e6),(170,2.331e6),(180,3.269e6),(190,4.36e6)],
 'CO2': [(216.58,517.95e3),(233.15,0.8006e6),(253.15,1.786e6),(273.15,3.485e6),(293.15,5.729e6),(300.0,6.71e6)],
 'C2H6':[(184.55,101325),(233.15,0.5538e6),(253.15,1.301e6),(273.15,2.677e6),(293.15,4.437e6)],
 'C3H8':[(231.05,101325),(253.15,0.2312e6),(273.15,0.4745e6),(293.15,0.8336e6),(313.15,1.383e6)],
 'N2':  [(77.35,101325),(90,0.3602e6),(100,0.7783e6),(110,1.441e6),(120,2.361e6)],
}
def kappa(T, Tc, w, k1, k2=0.0, k3=0.0):
    k0 = 0.378893 + 1.4897153*w - 0.17131848*w**2 + 0.0196554*w**3
    Tr = T/Tc; s = 1-np.sqrt(Tr)
    return k0 + (k1 + k2*(k3-Tr)*s)*s*(0.7-Tr)

def psat_eos(name, T, k1, k2=0.0, k3=0.0):
    Tc,Pc,w = COMPS[name]
    a0 = 0.457235*Rg*Rg*Tc*Tc/Pc
    b  = 0.077796*Rg*Tc/Pc
    kap = kappa(T,Tc,w,k1,k2,k3)
    alpha = (1+kap*(1-np.sqrt(T/Tc)))**2
    a = a0*alpha
    def P_of_v(v): return Rg*T/(v-b) - a/(v*v+2*b*v-b*b)
    # find Psat by equal fugacity: solve via bisection on pressure with cubic roots
    def lnphi(P, root):
        A = a*P/(Rg*T)**2; B = b*P/(Rg*T)
        # cubic Z^3 -(1-B)Z^2 + (A-3B^2-2B)Z -(AB -B^2 -B^3)=0
        c2 = -(1-B); c1 = A-3*B*B-2*B; c0 = -(A*B - B*B - B**3)
        roots = np.roots([1,c2,c1,c0])
        roots = sorted(z.real for z in roots if abs(z.imag)<1e-9 and z.real>B)
        if not roots: return None, None
        Z = roots[-1] if root=='v' else roots[0]
        s2 = np.sqrt(2)
        lp = Z-1-np.log(Z-B) - A/(2*s2*B)*np.log((Z+(1+s2)*B)/(Z+(1-s2)*B))
        return lp, Z
    lo, hi = 10.0, min(Pc*0.999, 1e8)
    # ensure bracket
    try:
        f = lambda lnP: (lnphi(np.exp(lnP),'l')[0] or 0) - (lnphi(np.exp(lnP),'v')[0] or 0)
        # scan
        lns = np.linspace(np.log(lo), np.log(hi), 200)
        vals = [f(x) for x in lns]
        for i in range(len(vals)-1):
            if vals[i] is not None and vals[i+1] is not None and vals[i]*vals[i+1] < 0:
                r = brentq(f, lns[i], lns[i+1], xtol=1e-10)
                return np.exp(r)
    except Exception as e:
        return None
    return None

print("=== fit k1 (k2=k3=0) ===")
fits={}
for name, data in PSAT.items():
    def resid(p):
        out=[]
        for T,Pt in data:
            P = psat_eos(name,T,p[0])
            out.append(np.log(P/Pt) if P else 10.0)
        return out
    sol = least_squares(resid,[0.1],method='lm',xtol=1e-12)
    r=np.array(sol.fun)
    fits[name]=sol.x[0]
    print(f"{name}: k1={sol.x[0]:.5f}  max|err|={np.max(np.abs(r))*100:.2f}%  rms={np.sqrt(np.mean(r**2))*100:.2f}%")
print(fits)
print("=== fit k1,k2 (k3=0) ===")
for name, data in PSAT.items():
    def resid(p):
        out=[]
        for T,Pt in data:
            P = psat_eos(name,T,p[0],p[1],0.0)
            out.append(np.log(P/Pt) if P else 10.0)
        return out
    sol = least_squares(resid,[fits[name],0.0],method='lm',xtol=1e-12)
    r=np.array(sol.fun)
    print(f"{name}: k1={sol.x[0]:.5f} k2={sol.x[1]:.5f}  max|err|={np.max(np.abs(r))*100:.2f}%  rms={np.sqrt(np.mean(r**2))*100:.2f}%")
