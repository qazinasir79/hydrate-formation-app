import numpy as np
kB = 1.380649e-23
xs, wg = np.polynomial.legendre.leggauss(400)

def langmuir_var(T, Rc, z, eps_k, sig, a, variant):
    # variant: 'V1' delta uses a/R, cross a/R  (Sloan printed form)
    #          'V2' delta uses 2a/R, cross 2a/R (exact shell avg of eq-21-type potential, offset 2a)
    #          'V3' delta uses 2a/R, cross a/R
    d = a/Rc if variant=='V1' else 2*a/Rc
    rmax = Rc*(1-d)*0.999999
    r = 0.5*rmax*(xs+1); jac=0.5*rmax
    u = 1 - r/Rc - d; v = 1 + r/Rc - d
    d10=(u**-10 - v**-10)/10; d11=(u**-11 - v**-11)/11
    d4=(u**-4 - v**-4)/4; d5=(u**-5 - v**-5)/5
    cross = a/Rc if variant=='V3' else d
    eps=eps_k*kB
    w = 2*z*eps*(sig**12/(Rc**11*r)*(d10+cross*d11) - sig**6/(Rc**5*r)*(d4+cross*d5))
    arg=np.clip(w/(kB*T),-700,700)
    return 4*np.pi/(kB*T)*np.sum(wg*np.exp(-arg)*r**2)*jac

sI_s=(3.917e-10,20); sI_l=(4.33e-10,24)
def th(C,f): x=C*f; return x/(1+x)

sets = {
 'paper_T2': dict(CH4=(232.2,3.505,0.3317), CO2=(491.78,3.418,0.6786), C2H6=(414.8,4.0,0.5695), C3H8=(493.7,4.519,0.6502), N2=(142.1,3.469,0.341)),
 'PP_mem':   dict(CH4=(159.95,3.464,0.3317), CO2=(171.44,3.418,0.6786), C2H6=(240.97,4.038,0.4799), C3H8=(212.57,4.519,0.6502), N2=(76.4,3.469,0.341)),
}
targets = {'CH4':(273.2, 0.87*2.6e6, 0.935, 0.969), 'CO2':(273.2, 0.90*1.28e6, 0.581, 0.986)}
for sname,S in sets.items():
    for var in ['V1','V2','V3']:
        for gas,(T,f,tS,tL) in targets.items():
            e,sg,aa = S[gas]
            Cs=langmuir_var(T,*sI_s,e,sg*1e-10,aa*1e-10,var); Cl=langmuir_var(T,*sI_l,e,sg*1e-10,aa*1e-10,var)
            print(f"{sname:9s} {var} {gas}: thetaS={th(Cs,f):.4f} (target {tS})  thetaL={th(Cl,f):.4f} (target {tL})")
