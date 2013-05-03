import CoolProp
from CoolProp.CoolProp import Props
import numpy as np

def test_superheated():
    for Fluid in CoolProp.__fluids__:
        for T in np.linspace(Props(Fluid,'Tmin')+1e-5,Props(Fluid,'Tcrit')-1e-5,5):
            p=Props('P','T',T,'Q',1.0,Fluid)
            for DTsup in np.linspace(0.01,100,5):
                if p > 1e-8:
                    yield check_rho,Fluid,T+DTsup,p
    
def test_subcooled():
    for Fluid in CoolProp.__fluids__:
        Tmin = Props(Fluid,'Tmin')+1e-5
        rhomin = Props('D', 'T', Tmin, 'Q', 0, Fluid)
        for Tsat in np.linspace(Tmin, Props(Fluid,'Tcrit')-1e-5, 5):
            p = Props('P', 'T', Tsat, 'Q', 0, Fluid)
            for T in np.linspace(max(Tmin,Tsat-10), Tsat-1e-5, 5):
                if T > Tmin and T < Tsat:
                    yield check_rho,Fluid,T,p
#    
def test_supercritical():
    for Fluid in CoolProp.__fluids__:
        for p in np.linspace(Props(Fluid,'pcrit')*1.01,Props(Fluid,'pcrit')*2,5):
            for T in np.linspace(Props(Fluid,'Tmin')+1e-5,2*Props(Fluid,'Tcrit'),5):
                yield check_rho,Fluid,T,p

def check_rho(Fluid,T,p):
    rhoEOS = Props('D','T',T,'P',p,Fluid)
    pEOS = Props('P','T',T,'D',rhoEOS,Fluid)
    if abs(pEOS/p-1) > 1e-1 and p>1e-6:
        raise AssertionError('{pEOS:g} {p:g}'.format(pEOS = pEOS, p = p))
    
if __name__=='__main__':
    import nose
    nose.runmodule()