import CoolProp
from CoolProp.CoolProp import Props
import numpy as np

def test_superheated():
    for Fluid in CoolProp.__fluids__:
        for T in np.linspace(Props(Fluid,'Tmin')+1e-5,Props(Fluid,'Tcrit')-1e-5,20):
            p=Props('P','T',T,'Q',1.0,Fluid)
            for DTsup in np.linspace(0.01,100,20):
                yield check_rho,Fluid,T+DTsup,p
    
def test_subcooled():
    for Fluid in CoolProp.__fluids__:
        for Tsat in np.linspace(Props(Fluid,'Tmin')+1e-5,Props(Fluid,'Tcrit')-1e-5,20):
            p = Props('P','T',Tsat,'Q',1.0,Fluid)
            for T in np.linspace(Props(Fluid,'Tmin')+1e-5,Tsat-1e-5,20):
                yield check_rho,Fluid,T,p
    
## def test_supercritical():
##     for Fluid in CoolProp.__fluids__:
##         for p in np.linspace(Props(Fluid,'pcrit'),Props(Fluid,'pcrit')*2,20):
##             for T in np.linspace(Props(Fluid,'Tcrit')+0.1,0.95*Props(Fluid,'Tcrit'),20):
##                 yield check_rho,Fluid,T,p

def check_rho(Fluid,T,p):
    rhoEOS = Props('D','T',T,'P',p,Fluid)
    pEOS = Props('P','T',T,'D',rhoEOS,Fluid)
    assert (abs(pEOS/p-1)<1e-1)
    
if __name__=='__main__':
    import nose
    nose.runmodule()
    