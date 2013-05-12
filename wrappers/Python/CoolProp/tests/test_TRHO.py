import CoolProp
from CoolProp.CoolProp import Props
import numpy as np

def test_Trho():
    for Fluid in CoolProp.__fluids__:
        for T in np.linspace(Props(Fluid,'Tmin')+1e-5,Props(Fluid,'Tcrit')+100,20):
            for rho in np.linspace(1e-10,Props(Fluid,'rhocrit')*3,20):
                yield check_rho,Fluid,T,rho

def check_rho(Fluid,T,rho):
    p = Props('P','T',T,'D',rho,Fluid)
    
if __name__=='__main__':
    import nose
    nose.runmodule()