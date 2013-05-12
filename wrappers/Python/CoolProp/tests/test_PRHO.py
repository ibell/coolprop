import CoolProp
from CoolProp.CoolProp import Props
import numpy as np

def test_Prho():
    for Fluid in CoolProp.__fluids__:
        for p in np.linspace(Props(Fluid,'ptriple')+1e-5,Props(Fluid,'pcrit')*10,20):
            for rho in np.linspace(1e-10,Props(Fluid,'rhocrit')*3,20):
                yield check_rho,Fluid,p,rho

def check_rho(Fluid,p,rho):
    T = Props('T','P',p,'D',rho,Fluid)
    
if __name__=='__main__':
    import nose
    nose.runmodule()