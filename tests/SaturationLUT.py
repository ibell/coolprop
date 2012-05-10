import CoolProp
from CoolProp.CoolProp import Props,UseSaturationLUT
import numpy as np

def test_saturation_pure():
    for Fluid in CoolProp.__pure_fluids__.split(','):
        for T in np.linspace(Props(Fluid,'Ttriple'),Props(Fluid,'Tcrit')-1e-3b,300):
            for state in ['L','V']:
                yield check_LUT_rho,Fluid,T,state

def check_LUT_rho(Fluid,T,state):
    if state=='L':
        rho=Props('D','T',T,'Q',0.0,Fluid)
    elif state =='V':
        rho=Props('D','T',T,'Q',1.0,Fluid)
    

## def test_saturation_pure():
##     for Fluid in CoolProp.__pure_fluids__.split(','):
##         for T in np.linspace(Props(Fluid,'Ttriple'),Props(Fluid,'Tcrit'),20):
##             yield check_LUT_rho,Fluid,T

## def check_LUT_rho(Fluid,T):
##     UseSaturationLUT(0)
##     rhoL=Props('D','T',T,'Q',0.0,Fluid)
##     UseSaturationLUT(1)
##     rhoL_LUT=Props('D','T',T,'Q',0.0,Fluid)
##     
##     assert abs(rhoL/rhoL_LUT-1)<1e-6