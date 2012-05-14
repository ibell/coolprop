import CoolProp
from CoolProp.CoolProp import Props,UseSaturationLUT
import numpy as np

def test_saturation_pure():
    for Fluid in CoolProp.__pure_fluids__.split(','):
        for T in np.linspace(Props(Fluid,'Ttriple'),0.995*Props(Fluid,'Tcrit'),300):
            for state in ['L','V']:
                yield check_LUT_rho,Fluid,T,state

def check_LUT_rho(Fluid,T,state):
    if state=='L':
        rho=Props('D','T',T,'Q',0.0,Fluid)
    elif state =='V':
        rho=Props('D','T',T,'Q',1.0,Fluid)
    
def test_saturation_LUT():
    for Fluid in CoolProp.__pure_fluids__.split(','):
        for T in np.linspace(Props(Fluid,'Ttriple'),0.995*Props(Fluid,'Tcrit'),300):
            for qual in ['L','V']:
                yield check_saturation_LUT,Fluid,T,qual

def check_saturation_LUT(Fluid,T,qual):
    if qual=='L':
        UseSaturationLUT(0)
        rhoL=Props('D','T',T,'Q',0.0,Fluid)
        UseSaturationLUT(1)
        rhoL_LUT=Props('D','T',T,'Q',0.0,Fluid)
        assert abs(rhoL/rhoL_LUT-1)<1e-3
    elif qual=='V':
        UseSaturationLUT(0)
        rhoV=Props('D','T',T,'Q',1.0,Fluid)
        UseSaturationLUT(1)
        rhoV_LUT=Props('D','T',T,'Q',1.0,Fluid)
        assert abs(rhoV/rhoV_LUT-1)<1e-3

if __name__=='__main__':
    test_saturation_pure()