from __future__ import division
import CoolProp
from CoolProp.CoolProp import Props, get_REFPROPname, rhosatL_anc, rhosatV_anc, psatL_anc, psatV_anc
import numpy as np
                    
def test_ancillaries():
    for Fluid in reversed(sorted(CoolProp.__fluids__)):
        T = np.linspace(Props(Fluid,'Tmin'),Props(Fluid,'Tcrit')-1,5)
        for T in T:
            yield check_ancillaries,Fluid,T

def check_ancillaries(Fluid,T):
        
    RPFluid = Fluid
##     if get_REFPROPname(Fluid) == 'N/A':
##         RPFluid = Fluid
##     else:
##         RPFluid = 'REFPROP-' + get_REFPROPname(Fluid)
        
    rhoL = Props('D','T',T,'Q',0,RPFluid)
    rhoV = Props('D','T',T,'Q',1,RPFluid)
    pL = Props('P','T',T,'Q',0,RPFluid)*1000
    pV = Props('P','T',T,'Q',1,RPFluid)*1000
            
    rhoL_anc = rhosatL_anc(Fluid, T)
    rhoV_anc = rhosatV_anc(Fluid, T)
    pL_anc = psatL_anc(Fluid, T)
    pV_anc = psatL_anc(Fluid, T)
    
    #  Check they are consistent
#     if abs(rhoL/rhoL_anc-1) > 0.02:
#         raise AssertionError('err: {err:g} % fit: {fit:g} eos: {eos:g}'.format(err = (rhoL/rhoL_anc-1)*100, fit = rhoL_anc, eos = rhoL))
    if abs(rhoV/rhoV_anc-1) > 0.02:
        raise AssertionError('err: {err:g} % fit: {fit:g} eos: {eos:g}'.format(err = (rhoV/rhoV_anc-1)*100, fit = rhoV_anc, eos = rhoV))
#     if abs(pL/pL_anc-1) > 0.02:
#         raise AssertionError('err: {err:g} % fit: {fit:g} eos: {eos:g}'.format(err = (pL/pL_anc-1)*100, fit = pL_anc, eos = pL))
    
if __name__=='__main__':
    import nose
    nose.runmodule()