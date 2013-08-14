import CoolProp
from CoolProp.CoolProp import Props, get_REFPROPname
import numpy as np

modes = ['normal']
# Check if REFPROP is supported, the Props call should work without throwing exception if it is supported
try:
    Props('D','T',300,'Q',1,'REFPROP-Propane')
    modes.append('REFPROP')
except ValueError:
    pass
            
twophase_inputs = [('T','D'),('T','Q'),('P','Q'),('P','H'),('P','S'),('P','D'),('H','S')]
singlephase_inputs = [('T','D'),('T','P'),('P','H'),('P','S'),('P','D'),('H','S')]

singlephase_outputs = ['T','P','H','S','A','O','C','G','V','L','C0','U']

def test_subcrit_singlephase_consistency():
    for Fluid in reversed(sorted(CoolProp.__fluids__)):
        T = (Props(Fluid,'Tmin')+Props(Fluid,'Tcrit'))/2.0
        for mode in modes:
            rhoL = Props('D','T',T,'Q',0,Fluid)
            rhoV = Props('D','T',T,'Q',1,Fluid)
            for rho in [rhoL*1.05,rhoV*0.95]:
                for inputs in singlephase_inputs:
                    yield check_consistency,Fluid,mode,T,rho,inputs
                    
def test_subcrit_twophase_consistency():
    for Fluid in reversed(sorted(CoolProp.__fluids__)):
        T = (Props(Fluid,'Tmin')+Props(Fluid,'Tcrit'))/2.0
        for mode in modes:
            rhoL = Props('D','T',T,'Q',0,Fluid)
            rhoV = Props('D','T',T,'Q',1,Fluid)
            for rho in [1/(0.5/rhoL+0.5/rhoV)]: # vapor mass quality of 0.5
                for inputs in twophase_inputs:
                    yield check_consistency,Fluid,mode,T,rho,inputs

def check_consistency(Fluid,mode,T,rho,inputs):
        
    if get_REFPROPname(Fluid) == 'N/A':
        return
        
    if mode == 'REFPROP':
        Fluid = 'REFPROP-' + get_REFPROPname(Fluid)
        
    #  Evaluate the inputs; if inputs is ('T','P'), calculate the temperature and the pressure
    Input1 = Props(inputs[0],'T',T,'D',rho,Fluid)
    Input2 = Props(inputs[1],'T',T,'D',rho,Fluid)
    
    #  Evaluate using the inputs given --> back to T,rho
    TEOS = Props('T',inputs[0],Input1,inputs[1],Input2,Fluid)
    DEOS = Props('D',inputs[0],Input1,inputs[1],Input2,Fluid)
    
    #  Check they are consistent
    if abs(TEOS/T-1) > 1e-1 or abs(DEOS/rho-1) > 1e-1:
        raise AssertionError('{T:g} K {D:g} kg/m^3 {inputs:s} = in1: {in1:g} in2: {in2:g}  T: {TEOS:g} D: {DEOS:g}'.format(T = T, 
                                                                                                                           D = rho, 
                                                                                                                           TEOS = TEOS, 
                                                                                                                           DEOS = DEOS, 
                                                                                                                           inputs = str(inputs), 
                                                                                                                           in1 = Input1, 
                                                                                                                           in2 = Input2)
                            )
    
if __name__=='__main__':
    import nose
    nose.runmodule()