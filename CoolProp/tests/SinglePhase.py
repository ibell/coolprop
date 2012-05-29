import CoolProp
from CoolProp.CoolProp import Props,UseSaturationLUT
import numpy as np

fp=open('log.txt','w')

def test_superheated():
    for Fluid in CoolProp.__fluids__:
        fp.write(Fluid+'\n')
        for T in np.linspace(Props(Fluid,'Ttriple')+0.1,0.95*Props(Fluid,'Tcrit'),20):
            p=Props('P','T',T,'Q',1.0,Fluid)
            for DTsup in np.linspace(0.01,100,20):
                yield check_rho,Fluid,T+DTsup,p
    
def test_subcooled():
    for Fluid in CoolProp.__fluids__:
        fp.write(Fluid+'\n')
        for T in np.linspace(Props(Fluid,'Ttriple')+0.1,0.95*Props(Fluid,'Tcrit'),20):
            p=Props('P','T',T,'Q',1.0,Fluid)
            for DTsub in np.linspace(0.01,10,20):
                yield check_rho,Fluid,T-DTsub,p
    
def test_supercritical():
    for Fluid in CoolProp.__fluids__:
        fp.write(Fluid+'\n')
        for p in np.linspace(Props(Fluid,'pcrit'),Props(Fluid,'pcrit')*2,20):
            for T in np.linspace(Props(Fluid,'Tcrit')+0.1,0.95*Props(Fluid,'Tcrit'),20):
                yield check_rho,Fluid,T,p

def check_rho(Fluid,T,p):
    rhoEOS = Props('D','T',T,'P',p,Fluid)
    pEOS = Props('P','T',T,'D',rhoEOS,Fluid)
    if not (abs(pEOS/p-1)*100<1e-3):
        fp.write('%g,%g,%g,%g\n'%(T,p,rhoEOS,pEOS))
    assert (abs(pEOS/p-1)*100<1e-3)
    
if __name__=='__main__':
    for args in test_superheated():
        try:
            args[0](*args[1::])
        except AssertionError,a:
            pass
    