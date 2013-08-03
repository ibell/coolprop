import CoolProp
from CoolProp.CoolProp import Props
from math import log10

###############################################################################
###############################################################################

def test_Tmin():
    for Fluid in CoolProp.__fluids__:
        yield check_Tmin,Fluid
        
def check_Tmin(Fluid):
    Tmin = Props(Fluid,'Tmin')
    Ttriple = Props(Fluid,'Ttriple')
    assert Tmin >= Ttriple

###############################################################################
###############################################################################
    
def test_ptriple():
    for Fluid in CoolProp.__fluids__:
        yield check_ptriple,Fluid
        
def check_ptriple(Fluid):
    Tmin = max(Props(Fluid,'Tmin'), Props(Fluid,'Ttriple'))
    ptriple = Props(Fluid,'ptriple')
    ptripleEOS = Props('P','T',Tmin,'Q',1,Fluid)
    if not abs(ptriple/ptripleEOS-1) < 0.01:
        raise ValueError('ptriple should be: '+str(ptripleEOS))

###############################################################################
###############################################################################

def test_accentric():
    for Fluid in CoolProp.__fluids__:
        yield check_accentric,Fluid
        
def check_accentric(Fluid):
    if Props(Fluid,"Tmin") < 0.7*Props(Fluid,'Tcrit'):
        accentricEOS = -log10(Props("P",'Q',1,'T',Props(Fluid,"Tcrit")*0.7,Fluid)/Props(Fluid,"pcrit"))-1
    else:
        return
    try:
        accentric = Props(Fluid,'accentric')
    except ValueError:
        raise ValueError('accentric should be: '+str(accentricEOS))
    
    if not (accentric/accentricEOS-1) < 0.01:
        raise ValueError('accentric should be: '+str(accentricEOS))
    
if __name__=='__main__':
    import nose
    nose.runmodule()