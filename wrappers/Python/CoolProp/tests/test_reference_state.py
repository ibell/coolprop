import CoolProp.CoolProp as CP
# Check reference state setting
def test_IIR():
    CP.set_reference_state('Propane','IIR')
    assert(abs(CP.Props('H','T',273.15,'Q',0,'Propane')-200) < 1e-10)
    assert(abs(CP.Props('S','T',273.15,'Q',0,'Propane')-1) < 1e-10)

def test_NBP():
    CP.set_reference_state('Propane','NBP')
    assert(abs(CP.Props('H','P',101.325,'Q',0,'Propane')-0) < 1e-10)
    assert(abs(CP.Props('S','P',101.325,'Q',0,'Propane')-0) < 1e-10)

def test_ASHRAE():
    CP.set_reference_state('Propane','ASHRAE')
    assert(abs(CP.Props('H','T',233.15,'Q',0,'Propane')-0) < 1e-10)
    assert(abs(CP.Props('S','T',233.15,'Q',0,'Propane')-0) < 1e-10)
    
if __name__=='__main__':
    import nose
    nose.runmodule()