import CoolProp.CoolProp as CP
# Check reference state setting
def test_IIR():
    CP.set_reference_state('Propane','IIR')
    h_target = 200
    s_target = 1
    h_EOS = CP.Props('H','T',273.15,'Q',0,'Propane')
    s_EOS = CP.Props('S','T',273.15,'Q',0,'Propane')
    if abs(h_EOS-h_target) > 1e-10:
        raise AssertionError('h_target {h_target:g} h_EOS: {h_EOS:g}'.format(**locals()))
    if abs(s_EOS-s_target) > 1e-10:
        raise AssertionError('s_target {s_target:g} s_EOS: {s_EOS:g}'.format(**locals()))

def test_NBP():
    CP.set_reference_state('Propane','NBP')
    h_target = 0
    s_target = 0
    h_EOS = CP.Props('H','P',101.325,'Q',0,'Propane')
    s_EOS = CP.Props('S','P',101.325,'Q',0,'Propane')
    if abs(h_EOS-h_target) > 1e-10:
        raise AssertionError('h_target {h_target:g} h_EOS: {h_EOS:g}'.format(**locals()))
    if abs(s_EOS-s_target) > 1e-10:
        raise AssertionError('s_target {s_target:g} s_EOS: {s_EOS:g}'.format(**locals()))

def test_ASHRAE():
    CP.set_reference_state('Propane','ASHRAE')
    h_target = 0
    s_target = 0
    h_EOS = CP.Props('H','T',233.15,'Q',0,'Propane')
    s_EOS = CP.Props('S','T',233.15,'Q',0,'Propane')
    if abs(h_EOS-h_target) > 1e-10:
        raise AssertionError('h_target {h_target:g} h_EOS: {h_EOS:g}'.format(**locals()))
    if abs(s_EOS-s_target) > 1e-10:
        raise AssertionError('s_target {s_target:g} s_EOS: {s_EOS:g}'.format(**locals()))
    
if __name__=='__main__':
    import nose
    nose.runmodule()