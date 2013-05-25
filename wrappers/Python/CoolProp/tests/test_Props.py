import unittest
from CoolProp.CoolProp import Props
import CoolProp
import numpy as np
       
def test_input_types():
    for Fluid in CoolProp.__fluids__:
        for Tvals in [0.5*Props(Fluid,'Tmin')+0.5*Props(Fluid,'Tcrit'),
                      [Props(Fluid,'Tmin')+1e-5,Props(Fluid,'Tcrit')-1e-5],
                      np.linspace(Props(Fluid,'Tmin')+1e-5, Props(Fluid,'Tcrit')-1e-5,30)
                      ]:
            yield check_type, Fluid, Tvals
        
def check_type(fluid, Tvals):
    Props('P','T',Tvals,'Q',0,fluid)

if __name__=='__main__':
    import nose
    nose.runmodule()