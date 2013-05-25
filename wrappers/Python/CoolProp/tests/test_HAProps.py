import CoolProp
from CoolProp.HumidAirProp import HAProps
import numpy as np

def test_TRP():
    for R in np.linspace(0, 1, 11):
        for p in [101.325]:#np.linspace(0.1, 1000, 10):
            for T in np.linspace(220,373.15,1000):
                for o in ['W','H','S','V']:
                    yield check_HAProps,o,'T',T,'R',R,'P',p


def check_HAProps(*args):
    HAProps(*args)
    
if __name__=='__main__':
    import nose
    nose.runmodule()