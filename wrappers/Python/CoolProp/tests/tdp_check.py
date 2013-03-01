import CoolProp.HumidAirProp as HAP
import numpy as np
import pylab

pylab.figure()
def Tdp(Tmin = 260,Tmax = 300, wmin = 0, wmax = 0.01):
    for T in np.linspace(Tmin,Tmax,100):
        TT = []
        WW = []
        for w in np.linspace(wmin,wmax,100):
            v = HAP.HAProps('D','T',T,'P',101.325,'W',w)
            print T,w,v
            TT.append(v)
            WW.append(w)
        
        pylab.plot(WW,TT)
            
    pylab.show()
if __name__=='__main__':
    Tdp()