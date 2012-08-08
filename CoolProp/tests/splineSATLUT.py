from CoolProp.CoolProp import Props,UseSaturationLUT
import numpy as np
from scipy.interpolate import UnivariateSpline,interp1d
from matplotlib import pyplot as plt
import CoolProp

print CoolProp.__pure_fluids__.split(',')
for Fluid in CoolProp.__pure_fluids__.split(','):
    if Fluid=='R32':
        continue
    if Fluid=='R744':
        continue
    
    TLUT = np.linspace(1.01*Props(Fluid,'Ttriple'),0.99*Props(Fluid,'Tcrit'),5000)
    rho_fit = np.zeros_like(TLUT)
    rho_EOS = np.zeros_like(TLUT)
    for i,T in enumerate(TLUT):
        UseSaturationLUT(0)
        rho_EOS[i] = Props('D','T',T,'Q',0.0,Fluid)
        UseSaturationLUT(1)
        rho_fit[i] = Props('D','T',T,'Q',0.0,Fluid)
    errperc=np.max(np.abs(rho_fit/rho_EOS*100-100))
    print 'Worst error with CoolProp is {0} % for {1}'.format(errperc,Fluid)
    
    if errperc>1.0:
        plt.plot(TLUT,rho_EOS,TLUT,rho_fit)
        plt.show()