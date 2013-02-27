from CoolProp import CoolProp as CP
from PDSim.misc.datatypes import Collector
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

fluid = 'R245fa'
Ttriple = CP.Props(fluid,'Ttriple')
Tcrit = CP.Props(fluid,'Tcrit')

RHO,TTT,RHO0,TTT0,ERR = Collector(),Collector(),Collector(),Collector(),Collector()

rhomax = CP.Props('D','T',Ttriple,'Q',0,'R245fa')
#Build a database of "experimental" data
for T in np.linspace(Ttriple,Tcrit+50,80):
    for rho in np.linspace(1e-10,rhomax,80):
        
        if (T > Tcrit or rho > CP.rhosatL_anc(fluid,T) or rho < CP.rhosatV_anc(fluid,T) ):
            mu = CP.Props('V','T',T,'D',rho,'R245fa')
            mu_ref = CP.Props('V','T',T,'D',rho,'REFPROP-R245fa') 
            
            RHO << rho
            TTT << T
            ERR << mu/mu_ref-1
        
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(np.array(RHO.vec),np.array(TTT.vec),ERR.vec)
plt.show()
