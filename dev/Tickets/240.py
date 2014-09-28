 
from __future__ import print_function
import sys
import numpy as np
import CoolProp
import CoolProp.CoolProp as CP


fluids = [ "DEB", "HCM", "HFE", "PMS1", "PMS2", "SAB", "HCB", "TCO" ]

if int(CoolProp.__version__[0])>4:
    for i in range(len(fluids)):
	fluids[i] = "INCOMP::{0}".format(fluids[i])
    
Tmin = [CP.PropsSI("Tmin","T",0,"P",0,fluid) for fluid in fluids]
Tmax = [CP.PropsSI("Tmax","T",0,"P",0,fluid) for fluid in fluids]
T = [np.linspace(Tmin[i]+1,Tmax[i]-1) for i in range(len(fluids))]
P = 100e5

print("Thermal conductivity in W/m/K")
for i in range(len(fluids)):
    print("{0:6s} ".format(fluids[i]),end="")
print("")

for i in range(len(T[0])):
    for j in range(len(fluids)):
	print("{0:6.4f} ".format(CP.PropsSI("L","T",T[j][i],"P",P,fluids[j])), end="")
    print("")

sys.exit(0)