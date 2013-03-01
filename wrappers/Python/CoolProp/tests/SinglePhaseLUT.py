import numpy as np
from CoolProp import State as ST
from CoolProp import CoolProp as CP
N=3
TT=np.zeros((2*N+1,2*N+1))
PP=np.zeros((2*N+1,2*N+1))
YY1=np.zeros((2*N+1,2*N+1))
YY2=np.zeros((2*N+1,2*N+1))

CP.debug(0)
CP.set_1phase_LUT_params("Air",N,N,200,500,200,1200)
ST.set_1phase_LUT_params("Air",N,N,200,500,200,1200)

for key in ['D','H','U','C','O','dpdT']:
    CP.UseSinglePhaseLUT(False)
    for i,T in enumerate(np.linspace(291,355,2*N+1)):
        for j,p in enumerate(np.linspace(201,299,2*N+1)):
            TT[i,j]=T
            PP[i,j]=p
            YY2[i,j]=CP.Props(key,'T',T,'P',p,"Air")
    #Using the LUT
    CP.UseSinglePhaseLUT(True)
    for i,T in enumerate(np.linspace(291,355,2*N+1)):
        for j,p in enumerate(np.linspace(201,299,2*N+1)):
            TT[i,j]=T
            PP[i,j]=p
            YY1[i,j]=CP.Props(key,'T',T,'P',p,"Air")
    print key,'Worst Error is ', np.max(np.abs(YY1/YY2-1)),' %'

ST.LUT(True)
ST.debug(0)
S=ST.State('Air',dict(T=300,P=300))
for i,T in enumerate(np.linspace(291,355,2*N+1)):
    for j,p in enumerate(np.linspace(201,299,2*N+1)):
        TT[i,j]=T
        PP[i,j]=p
        S.update(dict(T=T,P=p))
        print S
        #YY2[i,j]=CP.Props(key,'T',T,'P',p,"Air")
    