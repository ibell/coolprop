
from CoolProp import Props, Tcrit
import matplotlib as mpl
import matplotlib.pyplot as pyplot

import numpy as np

errCode=0
errString=''

print Props('H','T',300,'P',500,"R134a")
print Props('S','T',300,'P',500,"R134a")
print Props('D','T',300,'P',500,"R134a")

N=2000
Tsat=np.linspace(240,350,N)
hsatL=0*Tsat
hsatV=0*Tsat
psat=0*Tsat

for i in range(N):
    hsatL[i]=Props('H','T',Tsat[i],'Q',0,'R290')
    hsatV[i]=Props('H','T',Tsat[i],'Q',1,'R290')
    psat[i]=Props('P','T',Tsat[i],'Q',1,'R290')

pyplot.plot(hsatL*1000,psat,'k',hsatV*1000,psat,'k')
pyplot.xlabel('Enthalpy [J/kg]')
pyplot.ylabel('Pressure [kPa]')

##
##Ps=Props('P','T',277,'Q',1,'R290')
####hsL=Props('H','T',277,'P',Ps,'R290')*1000
####hsV=Props('H','T',277,'Q',1,'R290')*1000
##hs1=Props('H','T',277+5,'P',Ps,'R290')*1000
##hs2=208144
##
##Pd=Props('P','T',312.15,'Q',1,'R290')
##hd3=Props('H','T',312.15-3,'P',Pd,'R290')*1000
##hd4=481785

pts=np.loadtxt('Z:\\Documents\\Code\\Emerson\\LumpedCycle\\LumpedCycle\\pts.txt',delimiter=',')
print pts

print np.size(pts)

h_IHX=np.r_[pts[0,0],pts[1,0]]
p_IHX=np.r_[pts[0,1],pts[1,1]]

h_Comp=np.r_[pts[2,0],pts[3,0]]
p_Comp=np.r_[pts[2,1],pts[3,1]]

h_Cond=np.r_[pts[4,0],pts[5,0]]
p_Cond=np.r_[pts[4,1],pts[5,1]]

h_TXV=np.r_[pts[6,0],pts[7,0]]
p_TXV=np.r_[pts[6,1],pts[7,1]]

##print(Tcrit('R290'))

####pyplot.plot(np.r_[hsL,hsV],np.r_[Ps,Ps])
pyplot.plot(h_IHX,p_IHX,'b',lw=2)
pyplot.plot(h_Comp,p_Comp,'b',lw=2)
pyplot.plot(h_Cond,p_Cond,'b',lw=2)
pyplot.plot(h_TXV,p_TXV,'b',lw=2)
pyplot.show()






