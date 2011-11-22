from CoolProp.CoolProp import Props,Help as CPHelp,UseSaturationLUT,Tsat
from CoolProp.HumidAirProp import HAHelp,HumAir

print Props('D','T',300,'P',100,'REFPROP-MIX:R32[0.697615]&R125[0.302385]')
print Props('D','T',300,'P',100,'R410A')

Q=0.5
print Props('H','T',280,'Q',Q,'R290')
D=Props('D','T',280,'Q',Q,'R290')
print D
print Props('H','T',280,'D',D,'R290')

print Props('H','T',280,'Q',Q,'REFPROP-Propane')
D=Props('D','T',280,'Q',Q,'REFPROP-Propane')
print D
print Props('H','T',280,'D',D,'REFPROP-Propane')

CPHelp()

UseSaturationLUT(0)
print Props('P','T',275.2,'Q',1,'Water')
UseSaturationLUT(1)
print Props('P','T',275.2,'Q',1,'Water')

UseSaturationLUT(1)
print str(Tsat('Water',101.325,1.0,300)-273.15) +'C'
UseSaturationLUT(0)
print str(Tsat('Water',101.325,1.0,300)-273.15) +'C'

HAHelp()
Tdp,W,h,RH,v=HumAir(290,101.325,2,0.005)
## Tdp,W,h,RH,v=HumAir(290,101.325,5,h)

## print 'a'
## for i in range(10000):
##     Tdp,W,h,RH,v=HumAir(290,101.325,2,0.005)
##     Props('D','T',298.15,'P',101.325,'Air')
## print 'a'
## print 'b'
## for i in range(10000):
##     Tdp,W,h,RH,v=HumAir(290,101.325,5,h)
## print 'b'