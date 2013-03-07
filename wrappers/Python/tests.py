
import CoolProp as CP

print CP.__version__
print CP.__svnrevision__
from CoolProp.CoolProp import Props

import timeit

#Specific heat (kJ/kg/K) of 20% ethylene glycol as a function of T
print Props('C','T',298.15,'P',101.325,'EG-20%')

#Density of Air at standard atmosphere in kg/m^3
print Props('D','T',298.15,'P',101.325,'Air')

#Saturation temperature of Water at 1 atm
print Props('T','P',101.325,'Q',0,'Water')

#Saturated vapor density of R134a at 0C
print Props('H','T',273.15,'Q',1,'R134a')

#Using properties from REFPROP to get R410A density
print Props('D','T',300,'P',100,'REFPROP-MIX:R32[0.697615]&R125[0.302385]')

#Check that the same as using pseudo-pure
print Props('D','T',300,'P',100,'R410A')

