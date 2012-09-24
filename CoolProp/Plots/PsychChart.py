"""
This file implements a psychrometric chart for air at 1 atm
"""

from CoolProp.HumidAirProp import HAProps
import matplotlib.pyplot as plt
import numpy as np
import textwrap

import_template=(
"""
import numpy as np
import matplotlib.pyplot as plt
from CoolProp.HumidAirProp import HAProps

Tdb = np.linspace(-10,55,100)+273.15

#Make the figure and the axes
fig=plt.figure(figsize=(10,8))
ax=fig.add_axes((0.15,0.15,0.8,0.8))
ax.set_xlim(Tdb[0]-273.15,Tdb[-1]-273.15)
ax.set_ylim(0,0.03)
ax.set_xlabel(r"Dry bulb temperature [$^{\circ}$C]")
ax.set_ylabel(r"Humidity ratio ($m_{water}/m_{dry\ air}$) [-]")
"""
)

closure_template=(
"""
plt.show()
"""                  
)
Tdb = np.linspace(-10,55,100)+273.15

class SaturationLine(object):
    
    def plot(self,ax):
        w = [HAProps('W','T',T,'P',101.325,'R',1.0) for T in Tdb]
        ax.plot(Tdb-273.15,w,lw=2)
        
    def __str__(self):
        return textwrap.dedent("""
               # Saturation line
               w = [HAProps('W','T',T,'P',101.325,'R',1.0) for T in Tdb]
               ax.plot(Tdb-273.15,w,lw=2)
               """
               )
        
class HumidityLines(object):
    
    def __init__(self,RH_values):
        self.RH_values = RH_values
        
    def plot(self,ax):
        for RH in self.RH_values:
            w = [HAProps('W','T',T,'P',101.325,'R',RH) for T in Tdb]
            ax.plot(Tdb-273.15,w,'r',lw=1)
        
    def __str__(self):
        return textwrap.dedent("""
               # Humidity lines
               RHValues = {RHValues:s}
               for RH in RHValues:
                   w = [HAProps('W','T',T,'P',101.325,'R',RH) for T in Tdb]
                   ax.plot(Tdb-273.15,w,'r',lw=1)
               """.format(RHValues=str(self.RH_values))
               )

class EnthalpyLines(object):
    
    def __init__(self,H_values):
        self.H_values = H_values
        
    def plot(self,ax):
        for H in self.H_values:
            #Line goes from saturation to zero humidity ratio for this enthalpy
            T1 = HAProps('T','H',H,'P',101.325,'R',1.0)-273.15
            T0 = HAProps('T','H',H,'P',101.325,'R',0.0)-273.15
            w1 = HAProps('W','H',H,'P',101.325,'R',1.0)
            w0 = HAProps('W','H',H,'P',101.325,'R',0.0)
            ax.plot(np.r_[T1,T0],np.r_[w1,w0],'r',lw=1)
        
    def __str__(self):
        return textwrap.dedent("""
               # Humidity lines
               for H in {HValues:s}:
                   #Line goes from saturation to zero humidity ratio for this enthalpy
                   T1 = HAProps('T','H',H,'P',101.325,'R',1.0)-273.15
                   T0 = HAProps('T','H',H,'P',101.325,'R',0.0)-273.15
                   w1 = HAProps('W','H',H,'P',101.325,'R',1.0)
                   w0 = HAProps('W','H',H,'P',101.325,'R',0.0)
                   ax.plot(np.r_[T1,T0],np.r_[w1,w0],'r',lw=1)
               """.format(HValues=str(self.H_values))
               )
    
fig=plt.figure(figsize=(10,8))
ax=fig.add_axes((0.15,0.15,0.8,0.8))
ax.set_xlim(Tdb[0]-273.15,Tdb[-1]-273.15)
ax.set_ylim(0,0.03)
ax.set_xlabel(r"Dry bulb temperature [$^{\circ}$C]")
ax.set_ylabel(r"Humidity ratio ($m_{water}/m_{dry\ air}$) [-]")

SL = SaturationLine()
SL.plot(ax)

RHL = HumidityLines([0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])
RHL.plot(ax)

HL = EnthalpyLines(range(-20,100,10))
HL.plot(ax)

plt.show()

fp = open('PsychScript.py','w')
for chunk in [import_template,SL,RHL,HL,closure_template]:
    fp.write(str(chunk))
fp.close()
execfile('PsychScript.py')
