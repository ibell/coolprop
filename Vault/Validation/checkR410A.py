import sys
sys.path.append('../')
from numpy import linspace
from CoolProp.CoolProp import Props, PrintSaturationTable
from scipy import optimize
import numpy as np


PrintSaturationTable('R134asat.csv','R134a',290,315)
PrintSaturationTable('R290sat.csv','R290',290,315)
PrintSaturationTable('REFPROP-R290sat.csv','REFPROP-Propane',290,315)
Tv=linspace(290,315,26)

for i in range(len(Tv)):
    p1=Props('P','T',Tv[i],'Q',1.0,'REFPROP-Propane')
    p2=Props('P','T',Tv[i],'Q',1.0,'R290')
    d1=Props('U','T',Tv[i]+1,'P',p1,'REFPROP-Propane')
    d2=Props('U','T',Tv[i]+1,'P',p2,'R290')
    print d1,d2
    
