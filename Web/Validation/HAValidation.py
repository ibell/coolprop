


##############################
#### Virial Coefficients #####
##############################

def Virials(variables):
    from CoolProp.HumidAirProp import HAProps_Aux
    import numpy as np
    
    varString="%-10s"%('T')
    units="%-10s"%('C')
    #Build the header
    for var in variables:
        varString+="%-20s"%(var)
        units+="%-20s" %(HAProps_Aux(var,300,100,0.0)[1])
    print varString
    print units

    #Build the table
    for T in np.linspace(-60,200,27)+273.15:
        values="%-10.1f" %(T-273.15)
        for var in variables:
            values+="%-20.10e" %(HAProps_Aux(var,T,100,0.0)[0])
        print values
print ""
print "Pure fluid Virial Coefficients"
print "------------------------------"
Virials(['Baa','Caaa','Bww','Cwww'])
Virials(['Baw','Caaw','Caww'])

##############################
####### Henry Constant #######
##############################

print ""
print "Water saturation pressure p_ws [kPa]"
from CoolProp.HumidAirProp import HAProps_Aux
import numpy as np
Tv=np.linspace(-60,300,13)+273.15
print "%-10s %-20s"%('T','p_ws')
print "%-10s %-20s"%('C',HAProps_Aux('p_ws',Tv[-1],100,0.0)[1])
#Build the table
for T in Tv:
    values="%-10.2f" %(T-273.15)
    values+="%-20.10e" %(HAProps_Aux('p_ws',T,100,0.0)[0])
    print values
    
##############################
####### Henry Constant #######
##############################

print ""
print "Henry Constant (zero for T < 273.15 K)"
from CoolProp.HumidAirProp import HAProps_Aux
import numpy as np
Tv=np.linspace(0,300,11)+273.15
print "%-10s %-20s"%('T','beta_H')
print "%-10s %-20s"%('C',HAProps_Aux('beta_H',Tv[-1],100,0.0)[1])
#Build the table
for T in Tv:
    values="%-10.2f" %(T-273.15)
    values+="%-20.10e" %(HAProps_Aux('beta_H',T,100,0.0)[0])
    print values

##########################################
####### Isothermal Compressibility #######
##########################################

print ""
print "Isothermal Compressibility of water (kT) [1/Pa]"
from CoolProp.HumidAirProp import HAProps_Aux
import numpy as np
Tv=np.linspace(-60,300,13)+273.15
Pv=[101.325,200,500,1000]
variables="%-10s"%('T')
for p in Pv:
    variables+="%-20s"%("p = %-0.3f kPa "%(p))
print variables
#Build the actual table
for T in Tv:
    values="%-10.2f" %(T-273.15)
    for p in Pv:
        values+="%-20.10e" %(HAProps_Aux('kT',T,p,0.0)[0])
    print values
    
##########################################
####### Saturated Molar Volume Water #####
##########################################

print ""
print "Molar volume of saturated liquid water or ice (vbar_ws) [m^3/mol]"
from CoolProp.HumidAirProp import HAProps_Aux
import numpy as np
Tv=np.linspace(-60,300,13)+273.15
Pv=[101.325,200,500,1000]
variables="%-10s"%('T')
for p in Pv:
    variables+="%-20s"%("p = %-0.3f kPa "%(p))
print variables
#Build the actual table
for T in Tv:
    values="%-10.2f" %(T-273.15)
    for p in Pv:
        values+="%-20.10e" %(HAProps_Aux('vbar_ws',T,p,0.0)[0])
    print values
    
##########################################
########### Enhancement Factor ###########
##########################################

print ""
print "Enhancement factor (f) [no units]"
from CoolProp.HumidAirProp import HAProps_Aux
import numpy as np
Tv=np.linspace(-60,300,13)+273.15
Pv=[101.325,200,500,1000]
variables="%-10s"%('T')
for p in Pv:
    variables+="%-20s"%("p = %-0.3f kPa "%(p))
print variables
#Build the actual table
for T in Tv:
    values="%-10.2f" %(T-273.15)
    for p in Pv:
        values+="%-20.10e" %(HAProps_Aux('f',T,p,0.0)[0])
    print values