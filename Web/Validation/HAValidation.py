from CoolProp.HumidAirProp import HAProps


s5=' '*5
print "{W:10s}{Twb:10s}{v:10s}{h:10s}{s:10s}{R:10s}".format(W=s5+' W',Twb=s5+'Twb',v=s5+'v',h=s5+'h',s=s5+'s',R=s5+'RH')
for W in [0.0,0.05,0.1,0.20,0.30,0.40,0.50,0.60,0.70,0.80,0.90,1.0]:
    h=HAProps('H','T',473.13,'W',W,'P',101.325)
    Twb=HAProps('Twb','T',473.13,'W',W,'P',101.325)-273.15
    R=HAProps('R','T',473.13,'W',W,'P',101.325)*100
    v=HAProps('V','T',473.13,'W',W,'P',101.325)
    s=0
    print "{W:10.2f}{Twb:10.2f}{v:10.3f}{h:10.2f}{s:10.4f}{R:10.4f}".format(W=W,Twb=Twb,v=v,h=h,s=s,R=R)
    
s5=' '*5
print "{W:10s}{Twb:10s}{v:10s}{h:10s}{s:10s}{R:10s}".format(W=s5+' W',Twb=s5+'Twb',v=s5+'v',h=s5+'h',s=s5+'s',R=s5+'RH')
for W in [0.0,0.05,0.1,0.20,0.30,0.40,0.50,0.60,0.70,0.80,0.90,1.0]:
    h=HAProps('H','T',473.13,'W',W,'P',2000)
    Twb=HAProps('Twb','T',473.13,'W',W,'P',2000)-273.15
    R=HAProps('R','T',473.13,'W',W,'P',2000)*100
    v=HAProps('V','T',473.13,'W',W,'P',2000)
    s=HAProps('S','T',473.13,'W',W,'P',2000)
    print "{W:10.2f}{Twb:10.2f}{v:10.3f}{h:10.2f}{s:10.4f}{R:10.4f}".format(W=W,Twb=Twb,v=v,h=h,s=s,R=R)

print '--'
s5=' '*5
print "{W:10s}{Twb:10s}{v:10s}{h:10s}{s:10s}{R:10s}".format(W=s5+' W',Twb=s5+'Twb',v=s5+'v',h=s5+'h',s=s5+'s',R=s5+'RH')
for W in [0.0,0.05,0.1,0.15,0.20,0.25,0.30]:
    h=HAProps('H','T',473.13,'W',W,'P',5000)
    Twb=HAProps('Twb','T',473.13,'W',W,'P',5000)-273.15
    R=HAProps('R','T',473.13,'W',W,'P',5000)*100
    v=HAProps('V','T',473.13,'W',W,'P',5000)
    s=HAProps('S','T',473.13,'W',W,'P',5000)
    print "{W:10.2f}{Twb:10.2f}{v:10.3f}{h:10.2f}{s:10.4f}{R:10.4f}".format(W=W,Twb=Twb,v=v,h=h,s=s,R=R)
    
print '--'
s5=' '*5
print "{W:10s}{Twb:10s}{v:10s}{h:10s}{s:10s}{R:10s}".format(W=s5+' W',Twb=s5+'Twb',v=s5+'v',h=s5+'h',s=s5+'s',R=s5+'RH')
for W in [0.0,0.05,0.1]:
    h=HAProps('H','T',473.13,'W',W,'P',10000)
    Twb=HAProps('Twb','T',473.13,'W',W,'P',10000)-273.15
    R=HAProps('R','T',473.13,'W',W,'P',10000)*100
    v=HAProps('V','T',473.13,'W',W,'P',10000)
    s=HAProps('S','T',473.13,'W',W,'P',10000)
    print "{W:10.2f}{Twb:10.2f}{v:10.3f}{h:10.2f}{s:10.4f}{R:10.4f}".format(W=W,Twb=Twb,v=v,h=h,s=s,R=R)

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
####### Water Saturation #####
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
Tv=np.linspace(0,300,11)+273.16
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
print "Molar volume of saturated liquid water or ice (vbar_ws) [m^3/mol_H2O]"
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
variables="%-10s"%(u'T')
for p in Pv:
    variables+="%-20s"%("p = %-0.3f kPa "%(p))
print variables
#Build the actual table
for T in Tv:
    values="%-10.2f" %(T-273.15)
    for p in Pv:
        values+="%-20.10e" %(HAProps_Aux('f',T,p,0.0)[0])
    print values