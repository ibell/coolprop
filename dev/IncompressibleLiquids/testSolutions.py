#!/usr/bin/python
import sys
import CoolProp.CoolProp as CP
#print "{0:14.8f}".format(CP.Props('V','D',13,'P',500,'n-Pentane'))
#print "{0:14.8f}".format(CP.Props('V','H',158,'P',1000,'TX22'))
#T = 300
T = float(sys.argv[1])+273.15
P = float(sys.argv[2])*1e5
print "Temperature: "+str(T-273.15)+" C"
print "Pressure:    "+str(P/1e5)+" bar"
print 
Melinder = ["MEG", "MPG", "MEA", "MMA", "MGL", "MAM", "MKC", "MCA", "MMG", "MNA", "MKA", "MKF", "MLI"]
SecCool = ["ZitrecAC", "IceSlurryEA", "IceSlurryPG", "IceSlurryNA"]

fluids = []
fluids.extend(Melinder)
fluids.extend(SecCool)

for fluid in fluids:
    print "Fluid: "+str(fluid)
    try: 
        print "Density:    "+"{0:14.8f} kg/m3  ".format(CP.PropsU('D','T',T,'P',P,fluid+'-20%','SI'))
        print "Heat cap.:  "+"{0:14.8f} kJ/kg/K".format(CP.PropsU('C','T',T,'P',P,fluid+'-20%','SI')/1e3)
        print "Th. cond.:  "+"{0:14.8f} W/m/K  ".format(CP.PropsU('L','T',T,'P',P,fluid+'-20%','SI'))
        print "Dyn. visc.: "+"{0:14.8f} mPas   ".format(CP.PropsU('V','T',T,'P',P,fluid+'-20%','SI')*1e3)
        print "Entropy:    "+"{0:14.8f} kJ/kg/K".format(CP.PropsU('S','T',T,'P',P,fluid+'-20%','SI')/1e3)
        print "Freezing:   "+"{0:14.8f} C      ".format(CP.PropsU('Tfreeze','T',T,'P',P,fluid+'-20%','SI')-273.15)
    except ValueError as ve:
        print "Error in CoolProp, try adjusting T and p:"
        print ve
    print 
