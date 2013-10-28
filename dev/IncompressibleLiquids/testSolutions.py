#!/usr/bin/python
import sys
import CoolProp.CoolProp as CP
#print "{0:14.8f}".format(CP.Props('V','D',13,'P',500,'n-Pentane'))
#print "{0:14.8f}".format(CP.Props('V','H',158,'P',1000,'TX22'))
#T = 300
T = float(sys.argv[1])
P = float(sys.argv[2])
print "Temperature: "+str(T)
print "Pressure:    "+str(P)
print 
Melinder = ["MEG", "MPG", "MEA", "MMA", "MGL", "MAM", "MKC", "MCA", "MMG", "MNA", "MKA", "MKF", "MLI"]

for fluid in Melinder:
    print "Entropy:  "+"{0:14.8f}".format(CP.Props('S','T',T,'P',P,fluid+'-20%'))
    print "Freezing: "+"{0:14.8f}".format(CP.Props('Tfreeze','T',T,'P',P,fluid+'-20%')-273.15)
    print 
