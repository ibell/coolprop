from CoolProp import CoolProp as cp

#Print out the help files
cp.Help()

print 'The functions available are:'
print dir(cp)

print 'The critical temperature of R410A is: %g K' %(cp.Tcrit('R410A'),)
print 'The saturated vapor enthalpy of Propane at 275 K is: %g kJ/kg' %(cp.Props('H','T',275.0,'Q',1.0,'R290'))

rho_N2=cp.Props('H','T',298.0,'P',101.325,'Nitrogen')
print 'The density of nitrogen at STP is: %g kg/m^3' %(rho_N2)
