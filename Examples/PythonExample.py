import sys,os
sys.path.append('..')
if 'MYPYTHONHOME' in os.environ and os.environ['MYPYTHONHOME'] in sys.path:
    from CoolProp import CoolProp as cp
else:
    import CoolProp as cp

print 'The functions available are:'
print dir(cp)

print 'The critical temperature of R410A is: %g' %(cp.Tcrit('R410A'),)
print 'The saturated vapor enthalpy of Propane at 275 K is: %g kJ/kg' %(cp.Props('H','T',275.0,'Q',1.0,'R290'),)
