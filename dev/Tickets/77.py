import matplotlib
matplotlib.use('Qt4Agg')

import numpy 
import CoolProp.CoolProp as CP


print "Testing the derivatives and store results:"

keys = ["dpdT","dpdrho","Z","dZ_dDelta","dZ_dTau","VB","dBdT","VC","dCdT","phir","dphir_dTau","d2phir_dTau2","dphir_dDelta",
  "d2phir_dDelta2","d2phir_dDelta_dTau","d3phir_dDelta2_dTau","phi0","dphi0_dTau","d2phi0_dTau2","dphi0_dDelta","d2phi0_dDelta2",
  "IsothermalCompressibility"]

fluid = "n-Pentane"

T = numpy.array([ 30,100,150,210])+273.15
D = numpy.array([710, 20, 20,210])

for key in keys:
    liquid    = CP.DerivTerms(key,T[0],D[0],fluid)
    twophase  = CP.DerivTerms(key,T[1],D[1],fluid)
    gaseous   = CP.DerivTerms(key,T[2],D[2],fluid)
    supercrit = CP.DerivTerms(key,T[3],D[3],fluid)
    print '{:<25}: {:>10.5f}; {:>10.5f}; {:>10.5f}; {:>10.5f}'.format(key,liquid,twophase,gaseous,supercrit)
    
