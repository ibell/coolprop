import numpy
from scipy.optimize._minimize import minimize

class IncompLiquidFit(object):
    """ 
    A class for fitting data sheet data to predefined functions. 
    Some functions are only used during the fitting procedure.
    Note that the order in which you fit the different properties 
    might impact the coefficients. Usually, the fitting order should be:
    1) Density
    2) Heat capacity
    3) Thermal conductivity
    4) Viscosity
    """
    
    def __init__(self):
        self.DEBUG = False
        
        # parameters for the different fits
        self._cDensity = numpy.ones(4)       # Typically 4 parameters
        self._cHeatCapacity = numpy.ones(5)  # Typically 5 parameters
        self._cTConductivity = numpy.ones(3) # Typically 3 parameters
        self._cViscosity = numpy.ones(3)     # Typically 3 parameters
        self._cPsat = numpy.ones(3)          # Typically 3 parameters
        
        # bounds for fit
        self._Tmin = None
        self._TminPsat = None 
        self._Tmax = None 
        
        # some flags to set
        self._TinC = False   # Temperature in Celsius
        self._DynVisc = True # Data for dynamic viscosity
        self._minPoints = 5 
        
        
    def setParams(self,fluid):
        if fluid=='init':
            # initial parameters for the different fits
            self._cDensity =        [+7.5e+2, -0.7e+0, -1.3e-4, -2.1e-6]
            self._cHeatCapacity =   [+2.0e+0, +3.8e-3, +2.1e-6, -1.1e-8, +3.9e-11]
            self._cTConductivity =  [+1.1e-1, +1.5e-4, -1.6e-7]
            self._cViscosity =      [+5.3e+2, +1.5e+2, -2.7e+0]
            self._cPsat =           [-3.6e+3, +1.9e+2, +1.4e+1]
            
        elif fluid=='TherminolD12inCelsius':
            self._cDensity =        [776.257 ,  -0.696982, -0.000131384, -0.00000209079]
            self._cHeatCapacity =   [2.01422 , 0.00386884,   2.05029e-6, -1.12621e-8, 3.86282e-11]
            self._cTConductivity =  [0.112994, 0.00014781,  -1.61429e-7]
            self._cViscosity =      [530.944, 146.4, -2.68168]
            self._cPsat =           [-3562.69, 194, 13.8526]
            self._Tmin =     -85.0 + 273.15
            self._TminPsat =  40.0 + 273.15
            self._Tmax =     260.0 + 273.15
        elif fluid=='TherminolD12':
            self._cDensity =        [1.08315084e+04,-8.21176568e+01,2.23399244e-01,  -2.03753274e-04]
            self._cHeatCapacity =   [2.01422 , 0.00386884,   2.05029e-6, -1.12621e-8, 3.86282e-11]
            self._cTConductivity =  [0.112994, 0.00014781,  -1.61429e-7]
            self._cViscosity =      [530.944, 146.4, -2.68168]
            self._cPsat =           [-3562.69, 194, 13.8526]
            self._Tmin =     -85.0 + 273.15
            self._TminPsat =  40.0 + 273.15
            self._Tmax =     260.0 + 273.15

        else:
            raise (ValueError("No coefficients available for "+str(fluid)))
        
        
    def _checkT(self,T=0):
        Tmin = self.Props('Tmin')
        Tmax = self.Props('Tmax')
        if Tmin is None:
            raise (ValueError("Please specify the minimum temperature."))
        if Tmax is None:
            raise (ValueError("Please specify the maximum temperature."))
        if not (Tmin<=T<=Tmax):
            raise (ValueError("Temperature out of range: "+str(T)+" not in "+str(Tmin)+"-"+str(Tmax)+". "))
        
    def _checkP(self,T=0,P=0):
        Psat = self.Props('Psat',T=T)
        if P<Psat:
            raise (ValueError("Equations are valid for liquid phase only: "+str(P)+" < "+str(Psat)+". "))
    
    def _checkTP(self,T=0,P=0):
        self._checkT(T=T)
        self._checkP(T=T, P=P)
        
        
    def _basePolynomial(self,coefficients,x):
        """ Base function to produce polynomials of 
        order len(coefficients) with the coefficients
        """
        result = 0.
        for i in range(len(coefficients)): 
            result = result + coefficients[i] * x**i
        return result 

    def _baseExponential(self,coefficients,x,num=3):
        """ Base function to produce exponential 
        with defined coefficients
        """
        if len(coefficients)==num:
            if num==3: return numpy.exp(coefficients[0]/(x+coefficients[1]) - coefficients[2])
        else:
            print "Error!"


    def Props(self,out,T=0,P=0):
        if out=='D':
            self._checkTP(T=T,P=P)
            return self._basePolynomial(self._cDensity,T)
        elif out=='C':
            self._checkTP(T=T,P=P)
            return self._basePolynomial(self._cHeatCapacity,T)
        elif out=='L':
            self._checkTP(T=T,P=P)
            return self._basePolynomial(self._cThermalConductivity,T)
        elif out=='V':
            self._checkTP(T=T,P=P)
            return self._baseExponential(self._cViscosity,T)
#            visc = self._baseExponential(self._cViscosity,T)
#            if not self._DynVisc: return visc * self.Props('D',T=T,P=P)
#            else: return visc
        elif out=='Psat':
            self._checkT(T=T)
            if T<self._TminPsat:
                return 1e-14
            return self._baseExponential(self._cPsat,T)
        elif out=='Tmin':
            return self._Tmin
        elif out=='Tmax':
            return self._Tmax
        else:
            raise (ValueError("Error: You used an unknown output qualifier."))
        
        
    def _PropsFit(self,coefficients,inVal,T=0):
        """
        Calculates a property from a given set of 
        coefficents for a certain temperature. Is used
        to obtain data to feed to the optimisation
        procedures.
        """
        if inVal=='D':
            self._checkT(T=T)
            return self._basePolynomial(coefficients,T)
        elif inVal=='C':
            self._checkT(T=T)
            return self._basePolynomial(coefficients,T)
        elif inVal=='L':
            self._checkT(T=T)
            return self._basePolynomial(coefficients,T)
        elif inVal=='V':
            self._checkT(T=T)
            return self._baseExponential(coefficients,T)
        elif inVal=='Psat':
            self._checkT(T=T)
            if T<self._TminPsat:
                return 1e-14
            return self._baseExponential(coefficients,T)
        else:
            raise (ValueError("Error: You used an unknown property qualifier."))
        
        
    def getCoefficients(self,inVal):
        """
        Get the array with coefficients.
        """
        if inVal=='D':
            return self._cDensity
        elif inVal=='C':
            return self._cHeatCapacity
        elif inVal=='L':
            return self._cTConductivity
        elif inVal=='V':
            return self._cViscosity
        elif inVal=='Psat':
            return self._cPsat
        else:
            raise (ValueError("Error: You used an unknown property qualifier."))
        
    def setCoefficients(self,inVal,coeffs):
        """
        Set the array of coefficients.
        """
        if inVal=='D':
            self._cDensity = coeffs
        elif inVal=='C':
            self._cHeatCapacity = coeffs
        elif inVal=='L':
            self._cTConductivity = coeffs
        elif inVal=='V':
            self._cViscosity = coeffs
        elif inVal=='Psat':
            self._cPsat = coeffs
        else:
            raise (ValueError("Error: You used an unknown property qualifier."))
        
    def setTmin(self,T):
        self._Tmin = T
        
    def setTmax(self,T):
        self._Tmax = T
    
    def setTminPsat(self,T):
        self._TminPsat = T
        
        
    def fitCoefficients(self,xName,T=[],xData=[]):
        
        if (len(T)!=len(xData)):
            raise (ValueError("Error: There has to be the same number of temperature and data points."))
        if len(T)<self._minPoints:
            raise (ValueError("Error: You should use at least "+str(self._minPoints)+" points."))
        
        def fun(coefficients,xName,T,xData):
            calculated = numpy.array([self._PropsFit(coefficients,xName,T=Ti) for Ti in T])
            data = numpy.array(xData)
            res = numpy.sum((calculated-data)**2.)
            return res 
        
        initValues = self.getCoefficients(xName)[:]
        arguments  = (xName,T,xData)
        options    = {'maxiter': 1e4, 'maxfev': 1e6}
        res = minimize(fun, initValues, method='Powell', args=arguments, options=options)
        if res.success:
            return res.x
        else:
            print res
            return False 
        
    
### Some test case 

liqObj = IncompLiquidFit()
liqObj.setParams("init")
liqObj.setTmin(    -85 + 273.15)
liqObj.setTminPsat( 40 + 273.15)
liqObj.setTmax(    260 + 273.15)    

T      = numpy.array([   50,    60,    70,    80,    90,   100,   110,   120,   130,   140,   150])+273.15
rho    = numpy.array([  740,   733,   726,   717,   710,   702,   695,   687,   679,   670,   662])
c_p    = numpy.array([ 2235,  2280,  2326,  2361,  2406,  2445,  2485,  2528,  2571,  2607,  2645])
lam    = numpy.array([0.105, 0.104, 0.102, 0.100, 0.098, 0.096, 0.095, 0.093, 0.091, 0.089, 0.087])
mu_dyn = numpy.array([0.804, 0.704, 0.623, 0.556, 0.498, 0.451, 0.410, 0.374, 0.346, 0.317, 0.289])/1000. 
psat   = numpy.array([  0.5,   0.9,   1.4,   2.3,   3.9,   6.0,   8.7,  12.4,  17.6,  24.4,  33.2])*1000.  

inVal = 'D'
xData = rho
oldCoeffs = liqObj.getCoefficients(inVal)
newCoeffs = liqObj.fitCoefficients(inVal,T=T,xData=xData)
print "Density, old: "+str(oldCoeffs)
print "Density, new: "+str(newCoeffs)
print 

inVal = 'C'
xData = c_p
oldCoeffs = liqObj.getCoefficients(inVal)
newCoeffs = liqObj.fitCoefficients(inVal,T=T,xData=xData)
print "Heat c., old: "+str(oldCoeffs)
print "Heat c., new: "+str(newCoeffs)
print 

inVal = 'L'
xData = lam
oldCoeffs = liqObj.getCoefficients(inVal)
newCoeffs = liqObj.fitCoefficients(inVal,T=T,xData=xData)
print "Th. Co., old: "+str(oldCoeffs)
print "Th. Co., new: "+str(newCoeffs)
print 

inVal = 'V'
xData = mu_dyn
oldCoeffs = liqObj.getCoefficients(inVal)
newCoeffs = liqObj.fitCoefficients(inVal,T=T,xData=xData)
print "Viscos., old: "+str(oldCoeffs)
print "Viscos., new: "+str(newCoeffs)
print 