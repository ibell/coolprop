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
            self._cDensity =        [+9.2e+2, -0.5e+0, +2.8e-4, -1.1e-6]
            self._cHeatCapacity =   [+1.0e+3, +2.9e+0, +2.5e-3, +3.2e-6, -9.8e-9]
            self._cTConductivity =  [+1.1e-1, +7.8e-5, -3.5e-7]
            self._cViscosity =      [+7.1e+4, +2.3e+3, +3.4e+1]
            self._cPsat =           [-5.3e+5, +3.2e+3, -1.6e+2]
            
#        elif fluid=='TherminolD12inCelsius':
#            self._cDensity =        [776.257 ,  -0.696982, -0.000131384, -0.00000209079]
#            self._cHeatCapacity =   [2.01422 , 0.00386884,   2.05029e-6, -1.12621e-8, 3.86282e-11]
#            self._cTConductivity =  [0.112994, 0.00014781,  -1.61429e-7]
#            self._cViscosity =      [530.944, 146.4, -2.68168]
#            self._cPsat =           [-3562.69, 194, 13.8526]
#            self._Tmin =     -85.0 + 273.15
#            self._TminPsat =  40.0 + 273.15
#            self._Tmax =     260.0 + 273.15
#        elif fluid=='TherminolD12':
#            self._cDensity =        [1.08315084e+04,-8.21176568e+01,2.23399244e-01,  -2.03753274e-04]
#            self._cHeatCapacity =   [2.01422 , 0.00386884,   2.05029e-6, -1.12621e-8, 3.86282e-11]
#            self._cTConductivity =  [0.112994, 0.00014781,  -1.61429e-7]
#            self._cViscosity =      [530.944, 146.4, -2.68168]
#            self._cPsat =           [-3562.69, 194, 13.8526]
#            self._Tmin =     -85.0 + 273.15
#            self._TminPsat =  40.0 + 273.15
#            self._Tmax =     260.0 + 273.15

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
            result += coefficients[i] * x**i
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
            # Values for conductivity are very small,
            # algorithms prefer larger values
            if xName=='L':
                calculated = numpy.array([self._PropsFit(coefficients,xName,T=Ti) for Ti in T])*1e4
                data       = numpy.array(xData)*1e6
            # Fit logarithms for viscosity and saturation pressure
            elif xName=='V' or xName=='Psat':
                calculated = numpy.log([self._PropsFit(coefficients,xName,T=Ti) for Ti in T])
                data       = numpy.log(xData)
            else:
                calculated = numpy.array([self._PropsFit(coefficients,xName,T=Ti) for Ti in T])
                data       = numpy.array(xData)
            
            res = numpy.sum((calculated-data)**2.)
            return res 
        
        initValues = self.getCoefficients(xName)[:]
        arguments  = (xName,T,xData)
        options    = {'maxiter': 1e4, 'maxfev': 1e6}
        res = minimize(fun, initValues, method='Powell', args=arguments, tol=1e-15, options=options)
        if res.success:
            return res.x
        else:
            print res
            return False 


### Load the data 
from data_incompressible import TherminolD12 as DataContainer
data = DataContainer()


### Some test case 
liqObj = IncompLiquidFit()
liqObj.setParams("init")
liqObj.setTmin(data.Tmin)
liqObj.setTminPsat(data.TminPsat)
liqObj.setTmax(data.Tmax)    

### This is the actual fitting
tData = data.T
inVal = 'D'
xData = data.rho
oldCoeffs = liqObj.getCoefficients(inVal)
newCoeffs = liqObj.fitCoefficients(inVal,T=tData,xData=xData)
print "Density, old: "+str(oldCoeffs)
print "Density, new: "+str(newCoeffs)
print 

inVal = 'C'
xData = data.c_p
oldCoeffs = liqObj.getCoefficients(inVal)
newCoeffs = liqObj.fitCoefficients(inVal,T=tData,xData=xData)
print "Heat c., old: "+str(oldCoeffs)
print "Heat c., new: "+str(newCoeffs)
print 

inVal = 'L'
xData = data.lam
oldCoeffs = liqObj.getCoefficients(inVal)
newCoeffs = liqObj.fitCoefficients(inVal,T=tData,xData=xData)
print "Th. Co., old: "+str(oldCoeffs)
print "Th. Co., new: "+str(newCoeffs)
print 

inVal = 'V'
xData = data.mu_dyn
oldCoeffs = liqObj.getCoefficients(inVal)
newCoeffs = liqObj.fitCoefficients(inVal,T=tData,xData=xData)
print "Viscos., old: "+str(oldCoeffs)
print "Viscos., new: "+str(newCoeffs)
print 

inVal = 'Psat'
xData = data.psat
oldCoeffs = liqObj.getCoefficients(inVal)
newCoeffs = liqObj.fitCoefficients(inVal,T=tData,xData=xData)
print "P sat. , old: "+str(oldCoeffs)
print "P sat. , new: "+str(newCoeffs)
print 

# TODO: We should plot some graphs to control fit quality

