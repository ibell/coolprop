import numpy, matplotlib.pyplot
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
    5) Vapour pressure
    """
    
    def __init__(self):
        self.DEBUG = False
        
        # parameters for the different fits
        self._cDensity = numpy.ones(4)       # Typically 4 parameters
        self._cHeatCapacity = numpy.ones(4)  # Typically 4 parameters
        self._cTConductivity = numpy.ones(3) # Typically 3 parameters
        self._cViscosity = numpy.ones(3)     # Typically 3 parameters
        self._cPsat = numpy.ones(3)          # Typically 3 parameters
        
        # bounds for fit
        self._Tmin = None
        self._TminPsat = None 
        self._Tmax = None 
        self._Tref = 273.15 + 25. 
        
        # some flags to set
        self._TinC = False   # Temperature in Celsius
        self._DynVisc = True # Data for dynamic viscosity
        self._minPoints = 3 
        
        
    def setParams(self,fluid):
        if fluid=='init':
            # initial parameters for the different fits
            self._cDensity =        [+9.2e+2, -0.5e+0, +2.8e-4, -1.1e-6]
            self._cHeatCapacity =   [+1.0e+0, +3.6e-3, -2.9e-7, +1.7e-9]
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
    
    def _basePolynomialInt(self,coefficients,x1,x0=-1):
        """ Base function to produce the integral of
        order len(coefficients) with coefficients from
        x0 to x1.
        """
        if x0==-1: x0 = self._Tref
        result = 0. 
        for i in range(len(coefficients)): 
            result += 1./(i+1.) * coefficients[i] * (x1**(i+1.) - x0**(i+1.))
        return result 

    def _baseExponential(self,coefficients,x,num):
        """ Base function to produce exponential 
        with defined coefficients
        """
        #if len(coefficients)==num:
        if num==1: return numpy.exp(coefficients[0]/(x+coefficients[1]) - coefficients[2])
        if num==2: return numpy.exp(self._basePolynomial(coefficients, x))
        #else:
        #    print "Error!"


    def Props(self,out,T=0,P=0):
        if out=='D':
            self._checkTP(T=T,P=P)
            return self._basePolynomial(self._cDensity,T)
        elif out=='C':
            self._checkTP(T=T,P=P)
            return self._basePolynomial(self._cHeatCapacity,T)
        elif out=='L':
            self._checkTP(T=T,P=P)
            return self._basePolynomial(self._cTConductivity,T)
        elif out=='V':
            self._checkTP(T=T,P=P)
            return self._baseExponential(self._cViscosity,T,1)
#            visc = self._baseExponential(self._cViscosity,T)
#            if not self._DynVisc: return visc * self.Props('D',T=T,P=P)
#            else: return visc
        elif out=='Psat':
            self._checkT(T=T)
            if T<self._TminPsat:
                return 1e-14
            return self._baseExponential(self._cPsat,T,1)
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
            return self._baseExponential(coefficients,T,1)
        elif inVal=='Psat':
            self._checkT(T=T)
            if T<self._TminPsat:
                return 1e-14
            return self._baseExponential(coefficients,T,1)
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
    
    def setTref(self,T):
        self._Tref = T    
        
    def fitCoefficients(self,xName,T=[],xData=[]):
        
        if (len(T)!=len(xData)):
            raise (ValueError("Error: There has to be the same number of temperature and data points."))
        if len(T)<self._minPoints:
            raise (ValueError("Error: You should use at least "+str(self._minPoints)+" points."))
        
        def fun(coefficients,xName,T,xData):
            # Values for conductivity are very small,
            # algorithms prefer larger values
            if xName=='L':
                calculated = numpy.array([self._PropsFit(coefficients,xName,T=Ti) for Ti in T])*1e6
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
        options    = {'maxiter': 1e3, 'maxfev': 1e5}
        res = minimize(fun, initValues, method='Powell', args=arguments, tol=1e-15, options=options)
        if res.success:
            return res.x
        else:
            print res
            return False 


### Load the data 
from data_incompressible import SylthermXLT as DataContainer
data = DataContainer()


### Some test case 
liqObj = IncompLiquidFit()
liqObj.setParams("init")
liqObj.setTmin(data.Tmin)
liqObj.setTminPsat(data.TminPsat)
liqObj.setTmax(data.Tmax)    

numpy.set_printoptions(formatter={'float': lambda x: format(x, '+1.10E')})

print "minimum T: "+str(data.Tmin)
print "maximum T: "+str(data.Tmax)
print "min T pSat:"+str(data.TminPsat)
print 

# row and column sharing for test plots
f, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = matplotlib.pyplot.subplots(3, 2, sharex='col')
f.set_size_inches(matplotlib.pyplot.figaspect(1.2)*1.5)

### This is the actual fitting
tData = data.T
tDat1 = numpy.linspace(numpy.min(tData), numpy.max(tData), 100)
Pin = 1e20 # Dummy pressure

inVal = 'D'
xData = data.rho
oldCoeffs = liqObj.getCoefficients(inVal)
newCoeffs = liqObj.fitCoefficients(inVal,T=tData,xData=xData)
print "Density, old: "+str(oldCoeffs)
print "Density, new: "+str(newCoeffs)
print
liqObj.setCoefficients(inVal,newCoeffs)
fData = numpy.array([liqObj.Props(inVal, T=Tin, P=Pin) for Tin in tDat1])
ax1.plot(tData-273.15, xData, 'o', label="Data Sheet")
ax1.plot(tDat1-273.15, fData, label="CoolProp")
ax1.set_ylabel(r'$\mathregular{Density\/(kg\/m^{-3})}$')


inVal = 'C'
xData = data.c_p
oldCoeffs = liqObj.getCoefficients(inVal)
newCoeffs = liqObj.fitCoefficients(inVal,T=tData,xData=xData)
print "Heat c., old: "+str(oldCoeffs)
print "Heat c., new: "+str(newCoeffs)
print 
liqObj.setCoefficients(inVal,newCoeffs)
fData = numpy.array([liqObj.Props(inVal, T=Tin, P=Pin) for Tin in tDat1])
ax2.plot(tData-273.15, xData, 'o', label="Data Sheet")
ax2.plot(tDat1-273.15, fData, label="CoolProp")
ax2.set_ylabel(r'$\mathregular{Heat\/Cap.\/(kJ\/kg^{-1}\/K^{-1})}$')


inVal = 'L'
xData = data.lam
oldCoeffs = liqObj.getCoefficients(inVal)
newCoeffs = liqObj.fitCoefficients(inVal,T=tData,xData=xData)
print "Th. Co., old: "+str(oldCoeffs)
print "Th. Co., new: "+str(newCoeffs)
print 
liqObj.setCoefficients(inVal,newCoeffs)
fData = numpy.array([liqObj.Props(inVal, T=Tin, P=Pin) for Tin in tDat1])
ax3.plot(tData-273.15, xData*1e6 , 'o', label="Data Sheet")
ax3.plot(tDat1-273.15, fData*1e6 , label="CoolProp")
ax3.set_ylabel(r'$\mathregular{Th.\/Cond.\/(mW\/m^{-1}\/K^{-1})}$')


inVal = 'V'
tData = data.T[data.mu_dyn > 0]
if len(tData>1):
    tDat1 = numpy.linspace(numpy.min(tData), numpy.max(tData), 100)
    xData = data.mu_dyn[data.mu_dyn > 0]
    oldCoeffs = liqObj.getCoefficients(inVal)
    newCoeffs = liqObj.fitCoefficients(inVal,T=tData,xData=xData)
    print "Viscos., old: "+str(oldCoeffs)
    print "Viscos., new: "+str(newCoeffs)
    print 
    liqObj.setCoefficients(inVal,newCoeffs)
    fData = numpy.array([liqObj.Props(inVal, T=Tin, P=Pin) for Tin in tDat1])
    ax4.plot(tData-273.15, xData*1e3, 'o', label="Data Sheet")
    ax4.plot(tDat1-273.15, fData*1e3, label="CoolProp")

ax4.set_ylabel(r'$\mathregular{Dyn.\/Viscosity\/(mPa\/s)}$')
ax4.set_yscale('log')


inVal = 'Psat'
tData = data.T[data.psat > 0]
if len(tData>1):
    tDat1 = numpy.linspace(numpy.min(tData), numpy.max(tData), 100)
    xData = data.psat[data.psat > 0]
    oldCoeffs = liqObj.getCoefficients(inVal)
    newCoeffs = liqObj.fitCoefficients(inVal,T=tData,xData=xData)
    print "P sat. , old: "+str(oldCoeffs)
    print "P sat. , new: "+str(newCoeffs)
    print 
    liqObj.setCoefficients(inVal,newCoeffs)
    fData = numpy.array([liqObj.Props(inVal, T=Tin, P=Pin) for Tin in tDat1])
    ax5.plot(tData-273.15, xData, 'o', label="Data Sheet")
    ax5.plot(tDat1-273.15, fData, label="CoolProp")
    
ax5.set_ylabel(r'$\mathregular{Vap.\/Pressure\/(kPa)}$')
ax5.set_yscale('log')

ax5.set_xlabel(ur'$\mathregular{Temperature\/(\u00B0C)}$')
ax6.set_xlabel(ur'$\mathregular{Temperature\/(\u00B0C)}$')

ax1.legend(loc=3)
matplotlib.pyplot.tight_layout()
matplotlib.pyplot.savefig(data.Name+"_std.pdf")

### Print the output for the C++ file
print "name = std::string(\"\");"
print "description = std::string(\""+data.Name+"\");"
print "reference = std::string(\"\");"
print ""
print "Tmin     = "+str(data.Tmin)+";"
print "Tmax     = "+str(data.Tmax)+";"
print "TminPsat = "+str(data.TminPsat)+";"
print ""
print "cRho.clear();"
C = liqObj.getCoefficients('D')
for Ci in C:
    print "cRho.push_back(%+1.10E);" %(Ci)
    
print ""
print "cHeat.clear();"
C = liqObj.getCoefficients('C')
for Ci in C:
    print "cHeat.push_back(%+1.10E);" %(Ci)

print ""
print "cCond.clear();"
C = liqObj.getCoefficients('L')
for Ci in C:
    print "cCond.push_back(%+1.10E);" %(Ci)

print ""
print "cVisc.clear();"
C = liqObj.getCoefficients('V')
for Ci in C:
    print "cVisc.push_back(%+1.10E);" %(Ci)

print ""
print "cPsat.clear();"
C = liqObj.getCoefficients('Psat')
for Ci in C:
    print "cPsat.push_back(%+1.10E);" %(Ci)


