import numpy


class LiquidData(object):
    """ 
    A base class that defines all the variables needed 
    in order to make a proper fit. You can copy this code
    put in your data and some documentation for where the
    information came from. 
    """
    Tmin     = None # Minimum temperature in K
    TminPsat = None # Minimum saturation temperature in K
    Tmax     = None # Maximum temperature in K
    T        = None # Temperature for data points in K
    rho      = None # Density in kg/m3
    c_p      = None # Heat capacity in J/(kg.K)
    lam      = None # Thermal conductivity in W/(m.K)
    mu_dyn   = None # Dynamic viscosity in Pa.s
    psat     = None # Saturation pressure in Pa
    


class Example(LiquidData):
    """ 
    Heat transfer fluid TherminolD12 by Solutia 
    """
    Tmin    = -85 + 273.15
    TminPsat=  40 + 273.15
    Tmax    = 260 + 273.15    
    T       = numpy.array([   50,    60,    70,    80,    90,   100,   110,   120,   130,   140,   150])+273.15
    rho     = numpy.array([  740,   733,   726,   717,   710,   702,   695,   687,   679,   670,   662])
    c_p     = numpy.array([ 2235,  2280,  2326,  2361,  2406,  2445,  2485,  2528,  2571,  2607,  2645])
    lam     = numpy.array([0.105, 0.104, 0.102, 0.100, 0.098, 0.096, 0.095, 0.093, 0.091, 0.089, 0.087])
    mu_dyn  = numpy.array([0.804, 0.704, 0.623, 0.556, 0.498, 0.451, 0.410, 0.374, 0.346, 0.317, 0.289])/1000. 
    psat    = numpy.array([  0.5,   0.9,   1.4,   2.3,   3.9,   6.0,   8.7,  12.4,  17.6,  24.4,  33.2])*1000.
    
    
class TherminolD12(LiquidData):
    """ 
    Heat transfer fluid Therminol D12 by Solutia
    Source: http://twt.mpei.ac.ru/MCS/Worksheets/HEDH/.%5C..%5C..%5C..%5CTTHB%5CHEDH%5CHTF-D12.PDF
    """
    Tmin    = -85 + 273.15
    TminPsat=  40 + 273.15
    Tmax    = 260 + 273.15    
    T       = numpy.array([   50,    60,    70,    80,    90,   100,   110,   120,   130,   140,   150])+273.15
    rho     = numpy.array([  740,   733,   726,   717,   710,   702,   695,   687,   679,   670,   662])
    c_p     = numpy.array([ 2235,  2280,  2326,  2361,  2406,  2445,  2485,  2528,  2571,  2607,  2645])
    lam     = numpy.array([0.105, 0.104, 0.102, 0.100, 0.098, 0.096, 0.095, 0.093, 0.091, 0.089, 0.087])
    mu_dyn  = numpy.array([0.804, 0.704, 0.623, 0.556, 0.498, 0.451, 0.410, 0.374, 0.346, 0.317, 0.289])/1000. 
    psat    = numpy.array([  0.5,   0.9,   1.4,   2.3,   3.9,   6.0,   8.7,  12.4,  17.6,  24.4,  33.2])*1000.
    
    
    
    
    
    