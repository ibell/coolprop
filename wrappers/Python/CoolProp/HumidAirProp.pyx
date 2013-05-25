#This file gets directly included in CoolProp.pyx, separate here for cleanness of code

cpdef double HAProps(str OutputName, str Input1Name, double Input1, str Input2Name, double Input2, str Input3Name, double Input3):
    """
    Copyright Ian Bell, 2011 email: ian.h.bell@gmail.com

    The function is called like

    HAProps('H','T',298.15,'P',101.325,'R',0.5)

    which will return the enthalpy of the air for a set of inputs of dry bulb temperature of 25C, atmospheric pressure, and a relative humidity of 50%.

    This function implements humid air properties based on the analysis in ASHRAE RP-1845 which is available online: http://rp.ashrae.biz/page/ASHRAE-D-RP-1485-20091216.pdf

    It employs real gas properties for both air and water, as well as the most accurate interaction parameters and enhancement factors.  The IAPWS-95 formulation for the properties of water is used throughout in preference to the industrial formulation.  It is unclear why the industrial formulation is used in the first place.

    Since humid air is nominally a binary mixture, three variables are needed to fix the state.  At least one of the input parameters must be dry-bulb temperature, relative humidity, dew-point temperature, or humidity ratio.  The others will be calculated.  If the output variable is a transport property (conductivity or viscosity), the state must be able to be found directly - i.e. make sure you give temperature and relative humidity or humidity ratio.  The list of possible input variables are

    ========  ========    ========================================
    String    Aliases     Description
    ========  ========    ========================================
    T         Tdb         Dry-Bulb Temperature [K]
    B         Twb         Wet-Bulb Temperature [K]
    D         Tdp         Dew-Point Temperature [K]
    P                     Pressure [kPa]
    V         Vda         Mixture volume [m3/kg dry air]
    R         RH          Relative humidity in (0,1) [-]
    W         Omega       Humidity Ratio [kg water/kg dry air]
    H         Hda         Mixture enthalpy [kJ/kg dry air]
    S         Sda         Mixture entropy [kJ/kg dry air/K]
    C         cp          Mixture specific heat [kJ/kg dry air/K]
    M         Visc        Mixture viscosity [Pa-s]
    K                     Mixture thermal conductivity [W/m/K]
    ========  ========    ========================================

    There are also strings for the mixture volume and mixture enthalpy that will return the properties on a total humid air flow rate basis, they are given by 'Vha' [units of m^3/kg humid air] and 'Cha' [units of kJ/kg humid air/K] and 'Hha' [units of kJ/kg humid air] respectively.

    For more information, go to http://coolprop.sourceforge.net
    """
    #Convert all strings to byte-strings
    cdef bytes _OutputName = OutputName.encode('ascii')
    cdef bytes _Input1Name = Input1Name.encode('ascii')
    cdef bytes _Input2Name = Input2Name.encode('ascii')
    cdef bytes _Input3Name = Input3Name.encode('ascii')
    return _HAProps(_OutputName,_Input1Name,Input1,_Input2Name,Input2,_Input3Name,Input3)

cpdef tuple HAProps_Aux(str OutputName, double T, double p, double w):
    """
    Allows low-level access to some of the routines employed in HumidAirProps

    Returns tuples of the form ``(Value, Units)`` where ``Value`` is the actual value and ``Units`` is a string that describes the units

    The list of possible inputs is

    * Baa [First virial air-air coefficient]
    * Caaa [Second virial air coefficient]
    * Bww [First virial water-water coefficient]
    * Cwww [Second virial water coefficient]
    * Baw [First cross virial coefficient]
    * Caww [Second air-water-water virial coefficient]
    * Caaw [Second air-air-water virial coefficient]
    * beta_H 
    * kT
    * vbar_ws [Molar saturated volume of water vapor]
    * p_ws [Saturated vapor pressure of pure water (>=0.01C) or ice (<0.01 C)]
    * f [Enhancement factor]
    """
    #Convert all strings to byte-strings
    cdef bytes units = (' '*100).encode('ascii')
    cdef bytes _OutputName = OutputName.encode('ascii')
    
    output = _HAProps_Aux(_OutputName,T,p,w,units)
    units = units.strip()
    units = units[0:len(units)-1] #Leave off the null character
    return output, units