#cython: embedsignature = True
#
#
# This file provides wrapper functions of all the CoolProp functions
#
#
# Each of the functions from the CoolProp header are renamed in cython code to
# an underscored name so that the same name can be used in the exposed functions below


cdef extern from "CoolProp.h":
    
    double _Props "Props"(string Output, char Name1, double Prop1, char Name2, double Prop2, string Ref)
    double _IProps "IProps"(long Output, long Name1, double Prop1, long Name2, double Prop2, long Ref)
    double _Props1 "Props"(string Ref, string Output)
    string _Phase "Phase"(char *Fluid, double T, double p)
    double _DerivTerms "DerivTerms"(char *Term, double T, double rho, char * Ref)
    
    # LUT functions
    void _UseSaturationLUT "UseSaturationLUT"(bint OnOff)
    bint _SaturationLUTStatus "SaturationLUTStatus" ()
    void _UseSinglePhaseLUT "UseSinglePhaseLUT"(bint OnOff)
    bint _SinglePhaseLUTStatus "SinglePhaseLUTStatus"()
    
    # Conversion functions
    double _F2K "F2K"(double T_F)
    double _K2F "K2F"(double T)
    
    string _FluidsList "FluidsList"()
    string _get_REFPROPname "get_REFPROPname"(string Ref)
    string _get_errstring "get_errstring"()
    string _get_aliases "get_aliases"(string Ref)
    long _set_phase "set_phase" (string phase)
    long _get_Fluid_index "get_Fluid_index" (string Fluid)
    long _get_param_index "get_param_index" (string param)
    string _get_index_units "get_index_units" (long index)
    char * get_errstringc()
    
    #Ancillary equations
    double _psatL_anc "psatL_anc"(char*Fluid, double T)
    double _psatV_anc "psatV_anc"(char*Fluid, double T)
    double _rhosatL_anc "rhosatL_anc"(char*Fluid, double T)
    double _rhosatV_anc "rhosatV_anc"(char*Fluid, double T)
    
    int _get_debug "get_debug"()
    void _debug "debug"(int level)
    
    void _PrintSaturationTable "PrintSaturationTable"(char *FileName, char * Ref, double Tmin, double Tmax)
    
    string _get_EOSReference "get_EOSReference"(string Ref)
    string _get_TransportReference "get_TransportReference"(string Ref)

    int _set_1phase_LUT_params "set_1phase_LUT_params"(char *Ref, int nT, int np, double Tmin, double Tmax, double pmin, double pmax, bint rebuild)
    void _get_1phase_LUT_params "get_1phase_LUT_params"(int *nT, int *np, double *Tmin, double *Tmax, double *pmin, double *pmax)
    
    # Convenience functions
    int _IsFluidType "IsFluidType"(char *Ref, char *Type)

cdef extern from "HumidAirProp.h":
    double _HAProps "HAProps"(char *OutputName, char *Input1Name, double Input1, char *Input2Name, double Input2, char *Input3Name, double Input3)
    double _HAProps_Aux "HAProps_Aux"(char* Name,double T, double p, double W, char *units)
    
    
#Check for the existence of quantities
cdef bint _quantities_supported
try:
    import quantities as pq
    _quantities_supported = True
except ImportError:
    _quantities_supported = False

import cython
import math

## Variables for each of the possible variables
cdef long iMM = _get_param_index('M')
cdef long iT = _get_param_index('T')
cdef long iD = _get_param_index('D')
cdef long iH = _get_param_index('H')
cdef long iP = _get_param_index('P')
cdef long iC = _get_param_index('C')
cdef long iC0 = _get_param_index('C0')
cdef long iO = _get_param_index('O')
cdef long iV = _get_param_index('V')
cdef long iL = _get_param_index('L')
cdef long iS = _get_param_index('S')
cdef long iU = _get_param_index('U')
cdef long iDpdT = _get_param_index('dpdT')

cpdef double _convert_to_desired_units(double value, bytes parameter_type, bytes desired_units) except *:
    """
    A convenience function to convert a double value back to the units that
    people
    """
    #Convert parameter string to index
    cdef long index = _get_param_index(parameter_type)
    #Get the units string by the index
    cdef bytes default_units = _get_index_units(index)
    #Create the Quantity instance
    old = pq.Quantity(value, default_units)
    #Convert the units
    old.units = desired_units
    #Return the scaled units
    return old.magnitude

cpdef _convert_to_default_units(bytes parameter_type, object parameter):
    """
    A convenience function to convert a quantities instance to the default 
    units required for CoolProp
    """
    #Convert parameter string to index
    cdef long index = _get_param_index(parameter_type)
    #Get the units string by the index
    cdef bytes default_units = _get_index_units(index)
    #Rescale the units of the paramter to the default units
    parameter.units = default_units
    #Return the scaled units
    return parameter
     
cpdef double HAProps(bytes OutputName, bytes InputName1, double Input1, bytes Input2Name, double Input2, bytes Input3Name, double Input3):
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
    C         cp          Mixture specific heat [kJ/kg dry air/K]
    M         Visc        Mixture viscosity [Pa-s]
    K                     Mixture thermal conductivity [W/m/K]
    ========  ========    ========================================

    There are also strings for the mixture volume and mixture enthalpy that will return the properties on a total humid air flow rate basis, they are given by 'Vha' [units of m^3/kg humid air] and 'Cha' [units of kJ/kg humid air/K] and 'Hha' [units of kJ/kg humid air] respectively.

    For more information, go to http://coolprop.sourceforge.net
    """
    return _HAProps(OutputName,InputName1,Input1,Input2Name,Input2,Input3Name,Input3)

cpdef tuple HAProps_Aux(bytes OutputName, double T, double p, double w, bytes units):
    """
    Allows low-level access to some of the routines employed in HumidAirProps

    Returns tuples of the form (Value, Units) where value is the actual value and Units is a string that describes the units

    The list of possible inputs is

    * Baa
    * Caaa
    * Bww
    * Cwww
    * Baw
    * Caww
    * Caaw
    * beta_H
    * kT
    * vbar_ws
    * p_ws
    * f
    """
    output = _HAProps_Aux(OutputName,T,p,w,units)
    return output,units
    
cpdef double Props(bytes in1, bytes in2, in3=None, in4=None,in5=None,in6=None,in7=None) except *:
    """
    Call Type #1::

        Props(Fluid,PropName) --> float

    Where ``Fluid`` is a string with a valid CoolProp fluid name, and ``PropName`` is one of the following strings:
    
    =============  ============================
    ``Tcrit``      Critical temperature [K]
    ``pcrit``      Critical pressure [kPa]
    ``rhocrit``    Critical density [kg/m3]
    ``molemass``   Molecular mass [kg/kmol]
    ``Ttriple``    Triple-point temperature [K]
    ``accentric``  Accentric factor [-]
    =============  ============================
   
    This type of call is used to get fluid-specific parameters that are not 
    dependent on the state 
     
    Call Type #2:
    
    Alternatively, Props can be called in the form::
    
        Props(OutputName,InputName1,InputProp1,InputName2,InputProp2,Fluid) --> float
    
    where ``Fluid`` is a string with a valid CoolProp fluid name.  The value 
    ``OutputName`` is either a single-character or a string alias.  This list 
    shows the possible values
    
    ==========================  ======================================================
    ``OutputName``              Description
    ==========================  ======================================================
    ``Q``                       Quality [-]
    ``T``                       Temperature [K]
    ``P``                       Pressure [kPa]
    ``D``                       Density [kg/m3]
    ``C0``                      Ideal-gas specific heat at constant pressure [kJ/kg]
    ``C``                       Specific heat at constant pressure [kJ/kg]
    ``O``                       Specific heat at constant volume [kJ/kg]
    ``U``                       Internal energy [kJ/kg]
    ``H``                       Enthalpy [kJ/kg]
    ``S``                       Entropy [kJ/kg/K]
    ``A``                       Speed of sound [m/s]
    ``G``                       Gibbs function [kJ/kg]
    ``V``                       Viscosity [Pa-s]
    ``L``                       Thermal conductivity [kW/m/K]
    ``I`` or `SurfaceTension`   Surface Tension [N/m]
    ``w`` or `accentric`        Accentric Factor [-]
    ==========================  ======================================================
    
    The following sets of input values are valid (order doesn't matter):
    
    =========================  ======================================
    ``InputName1``             ``InputName2``
    =========================  ======================================
    ``T``                      ``P``
    ``T``                      ``D``
    ``T``                      ``Q``
    ``P``                      ``Q``
    ``H``                      ``P``
    ``S``                      ``P``
    =========================  ======================================
    
    
    If `InputName1` is `T` and `OutputName` is ``I`` or ``SurfaceTension``, the second input is neglected
    since surface tension is only a function of temperature
    
    Call Type #3:
    New in 2.2
    If you provide InputName1 or InputName2 as a derived class of Quantity, the value will be internally
    converted to the required units as long as it is dimensionally correct.  Otherwise a ValueError will 
    be raised by the conversion
    """
    cdef char _in2
    cdef char _in4
    cdef bytes errs

    if (in3 is None
        and in4 is None
        and in5 is None
        and in6 is None
        and in7 is None):
        
        val = _Props1(in1,in2)
        if math.isnan(val) or abs(val)>1e20:
            raise ValueError(_get_errstring())
        else:
            return val
    else:
        if _quantities_supported:
            if isinstance(in3,pq.Quantity):
                in3 = _convert_to_default_units(in2,in3).magnitude
            if isinstance(in5,pq.Quantity):
                in5 = _convert_to_default_units(in4,in5).magnitude
                
        _in2 = <char>((<bytes>in2)[0])
        _in4 = <char>((<bytes>in4)[0])
        val = _Props(in1, _in2, in3, _in4, in5, in6)
        if math.isnan(val) or abs(val)>1e20:
            err_string = _get_errstring()
            if not len(err_string) == 0:
                raise ValueError(err_string)
            else:
                raise ValueError("Props failed ungracefully with inputs:"+str(in1)+','+str(in2)+','+str(in3)+','+str(in4)+','+str(in5)+','+str(in6)+'; please file a ticket at https://sourceforge.net/p/coolprop/tickets/')
        else:
            if not _quantities_supported:
                return val
            else:
                if in7 is not None:
                    #Convert the units to the units given by in7
                    return _convert_to_desired_units(val,in1,in7)
                else:
                    return val
    
cpdef double DerivTerms(bytes Output, double T, double rho, bytes Fluid):
    """

    .. |cubed| replace:: \ :sup:`3`\ 
    .. |squared| replace:: \ :sup:`2`\ 
    .. |IC| replace:: ``IsothermalCompressibility``
    
    Call signature::
    
        DerivTerms(OutputName, T, rho, Fluid) --> float
    
    where ``Fluid`` is a string with a valid CoolProp fluid name, and ``T`` and ``rho`` are the temperature in K and density in kg/m |cubed| .  The value 
    ``OutputName`` is one of the strings in the table below:
    
    ========================  =====================================================================================================================================
    OutputName                Description
    ========================  =====================================================================================================================================
    ``dpdT``                  Derivative of pressure with respect to temperature at constant density [kPa/K]
    ``dpdrho``                Derivative of pressure with respect to density at constant temperature [kPa/(kg/m\ |cubed|\ )]
    ``Z``                     Compressibility factor [-]
    ``dZ_dDelta``             Derivative of Z with respect to reduced density [-]
    ``dZ_dTau``               Derivative of Z with respect to inverse reduced temperature [-]
    ``B``                     Second virial coefficient [m\ |cubed|\ /kg]
    ``dBdT``                  Derivative of second virial coefficient with respect to temperature [m\ |cubed|\ /kg/K]
    ``C``                     Third virial coefficient [m\ :sup:`6`\ /kg\ |squared|\ ]
    ``dCdT``                  Derivative of third virial coefficient with respect to temperature [m\ :sup:`6`\ /kg\ |squared|\ /K]
    ``phir``                  Residual non-dimensionalized Helmholtz energy [-]
    ``dphir_dTau``            Partial of residual non-dimensionalized Helmholtz energy with respect to inverse reduced temperature [-]
    ``d2phir_dTau2``          Second partial of residual non-dimensionalized Helmholtz energy with respect to inverse reduced temperature [-]
    ``dphir_dDelta``          Partial of residual non-dimensionalized Helmholtz energy with respect to reduced density [-]
    ``d2phir_dDelta2``        Second partial of residual non-dimensionalized Helmholtz energy with respect to reduced density [-]
    ``d2phir_dDelta_dTau``    First cross-partial of residual non-dimensionalized Helmholtz energy [-]
    ``d3phir_dDelta2_dTau``   Second cross-partial of residual non-dimensionalized Helmholtz energy [-]
    ``phi0``                  Ideal-gas non-dimensionalized Helmholtz energy [-]
    ``dphi0_dTau``            Partial of ideal-gas non-dimensionalized Helmholtz energy with respect to inverse reduced temperature [-]
    ``d2phi0_dTau2``          Second partial of ideal-gas non-dimensionalized Helmholtz energy with respect to inverse reduced temperature [-]
    ``dphi0_dDelta``          Partial of ideal-gas non-dimensionalized Helmholtz energy with respect to reduced density [-]
    ``d2phi0_dDelta2``        Second partial of ideal-gas non-dimensionalized Helmholtz energy with respect to reduced density [-]
    |IC|                      Isothermal compressibility [1/kPa]
    ========================  =====================================================================================================================================
    """
    return _DerivTerms(Output,T,rho,Fluid)

cpdef string Phase(bytes Fluid, double T, double p):
    """
    Given a set of temperature and pressure, returns one of the following strings

    * Gas
    * Liquid
    * Supercritical
    * Two-Phase
    
    Phase diagram::
    
            |         |     
            |         |    Supercritical
            |         |
        p   | Liquid (b)------------
            |        /
            |       / 
            |      /       Gas
            |     / 
            |   (a)
            |  /
            |------------------------
    
                       T
    
           a: triple point
           b: critical point
           a-b: Saturation line
    """
    return _Phase(Fluid,T,p)

cpdef UseSaturationLUT(bint OnOff):
    """
    Turn the saturation lookup table on or off
    
    Parameters
    ----------
    OnOff : boolean
        If ``True``, turn on saturation lookup table
    """
    _UseSaturationLUT(OnOff)
    
cpdef bint SaturationLUTStatus():
    """
    Get the saturation lookup table status
    
    Returns
    -------
    Status : boolean
        ``True`` if saturation LUT is enabled, ``False`` otherwise
    
    """
    return _SaturationLUTStatus()

cpdef UseSinglePhaseLUT(bint OnOff):
    """
    Turn the SinglePhase lookup table on or off
    
    Parameters
    ----------
    OnOff : boolean
        If ``True``, turn on SinglePhase lookup table
    """
    _UseSinglePhaseLUT(OnOff)
    
cpdef bint SinglePhaseLUTStatus():
    """
    Get the SinglePhase lookup table status
    
    Returns
    -------
    Status : boolean
        ``True`` if SinglePhase LUT is enabled, ``False`` otherwise
    
    """
    return _SinglePhaseLUTStatus()

cpdef F2K(double T_F):
    """
    Convert temperature in degrees Fahrenheit to Kelvin
    """
    return _F2K(T_F)

cpdef K2F(double T_K):
    """
    Convert temperature in Kelvin to degrees Fahrenheit
    """
    return _K2F(T_K)

cpdef list FluidsList():
    """
    Return a list of strings of all fluid names
    
    Returns
    -------
    FluidsList : list of strings of fluid names
        All the fluids that are included in CoolProp
    
    Notes
    -----
    
    Here is an example::
        
       In [0]: from CoolProp.CoolProp import FluidsList
    
       In [1]: FluidsList()
       
    """ 
    return _FluidsList().split(',')

cpdef get_aliases(bytes Fluid):
    """
    Return a comma separated string of aliases for the given fluid
    """
    return _get_aliases(Fluid)
    
cpdef string get_REFPROPname(bytes Fluid):
    """
    Return the REFPROP compatible name for the fluid (only useful on windows)
    
    Some fluids do not use the REFPROP name.  For instance, 
    ammonia is R717, and propane is R290.  You can still can still call CoolProp
    using the name ammonia or R717, but REFPROP requires that you use a limited
    subset of names.  Therefore, this function that returns the REFPROP compatible
    name.  To then use this to call REFPROP, you would do something like::
    
       In [0]: from CoolProp.CoolProp import get_REFPROPname, Props
    
       In [1]: Fluid = 'REFPROP-' + get_REFPROPname('R290')
       
       In [2]: Props('D', 'T', 300, 'P', 300, Fluid)
    """
    return _get_REFPROPname(Fluid)

cpdef string get_errstr():
    """
    Return the current error string
    """
    return _get_errstring()

cpdef set_debug(int level):
    """
    Set the current debug level as integer in the range [0,10]
    
    Parameters
    ----------
    level : int
        If level is 0, no output will be written to screen, if >0, 
        some output will be written to screen.  The larger level is, 
        the more verbose the output will be
    """
    _debug(level)

cpdef get_debug():
    """
    Return the current debug level as integer
    """
    return _get_debug()

cpdef PrintSaturationTable(bytes FileName, bytes Fluid, double Tmin, double Tmax):
    """
    Write a saturation table to a file for the given fluid
    
    Parameters
    ----------
    FileName : string
    Fluid : string
    Tmin : float
        Minimum temp [K]
    Tmax : float
        Maximum temp [K]
    
    """
    _PrintSaturationTable(FileName, Fluid, Tmin, Tmax)
    
cpdef string get_EOSReference(bytes Fluid):
    """
    Return a string with the reference for the equation of state
    """
    return _get_EOSReference(Fluid)

cpdef string get_TransportReference(bytes Fluid):
    """
    Return a string with the reference for the transport properties (thermal conductivity, viscosity, surface tension)
    """
    return _get_TransportReference(Fluid)

cpdef int set_1phase_LUT_params(bytes Fluid, int nT, int np, double Tmin, double Tmax, double pmin, double pmax, bint rebuild=True):
    """
    Set the parameters for the lookup table in the single-phase region and optionally build the LUT
    
    Parameters
    ----------
    Fluid : string
        Fluid name
    nT : int
        Number of points for T linearly spaced
    np : int
        Number of points for p linearly spaced
    Tmin : float
        Minimum temperature [K]
    tmax : float
        Maximum temperature [K]
    pmin : float
        Minimum pressure [kPa]
    pmax : float
        Maximum pressure [kPa]
    rebuild : boolean
        If ``True``, build the LUT right when the function is called
    
    """
    return _set_1phase_LUT_params(Fluid,nT,np,Tmin,Tmax,pmin,pmax,rebuild)
    
cpdef dict get_1phase_LUT_params():
    cdef int *nT = NULL, *np = NULL
    cdef double *Tmin = NULL, *Tmax = NULL, *pmin = NULL, *pmax = NULL
    _get_1phase_LUT_params(nT,np,Tmin,Tmax,pmin,pmax)
    #In cython, nT[0] to dereference rather than *nT
    return dict(nT = nT[0],
                np = np[0],
                Tmin = Tmin[0],
                Tmax = Tmax[0],
                pmin = pmin[0],
                pmax = pmax[0]
                )
    
cpdef bint IsFluidType(bytes Ref, bytes Type):
    """
    Check if a fluid is of a given type
    
    Valid types are:
    - Brine
    - PseudoPure (or equivalently PseudoPureFluid)
    - PureFluid
    """
    if _IsFluidType(Ref,Type):
        return True
    else:
        return False
    
cpdef rhosatL_anc(bytes Fluid, double T):
    return _rhosatL_anc(Fluid,T)

cpdef rhosatV_anc(bytes Fluid, double T):
    return _rhosatV_anc(Fluid,T)

cpdef psatL_anc(bytes Fluid, double T):
    return _psatL_anc(Fluid,T)

cpdef psatV_anc(bytes Fluid, double T):
    return _psatV_anc(Fluid,T)
    

cpdef int debug(int level):
    """
    Sets the debug level
    
    Parameters
    ----------
    level : int
        Flag indicating how verbose the debugging should be.
            0 : no debugging output
            ...
            ...
            10 : very annoying debugging output - every function call debugged
    """
    _debug(level)
   
cpdef int LUT(bint LUTval):
    """
    
    LUTval : boolean
        If ``True``, turn on the use of lookup table.  Parameters must have 
        already been set through the use of set_1phase_LUT_params

    """
    if LUTval:
        _LUT_Enabled = True
        print 'Turning on singlephase LUT'
        UseSinglePhaseLUT(True)
    else:
        _LUT_Enabled = False
        UseSinglePhaseLUT(False)
        
from math import pow as pow_
cdef bint _LUT_Enabled



#A dictionary mapping parameter index to string for use with non-CoolProp fluids
cdef dict paras = {iMM : 'M',
                   iT : 'T',
                   iH : 'H',
                   iP : 'P',
                   iC : 'C',
                   iC0 : 'C0',
                   iO : 'O',
                   iV : 'V',
                   iL : 'L',
                   iS : 'S',
                   iU : 'U',
                   iDpdT : 'dpdT'}

@cython.final
cdef class State:
    """
    A class that contains all the code that represents a thermodynamic state
    """
    
    def __init__(self, bytes Fluid, dict StateDict, double xL=-1.0, bytes Liquid=str(''), phase = str('Gas')):
        self.Fluid = Fluid
        self.iFluid = _get_Fluid_index(Fluid)
        #Try to get the fluid from CoolProp
        if self.iFluid >= 0:
            #It is a CoolProp Fluid so we can use the faster integer passing function
            self.is_CPFluid = True
        else:
            self.is_CPFluid = False
        self.xL = xL
        self.Liquid = Liquid
        self.phase = phase
        #Parse the inputs provided
        self.update(StateDict)
        #Set the phase flag
        if self.phase == str('Gas') or self.phase == str('Liquid') or self.phase == str('Supercritical'):
            _set_phase(self.phase)
            
    def __reduce__(self):
        d={}
        d['xL']=self.xL
        d['Liquid']=self.Liquid
        d['Fluid']=self.Fluid
        d['T']=self.T_
        d['rho']=self.rho_
        return rebuildState,(d,)
          
          
    cpdef update_Trho(self, double T, double rho):
        """
        Just use the temperature and density directly
        """
        self.T_ = T
        self.rho_ = rho
        cdef double p
        
        if self.is_CPFluid:
            p = _IProps(iP,iT,T,iD,rho,self.iFluid)
        else:
            p = _Props('P','T',T,'D',rho,self.Fluid)
        
        if abs(p)<1e90:
            self.p_=p
        else:
            errstr = get_errstringc()
            raise ValueError(errstr)
        
    cpdef update(self,dict params, double xL=-1.0):
        """
        *params* is a list(or tuple) of strings that represent the parameters 
        that have been updated and will be used to fix the rest of the state. 
        ['T','P'] for temperature and pressure for instance
        """
            
        cdef double p
        cdef bytes errstr
        
        # If no value for xL is provided, it will have a value of -1 which is 
        # impossible, so don't update xL
        if xL > 0:
            #There is liquid
            self.xL=xL
            self.hasLiquid=True
        else:
            #There's no liquid
            self.xL=0.0
            self.hasLiquid=False
        
        #You passed in a dictionary, use the values to update the state
        if 'T' not in params:
            raise AttributeError('T must be provided in params dict in State.update')
            
        #Consume the 'T' key since it is required (TODO?)
        self.T_=float(params.pop('T'))
            
        #Given temperature and pressure, determine density of gas 
        # (or gas and oil if xL is provided)
        if abs(self.xL)<=1e-15:
            #Get the density if T,P provided, or pressure if T,rho provided
            if 'P' in params:
                self.p_=params['P']
                
                if self.is_CPFluid:
                    rho = _IProps(iD,iT,self.T_,iP,self.p_,self.iFluid)
                else:
                    rho = _Props('D','T',self.T_,'P',self.p_,self.Fluid)
                
                if abs(rho)<1e90:
                    self.rho_=rho
                else:
                    errstr = get_errstringc()
                    raise ValueError(errstr)
            elif 'D' in params:
                
                self.rho_=params['D']
                
                if self.is_CPFluid:
                    p = _IProps(iP,iT,self.T_,iD,self.rho_,self.iFluid)
                else:
                    p = _Props('P','T',self.T_,'D',self.rho_,self.Fluid)
                
                if abs(p)<1e90:
                    self.p_=p
                else:
                    errstr = get_errstringc()
                    raise ValueError(errstr+str(params))
            else:
                raise KeyError("Dictionary must contain the key 'T' and one of 'P' or 'D'")
            
        elif self.xL>0 and self.xL<=1:
            raise ValueError('xL is out of range - value for xL is [0,1]')
        else:
            raise ValueError('xL must be between 0 and 1')
        
    cpdef double Props(self, long iOutput):
        if iOutput<0:
            raise ValueError('Your output is invalid') 
        
        if self.is_CPFluid:
            return _IProps(iOutput,iT,self.T_,iD,self.rho_,self.iFluid)
        else:
            return _Props(paras[iOutput],'T',self.T_,'D',self.rho_,self.Fluid)
            
    cpdef double get_MM(self):
        return _Props1(self.Fluid,'molemass')
    
    cpdef double get_rho(self): 
        return self.rho_
    property rho:
        def __get__(self):
            return self.rho_
            
    cpdef double get_p(self): 
        return self.p_
    property p:
        def __get__(self):
            return self.p_
    
    cpdef double get_T(self): 
        return self.T_
    property T:
        def __get__(self):
            return self.T_
    
    cpdef double get_h(self): 
        return self.Props(iH)
    property h:
        def __get__(self):
            return self.get_h()
          
    cpdef double get_u(self): 
        return self.Props(iU)
    property u:
        def __get__(self):
            return self.get_u()
            
    cpdef double get_s(self): 
        return self.Props(iS)
    property s:
        def __get__(self):
            return self.get_s()
    
    cpdef double get_cp0(self):
        return self.Props(iC0)
    
    cpdef double get_cp(self): 
        return self.Props(iC)
    property cp:
        def __get__(self):
            return self.get_cp()
            
    cpdef double get_cv(self): 
        return self.Props(iO)
    property cv:
        def __get__(self):
            return self.get_cv()
            
    cpdef double get_visc(self):
        return self.Props(iV)
    property visc:
        def __get__(self):
            return self.get_visc()

    cpdef double get_cond(self):
        return self.Props(iL)
    property k:
        def __get__(self):
            return self.get_cond()
            
    property Prandtl:
        def __get__(self):
            return self.cp * self.visc / self.k
            
    cpdef double get_dpdT(self):
        return self.Props(iDpdT)
    property dpdT:
        def __get__(self):
            return self.get_dpdT()
        
    cpdef speed_test(self, int N):
        from time import clock
        cdef int i
        cdef char * k
        cdef string Fluid = self.Fluid
        cdef long IT = 'T'
        cdef long ID = 'D'
        import CoolProp as CP
        
        print 'Direct c++ call to CoolProp without the Python call layer'
        print "'M' involves basically no computational effort and is a good measure of the function call overhead"
        keys = ['H','P','S','U','C','O','V','L','M','C0']
        for key in keys:
            t1=clock()
            for i in range(N):
                _Props(key,'T',self.T_,'D',self.rho_,Fluid)
            t2=clock()
            print 'Elapsed time for {0:d} calls for "{1:s}" at {2:g} us/call'.format(N,key,(t2-t1)/N*1e6)
        
        print 'Call to the c++ layer using integers'
        keys = [iH,iP,iS,iU,iC,iO,iV,iL,iMM,iC0]
        for key in keys:
            t1=clock()
            for i in range(N):
                _IProps(key,iT,self.T_,iD,self.rho_,self.iFluid)
            t2=clock()
            print 'Elapsed time for {0:d} calls for "{1:s}" at {2:g} us/call'.format(N,paras[key],(t2-t1)/N*1e6)
            
        print 'Call to the Python call layer'
        print "'M' involves basically no computational effort and is a good measure of the function call overhead"
        keys = ['H','P','S','U','C','O','V','L','M','C0']
        for key in keys:
            t1=clock()
            for i in range(N):
                CP.Props(key,'T',self.T_,'D',self.rho_,Fluid)
            t2=clock()
            print 'Elapsed time for {0:d} calls for "{1:s}" at {2:g} us/call'.format(N,key,(t2-t1)/N*1e6)
    
    def __str__(self):
        """
        Return a string representation of the state
        """
        units={'T': 'K', 
               'p': 'kPa', 
               'rho': 'kg/m^3',
               'h':'kJ/kg',
               'u':'kJ/kg',
               's':'kJ/kg/K',
               'visc':'Pa-s',
               'k':'kW/m/K',
               'cp':'kJ/kg/K',
               'cv':'kJ/kg/K',
               'dpdT':'kPa/K'}
        s=''
        for k in ['T','p','rho','h','u','s','visc','k','cp','cv','dpdT','Prandtl']:
            if k in units:
                s+=k+' = '+str(getattr(self,k))+' '+units[k]+'\n'
            else:
                s+=k+' = '+str(getattr(self,k))+' NO UNITS'+'\n'
        return s.rstrip()
        
    cpdef copy(self):
        cdef double T = self.T_*(1.0+1e-20)
        cdef double rho = self.rho_*(1.0+1e-20)
        ST=State(self.Fluid,{'T':T,'D':rho})
        return ST
    
def rebuildState(d):
    S=State(d['Fluid'],{'T':d['T'],'D':d['rho']})
    S.xL = d['xL']
    S.Liquid=d['Liquid']
    return S