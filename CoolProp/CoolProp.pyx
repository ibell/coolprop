#cython: embedsignature = True

#
#
# This file provides wrapper functions of all the CoolProp functions
#
# The CoolProp header wrapper is provided in the file CoolProp.pyx
#
# Each of the functions from the CoolProp header are renamed in cython code to
# an underscored name so that the same name can be used in the exposed functions below

import cython
import math

cpdef double Props(bytes in1, bytes in2, in3=None, in4=None,in5=None,in6=None):
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
    =============  ============================
   
    Call Type #2:
    
    Alternatively, Props can be called in the form::
    
        Props(OutputName,InputName1,InputProp1,InputName2,InputProp2,Fluid) --> float
    
    where ``Fluid`` is a string with a valid CoolProp fluid name.  The value 
    ``OutputName`` is either a single-character or a string alias.  This list 
    shows the possible values
    
    ========================  ======================================================
    ``OutputName``            Description
    ========================  ======================================================
    `T`                       Temperature [K]
    `P`                       Pressure [kPa]
    `D`                       Density [kg/m3]
    `C0`                      Ideal-gas specific heat at constant pressure [kJ/kg]
    `C`                       Specific heat at constant pressure [kJ/kg]
    `O`                       Specific heat at constant volume [kJ/kg]
    `U`                       Internal energy [kJ/kg]
    `H`                       Enthalpy [kJ/kg]
    `S`                       Entropy [kJ/kg/K]
    `A`                       Speed of sound [m/s]
    `G`                       Gibbs function [kJ/kg]
    `V`                       Viscosity [Pa-s]
    `L`                       Thermal conductivity [kW/m/K]
    `I` or `SurfaceTension`   Surface Tension [N/m]
    ========================  ======================================================
    
    If `InputName1` is `T` and `InputName2` is `D`, any of the outputs are valid
    
    If `InputName1` is `T` and `InputName2` is `P`, any of the outputs are valid
    
    If `InputName1` is `T` and `InputName2` is `Q`, any of the outputs are valid
    
    If `InputName1` is `T` and `OutputName` is `I`, the second input is neglected
    since surface tension is only a function of temperature
    """
    cdef char _in2
    cdef char _in4
    
    if (in3 is None
        and in4 is None
        and in5 is None
        and in6 is None):
        
        val = _Props1(in1,in2)
        if math.isnan(val) or abs(val)>1e20:
            raise ValueError(_get_errstring())
        else:
            return val
    else:
        _in2 = <char>((<bytes>in2)[0])
        _in4 = <char>((<bytes>in4)[0])
        val = _Props(in1, _in2, in3, _in4, in5, in6)
        if math.isnan(val) or abs(val)>1e20:
            raise ValueError(_get_errstring())
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
    cdef int *nT, *np
    cdef double *Tmin, *Tmax, *pmin, *pmax
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