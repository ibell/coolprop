%module CoolProp

// This %include allows the use of std::string natively
%include "std_string.i"

// This stuff will get included verbatim in CoolProp_wrap.cpp
%{
#include "CoolProp.h"
%}

// Have SWIG generate python docstrings for all the functions automatically
%feature("autodoc", "1");

#ifdef SWIGPYTHON
/* 
Something went wrong because an exception was thrown, yielding an infinite output value.
get_errstring() returns a std::string string, which is converted to a c-string and returned in the ValueError call
*/
%exception Props {
  $action
  if (result>1e20 || -result>1e20) {
     PyErr_SetString(PyExc_ValueError,get_errstring().c_str());
     return NULL;
  }
}
#endif

//Add the docstring for the Props function
%define PROPSDOCSTRING
"
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

"
%enddef
%feature("autodoc",PROPSDOCSTRING) Props(char *, char, double, char, double, char*);

%define PROPSDOCSTRING2
"
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

"
%enddef
%feature("autodoc",PROPSDOCSTRING2) Props(char *, char*);


%define DERIVTERMS
"
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
"
%enddef
%feature("autodoc",DERIVTERMS) DerivTerms(char *, double T, double rho, char*);

%define PHASE
"
Given a set of temperature and pressure, returns one of the following strings

* Gas
* Liquid
* Supercritical
* Two-Phase

Phase diagram::

		|         |     
		|         |    Supercritical
		|         |
	p	| Liquid (b)------------
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
"
%enddef
%feature("autodoc",PHASE) Phase(std::string,double,double);

%ignore Props(char, char, double, char, double, char*);
%ignore Props(std::string, char, double, char, double, std::string);
%ignore Props(std::string, std::string);

%include "CoolProp.h"