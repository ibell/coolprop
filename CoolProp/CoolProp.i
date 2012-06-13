%module CoolProp

// This %include allows the use of std::string natively in Python
%include "std_string.i"

// This stuff will get included verbatim in CoolProp_wrap.cpp
%{
#include "CoolProp.h"
%}

// Have SWIG generate python docstrings for all the functions automatically
%feature("autodoc", "1");

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

========================  =================================================
``OutputName``            Description
========================  =================================================
`T`                       Temperature [K]
`P`                       Pressure [kPa]
`D`                       Density [kg/m3]
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
========================  =================================================

If `InputName1` is `T` and `InputName2` is `D`, any of the outputs are valid

If `InputName1` is `T` and `InputName2` is `P`, any of the outputs are valid

If `InputName1` is `T` and `InputName2` is `Q`, any of the outputs are valid

If `InputName1` is `T` and `OutputName` is `I`, the second input is neglected
since surface tension is only a function of temperature

"
%enddef
%feature("autodoc",PROPSDOCSTRING2) Props(char *, char*);

%ignore Props(char, char, double, char, double, char*);
%ignore Props(std::string, char, double, char, double, std::string);
%ignore Props(std::string, std::string);

%include "CoolProp.h"