%module CoolProp
//%include "typemaps.i"
%{
// This stuff will get included verbatim in CoolProp_wrap.c
#include "CoolProp.h"
%}

//All of the parameters that match the list below will be converted to output parameters
//%apply double *OUTPUT { double *rhoLout, double *rhoVout, double *pout };

// Have SWIG generate python docstrings
%feature("autodoc", "1");

%include "CoolProp.h"
