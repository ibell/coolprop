%module CoolProp
%include "typemaps.i"
%{
// This stuff will get included verbatim in CoolProp_wrap.c
#include "CoolProp.h"
%}

//All of the parameters that match the list below will be converted to output parameters
%apply double *OUTPUT { double *rhoLout, double *rhoVout, double *pout };

%include "CoolProp.h"


