%module HumidAirProp
%include "typemaps.i"

%{
#include "HumAir.h"
%}

//All of the parameters that match the list below will be converted to output parameters
%apply double *OUTPUT { double *Tdp_out, double *W_out, double *h_out, double *RH_out, double *v_out};

%include "HumAir.h"


