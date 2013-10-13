%module CoolProp

// This %include allows the use of std::string natively
%include "std_string.i"

#ifdef SWIGOCTAVE
    %ignore DerivTerms(char *Term, double T, double rho, Fluid * pFluid, bool SinglePhase, bool TwoPhase);
#endif

#ifdef SWIGCSHARP
    %ignore get_global_param_string(char *param, char * Output);
    %ignore get_fluid_param_string(char *fluid, char *param, char * Output);
#endif

// This stuff will get included verbatim in CoolProp_wrap.cpp
%{
#include "CoolProp.h"
#include "CoolPropDLL.h"
#include "GlobalConstants.h"
#include "HumidAirProp.h"
%}

%include "CoolPropDLL.h"
%include "GlobalConstants.h"
%include "HumidAirProp.h"

