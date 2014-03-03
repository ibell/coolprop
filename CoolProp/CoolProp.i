%module CoolProp

// This %include allows the use of std::string natively
%include "std_string.i"

#ifdef SWIGOCTAVE
    %ignore DerivTerms(char *Term, double T, double rho, Fluid * pFluid, bool SinglePhase, bool TwoPhase);
#endif

// This stuff will get included verbatim in CoolProp_wrap
%{
#include "CoolPropTools.h"
#include "CoolProp.h"
#include "CoolPropDLL.h"
#include "GlobalConstants.h"
#include "HumidAirProp.h"
%}

%include "CoolPropTools.h"
%include "CoolPropDLL.h"
%include "GlobalConstants.h"
%include "HumidAirProp.h"

