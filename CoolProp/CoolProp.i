%module CoolProp

// This %include allows the use of std::string natively
%include "std_string.i"

// For C#, need to remove the overloads so that SWIG knows what to do
#ifdef SWIGCSHARP
    %ignore Props(char*, char*);
    %ignore Props1(char*, char*);
    %ignore set_phase(char*);
    %ignore get_param_index(char*);
    %ignore get_Fluid_index(char*);
    %ignore set_1phase_LUT_params(char *,int,int,double,double,double,double,bool);
    %ignore set_1phase_LUT_params(char *,int,int,double,double,double,double);
    %ignore Props(char*, char, double, char, double, char*);
#endif

// This stuff will get included verbatim in CoolProp_wrap.cpp
%{
#include "CoolProp.h"
%}

%include "CoolProp.h"

