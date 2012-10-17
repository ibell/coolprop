%module CoolProp

// This %include allows the use of std::string natively
%include "std_string.i"

// This stuff will get included verbatim in CoolProp_wrap.cpp
%{
#include "CoolProp.h"
%}

%include "CoolProp.h"