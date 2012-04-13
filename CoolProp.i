%module CoolProp
//%include "typemaps.i"
%{
// This stuff will get included verbatim in CoolProp_wrap.c
#include "CoolProp.h"
#define SWIG_FILE_WITH_INIT
%}

// Have SWIG generate python docstrings
%feature("autodoc", "1");

%include "CoolProp.h"
