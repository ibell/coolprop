%module CoolProp

%include "std_string.i"

// This stuff will get included verbatim in CoolProp_wrap.cpp
%{
#include "CoolProp.h"
%}

// Have SWIG generate python docstrings
%feature("autodoc", "0");

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
The function Props has two ways it can be called.

For more information, go to http://coolprop.sourceforge.net
"
%enddef
%feature("autodoc", PROPSDOCSTRING) Props;

%include "CoolProp.h"