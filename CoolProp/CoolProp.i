// Some ideas came from http://www.scipy.org/Cookbook/SWIG_NumPy_examples

%module CoolProp

%include "std_string.i"

%{
// This stuff will get included verbatim in CoolProp_wrap.c
//#define SWIG_FILE_WITH_INIT
#include "CoolProp.h"
%}

//%include "numpy.i"

//%init %{
//    import_array();
//%}

%exception Props {
  $action
  if (result>1e20 || -result>1e20) {
     PyErr_SetString(PyExc_ValueError,"Props returned an invalid number");
     return NULL;
  }
}
%include "CoolProp.h"
    
// Have SWIG generate python docstrings
%feature("autodoc", "1");