// Some ideas came from http://www.scipy.org/Cookbook/SWIG_NumPy_examples

%module CoolProp

%{
// This stuff will get included verbatim in CoolProp_wrap.c
#define SWIG_FILE_WITH_INIT
#include "CoolProp.h"
%}

%include "numpy.i"

%init %{
    import_array();
%}

%apply (double* ARGOUT_ARRAY1, int DIM1) {(double* OutVec, int n)}
%apply (double* IN_ARRAY1, int DIM1) {(double* Prop1, int len1), (double* Prop2, int len2)}

%include "CoolProp.h"
    
// Have SWIG generate python docstrings
%feature("autodoc", "1");

//~ %rename (Props) myPropsV;

//~ %inline %{
    //~ void myPropsV(char *Output,char Name1, double *Prop1, int len1, char Name2, double *Prop2, int len2, char * Ref, double *OutVec, int n) 
    //~ {
        //~ if (len1 != len2) 
        //~ {
            //~ PyErr_Format(PyExc_ValueError, "Arrays of lengths (%d,%d) given.  Must be the same length", len1, len2);
            //~ return 0.0;
        //~ }
        //~ PropsV(Output,Name1, Prop1, Name2, Prop2, Ref, OutVec, n);
    //~ }
//~ %}
