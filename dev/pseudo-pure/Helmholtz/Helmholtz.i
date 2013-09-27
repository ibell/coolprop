%module helmholtz

//// *************************** EXCEPTION HANDLING ****************************
//// *************************** EXCEPTION HANDLING ****************************
//// *************************** EXCEPTION HANDLING ****************************
%include exception.i

// A generic exception handler.  Any exceptions thrown in C++ will be caught here
%exception {
	try {
		$action
	}
    catch(std::exception &e) {
		SWIG_exception(SWIG_RuntimeError,e.what());
	}
    catch(...) {
		SWIG_exception(SWIG_RuntimeError,"Unknown exception");
	}
}

// This allows for the use of STL vectors
%include "std_vector.i"
// This allows for the use of STL strings
%include "std_string.i"

namespace std {
   %template(vectord) vector<double>;
};

// This stuff will get included verbatim in CoolProp_wrap.cpp
%{
#include "../../../CoolProp/Helmholtz.h"
#include "../../../CoolProp/CoolPropTools.h"
%}

// This is where the parsing actually happens
%include "../../../CoolProp/Helmholtz.h"