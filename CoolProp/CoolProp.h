/*
Add some pre-processor directives to this file so that it can either be built as 
usual, or if the COOLPROP_LIB macro is defined, it will export the functions in 
this file for building a static or dynamic library.  

The __stdcall calling convention is used by default.  By providing the macro CONVENTION, the 
calling convention can be changed at build time.

Any functions that are exported to DLL must also add the "EXPORT_CODE function_name CONVENTION ..." code 
to the CoolProp.cpp implementation.  See the Props function for instance
*/

/*! \mainpage CoolProp Core Code Documentation

Welcome to the home page for the C++ sources of CoolProp.  This information may be useful for developers or just the plain inquisitive

You might want to start by looking at CoolProp.h
*/

#ifndef CoolProp_H
#define CoolProp_H

	#include "FluidClass.h"
	#include "CoolPropDLL.h"
	#include "Units.h"

	double Props(std::string Fluid,std::string Output);
	double Props(char *Fluid, char *Output);
	double Props1(std::string Fluid,std::string Output);
	double Props(char Output,char Name1, double Prop1, char Name2, double Prop2, char * Ref);
	double Props(std::string Output,char Name1, double Prop1, char Name2, double Prop2, std::string Ref);

	double DerivTerms(char *Term, double T, double rho, Fluid * pFluid);
	double DerivTerms(char *Term, double T, double rho, Fluid * pFluid, bool SinglePhase, bool TwoPhase);

	std::string get_global_param_string(std::string ParamName);
	std::string get_fluid_param_string(std::string FluidName, std::string ParamName);

	// Getter and setter for debug level
	// ---------------------------------
	/// Get the debug level
	/// @returns level The level of the verbosity for the debugging output (0-10) 0: no debgging output
	int get_debug_level();
	/// Set the debug level
	/// @param level The level of the verbosity for the debugging output (0-10) 0: no debgging output
	void set_debug_level(int level);

	void set_phase(std::string Phase_str);

	std::string Phase(std::string Fluid, double T, double p);
	std::string Phase_Trho(std::string Fluid, double T, double rho);
    std::string Phase_Tp(std::string Fluid, double T, double p);

	std::string get_BibTeXKey(std::string Ref, std::string item);

	bool add_REFPROP_fluid(std::string FluidName);

	long get_param_index(std::string param);
	long get_Fluid_index(std::string FluidName);
	std::string get_index_units(long index);
	Fluid * get_fluid(long iFluid);
	void set_err_string(std::string err_string);

	int get_standard_unit_system();
	void set_standard_unit_system(int);

	// Define some constants that will be used throughout
	#include "GlobalConstants.h"
#endif
