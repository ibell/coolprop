#ifndef COOLPROPDLL_H
#define COOLPROPDLL_H

    #include "CoolPropTools.h"
	
	// Functions with the call type like
	// EXPORT_CODE void CONVENTION AFunction(double, double);
	// will be exported to the DLL

    /*
    ####################################################################################
    Overloads for DLL wrapping puposes

    These functions take strings which are 0-terminated.  These functions pass directly to 
    equivalently named functions in CoolProp.h that take std::string
    ####################################################################################
    */

    /**
    \overload 
    \sa Props1SI(std::string, std::string)
    */
    EXPORT_CODE double CONVENTION Props1SI(const char *FluidName, const char* Output);
    /**
    \overload 
    \sa PropsSI(std::string, std::string, double, std::string, double, std::string)
    */
    EXPORT_CODE double CONVENTION PropsSI(const char *Output,const char *Name1, double  Prop1, const char *Name2, double  Prop2, const char *Ref);
	/**
    \overload 
    \sa Props(std::string, std::string, double, std::string, double, std::string)
    */
	EXPORT_CODE double CONVENTION PropsS(const char *Output,const char *Name1, double  Prop1, const char *Name2, double  Prop2, const char *Ref);
    /**
    \overload 
    \sa Props(std::string, std::string, double, std::string, double, std::string)
    */
	EXPORT_CODE double CONVENTION Props(const char *Output,const char  Name1, double  Prop1, const char  Name2, double  Prop2, const char *Ref);
    /**
    \overload 
    \sa Props1(std::string, std::string)
    */
	EXPORT_CODE double CONVENTION Props1(const char *FluidName, const char *Output);
    /**
    \overload
    \sa IsFluidType(std::string, std::string)
    */
	EXPORT_CODE int CONVENTION IsFluidType(const char *Ref, const char *Type);

	EXPORT_CODE double CONVENTION DerivTerms(const char *Term, double T, double rho, const char *Ref);
	EXPORT_CODE long CONVENTION Phase(const char *Fluid, double T, double p, char *Phase_str);
	EXPORT_CODE long CONVENTION Phase_Trho(const char *Fluid, double T, double p, char *Phase_str);
	EXPORT_CODE long CONVENTION Phase_Tp(const char *Fluid, double T, double rho, char *Phase_str);
	EXPORT_CODE void CONVENTION set_phase(const char *Phase_str);

    EXPORT_CODE double CONVENTION fromSI(const char *input, const double value, const char *new_system);
	EXPORT_CODE double CONVENTION   toSI(const char *input, const double value, const char *old_system);
    
    // When using SWIG, it is extremely difficult to deal with char* for output strings, so just use 
	// the std::string version since SWIG can handle std::string just fine
	#if defined(SWIG)
		std::string get_global_param_string(std::string ParamName);
		std::string get_fluid_param_string(std::string FluidName, std::string ParamName);
	#else
		EXPORT_CODE long CONVENTION get_global_param_string(const char *param, char *Output);
		EXPORT_CODE long CONVENTION get_fluid_param_string(const char *fluid, const char *param, char *Output);
	#endif

    /**
    \overload
    \sa set_reference_stateS(std::string, std::string)
    */
    EXPORT_CODE int CONVENTION set_reference_stateS(const char *Ref, const char *reference_state);
	/**
    \overload
    \sa set_reference_stateD(std::string, double, double, double, double)
    */
    EXPORT_CODE int CONVENTION set_reference_stateD(const char *Ref, double T, double rho, double h0, double s0);
    /*
    ####################################################################################
    Implemented functions

    These functions take inputs that are compatible with DLL passing and are 
    implemented in CoolPropDLL.cpp
    ####################################################################################
    */

	// This version uses the indices in place of the strings for speed.  Get the parameter indices
	// from get_param_index('D') for instance and the Fluid index from get_Fluid_index('Air') for instance
	EXPORT_CODE double CONVENTION IPropsSI(long iOutput, long iName1, double Prop1, long iName2, double Prop2, long iFluid);
    EXPORT_CODE double CONVENTION IProps(long iOutput, long iName1, double Prop1, long iName2, double Prop2, long iFluid);
	
    /// Convert from degrees Fahrenheit to Kelvin (useful primarily for testing)
	EXPORT_CODE double CONVENTION F2K(double T_F);
	/// Convert from Kelvin to degrees Fahrenheit (useful primarily for testing)
    EXPORT_CODE double CONVENTION K2F(double T_K);

    /// Returns the value for the integer flag corresponding to the current set of units
	/// @returns val The integer value for the current set of units, one of enumerated values UNIT_SYSTEM_SI, UNIT_SYSTEM_KSI
    /// \sa GlobalConstants.h
	EXPORT_CODE int CONVENTION get_standard_unit_system();
	/// Sets the flag for the integer flag corresponding to the current set of units
	/// @param val The integer value for the current set of units, one of enumerated values UNIT_SYSTEM_SI, UNIT_SYSTEM_KSI
    /// \sa GlobalConstants.h
	EXPORT_CODE void CONVENTION set_standard_unit_system(int val);
    
	EXPORT_CODE long CONVENTION get_param_index(const char *param);
	EXPORT_CODE long CONVENTION get_Fluid_index(const char *param);
	EXPORT_CODE long CONVENTION redirect_stdout(const char *file);

    // ---------------------------------
	// Getter and setter for debug level
	// ---------------------------------

	/// Get the debug level
	/// @returns level The level of the verbosity for the debugging output (0-10) 0: no debgging output
	EXPORT_CODE int CONVENTION get_debug_level();
	/// Set the debug level
	/// @param level The level of the verbosity for the debugging output (0-10) 0: no debgging output
	EXPORT_CODE void CONVENTION set_debug_level(int level);

	/// -------------------------------------------
	///     TTSE Tabular Taylor Series Expansion
	/// -------------------------------------------
	/// Enable the TTSE
	EXPORT_CODE bool CONVENTION enable_TTSE_LUT(const char *FluidName);
	/// Check if TTSE is enabled
	EXPORT_CODE bool CONVENTION isenabled_TTSE_LUT(const char *FluidName);
	/// Disable the TTSE
	EXPORT_CODE bool CONVENTION disable_TTSE_LUT(const char *FluidName);
	/// Enable the writing of TTSE tables to file for this fluid
	EXPORT_CODE bool CONVENTION enable_TTSE_LUT_writing(const char *FluidName);
	/// Check if the writing of TTSE tables to file is enabled
	EXPORT_CODE bool CONVENTION isenabled_TTSE_LUT_writing(const char *FluidName);
	/// Disable the writing of TTSE tables to file for this fluid
	EXPORT_CODE bool CONVENTION disable_TTSE_LUT_writing(const char *FluidName);
	/// Over-ride the default size of both of the saturation LUT
	EXPORT_CODE bool CONVENTION set_TTSESat_LUT_size(const char *FluidName, int);
	/// Over-ride the default size of the single-phase LUT
	EXPORT_CODE bool CONVENTION set_TTSESinglePhase_LUT_size(const char *FluidName, int Np, int Nh);
	/// Over-ride the default range of the single-phase LUT
	EXPORT_CODE bool CONVENTION set_TTSESinglePhase_LUT_range(const char *FluidName, double hmin, double hmax, double pmin, double pmax);
	/// Get the current range of the single-phase LUT
	EXPORT_CODE bool CONVENTION get_TTSESinglePhase_LUT_range(const char *FluidName, double *hmin, double *hmax, double *pmin, double *pmax);
	/// Set the TTSE mode (normal or bicubic)
	EXPORT_CODE int CONVENTION set_TTSE_mode(const char *FluidName, const char * Value);


#endif
