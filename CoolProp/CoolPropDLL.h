#ifndef COOLPROPDLL_H
#define COOLPROPDLL_H

    #include "CoolPropTools.h"
	
	// Functions with the call type like
	// EXPORT_CODE void CONVENTION AFunction(double, double);
	// will be exported to the DLL

	// When using SWIG, it is extremely difficult to deal with char* for output strings, so just use 
	// the std::string version since SWIG can handle std::string just fine
	#ifdef SWIG
		std::string get_global_param_string(std::string ParamName);
		std::string get_fluid_param_string(std::string FluidName, std::string ParamName);
	#else
		EXPORT_CODE long CONVENTION get_global_param_string(char *param, char * Output);
		EXPORT_CODE long CONVENTION get_fluid_param_string(char *fluid, char *param, char * Output);
	#endif

	// They can only use data types that play well with DLL wrapping (int, long, double, char*, void, etc.)
	EXPORT_CODE double CONVENTION PropsS(char *Output,char *Name1, double  Prop1, char *Name2, double  Prop2, char *Ref);
	EXPORT_CODE double CONVENTION Props(char *Output,char  Name1, double  Prop1, char  Name2, double  Prop2, char *Ref);
	EXPORT_CODE double CONVENTION PropsSI(char *Output,char *Name1, double  Prop1, char *Name2, double  Prop2, char *Ref);

	/// Return a fluid value that does not depend on the thermodynamic state
	/// @param FluidName The name of the fluid
	/// @param Output The name of the output parameter, some options are "Ttriple", "Tcrit", "pcrit", "Tmin", "molemass", "rhocrit", "accentric" (not all parameters are valid for all fluids)
	/// @returns val The value, or _HUGE if not valid
	EXPORT_CODE double CONVENTION Props1SI(char *FluidName, char* Output);

	EXPORT_CODE double CONVENTION Props1(char *Ref, char *Output);

	// This version uses the indices in place of the strings for speed.  Get the parameter indices
	// from get_param_index('D') for instance and the Fluid index from get_Fluid_index('Air') for instance
	EXPORT_CODE double CONVENTION IProps(long iOutput, long iName1, double Prop1, long iName2, double Prop2, long iFluid);

	// Convenience functions
	EXPORT_CODE int CONVENTION IsFluidType(char *Ref, char *Type);
	EXPORT_CODE double CONVENTION DerivTerms(char *Term, double T, double rho, char * Ref);
	EXPORT_CODE long CONVENTION Phase(char *Fluid, double T, double p, char *Phase_str);
	EXPORT_CODE long CONVENTION Phase_Trho(char *Fluid, double T, double p, char *Phase_str);
	EXPORT_CODE long CONVENTION Phase_Tp(char *Fluid, double T, double rho, char *Phase_str);
	EXPORT_CODE void CONVENTION set_phase(char *Phase_str);
	EXPORT_CODE double CONVENTION F2K(double T_F);
	EXPORT_CODE double CONVENTION K2F(double T);
	
	EXPORT_CODE double CONVENTION fromSI(char *input, double value, char *new_system);
	EXPORT_CODE double CONVENTION   toSI(char *input, double value, char *old_system);

	EXPORT_CODE long CONVENTION get_param_index(char * param);
	EXPORT_CODE long CONVENTION get_Fluid_index(char * param);
	EXPORT_CODE long CONVENTION get_index_units(long param, char * units);

	EXPORT_CODE long CONVENTION redirect_stdout(char* file);

	// Getter and setter for debug level
	// ---------------------------------
	/// Get the debug level
	/// @returns level The level of the verbosity for the debugging output (0-10) 0: no debgging output
	EXPORT_CODE int CONVENTION get_debug_level();
	/// Set the debug level
	/// @param level The level of the verbosity for the debugging output (0-10) 0: no debgging output
	EXPORT_CODE void CONVENTION set_debug_level(int level);

	EXPORT_CODE double CONVENTION rhosatL_anc(char* Fluid, double T);
	EXPORT_CODE double CONVENTION rhosatV_anc(char* Fluid, double T);
	EXPORT_CODE double CONVENTION psatL_anc(char* Fluid, double T);
	EXPORT_CODE double CONVENTION psatV_anc(char* Fluid, double T);

	/// -------------------------------------------
	///     TTSE Tabular Taylor Series Expansion
	/// -------------------------------------------
	/// Enable the TTSE
	EXPORT_CODE bool CONVENTION enable_TTSE_LUT(char *FluidName);
	/// Check if TTSE is enabled
	EXPORT_CODE bool CONVENTION isenabled_TTSE_LUT(char *FluidName);
	/// Disable the TTSE
	EXPORT_CODE bool CONVENTION disable_TTSE_LUT(char *FluidName);
	/// Enable the writing of TTSE tables to file for this fluid
	EXPORT_CODE bool CONVENTION enable_TTSE_LUT_writing(char *FluidName);
	/// Check if the writing of TTSE tables to file is enabled
	EXPORT_CODE bool CONVENTION isenabled_TTSE_LUT_writing(char *FluidName);
	/// Disable the writing of TTSE tables to file for this fluid
	EXPORT_CODE bool CONVENTION disable_TTSE_LUT_writing(char *FluidName);
	/// Over-ride the default size of both of the saturation LUT
	EXPORT_CODE bool CONVENTION set_TTSESat_LUT_size(char *FluidName, int);
	/// Over-ride the default size of the single-phase LUT
	EXPORT_CODE bool CONVENTION set_TTSESinglePhase_LUT_size(char *FluidName, int Np, int Nh);
	/// Over-ride the default range of the single-phase LUT
	EXPORT_CODE bool CONVENTION set_TTSESinglePhase_LUT_range(char *FluidName, double hmin, double hmax, double pmin, double pmax);
	/// Get the current range of the single-phase LUT
	EXPORT_CODE bool CONVENTION get_TTSESinglePhase_LUT_range(char *FluidName, double *hmin, double *hmax, double *pmin, double *pmax);

	/// Set the TTSE mode (normal or bicubic)
	EXPORT_CODE int CONVENTION set_TTSE_mode(char *FluidName, char * Value);

	EXPORT_CODE int CONVENTION set_reference_stateS(char *Ref, char *reference_state);
	EXPORT_CODE int CONVENTION set_reference_stateD(char *Ref, double T, double rho, double h0, double s0);

	/// Returns the value for the integer flag corresponding to the current set of units
	/// @returns val The integer value for the current set of units, one of enumerated values UNIT_SYSTEM_SI, UNIT_SYSTEM_KSI (see GlobalConstants.h)
	EXPORT_CODE int CONVENTION get_standard_unit_system();
	/// Sets the flag for the integer flag corresponding to the current set of units
	/// @param val The integer value for the current set of units, one of enumerated values UNIT_SYSTEM_SI, UNIT_SYSTEM_KSI (see GlobalConstants.h)
	EXPORT_CODE void CONVENTION set_standard_unit_system(int val);

	// Expose some functions that are useful for ECS debugging
	EXPORT_CODE double CONVENTION viscosity_dilute(char* FluidName, double T);
	EXPORT_CODE double CONVENTION viscosity_residual(char* FluidName, double T, double rho);
	EXPORT_CODE double CONVENTION conductivity_critical(char* FluidName, double T, double rho);
	EXPORT_CODE double CONVENTION conductivity_background(char* FluidName, double T, double rho);
	EXPORT_CODE double CONVENTION conformal_Trho(char* FluidName, char* ReferenceFluidName, double T, double rho, double *Tconform, double *rhoconform);

#endif
