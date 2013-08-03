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

	#if defined(COOLPROP_LIB)
	#define EXPORT_CODE extern "C" __declspec(dllexport)
	#ifndef CONVENTION
		#define CONVENTION __stdcall
	#endif
	#else
        #ifndef EXPORT_CODE
            #if defined(EXTERNC)
                #define EXPORT_CODE extern "C"
            #else
                #define EXPORT_CODE
            #endif
        #endif
		#ifndef CONVENTION
			#define CONVENTION
		#endif
	#endif

	// Hack for PowerPC compilation to only use extern "C"
	#if defined(__powerpc__) || defined(EXTERNC)
	#define EXPORT_CODE extern "C"
	#endif

	#include "FluidClass.h"
	
	// Functions with the call type like
    // EXPORT_CODE void CONVENTION AFunction(double, double);
	// will be exported to the DLL

	// They can only use data types that play well with DLL wrapping

	// This version uses the indices in place of the strings for speed.  Get the parameter indices
	// from get_param_index('D') for instance and the Fluid index from get_Fluid_index('Air') for instance
	EXPORT_CODE double CONVENTION IProps(long iOutput, long iName1, double Prop1, long iName2, double Prop2, long iFluid);

	EXPORT_CODE double CONVENTION Props(char *Output,char Name1, double Prop1, char Name2, double Prop2, char * Ref);
	EXPORT_CODE double CONVENTION Props1(char *Output, char * Ref);

	// Convenience functions
	EXPORT_CODE int CONVENTION IsFluidType(char *Ref, char *Type);
	EXPORT_CODE double CONVENTION DerivTerms(char *Term, double T, double rho, char * Ref);
	EXPORT_CODE void CONVENTION Phase(char *Fluid, double T, double p, char *Phase_str);
	EXPORT_CODE long CONVENTION Phase_Trho(char *Fluid, double T, double p, char *Phase_str);
	EXPORT_CODE long CONVENTION Phase_Tp(char *Fluid, double T, double rho, char *Phase_str);
	EXPORT_CODE void CONVENTION set_phase(char *Phase_str);
	EXPORT_CODE double CONVENTION F2K(double T_F);
	EXPORT_CODE double CONVENTION K2F(double T);
	EXPORT_CODE void CONVENTION PrintSaturationTable(char *FileName, char * Ref, double Tmin, double Tmax);
	
	EXPORT_CODE void CONVENTION FluidsList(char*);
	EXPORT_CODE void CONVENTION get_aliases(char* Ref, char *aliases);
	EXPORT_CODE void CONVENTION get_REFPROPname(char* Ref, char*);
	EXPORT_CODE void CONVENTION get_errstring(char*);
	EXPORT_CODE char* CONVENTION get_errstringc(void);
	EXPORT_CODE long CONVENTION get_errstring_copy(char *);
	EXPORT_CODE long CONVENTION get_svnrevision(void);
	EXPORT_CODE long CONVENTION get_version(char * pversion);
	EXPORT_CODE long CONVENTION get_param_index(char * param);
	EXPORT_CODE long CONVENTION get_Fluid_index(char * param);
	EXPORT_CODE long CONVENTION get_ASHRAE34(char * fluid, char *output);
	EXPORT_CODE long CONVENTION get_CAS_code(char * fluid, char *output);
	EXPORT_CODE void CONVENTION get_index_units(long param, char * units);

	

	EXPORT_CODE int CONVENTION get_debug();
	EXPORT_CODE void CONVENTION debug(int level);

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

	/// Get the TTSE mode (normal or bicubic)
	EXPORT_CODE int CONVENTION get_TTSE_mode(char *FluidName, char * Value);
	/// Set the TTSE mode (normal or bicubic)
	EXPORT_CODE int CONVENTION set_TTSE_mode(char *FluidName, char * Value);

	// Expose some functions that are useful for ECS debugging
	EXPORT_CODE double CONVENTION viscosity_dilute(char* FluidName, double T);
	EXPORT_CODE double CONVENTION viscosity_residual(char* FluidName, double T, double rho);
	EXPORT_CODE double CONVENTION conductivity_critical(char* FluidName, double T, double rho);
	EXPORT_CODE double CONVENTION conductivity_background(char* FluidName, double T, double rho);
	EXPORT_CODE double CONVENTION conformal_Trho(char* FluidName, char* ReferenceFluidName, double T, double rho, double *Tconform, double *rhoconform);

	// ------------------------------------------------------------------------------------------------
	// All the functions below this comment do NOT get exported to REFPROP DLL due to the fact that the 
	// DLL MUST use extern "C" for all exported functions, which does not allow for function overloads 
	// or the use of any c++ types like std::string or std::vector
	// ------------------------------------------------------------------------------------------------
	double Props(std::string Fluid,std::string Output);
	double Props(char *Fluid, char *Output);
	double Props1(std::string Fluid,std::string Output);
	double Props(char Output,char Name1, double Prop1, char Name2, double Prop2, char * Ref);
	double Props(std::string Output,char Name1, double Prop1, char Name2, double Prop2, std::string Ref);

	double DerivTerms(char *Term, double T, double rho, Fluid * pFluid);
	double DerivTerms(char *Term, double T, double rho, Fluid * pFluid, bool SinglePhase, bool TwoPhase);

	int debug();
	void set_debug(int level);
	void set_phase(std::string Phase_str);

	std::string Phase(std::string Fluid, double T, double p);
	std::string Phase_Trho(std::string Fluid, double T, double rho);
    std::string Phase_Tp(std::string Fluid, double T, double p);
	std::string get_BibTeXKey(std::string Ref, std::string item);
	std::string get_EOSReference(std::string Ref);
	std::string get_TransportReference(std::string Ref);
	std::string get_ASHRAE34(std::string Ref);
	std::string get_CAS_code(std::string Ref);
	std::string FluidsList(void);
	std::string get_aliases(std::string Ref);
	std::string get_REFPROPname(std::string Ref);
	std::string get_errstring(void);
	std::string get_version(void);
	std::string get_TTSE_mode(std::string FluidName);
	bool add_REFPROP_fluid(std::string FluidName);

	long get_param_index(std::string param);
	long get_Fluid_index(std::string param);
	std::string get_index_units(long index);

	Fluid * get_fluid(long iFluid);

	// Define some constants that will be used throughout
	enum params {iB,iT,iP,iD,iC,iC0,iO,iU,iH,iS,iA,iG,iQ,iV,iL,iI,iMM,iTcrit,iTtriple,iPtriple,iPcrit,iRhocrit,iAccentric,iDpdT,iDrhodT_p,iTmin,iDipole,iPhase,iPHASE_LIQUID,iPHASE_GAS,iPHASE_SUPERCRITICAL,iPHASE_TWOPHASE,iODP,iGWP20,iGWP100,iGWP500, iCritSplineT,iHcrit,iScrit};
	enum phases {iLiquid, iSupercritical, iGas, iTwoPhase};
#endif
