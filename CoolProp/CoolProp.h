/*
Add some pre-processor directives to this file so that it can either be built as 
usual, or if the COOLPROP_LIB macro is defined, it will export the functions in 
this file for building a static or dynamic library
*/
#if defined(COOLPROP_LIB)
	#define EXPORT_CODE __declspec(dllexport)
	#define CONVENTION __cdecl
#else
	#define EXPORT_CODE 
	#define CONVENTION 
#endif

	/*
	Following the naming conventions of MATLAB linked with REFPROP,
	each outputproperty is represented by one character:

	P   Pressure [kPa]
	T   Temperature [K]
	D   Density [kg/m3]
	H   Enthalpy [kJ/kg]
	S   Entropy [kJ/(kg/K)]
	U   Internal energy [kJ/kg]
	C   Cp [kJ/(kg K)]
	O   Cv [kJ/(kg K)]
	K   Ratio of specific heats (Cp/Cv) [-]
	A   Speed of sound [m/s]
	X   liquid phase and gas phase composition (mass fractions)
	V   Dynamic viscosity [Pa*s]
	L   Thermal conductivity [kW/(m K)]
	Q   Quality (vapor fraction) (kg/kg)
	I   Surface tension [N/m]
	F	Freezing point of secondary fluid [K] **NOT IN MATLAB-REFPROP **
	M	Maximum temperature for secondary fluid [K] **NOT IN MATLAB-REFPROP **
	B	Critical Temperature [K] **NOT IN MATLAB-REFPROP **
	E	Critical Pressure [K] **NOT IN MATLAB-REFPROP **
	R   

	******** To call **************
	To call the function Props, for instance for R410A at 300K, 400 kPa, you would do:
	Props("H","T",300,"P",400,"R410A")

	Or to call a pure fluid from REFPROP (for instance Propane).  
	The name of the refrigerant is "REPFROP-" plus the REFPROP defined name of the fluid, for instance
	"Propane" for propane (R290)

	See the folder C:\Program Files\REFPROP\fluids for the names of the fluids
	
	To call Propane:
	Props("H","T",300,"P",400,"REFPROP-Propane")

	**************** Inputs ***************
	The limited list of inputs that are allowed are:

	Prop1    ||    Prop2
	--------------------
	  T      ||      P
	  T      ||      Q

	*/

/*! \mainpage CoolProp Core Code Documentation

Welcome to the home page for the C++ sources of CoolProp.  This information may be useful for developers or just the plain inquisitive

You might want to start by looking at CoolProp.h
*/

#ifndef CoolProp_H
#define CoolProp_H

	#include "CoolPropTools.h"
	#include "PropErrorCodes.h"
	#include "PropMacros.h"

	// Always include Air and Water, and R134a since it is the default reference fluid for ECS
	#include "Water.h"
	#include "Air.h"
	#include "R134a.h"

	// If the preprocessor key ONLY_AIR_WATER is defined, only air and water will be included (and R134a since default for ECS)
	#if !defined(ONLY_AIR_WATER)
	#include "R410A.h"
	#include "R290.h"
	#include "R32.h"
	#include "R744.h"
	#include "R404A.h"
	#include "R507A.h"
	#include "R407C.h"
	#include "R717.h"
	#include "Argon.h"
	#include "R1234yf.h"
	#include "Nitrogen.h"
	#include "Brine.h"
	#include "IndustrialFluids.h"
	#include "R22.h"
	#endif

	
	// Functions within this extern "C" bracket will be export to the DLL
	extern "C" {
		EXPORT_CODE void CONVENTION Help(void);
		EXPORT_CODE void CONVENTION UseSaturationLUT(bool OnOff);
		EXPORT_CODE bool CONVENTION SaturationLUTStatus();
		EXPORT_CODE void CONVENTION UseSinglePhaseLUT(bool OnOff);
		EXPORT_CODE bool CONVENTION SinglePhaseLUTStatus(void);

		EXPORT_CODE double CONVENTION Props(char *Output,char Name1, double Prop1, char Name2, double Prop2, char * Ref);

		// Convenience functions
		EXPORT_CODE int CONVENTION IsFluidType(char *Ref, char *Type);
		EXPORT_CODE double CONVENTION T_hp(char *Ref, double h, double p, double T_guess);
		EXPORT_CODE double CONVENTION h_sp(char *Ref, double s, double p, double T_guess);
		EXPORT_CODE double CONVENTION DerivTerms(char *Term, double T, double rho, char * Ref);
		EXPORT_CODE void CONVENTION Phase(char *Fluid, double T, double p, char *Phase_str);

		EXPORT_CODE double CONVENTION F2K(double T_F);
		EXPORT_CODE double CONVENTION K2F(double T);
		EXPORT_CODE void CONVENTION PrintSaturationTable(char *FileName, char * Ref, double Tmin, double Tmax);
		
		EXPORT_CODE void CONVENTION FluidsList(char*);
		EXPORT_CODE void CONVENTION get_REFPROPname(char* Ref,char*);
		EXPORT_CODE void CONVENTION get_errstring(char*);
		EXPORT_CODE char* CONVENTION get_errstringc(void);

		/*
		returns 1 if parameters set properly
		*/
		EXPORT_CODE int CONVENTION set_1phase_LUT_params(char *Ref, int nT, int np, double Tmin, double Tmax, double pmin, double pmax, bool rebuild);
		EXPORT_CODE void CONVENTION get_1phase_LUT_params(int *nT, int *np, double *Tmin, double *Tmax, double *pmin, double *pmax);

		EXPORT_CODE int CONVENTION get_debug();
		EXPORT_CODE void CONVENTION debug(int level); 
	}

	// ------------------------------------------------------------------------------------------------
	// All the functions below this comment do NOT get exported to REFPROP DLL due to the fact that the 
	// DLL MUST use extern "C" for all exported functions, which does not allow for function overloads 
	// or the use of any c++ types like std::string or std::vector
	// ------------------------------------------------------------------------------------------------
	double Props(std::string Fluid,std::string Output);
	double Props(char *Fluid, char *Output);
	double Props(char Output,char Name1, double Prop1, char Name2, double Prop2, char * Ref);
	double Props(std::string Output,char Name1, double Prop1, char Name2, double Prop2, std::string Ref);

	int set_1phase_LUT_params(std::string Ref, int nT, int np, double Tmin, double Tmax, double pmin, double pmax, bool rebuild);
	// Without rebuild provided, default to not rebuild automatically (rebuild=False)
	int set_1phase_LUT_params(char *Ref,int nT, int np, double Tmin, double Tmax, double pmin, double pmax);
	int set_1phase_LUT_params(std::string Ref, int nT, int np, double Tmin, double Tmax, double pmin, double pmax);
	int debug();

	std::string Phase(std::string Fluid, double T, double p);
	std::string get_EOSReference(std::string Ref);
	std::string get_TransportReference(std::string Ref);
	std::string FluidsList(void);
	std::string get_REFPROPname(std::string Ref);
	std::string get_errstring(void);
#endif
