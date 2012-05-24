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

#include "CoolPropTools.h"
#include "PropErrorCodes.h"
#include "PropMacros.h"

// Always include Air and Water
#include "Water.h"
#include "Air.h"

// If the proprocessor key ONLY_AIR_WATER is defined, only air and water will be included
#if !defined(ONLY_AIR_WATER)
#include "R134a.h"
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
#endif

#define PHASE_SUPERCRITICAL 1
#define PHASE_SUPERHEATED 4
#define PHASE_SUBCOOLED 2
#define PHASE_TWOPHASE 3    

#define FLUIDTYPE_REFPROP 0
#define FLUIDTYPE_BRINE 1
#define FLUIDTYPE_REFRIGERANT_PURE 2
#define FLUIDTYPE_REFRIGERANT_PSEUDOPURE 3

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

#ifndef CoolProp_H
#define CoolProp_H

	EXPORT_CODE void CONVENTION Help(void);
	EXPORT_CODE void CONVENTION UseSaturationLUT(bool OnOff);
	EXPORT_CODE void CONVENTION UseSinglePhaseLUT(bool OnOff);
    EXPORT_CODE bool CONVENTION SinglePhaseLUTStatus(void);
	EXPORT_CODE double CONVENTION Props(std::string Fluid,std::string Output);
	EXPORT_CODE double CONVENTION Props(char *Fluid, char *Output);
	EXPORT_CODE double CONVENTION Props(char Output,char Name1, double Prop1, char Name2, double Prop2, char * Ref);
    EXPORT_CODE double CONVENTION Props(char Output,char Name1, double Prop1, char Name2, double Prop2, std::string Ref);
	EXPORT_CODE double CONVENTION Props(char *Output,char Name1, double Prop1, char Name2, double Prop2, char * Ref);

	// Convenience functions
	EXPORT_CODE int CONVENTION IsFluidType(char *Ref, char *Type);
	EXPORT_CODE double CONVENTION T_hp(char *Ref, double h, double p, double T_guess);
	EXPORT_CODE double CONVENTION h_sp(char *Ref, double s, double p, double T_guess);
	EXPORT_CODE double CONVENTION DerivTerms(char *Term, double T, double rho, char * Ref);

	EXPORT_CODE double CONVENTION F2K(double T_F);
	EXPORT_CODE double CONVENTION K2F(double T);
	EXPORT_CODE void CONVENTION PrintSaturationTable(char *FileName, char * Ref, double Tmin, double Tmax);
	EXPORT_CODE int CONVENTION Phase(double T, double rho, char * Ref);
	
	EXPORT_CODE std::string CONVENTION FluidsList(void);
#endif
