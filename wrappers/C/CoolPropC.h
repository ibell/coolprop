#ifndef CoolPropC_H
#define CoolPropC_H

#define COOLPROP_LIB

#include <stdbool.h>

#ifdef __cplusplus
#include <CoolProp/CoolPropTools.h>
//#include <CoolProp/CoolProp.h>
extern "C" {
#endif

#ifndef CONVENTION
#  define CONVENTION
#endif

/// Return a fluid value that does not depend on the thermodynamic state
/// @param FluidName The name of the fluid
/// @param Output The name of the output parameter, some options are "Ttriple", "Tcrit", "pcrit", "Tmin", "molemass", "rhocrit", "accentric" (not all parameters are valid for all fluids)
/// @returns val The value, or _HUGE if not valid
double CONVENTION props1_(char * FluidName, char * Output);

/// Return a value that depends on the thermodynamic state
/// @param Output The output parameter, one of "T","D","H",etc.
/// @param Name1 The first state variable name, one of "T","D","H",etc.
/// @param Prop1 The first state variable value
/// @param Name2 The second state variable name, one of "T","D","H",etc.
/// @param Prop2 The second state variable value
/// @param FluidName The fluid name
double CONVENTION props_(char * Output, char *Name1, double *Prop1, char *Name2, double *Prop2, char * FluidName);
//double Props(char Output,char Name1, double Prop1, char Name2, double Prop2, char * FluidName);

/// Return some low level derivative terms, see source for a complete list
/// @param Term String, some options are "phir" (residual Helmholtz energy),"dphir_dDelta", "dphir_dTau", etc.
/// @param T Temperature [K]
/// @param rho Density [kg/m^3]
/// @param FluidName String
double CONVENTION derivterms_(char * Term, double *T, double *rho, char * FluidName);
//double DerivTerms(std::string Term, double T, double rho, std::string FluidName);

/// Set the reference state based on a string representation (consistent naming with RFPROP)
/// @param FluidName The name of the fluid
/// @param reference_state The reference state to use, one of "IIR" (h=200 kJ/kg, s=1 kJ/kg/K at 0C sat. liq.) "ASHRAE" (h=0,s=0 @ -40C sat liq), "NBP" (h=0,s=0 @ 1.0 bar sat liq.)
int CONVENTION set_reference_states_(char * FluidName, char * reference_state);
//int set_reference_stateS(std::string FluidName, std::string reference_state);

/// Enable the TTSE
bool CONVENTION enable_ttse_lut_(char *FluidName);
//bool CONVENTION enable_TTSE_LUT(char *FluidName);

/// Check if TTSE is enabled
bool CONVENTION isenabled_ttse_lut_(char *FluidName);
//bool CONVENTION isenabled_TTSE_LUT(char *FluidName);

/// Disable the TTSE
bool CONVENTION disable_ttse_lut_(char *FluidName);
//bool CONVENTION disable_TTSE_LUT(char *FluidName);

/// Set the TTSE mode (normal or bicubic)
int CONVENTION set_ttse_mode_(char *FluidName, char * Value);
//int CONVENTION set_TTSE_mode(char *FluidName, char * Value);


#ifdef __cplusplus
}
#endif

#endif 