#if defined(_MSC_VER)
#define _CRTDBG_MAP_ALLOC
#define _CRT_SECURE_NO_WARNINGS
#include <crtdbg.h>
#endif

#include "CoolProp.h"
#include "CPState.h"
#include "FluidClass.h"

static int debug_level=0;

EXPORT_CODE long CONVENTION redirect_stdout(const char* file)
{
	freopen(file, "a+", stdout);
	return 0;
}
EXPORT_CODE int CONVENTION set_reference_stateS(const char *Ref, const char *reference_state)
{
	return set_reference_stateS(std::string(Ref), std::string(reference_state));
}
EXPORT_CODE int CONVENTION set_reference_stateD(const char *Ref, double T, double rho, double h0, double s0)
{
	return set_reference_stateD(std::string(Ref), T, rho, h0, s0);
}
// All the function interfaces that point to the single-input Props function
EXPORT_CODE double CONVENTION Props1(const char *FluidName, const char *Output)
{
    return Props1(std::string(FluidName), std::string(Output));
}
EXPORT_CODE double CONVENTION Props1SI(const char *FluidName, const char *Output)
{
    return Props1SI(std::string(FluidName),std::string(Output));
}
EXPORT_CODE double CONVENTION PropsS(const char *Output,const char* Name1, double Prop1, const char* Name2, double Prop2, const char * Ref)
{
	double val = Props(Output,Name1[0],Prop1,Name2[0],Prop2,Ref);
	return val;
}
EXPORT_CODE double CONVENTION Props(const char *Output,char Name1, double Prop1, char Name2, double Prop2, const char * Ref)
{
    // Get parameter indices
    long iOutput = get_param_index(Output);
    long iName1 = get_param_index(std::string(1,Name1));
    long iName2 = get_param_index(std::string(1,Name2));
    char n1[] = "\0", n2[] = "\0";
    n1[0] = Name1;
    n2[0] = Name2;

    // Convert inputs to SI
    Prop1 = convert_from_unit_system_to_SI(iName1, Prop1, get_standard_unit_system());
    Prop2 = convert_from_unit_system_to_SI(iName2, Prop2, get_standard_unit_system());
    
	// Call the SI function
    double val = PropsSI(Output,(const char*)n1,Prop1,(const char*)n2,Prop2,Ref);

	// Convert back to unit system
    return convert_from_SI_to_unit_system(iOutput,val,get_standard_unit_system());
}

EXPORT_CODE double CONVENTION PropsSI(const char *Output, const char *Name1, double Prop1, const char *Name2, double Prop2, const char * FluidName)
{
    return PropsSI(std::string(Output),Name1,Prop1,Name2,Prop2,FluidName);
}
EXPORT_CODE double CONVENTION K2F(double T)
{ return T * 9 / 5 - 459.67; }

EXPORT_CODE double CONVENTION F2K(double T_F)
{ return (T_F + 459.67) * 5 / 9;}

EXPORT_CODE double CONVENTION DerivTerms(const char *Term, double T, double rho, const char * Ref)
{
	// Go to the std::string, std::string version
	double val = DerivTerms(std::string(Term),T,rho,std::string(Ref));
	return val;
}

EXPORT_CODE double CONVENTION fromSI(const char *input, double value, const char *new_system)
{
	return convert_from_SI_to_unit_system(input, value, new_system);
}
EXPORT_CODE double CONVENTION toSI(const char *input, double value, const char *old_system)
{
	return convert_from_unit_system_to_SI(input, value, old_system);
}

EXPORT_CODE int CONVENTION get_debug_level(){return debug_level;}
EXPORT_CODE void CONVENTION set_debug_level(int level){debug_level=level;}
EXPORT_CODE long CONVENTION get_Fluid_index(const char * param)
{
	return get_Fluid_index(std::string(param));
}
EXPORT_CODE int CONVENTION set_TTSE_mode(const char* fluid, const char *value)
{
	long iFluid = get_Fluid_index(fluid);
	if (iFluid > -1)
	{
		Fluid *pFluid = get_fluid(iFluid);
		
		// Try to build the LUT; Nothing will happen if the tables are already built
		CoolPropStateClass C(fluid);
		pFluid->build_TTSE_LUT();

		if (!strcmp(value,"TTSE"))
		{
			pFluid->TTSESinglePhase.set_mode(TTSE_MODE_TTSE); return true;
		}
		else if (!strcmp(value,"BICUBIC"))
		{
			pFluid->TTSESinglePhase.set_mode(TTSE_MODE_BICUBIC); return true;
		}
		else
		{
			return false;
		}
	}
	else
	{
		return false;
	}
}
EXPORT_CODE long CONVENTION get_param_index(const char * param)
{
	return get_param_index(std::string(param));
}
#ifndef SWIG
EXPORT_CODE long CONVENTION get_global_param_string(const char *param, char * Output)
{
	strcpy(Output,get_global_param_string(std::string(param)).c_str());
	return 0;
}
EXPORT_CODE long CONVENTION get_fluid_param_string(const char *fluid, const char *param, char * Output)
{
	strcpy(Output, get_fluid_param_string(std::string(fluid), std::string(param)).c_str());
	return 0;
}
#endif
EXPORT_CODE long CONVENTION Phase(const char *Fluid, double T, double p, char *Phase_str)
{
	strcpy(Phase_str,(char*)Phase(std::string(Fluid),T,p).c_str());
	return 0;
}
EXPORT_CODE long CONVENTION Phase_Tp(const char *Fluid, double T, double p, char *Phase_str)
{
	strcpy(Phase_str,(char*)Phase(std::string(Fluid),T,p).c_str());
	return 0;
}
EXPORT_CODE long CONVENTION Phase_Trho(const char *Fluid, double T, double rho, char *Phase_str)
{
	strcpy(Phase_str,(char*)Phase_Trho(std::string(Fluid),T,rho).c_str());
	return 0;
}
// A function to enforce the state if known
EXPORT_CODE void CONVENTION set_phase(const char *Phase_str){
	set_phase(std::string(Phase_str));
}
/// Enable the TTSE for this fluid
EXPORT_CODE bool CONVENTION enable_TTSE_LUT(const char *FluidName){
	long iFluid = get_Fluid_index(FluidName); if (iFluid<0){ return false; };
	Fluid *pFluid = get_fluid(iFluid);
	pFluid->enable_TTSE_LUT();
	return true;
};
/// Check if TTSE is enabled
EXPORT_CODE bool CONVENTION isenabled_TTSE_LUT(const char *FluidName){
	long iFluid = get_Fluid_index(FluidName); if (iFluid<0){ return false; };
	Fluid *pFluid = get_fluid(iFluid);
	return pFluid->isenabled_TTSE_LUT();
}
/// Disable the TTSE for this fluid
EXPORT_CODE bool CONVENTION disable_TTSE_LUT(const char *FluidName){
	long iFluid = get_Fluid_index(FluidName); if (iFluid<0){ return false; };
	Fluid *pFluid = get_fluid(iFluid);
	pFluid->disable_TTSE_LUT();
	return true;
}
/// Enable the writing of TTSE tables to file for this fluid
EXPORT_CODE bool CONVENTION enable_TTSE_LUT_writing(const char *FluidName){
	long iFluid = get_Fluid_index(FluidName); if (iFluid<0){ return true;};
	Fluid *pFluid = get_fluid(iFluid);
	pFluid->enable_TTSE_LUT_writing();
	return true;
};
/// Check if the writing of TTSE tables to file is enabled
EXPORT_CODE bool CONVENTION isenabled_TTSE_LUT_writing(const char *FluidName){
	long iFluid = get_Fluid_index(FluidName); if (iFluid<0){ return false;};
	Fluid *pFluid = get_fluid(iFluid);
	return pFluid->isenabled_TTSE_LUT_writing();
}
/// Disable the writing of TTSE tables to file for this fluid
EXPORT_CODE bool CONVENTION disable_TTSE_LUT_writing(const char *FluidName){
	long iFluid = get_Fluid_index(FluidName); if (iFluid<0){ return false;};
	Fluid *pFluid = get_fluid(iFluid);
	pFluid->disable_TTSE_LUT_writing();
	return true;
}
/// Over-ride the default size of both of the saturation LUT
EXPORT_CODE bool CONVENTION set_TTSESat_LUT_size(const char *FluidName, int Nsat){
	long iFluid = get_Fluid_index(FluidName); if (iFluid<0){ return false; };
	Fluid *pFluid = get_fluid(iFluid);
	pFluid->set_TTSESat_LUT_size(Nsat);
	return true;
}
/// Over-ride the default size of the single-phase LUT
EXPORT_CODE bool CONVENTION set_TTSESinglePhase_LUT_size(const char *FluidName, int Np, int Nh){
	long iFluid = get_Fluid_index(FluidName); if (iFluid<0){ return false;};
	Fluid *pFluid = get_fluid(iFluid);
	pFluid->set_TTSESinglePhase_LUT_size(Np,Nh);
	return true;
}
/// Over-ride the default range of the single-phase LUT
EXPORT_CODE bool CONVENTION set_TTSESinglePhase_LUT_range(const char *FluidName, double hmin, double hmax, double pmin, double pmax){
	long iFluid = get_Fluid_index(FluidName); if (iFluid<0){ return false;};
	Fluid *pFluid = get_fluid(iFluid);
	pFluid->set_TTSESinglePhase_LUT_range(hmin,hmax,pmin,pmax);
	return true;
}
/// Get the current range of the single-phase LUT
EXPORT_CODE bool CONVENTION get_TTSESinglePhase_LUT_range(const char *FluidName, double *hmin, double *hmax, double *pmin, double *pmax){
	long iFluid = get_Fluid_index(FluidName); if (iFluid<0){ return false;};
	Fluid *pFluid = get_fluid(iFluid);
	pFluid->get_TTSESinglePhase_LUT_range(hmin,hmax,pmin,pmax);
	if (!ValidNumber(*hmin) && !ValidNumber(*hmax) && !ValidNumber(*pmin) && !ValidNumber(*pmax))
	{
		return false;
	}
	else
	{	
		return true;
	}
}

/// Returns the value for the integer flag corresponding to the current set of units
/// @returns val The integer value for the current set of units, one of enumerated values UNIT_SYSTEM_SI, UNIT_SYSTEM_KSI (see GlobalConstants.h)
EXPORT_CODE int CONVENTION get_standard_unit_system(void){
    return _get_standard_unit_system();
}
/// Sets the flag for the integer flag corresponding to the current set of units
/// @param val The integer value for the current set of units, one of enumerated values UNIT_SYSTEM_SI, UNIT_SYSTEM_KSI (see GlobalConstants.h)
EXPORT_CODE void CONVENTION set_standard_unit_system(int val){
    return _set_standard_unit_system(val);
}
