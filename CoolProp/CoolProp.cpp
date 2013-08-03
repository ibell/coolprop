#if defined(_MSC_VER)
#define _CRTDBG_MAP_ALLOC
#define _CRT_SECURE_NO_WARNINGS
#include <crtdbg.h>
#endif

#include "CoolProp.h"
#include "REFPROP.h"


#if defined(__ISWINDOWS__)
#include <windows.h>
#else
#ifndef DBL_EPSILON
	#include <limits>
	#define DBL_EPSILON std::numeric_limits<double>::epsilon()
#endif
#endif

#include <iostream>
#include <stdlib.h>
#include <vector>
#include <exception>
#include <stdio.h>
#include "string.h"
#include "FluidClass.h"
#include "CoolPropTools.h"
#include "CPExceptions.h"
#include "Brine.h"
#include "Solvers.h"
#include "CPState.h"
#include "IncompLiquid.h"
#include "TTSE.h"
#include "purefluids/R290.h"
#include "purefluids/R134a.h"

// Function prototypes
double rho_TP(double T, double p);
double _Props(std::string Output,std::string Name1, double Prop1, std::string Name2, double Prop2, std::string Ref);
double _CoolProp_Fluid_Props(long iOutput, long iName1, double Value1, long iName2, double Value2, Fluid *pFluid, bool SinglePhase = false);

static std::string err_string;
static int debug_level=0;
static Fluid * pFluid;

// This is very hacky, but pull the subversion revision from the file
#include "svnrevision.h" // Contents are like "long svnrevision = 286;"
#include "version.h" // Contents are like "char version [] ="2.5";"

int global_Phase = -1;
bool global_SinglePhase = false;
bool global_SaturatedL = false;
bool global_SaturatedV = false;

// This is a map of all possible strings to a unique identifier
std::pair<std::string, long> map_data[] = {
	std::make_pair(std::string("E"),iPcrit),
	std::make_pair(std::string("M"),iMM),
	std::make_pair(std::string("w"),iAccentric),
	std::make_pair(std::string("R"),iTtriple),
	std::make_pair(std::string("N"),iRhocrit),
	std::make_pair(std::string("B"),iTcrit),

	std::make_pair(std::string("pcrit"),iPcrit),
	std::make_pair(std::string("molemass"),iMM),
	std::make_pair(std::string("accentric"),iAccentric),
	std::make_pair(std::string("dipole"),iDipole),
	std::make_pair(std::string("Tmin"),iTmin),
	std::make_pair(std::string("t"),iTmin),
	std::make_pair(std::string("Ttriple"),iTtriple),
	std::make_pair(std::string("ptriple"),iPtriple),
	std::make_pair(std::string("rhocrit"),iRhocrit),
	std::make_pair(std::string("Tcrit"),iTcrit),
	std::make_pair(std::string("Hcrit"),iHcrit),
	std::make_pair(std::string("Scrit"),iScrit),

	std::make_pair(std::string("Q"),iQ),
	std::make_pair(std::string("T"),iT),
    std::make_pair(std::string("P"),iP),
	std::make_pair(std::string("D"),iD),
	std::make_pair(std::string("C"),iC),
	std::make_pair(std::string("C0"),iC0),
	std::make_pair(std::string("O"),iO),
	std::make_pair(std::string("U"),iU),
	std::make_pair(std::string("H"),iH),
	std::make_pair(std::string("S"),iS),
	std::make_pair(std::string("A"),iA),
	std::make_pair(std::string("G"),iG),
	std::make_pair(std::string("V"),iV),
	std::make_pair(std::string("L"),iL),
	std::make_pair(std::string("I"),iI),
	std::make_pair(std::string("SurfaceTension"),iI),
	std::make_pair(std::string("dpdT"),iDpdT),
	std::make_pair(std::string("drhodT|p"),iDrhodT_p),
	std::make_pair(std::string("M"),iMM),
	std::make_pair(std::string("Phase"),iPhase),
	std::make_pair(std::string("PHASE_LIQUID"),iPHASE_LIQUID),
	std::make_pair(std::string("PHASE_GAS"),iPHASE_GAS),
	std::make_pair(std::string("PHASE_SUPERCRITICAL"),iPHASE_SUPERCRITICAL),
	std::make_pair(std::string("PHASE_TWOPHASE"),iPHASE_TWOPHASE),
	std::make_pair(std::string("ODP"),iODP),
	std::make_pair(std::string("GWP20"),iGWP20),
	std::make_pair(std::string("GWP100"),iGWP100),
	std::make_pair(std::string("GWP500"),iGWP500),
	std::make_pair(std::string("CritSplineT"),iCritSplineT),
};
//Now actually construct the map
std::map<std::string, long> param_map(map_data,
    map_data + sizeof map_data / sizeof map_data[0]);

// This is a map of all unique identifiers to std::strings with the default units
// that are used internally
std::pair<long, std::string> units_data[] = {
	std::make_pair(iPcrit, std::string("kPa")),
	std::make_pair(iMM, std::string("kg/kmol")),
	std::make_pair(iAccentric, std::string("-")),
	std::make_pair(iTtriple, std::string("K")),
	std::make_pair(iRhocrit, std::string("kg/m^3")),
	std::make_pair(iTcrit, std::string("K")),

	std::make_pair(iQ, std::string("")),
	std::make_pair(iT, std::string("K")),
    std::make_pair(iP, std::string("kPa")),
	std::make_pair(iD, std::string("kg/m^3")),
	std::make_pair(iC, std::string("kJ/kg/K")),
	std::make_pair(iC0, std::string("kJ/kg/K")),
	std::make_pair(iO, std::string("kJ/kg/K")),
	std::make_pair(iU, std::string("kJ/kg")),
	std::make_pair(iH, std::string("kJ/kg")),
	std::make_pair(iS, std::string("kJ/kg/K")),
	std::make_pair(iA, std::string("m/s")),
	std::make_pair(iG, std::string("kJ/kg")),
	std::make_pair(iV, std::string("Pa*s")),
	std::make_pair(iL, std::string("kW/m/K")),
	std::make_pair(iI, std::string("N/m")),
	std::make_pair(iDpdT, std::string("kPa/K")),
	std::make_pair(iDrhodT_p, std::string("kg/K/m^3"))
};

//Now actually construct the map
std::map<long, std::string> units_map(units_data,
		units_data + sizeof units_data / sizeof units_data[0]);

FluidsContainer Fluids = FluidsContainer();

EXPORT_CODE long CONVENTION get_svnrevision(){return svnrevision;}
EXPORT_CODE long CONVENTION get_version(char * pversion){
	strcpy(pversion,version);
	return 0;
}
std::string get_version(){return std::string(version);}

void set_debug(int level){debug_level = level;}
EXPORT_CODE int CONVENTION get_debug(){return debug_level;}
int  debug(){return debug_level;}
EXPORT_CODE void CONVENTION debug(int level){debug_level=level;}

std::string get_errstring(void){
    std::string temp = err_string;
    err_string = std::string("");
    return temp;
    }
EXPORT_CODE void CONVENTION get_errstring(char* str){
    str=(char*) get_errstring().c_str();
	err_string.clear();
    };

EXPORT_CODE char * CONVENTION get_errstringc(void){
    std::string temp = err_string;
	err_string.clear();
    return (char*)temp.c_str();
    }

EXPORT_CODE long CONVENTION get_errstring_copy(char* str){
    strcpy(str, get_errstring().c_str());
	err_string.clear();
	return strlen(str);
    };

// A function to enforce the state if known
EXPORT_CODE void CONVENTION set_phase(char *Phase_str){
	set_phase(std::string(Phase_str));
}

/// Enable the TTSE for this fluid
EXPORT_CODE bool CONVENTION enable_TTSE_LUT(char *FluidName){
	long iFluid = get_Fluid_index(FluidName); if (iFluid<0){ return false; };
	pFluid = Fluids.get_fluid(iFluid);
	pFluid->enable_TTSE_LUT();
	return true;
};
/// Check if TTSE is enabled
EXPORT_CODE bool CONVENTION isenabled_TTSE_LUT(char *FluidName){
	long iFluid = get_Fluid_index(FluidName); if (iFluid<0){ return false; };
	pFluid = Fluids.get_fluid(iFluid);
	return pFluid->isenabled_TTSE_LUT();
}
/// Disable the TTSE for this fluid
EXPORT_CODE bool CONVENTION disable_TTSE_LUT(char *FluidName){
	long iFluid = get_Fluid_index(FluidName); if (iFluid<0){ return false; };
	pFluid = Fluids.get_fluid(iFluid);
	pFluid->disable_TTSE_LUT();
	return true;
}
/// Enable the writing of TTSE tables to file for this fluid
EXPORT_CODE bool CONVENTION enable_TTSE_LUT_writing(char *FluidName){
	long iFluid = get_Fluid_index(FluidName); if (iFluid<0){ return true;};
	pFluid = Fluids.get_fluid(iFluid);
	pFluid->enable_TTSE_LUT_writing();
	return true;
};
/// Check if the writing of TTSE tables to file is enabled
EXPORT_CODE bool CONVENTION isenabled_TTSE_LUT_writing(char *FluidName){
	long iFluid = get_Fluid_index(FluidName); if (iFluid<0){ return false;};
	pFluid = Fluids.get_fluid(iFluid);
	return pFluid->isenabled_TTSE_LUT_writing();
}
/// Disable the writing of TTSE tables to file for this fluid
EXPORT_CODE bool CONVENTION disable_TTSE_LUT_writing(char *FluidName){
	long iFluid = get_Fluid_index(FluidName); if (iFluid<0){ return false;};
	pFluid = Fluids.get_fluid(iFluid);
	pFluid->disable_TTSE_LUT_writing();
	return true;
}
/// Over-ride the default size of both of the saturation LUT
EXPORT_CODE bool CONVENTION set_TTSESat_LUT_size(char *FluidName, int Nsat){
	long iFluid = get_Fluid_index(FluidName); if (iFluid<0){ return false; };
	pFluid = Fluids.get_fluid(iFluid);
	pFluid->set_TTSESat_LUT_size(Nsat);
	return true;
}
/// Over-ride the default size of the single-phase LUT
EXPORT_CODE bool CONVENTION set_TTSESinglePhase_LUT_size(char *FluidName, int Np, int Nh){
	long iFluid = get_Fluid_index(FluidName); if (iFluid<0){ return false;};
	pFluid = Fluids.get_fluid(iFluid);
	pFluid->set_TTSESinglePhase_LUT_size(Np,Nh);
	return true;
}
/// Over-ride the default range of the single-phase LUT
EXPORT_CODE bool CONVENTION set_TTSESinglePhase_LUT_range(char *FluidName, double hmin, double hmax, double pmin, double pmax){
	long iFluid = get_Fluid_index(FluidName); if (iFluid<0){ return false;};
	pFluid = Fluids.get_fluid(iFluid);
	pFluid->set_TTSESinglePhase_LUT_range(hmin,hmax,pmin,pmax);
	return true;
}
/// Get the current range of the single-phase LUT
EXPORT_CODE bool CONVENTION get_TTSESinglePhase_LUT_range(char *FluidName, double *hmin, double *hmax, double *pmin, double *pmax){
	long iFluid = get_Fluid_index(FluidName); if (iFluid<0){ return false;};
	pFluid = Fluids.get_fluid(iFluid);
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

void set_phase(std::string Phase_str){
	if (!Phase_str.compare("Two-Phase")){
		global_SinglePhase = false;
		global_SaturatedL = false;
		global_SaturatedV = false;
		global_Phase = iTwoPhase;
	}
	else if (!Phase_str.compare("Liquid")){
		global_SinglePhase = true;
		global_SaturatedL = false;
		global_SaturatedV = false;
		global_Phase = iLiquid;
	}
	else if (!Phase_str.compare("Gas")){
		global_SinglePhase = true;
		global_SaturatedL = false;
		global_SaturatedV = false;
		global_Phase = iGas;
	}
	else if (!Phase_str.compare("Supercritical")){
		global_SinglePhase = true;
		global_SaturatedL = false;
		global_SaturatedV = false;
		global_Phase = iSupercritical;
	}
	else if (!Phase_str.compare("SaturatedL")){
		global_SinglePhase = false;
		global_SaturatedL = true;
		global_SaturatedV = false;
		global_Phase = iTwoPhase;
	}
	else if (!Phase_str.compare("SaturatedV")){
		global_SinglePhase = false;
		global_SaturatedL = false;
		global_SaturatedV = true;
		global_Phase = iTwoPhase;
	}
}

// Returns a pointer to the fluid class
Fluid* get_fluid(long iFluid){
	return Fluids.get_fluid(iFluid);
}
long get_Fluid_index(std::string FluidName)
{
	// Try to get the fluid from Fluids by name
	pFluid = Fluids.get_fluid(FluidName);
	// If NULL, didn't find it (or its alias)
	if (pFluid!=NULL)
	{
		// Find the fluid index
		return Fluids.get_fluid_index(pFluid);
	}
	else
		return -1;
}
EXPORT_CODE long CONVENTION get_Fluid_index(char * param)
{
	return get_Fluid_index(std::string(param));
}

bool add_REFPROP_fluid(std::string FluidName)
{
	double x[100];
	// Starts with REFPROP- keep going
	if (!(FluidName.find("REFPROP-") == 0)) return false;
	// Stop here if there is no REFPROP support
	if (!REFPROPFluidClass::refpropSupported()) return false;
	// Try to load this fluid, index >= 0 if already added
	long iFluid = get_Fluid_index(FluidName);
	// If not added yet, and a valid fluid, then continue
	if (iFluid < 0 && set_REFPROP_fluid(FluidName, x)) // If you can set the fluid, it's a valid fluid
	{
		Fluids.add_REFPROP_fluid(FluidName,std::vector<double>(1,1));
		return true;
		std::cout << format("Added the fluid %s\n",FluidName.c_str()).c_str();
	}
	return true;
}
std::string get_index_units(long index)
{
	std::map<long, std::string>::iterator it;
	// Try to find using the map
	it = units_map.find(index);
	// If it is found the iterator will not be equal to end
	if (it != units_map.end() )
	{
		// Return the index of the parameter
		return (*it).second;
	}
	else
	{
		return std::string("Didn't match parameter");
	}
}
EXPORT_CODE void CONVENTION get_index_units(long param, char * units)
{
	strcpy(units, (char*)get_index_units(param).c_str());
	return;
}

std::string get_ASHRAE34(std::string fluid)
{
	long iFluid = get_Fluid_index(fluid);
	if (iFluid > -1)
	{
		Fluid *pFluid = get_fluid(iFluid);
		return pFluid->environment.ASHRAE34;
	}
	else
	{
		return "Fluid name invalid";
	}
}

EXPORT_CODE int CONVENTION get_TTSE_mode(char* fluid, char *value)
{
	long iFluid = get_Fluid_index(fluid);
	if (iFluid > -1)
	{
		Fluid *pFluid = get_fluid(iFluid);
		int iMode = pFluid->TTSESinglePhase.get_mode();
		switch (iMode)
		{
		case TTSE_MODE_TTSE:
			strcpy(value,"TTSE"); 
			std::cout << "set TTSE mode to "<< value << std::endl;
			return true;
		case TTSE_MODE_BICUBIC:
			strcpy(value,"BICUBIC"); 
			std::cout << "set TTSE mode to "<< value << std::endl;
			return true;
		default:
			strcpy(value,"UNKNOWN"); return false;
		}
	}
	else
	{
		return false;
	}
}
std::string get_TTSE_mode(std::string fluid)
{
	long iFluid = get_Fluid_index(fluid);
	if (iFluid > -1)
	{
		Fluid *pFluid = get_fluid(iFluid);
		int iMode = pFluid->TTSESinglePhase.get_mode();
		switch (iMode)
		{
		case TTSE_MODE_TTSE:
			return "TTSE";
		case TTSE_MODE_BICUBIC:
			return "BICUBIC"; 
		default:
			return "UNKNOWN";
		}
	}
	else
	{
		return false;
	}
}

EXPORT_CODE int CONVENTION set_TTSE_mode(char* fluid, char *value)
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


EXPORT_CODE long CONVENTION get_ASHRAE34(char* fluid, char *value)
{
	strcpy(value, (char*)get_ASHRAE34(fluid).c_str());
	return 1;
}
std::string get_CAS_code(std::string fluid)
{
	long iFluid = get_Fluid_index(fluid);
	if (iFluid > -1)
	{
		Fluid *pFluid = get_fluid(iFluid);
		return pFluid->params.CAS;
	}
	else
	{
		return "Fluid name invalid";
	}
}
EXPORT_CODE long CONVENTION get_CAS_code(char* fluid, char *value)
{
	strcpy(value, (char*)get_CAS_code(fluid).c_str());
	return 1;
}

long get_param_index(std::string param)
{
	std::map<std::string,long>::iterator it;
	// Try to find using the map
	it = param_map.find(param);
	// If it is found the iterator will not be equal to end
	if (it != param_map.end() )
	{
		// Return the index of the parameter
		return (*it).second;
	}
	else
	{
		return -1;
	}
}
EXPORT_CODE long CONVENTION get_param_index(char * param)
{
	return get_param_index(std::string(param));
}
static int IsCoolPropFluid(std::string FluidName)
{
	// Try to get the fluid from Fluids by name
	try
	{
		pFluid = Fluids.get_fluid(FluidName);
	}
	catch (NotImplementedError)
	{
		return false;
	}
	// If NULL, didn't find it (or its alias)
	if (pFluid!=NULL)
	{
		return true;
	}
	else
		return false;
}

static int IsBrine(char* Ref)
{
	// First check whether it is one of the Brines that does
	// not have a pure-fluid equivalent in CoolProp
    if (
        strcmp(Ref,"HC-10")==0 || 
        strncmp(Ref,"PG",2)==0 || 
		strncmp(Ref,"EG",2)==0 || 
		strncmp(Ref,"EA",2)==0 ||
		strncmp(Ref,"MA",2)==0 || 
		strncmp(Ref,"Glycerol",8)==0 ||
		strncmp(Ref,"K2CO3",5)==0 || 
		strncmp(Ref,"CaCl2",5)==0 || 
		strncmp(Ref,"MgCl2",5)==0 || 
		strncmp(Ref,"NaCl",4)==0 || 
		strncmp(Ref,"KAC",3)==0 || 
		strncmp(Ref,"KFO",3)==0 || 
		strncmp(Ref,"LiCl",3)==0 || 
        strncmp(Ref,"NH3/H2O",7)==0
       )
    {
        return 1;
    }
	// Then check for diluants that are also pure fluids in CoolProp
	else if ( (strncmp(Ref,"Methanol",8)==0 && Ref[8] == '-') ||
		      (strncmp(Ref,"Ethanol",7)==0 && Ref[7] == '-') ||
			  (strncmp(Ref,"NH3",3)==0 && Ref[3] == '-')
		)
	{
		return 1;
	}
    else
    {
        return 0;
    }
}
static int IsREFPROP(std::string Ref)
{
    if (!Ref.compare(0,8,"REFPROP-"))
        return 1;
    else
        return 0;
}
EXPORT_CODE int CONVENTION IsFluidType(char *Ref, char *Type)
{
	pFluid = Fluids.get_fluid(Ref);

	if (IsBrine(Ref)){
		if (!strcmp(Type,"Brine")){
			return 1;
		}
		else{
			return 0;
		}
	}
	else if (IsIncompressibleLiquid(Ref)){
		if (!strcmp(Type,"Liquid")){
			return 1;
		}
		else{
			return 0;
		}
	}
	else if (IsREFPROP(Ref)){
		if (!strcmp(Type,"PureFluid")){
			return 1;
		}
		else{
			return 0;
		}
	}
	else if (!pFluid->pure()){
		if (!strcmp(Type,"PseudoPure") || !strcmp(Type,"PseudoPureFluid")){
			return 1;
		}
		else{
			return 0;
		}
	}
	else if (pFluid->pure()){
		if (!strcmp(Type,"PureFluid")){
			return 1;
		}
		else{
			return 0;
		}
	}
    else
    {
        return 0;
    }
}
EXPORT_CODE double CONVENTION conformal_Trho(char* FluidName, char* ReferenceFluidName, double T, double rho, double *Tconform, double *rhoconform)
{
	long iFluid = get_Fluid_index(FluidName);
	long iRefFluid = get_Fluid_index(ReferenceFluidName);
	if (iFluid < 0 || iRefFluid < 0)
	{
		return _HUGE;
	}
	else
	{
		try{
			Fluid *pFluid = get_fluid(iFluid);
			Fluid *pRefFluid = get_fluid(iRefFluid);

			double rhobar=rho/pFluid->params.molemass;
			double rhocbar=pFluid->reduce.rho/pFluid->params.molemass;
			double Tc=pFluid->reduce.T;

			double Tc0 = pRefFluid->reduce.T;
			double rhoc0bar=pRefFluid->reduce.rho/pRefFluid->params.molemass;
			double T0=T*Tc0/Tc;
			double rho0bar = rhobar*rhoc0bar/rhocbar;  // Must be the ratio of MOLAR densities!!
			double rho0 = rho0bar*pRefFluid->params.molemass;
			std::string errstring;
			std::vector<double> Trho = pFluid->ConformalTemperature(pFluid,pRefFluid,T,rho,T0,rho0,&errstring);
			if (errstring.size()>0){
				return _HUGE;
			}
			else{
				*Tconform = Trho[0];
				*rhoconform = Trho[1];
				return 0;
			}
		}
		catch (std::exception){
			return _HUGE;
		}
	}
}

EXPORT_CODE double CONVENTION viscosity_dilute(char* FluidName, double T)
{
	long iFluid = get_Fluid_index(FluidName);
	if (iFluid < 0)
	{
		return _HUGE;
	}
	else
	{
		double e_k, sigma;
		Fluid *pFluid = get_fluid(iFluid);
		pFluid->ECSParams(&e_k, &sigma);
		return pFluid->viscosity_dilute(T,e_k,sigma);
	}
}
EXPORT_CODE double CONVENTION viscosity_residual(char* FluidName, double T, double rho)
{
	long iFluid = get_Fluid_index(FluidName);
	if (iFluid < 0)
	{
		return _HUGE;
	}
	else
	{
		Fluid *pFluid = get_fluid(iFluid);
		try{
			return pFluid->viscosity_residual(T,rho);
		}
		catch (NotImplementedError)
		{
			return _HUGE;
		}
	}
}

EXPORT_CODE double CONVENTION conductivity_critical(char* FluidName, double T, double rho)
{
	long iFluid = get_Fluid_index(FluidName);
	if (iFluid < 0)
	{
		return _HUGE;
	}
	else
	{
		Fluid *pFluid = get_fluid(iFluid);
		return pFluid->conductivity_critical(T,rho);
	}
}
EXPORT_CODE double CONVENTION conductivity_background(char* FluidName, double T, double rho)
{
	long iFluid = get_Fluid_index(FluidName);
	if (iFluid < 0)
	{
		return _HUGE;
	}
	else
	{
		Fluid *pFluid = get_fluid(iFluid);
		return pFluid->conductivity_background(T,rho);
	}
}

EXPORT_CODE double CONVENTION rhosatL_anc(char* Fluid, double T)
{
	try{
		// Try to load the CoolProp Fluid
		pFluid = Fluids.get_fluid(std::string(Fluid));
		return pFluid->rhosatL(T);
	}
	catch(NotImplementedError){
		return -_HUGE;
	}
}
EXPORT_CODE double CONVENTION rhosatV_anc(char* Fluid, double T)
{
	try{
		// Try to load the CoolProp Fluid
		pFluid = Fluids.get_fluid(std::string(Fluid));
		return pFluid->rhosatV(T);
	}
	catch(NotImplementedError){
		return -_HUGE;
	}
}
EXPORT_CODE double CONVENTION psatL_anc(char* Fluid, double T)
{
	try{
		// Try to load the CoolProp Fluid
		pFluid = Fluids.get_fluid(std::string(Fluid));
		return pFluid->psatL_anc(T);
	}
	catch(NotImplementedError){
		return -_HUGE;
	}
}
EXPORT_CODE double CONVENTION psatV_anc(char* Fluid, double T)
{
	try{
		// Try to load the CoolProp Fluid
		pFluid = Fluids.get_fluid(std::string(Fluid));
		return pFluid->psatV_anc(T);
	}
	catch(NotImplementedError){
		return -_HUGE;
	}
}

double _T_hp_secant(std::string Ref, double h, double p, double T_guess)
{
    double x1=0,x2=0,x3=0,y1=0,y2=0,eps=1e-8,change=999,f=999,T=300;
    int iter=1;

    while ((iter<=3 || fabs(f)>eps) && iter<100)
    {
        if (iter==1){x1=T_guess; T=x1;}
        if (iter==2){x2=T_guess+0.1; T=x2;}
        if (iter>2) {T=x2;}
            f=Props("H",'T',T,'P',p,Ref)-h;
        if (iter==1){y1=f;}
        if (iter>1)
        {
            y2=f;
            x3=x2-y2/(y2-y1)*(x2-x1);
            change=fabs(y2/(y2-y1)*(x2-x1));
            y1=y2; x1=x2; x2=x3;
        }
        iter=iter+1;
        if (iter>100)
        {
			throw SolutionError(format("iter %d: T_hp not converging with inputs(%s,%g,%g,%g) value: %0.12g\n",iter,(char*)Ref.c_str(),h,p,T_guess,f));
        }
    }
    return T;
}

EXPORT_CODE void CONVENTION Phase(char *Fluid,double T, double p, char *Phase_str)
{
	strcpy(Phase_str,(char*)Phase(std::string(Fluid),T,p).c_str());
}
EXPORT_CODE long CONVENTION Phase_Tp(char *Fluid,double T, double p, char *Phase_str)
{
	strcpy(Phase_str,(char*)Phase(std::string(Fluid),T,p).c_str());
	return 0;
}
EXPORT_CODE long CONVENTION Phase_Trho(char *Fluid,double T, double rho, char *Phase_str)
{
	strcpy(Phase_str,(char*)Phase_Trho(std::string(Fluid),T,rho).c_str());
	return 0;
}

std::string Phase_Trho(std::string Fluid, double T, double rho)
{
	try{
		// Try to load the CoolProp Fluid
		pFluid = Fluids.get_fluid(Fluid);
		double pL,pV,rhoL,rhoV;
		return pFluid->phase_Trho(T,rho, &pL, &pV, &rhoL, &rhoV);
	}
	catch(NotImplementedError){
		return std::string("");
	}
}

std::string Phase(std::string Fluid, double T, double p)
{
	try{
		// Try to load the CoolProp Fluid
		pFluid = Fluids.get_fluid(Fluid);
		double pL,pV,rhoL,rhoV;
		return pFluid->phase_Tp(T, p, &pL, &pV, &rhoL, &rhoV);
	}
	catch(NotImplementedError){
		return std::string("");
	}
}

std::string Phase_Tp(std::string Fluid, double T, double p)
{
	return Phase(Fluid,T,p);
}

double _Props1(char *Fluid, char *Output)
{
	// Try to load the CoolProp Fluid
	pFluid = Fluids.get_fluid(Fluid);
	if (pFluid != NULL)
	{
		// It's a CoolProp fluid
		// Convert the parameter to integer
		long iOutput = get_param_index(Output);
		if (iOutput < 0){
			throw ValueError(format("Your output key [%s] is not valid. (names are case sensitive)",Output));
		}
		// Get the output using the conventional function
		return _CoolProp_Fluid_Props(iOutput,0,0,0,0,pFluid);
	}
	else if (IsREFPROP(Fluid))
	{
		// REFPROP fluid, or something else that is invalid
		try{
			if (!strcmp(Output,"Ttriple"))
				return Props('R','T',0,'P',0,Fluid);
			else if (!strcmp(Output,"Tcrit"))
				return Props('B','T',0,'P',0,Fluid);
			else if (!strcmp(Output,"pcrit"))
				return Props('E','T',0,'P',0,Fluid);
			else if (!strcmp(Output,"Tmin"))
				return Props('t','T',0,'P',0,Fluid);
			else if (!strcmp(Output,"molemass"))
				return Props('M','T',0,'P',0,Fluid);
			else if (!strcmp(Output,"rhocrit"))
				return Props('N','T',0,'P',0,Fluid);
			else if (!strcmp(Output,"accentric"))
				return Props('w','T',0,'P',0,Fluid);
			else 
				throw ValueError(format("Output parameter \"%s\" is invalid for REFPROP fluid",Output));
		}
		// Catch any error that subclasses the std::exception
		catch(std::exception &e){
			err_string = std::string("CoolProp error: ").append(e.what());
			return _HUGE;
		}
	}
	else
	{
		throw ValueError(format("Fluid \"%s\" is an invalid fluid",Fluid));
	}
}

// All the function interfaces that point to the single-input Props function
EXPORT_CODE double CONVENTION Props1(char* Ref, char * Output)
{
	//FILE *fp;
	//fp = fopen("c:\\CoolProp\\log_Props1.txt", "a");
	//fprintf(fp,"%s %s\n",Ref,Output);
	//fclose(fp);

	try{
	// Redirect to the Props function - should have called it Props1 from the outset
	return _Props1(Ref, Output);
	}
	catch(std::exception &e)
	{
		err_string = std::string("CoolProp error: ").append(e.what());
		return _HUGE;
	}
	catch(...){
		err_string = std::string("CoolProp error: Indeterminate error");
		return _HUGE;
	}
}
double Props1(std::string Ref, std::string Output)
{
	// Redirect to the Props function - should have called it Props1 from the outset
	return Props1((char*)Ref.c_str(), (char*)Output.c_str());
}
double Props(std::string Ref, std::string Output)
{
	return Props1((char*)Ref.c_str(), (char*)Output.c_str());
}
double Props(char* Ref, char* Output)
{
	return Props1(Ref, Output);
}

EXPORT_CODE double CONVENTION Props(char *Output,char Name1, double Prop1, char Name2, double Prop2, char * Ref)
{
	double val = Props(std::string(Output),Name1,Prop1,Name2,Prop2,std::string(Ref));

	//FILE *fp;
	//fp = fopen("c:\\CoolProp\\log_Props.txt", "a");
	//fprintf(fp,"%s,%c,%g,%c,%g,%s-->%g\n",Output,Name1,Prop1,Name2,Prop2,Ref,val);
	//fclose(fp);

	// Go to the std::string, std::string version
	return val;
}
double Props(char Output,char Name1, double Prop1, char Name2, double Prop2, char* Ref)
{
	// Go to the std::string, std::string version
	return Props(std::string(1,Output),Name1,Prop1,Name2,Prop2,std::string(Ref));
}
double Props(std::string Output,char Name1, double Prop1, char Name2, double Prop2, std::string Ref)
{
	// In this function the error catching happens;
	try{
		return _Props(Output,std::string(1,Name1),Prop1,std::string(1,Name2),Prop2,Ref);
	}
	catch(const CoolPropBaseError& e){
		err_string = std::string("CoolProp error: ").append(e.what());
		return _HUGE;
	}
	catch(const std::exception& e){
			err_string = std::string("CoolProp error: ").append(e.what());
			return _HUGE;
		}
	catch(...){
		err_string = std::string("CoolProp error: Indeterminate error");
		return _HUGE;
	}
}
// Make this a wrapped function so that error bubbling can be done properly
double _Props(std::string Output,std::string Name1, double Prop1, std::string Name2, double Prop2, std::string Ref)
{
	if (debug()>5){
		std::cout<<__FILE__<<": "<<Output.c_str()<<","<<Name1.c_str()<<","<<Prop1<<","<<Name2.c_str()<<","<<Prop2<<","<<Ref.c_str()<<std::endl;
	}
	/* 
    If the fluid name is not actually a refrigerant name, but a string beginning with "REFPROP-",
    then REFPROP is used to calculate the desired property.
    */
    if (IsREFPROP(Ref))  // First eight characters match "REFPROP-"
    {
        // Stop here if there is no REFPROP support
    	if (REFPROPFluidClass::refpropSupported()) {
			return REFPROP(Output,Name1,Prop1,Name2,Prop2,Ref);
    	} else {
    		throw AttributeError(format("Your refrigerant [%s] is from REFPROP, but CoolProp does not support REFPROP on this platform, yet.",Ref.c_str()));
    		return -_HUGE;
    	}
    }
	else if (IsCoolPropFluid(Ref))
	{
		pFluid = Fluids.get_fluid(Ref);
		// Convert all the parameters to integers
		long iOutput = get_param_index(Output);
		if (iOutput<0) 
			throw ValueError(format("Your output key [%s] is not valid. (names are case sensitive)",Output.c_str()));
		long iName1 = get_param_index(std::string(Name1));  
		if (iName1<0) 
			throw ValueError(format("Your input key #1 [%s] is not valid. (names are case sensitive)",Name1.c_str()));
		long iName2 = get_param_index(std::string(Name2));  
		if (iName2<0) 
			throw ValueError(format("Your input key #2 [%s] is not valid. (names are case sensitive)",Name2.c_str()));
		// Call the internal method that uses the parameters converted to longs
		return _CoolProp_Fluid_Props(iOutput,iName1,Prop1,iName2,Prop2,pFluid);
	}
    

    // It's a brine, call the brine routine
	else if (IsBrine((char*)Ref.c_str()))
    {
		//Enthalpy and pressure are the inputs
		if ((Name1.c_str()[0]=='H' && Name2.c_str()[0]=='P') || (Name2.c_str()[0]=='H' && Name1.c_str()[0]=='P'))
        {
			if (Name2.c_str()[0]=='H' && Name1.c_str()[0]=='P')
			{
				std::swap(Prop1,Prop2);
				std::swap(Name1,Name2);
			}
			// Start with a guess of 10 K below max temp of fluid
			double Tguess = SecFluids('M',Prop1,Prop2,(char*)Ref.c_str())-10;
			// Solve for the temperature
			double T =_T_hp_secant(Ref,Prop1,Prop2,Tguess);
			// Return whatever property is desired
			return SecFluids(Output[0],T,Prop2,(char*)Ref.c_str());
		}
		else if ((Name1.c_str()[0] == 'T' && Name2.c_str()[0] =='P') || (Name1.c_str()[0] == 'P' && Name2.c_str()[0] == 'T'))
        {
			if (Name1.c_str()[0] =='P' && Name2.c_str()[0] =='T'){
				std::swap(Prop1,Prop2);
			}
			return SecFluids(Output[0],Prop1,Prop2,(char*)Ref.c_str());
		}
		else
		{
			throw ValueError("For brine, inputs must be (order doesnt matter) 'T' and 'P', or 'H' and 'P'");
		}
    }
	// It's an incompressible liquid, call the routine
	else if (IsIncompressibleLiquid((char*)Ref.c_str()))
    {
		//Enthalpy and pressure are the inputs
		if ((Name1.c_str()[0]=='H' && Name2.c_str()[0]=='P') || (Name2.c_str()[0]=='H' && Name1.c_str()[0]=='P'))
        {
			if (Name2.c_str()[0]=='H' && Name1.c_str()[0]=='P')
			{
				std::swap(Prop1,Prop2);
				std::swap(Name1,Name2);
			}
			
			// Solve for the temperature
			double T =_T_hp_secant(Ref,Prop1,Prop2,300);
			// Return whatever property is desired
			return IncompLiquid(get_param_index(Output),T,Prop2,(char*)Ref.c_str());
		}
		else if ((Name1.c_str()[0] == 'T' && Name2.c_str()[0] =='P') || (Name1.c_str()[0] == 'P' && Name2.c_str()[0] == 'T'))
        {
			if (Name1.c_str()[0] =='P' && Name2.c_str()[0] =='T'){
				std::swap(Prop1,Prop2);
			}
			return IncompLiquid(get_param_index(Output),Prop1,Prop2,Ref);
		}
		else
		{
			throw ValueError("For brine, inputs must be (order doesnt matter) 'T' and 'P', or 'H' and 'P'");
		}
    }
	else
	{
		throw ValueError(format("Your fluid name [%s] is not a CoolProp fluid, a REFPROP fluid, a brine or a liquid",Ref.c_str()));
	}
}
double _CoolProp_Fluid_Props(long iOutput, long iName1, double Prop1, long iName2, double Prop2, Fluid *pFluid, bool SinglePhase)
{
	double val = _HUGE;
	// This private method uses the indices directly for speed

	if (debug()>3){
		std::cout<<__FILE__<<" _CoolProp_Fluid_Props: "<<iOutput<<","<<iName1<<","<<Prop1<<","<<iName2<<","<<Prop2<<","<<pFluid->get_name().c_str()<<std::endl;
	}

	// Generate a State instance wrapped around the Fluid instance
	CoolPropStateClass CPS = CoolPropStateClass(pFluid);

	// Check if it is an output that doesn't require a state input
	// Deal with it and return

	switch (iOutput)
	{
		case iMM:
		case iPcrit:
		case iTcrit:
		case iTtriple:
		case iPtriple:
		case iRhocrit:
		case iTmin:
		case iAccentric:
		case iPHASE_LIQUID:
		case iPHASE_GAS:
		case iPHASE_SUPERCRITICAL:
		case iPHASE_TWOPHASE:
		case iGWP20:
		case iGWP100:
		case iGWP500:
		case iODP:
		case iCritSplineT:
		case iScrit:
		case iHcrit:
			return CPS.keyed_output(iOutput);
	}

	// Update the class
	CPS.update(iName1,Prop1,iName2,Prop2);
	// Get the output
	val = CPS.keyed_output(iOutput);
	
	if (debug()>5){
		std::cout<<__FILE__<<" _CoolProp_Fluid_Props return: "<<val<<std::endl;
	}
	// Return the value
	return val;
}
EXPORT_CODE double CONVENTION IProps(long iOutput, long iName1, double Prop1, long iName2, double Prop2, long iFluid)
{
	pFluid = Fluids.get_fluid(iFluid);
	// Didn't work
	if (pFluid == NULL){
		err_string=std::string("CoolProp error: ").append(format("%d is an invalid fluid index to IProps",iFluid));
		return _HUGE;
	}
	else{
		// In this function the error catching happens;
		try{
			return _CoolProp_Fluid_Props(iOutput,iName1,Prop1,iName2,Prop2,pFluid);
		}
		catch(std::exception &e){
			err_string=std::string("CoolProp error: ").append(e.what());
			return _HUGE;
		}
		catch(...){
			err_string=std::string("CoolProp error: Indeterminate error");
			return _HUGE;
		}
	}
}

EXPORT_CODE double CONVENTION K2F(double T)
{ return T * 9 / 5 - 459.67; }

EXPORT_CODE double CONVENTION F2K(double T_F)
{ return (T_F + 459.67) * 5 / 9;}

EXPORT_CODE void CONVENTION PrintSaturationTable(char *FileName, char * Ref, double Tmin, double Tmax)
{
    double T;
    FILE *f;
    f=fopen(FileName,"w");
    fprintf(f,"T,pL,pV,rhoL,rhoV,uL,uV,hL,hV,sL,sV,cvL,cvV,cpL,cpV,kL,kV,muL,muV\n");
    fprintf(f,"[K],[kPa],[kPa],[kg/m^3],[kg/m^3],[kJ/kg],[kJ/kg],[kJ/kg],[kJ/kg],[kJ/kg-K],[kJ/kg-K],[kJ/kg-K],[kJ/kg-K],[kJ/kg-K],[kJ/kg-K],[kW/m-K],[kW/m-K],[Pa-s],[Pa-s]\n");
    fprintf(f,"--------------------------------------------------------------------------\n");

    for (T=Tmin;T<Tmax;T+=1.0)
    {
    fprintf(f,"%0.3f,",T);
    fprintf(f,"%0.6f,",Props('P','T',T,'Q',0,Ref));
    fprintf(f,"%0.6f,",Props('P','T',T,'Q',1,Ref));
    fprintf(f,"%0.6f,",Props('D','T',T,'Q',0,Ref));
    fprintf(f,"%0.6f,",Props('D','T',T,'Q',1,Ref));
    fprintf(f,"%0.6f,",Props('U','T',T,'Q',0,Ref));
    fprintf(f,"%0.6f,",Props('U','T',T,'Q',1,Ref));
    fprintf(f,"%0.6f,",Props('H','T',T,'Q',0,Ref));
    fprintf(f,"%0.6f,",Props('H','T',T,'Q',1,Ref));
    fprintf(f,"%0.6f,",Props('S','T',T,'Q',0,Ref));
    fprintf(f,"%0.6f,",Props('S','T',T,'Q',1,Ref));
    fprintf(f,"%0.6f,",Props('O','T',T,'Q',0,Ref));
    fprintf(f,"%0.6f,",Props('O','T',T,'Q',1,Ref));
    fprintf(f,"%0.6f,",Props('C','T',T,'Q',0,Ref));
    fprintf(f,"%0.6f,",Props('C','T',T,'Q',1,Ref));
    fprintf(f,"%0.9f,",Props('L','T',T,'Q',0,Ref));
    fprintf(f,"%0.9f,",Props('L','T',T,'Q',1,Ref));
    fprintf(f,"%0.6g,",Props('V','T',T,'Q',0,Ref));
    fprintf(f,"%0.6g,",Props('V','T',T,'Q',1,Ref));
    fprintf(f,"\n");
    }
    fclose(f);
}

EXPORT_CODE void CONVENTION FluidsList(char* str)
{
	str=(char*)FluidsList().c_str();
	return;
}
std::string FluidsList()
{
	return Fluids.FluidList();
}

EXPORT_CODE double CONVENTION DerivTerms(char *Term, double T, double rho, char * Ref)
{
	pFluid=Fluids.get_fluid(Ref);
	return DerivTerms(Term,T,rho,pFluid);
}

/// Calculate some interesting derivatives
double DerivTerms(char *Term, double T, double rho, Fluid * pFluid)
{
    double rhoc =pFluid->reduce.rho;
	double delta=rho/rhoc;
	double tau=pFluid->reduce.T/T;
	double R=pFluid->R();
	double dtau_dT = -pFluid->reduce.T/T/T;

	if (!strcmp(Term,"dhdp|rho") || !strcmp(Term,"dhdp|v"))
	{
		// All the terms are re-calculated here in order to cut down on the number of calls
		double dphir_dDelta = pFluid->dphir_dDelta(tau,delta);	
		double d2phir_dDelta_dTau = pFluid->d2phir_dDelta_dTau(tau,delta);
		double d2phir_dDelta2 = pFluid->d2phir_dDelta2(tau,delta);
		
		double dpdrho = R*T*(1+2*delta*dphir_dDelta+delta*delta*d2phir_dDelta2);
		double dpdT = R*rho*(1+delta*dphir_dDelta-delta*tau*d2phir_dDelta_dTau);
		double cp = -tau*tau*R*(pFluid->d2phi0_dTau2(tau,delta)+pFluid->d2phir_dTau2(tau,delta))+T/rho/rho*(dpdT*dpdT)/dpdrho;
		double drhodT = -dpdT/dpdrho;
		return -cp/dpdrho/drhodT-T*drhodT*(-1/pow(rho,2))+1/rho;
	}
	else if (!strcmp(Term,"Z"))
	{
		return 1+delta*pFluid->dphir_dDelta(tau,delta);
	}
	else if (!strcmp(Term,"dZ_dDelta"))
    {
        return delta*pFluid->d2phir_dDelta2(tau,delta)+pFluid->dphir_dDelta(tau,delta);
    }
	else if (!strcmp(Term,"dZ_dTau"))
    {
        return delta*pFluid->d2phir_dDelta_dTau(tau,delta);
    }
	else if (!strcmp(Term,"B"))
	{
		// given by B*rhoc=lim(delta --> 0) [dphir_ddelta(tau)]
		return 1.0/rhoc*pFluid->dphir_dDelta(tau,1e-12);
	}
	else if (!strcmp(Term,"dBdT"))
	{
		return 1.0/rhoc*pFluid->d2phir_dDelta_dTau(tau,1e-12)*dtau_dT;
	}
	else if (!strcmp(Term,"C"))
	{
		// given by C*rhoc^2=lim(delta --> 0) [d2phir_dDelta2(tau)]
		return 1.0/(rhoc*rhoc)*pFluid->d2phir_dDelta2(tau,1e-12);
	}
	else if (!strcmp(Term,"dCdT"))
	{
		return 1.0/(rhoc*rhoc)*pFluid->d3phir_dDelta2_dTau(tau,1e-12)*dtau_dT;
    }
	else if (!strcmp(Term,"phir"))
    {
        return pFluid->phir(tau,delta);
    }
	else if (!strcmp(Term,"dphir_dTau"))
    {
        return pFluid->dphir_dTau(tau,delta);
    }
	else if (!strcmp(Term,"dphir_dDelta"))
    {
        return pFluid->dphir_dDelta(tau,delta);
    }
	else if (!strcmp(Term,"d2phir_dTau2"))
    {
        return pFluid->d2phir_dTau2(tau,delta);
    }
	else if (!strcmp(Term,"d2phir_dDelta2"))
    {
        return pFluid->d2phir_dDelta2(tau,delta);
    }
	else if (!strcmp(Term,"d2phir_dDelta_dTau"))
    {
        return pFluid->d2phir_dDelta_dTau(tau,delta);
    }
	else if (!strcmp(Term,"d3phir_dDelta3"))
    {
        return pFluid->d3phir_dDelta3(tau,delta);
    }
	else if (!strcmp(Term,"d3phir_dDelta2_dTau"))
    {
        return pFluid->d3phir_dDelta2_dTau(tau,delta);
    }
	else if (!strcmp(Term,"d3phir_dDelta_dTau2"))
    {
        return pFluid->d3phir_dDelta_dTau2(tau,delta);
    }
	else if (!strcmp(Term,"d3phir_dTau3"))
    {
        return pFluid->d3phir_dTau3(tau,delta);
    }
	else if (!strcmp(Term,"phi0"))
    {
        return pFluid->phi0(tau,delta);
    }
    else if (!strcmp(Term,"dphi0_dTau"))
    {
        return pFluid->dphi0_dTau(tau,delta);
    }
	else if (!strcmp(Term,"d2phi0_dTau2"))
    {
        return pFluid->d2phi0_dTau2(tau,delta);
    }
	else if (!strcmp(Term,"dphi0_dDelta"))
    {
        return pFluid->dphi0_dDelta(tau,delta);
    }
	else if (!strcmp(Term,"d2phi0_dDelta2"))
    {
        return pFluid->d2phi0_dDelta2(tau,delta);
    }
	else if (!strcmp(Term,"d2phi0_dDelta_dTau"))
    {
        return pFluid->d2phi0_dDelta_dTau(tau,delta);
    }
	else if (!strcmp(Term,"d3phi0_dTau3"))
    {
        return pFluid->d3phi0_dTau3(tau,delta);
    }

	CoolPropStateClass CPS = CoolPropStateClass(pFluid);
	CPS.update(iT,T,iD,rho);

	if (!strcmp(Term,"dpdT") || !strcmp(Term,"dpdT|rho")){
		return CPS.dpdT_constrho();
	}
    else if (!strcmp(Term,"dpdrho") || !strcmp(Term,"dpdrho|T")){
		return CPS.dpdrho_constT();
	}
	else if (!strcmp(Term,"dhdT") || !strcmp(Term,"dhdT|rho")){
		return CPS.dhdT_constrho();
	}
	else if (!strcmp(Term,"dhdrho") || !strcmp(Term,"dhdrho|T")){
		return CPS.dhdrho_constT();
	}
	else if (!strcmp(Term,"drhodT|p")){
		return CPS.drhodT_constp();
	}
	else if (!strcmp(Term,"drhodh|p")){
		return CPS.drhodh_constp();
	}
	else if (!strcmp(Term,"drhodp|h")){
		return CPS.drhodp_consth();
	}
	else if (!strcmp(Term,"IsothermalCompressibility")){
		return CPS.isothermal_compressibility();
	}
	else
	{
		printf("Sorry DerivTerms is a work in progress, your derivative term [%s] is not available!!",Term);
		return _HUGE;
	}
}

std::string get_BibTeXKey(std::string Ref, std::string item)
{
	pFluid=Fluids.get_fluid(Ref);
	if (pFluid!=NULL)
	{
		if (!item.compare("EOS")){ return pFluid->BibTeXKeys.EOS; }
		else if (!item.compare("CP0")){ return pFluid->BibTeXKeys.CP0; }
		else if (!item.compare("VISCOSITY")){ return pFluid->BibTeXKeys.VISCOSITY; }
		else if (!item.compare("CONDUCTIVITY")){ return pFluid->BibTeXKeys.CONDUCTIVITY; }
		else if (!item.compare("ECS_LENNARD_JONES")){ return pFluid->BibTeXKeys.ECS_LENNARD_JONES; }
		else if (!item.compare("ECS_FITS")){ return pFluid->BibTeXKeys.ECS_FITS; }
		else if (!item.compare("SURFACE_TENSION")){ return pFluid->BibTeXKeys.SURFACE_TENSION; }
		else{ return "Bad key";}
	}
	else{
		return std::string("");
	}
}

std::string get_EOSReference(std::string Ref)
{
	pFluid=Fluids.get_fluid(Ref);
	if (pFluid!=NULL)
	{
		return pFluid->get_EOSReference();
	}
	else
		return std::string("");
}
std::string get_TransportReference(std::string Ref)
{
	pFluid=Fluids.get_fluid(Ref);
	if (pFluid!=NULL)
	{
		return pFluid->get_TransportReference();
	}
	else
		return std::string("");
}
EXPORT_CODE void CONVENTION get_aliases(char* Ref, char *aliases)
{
	strcpy(aliases, get_aliases(std::string(Ref)).c_str());
}
std::string get_aliases(std::string Ref)
{
	pFluid=Fluids.get_fluid(Ref);
	if (pFluid == NULL){
		return std::string("Fluid not found");
	}
	std::string s;
	
	std::vector<std::string> v = pFluid->get_aliases();
	for (unsigned long i = 0; i< v.size(); i++){
		if (i==0)
		{
			s = v[i];
		}
		else
		{
			s += ", " + v[i];
		}
	}
	return s;
}
EXPORT_CODE void CONVENTION get_REFPROPname(char* Ref, char * str)
{
	str= (char*)get_REFPROPname(std::string(Ref)).c_str();
}
std::string get_REFPROPname(std::string Ref)
{
	pFluid=Fluids.get_fluid(Ref);

	if ( (pFluid->get_REFPROPname()).size()!=0 ){
		return pFluid->get_REFPROPname();
	}
	else{
		return pFluid->get_name();
	}
}
