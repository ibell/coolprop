#if defined(_MSC_VER)
#define _CRTDBG_MAP_ALLOC
#define _CRT_SECURE_NO_WARNINGS
#include <crtdbg.h>
#endif

#include "CoolProp.h"
#include "CPState.h"
#include "FluidClass.h"

static int debug_level=0;

EXPORT_CODE long CONVENTION redirect_stdout(char* file)
{
	freopen(file, "a+", stdout);
	return 0;
}

//double _Props1(char *Fluid, char *Output);
//double _Props(std::string Output,std::string Name1, double Prop1, std::string Name2, double Prop2, std::string Ref);

EXPORT_CODE int CONVENTION set_reference_stateS(char *Ref, char *reference_state)
{
	return set_reference_stateS(std::string(Ref), std::string(reference_state));
}
EXPORT_CODE int CONVENTION set_reference_stateD(char *Ref, double T, double rho, double h0, double s0)
{
	return set_reference_stateD(std::string(Ref), T, rho, h0, s0);
}
EXPORT_CODE double CONVENTION PropsS(char *Output,char* Name1, double Prop1, char* Name2, double Prop2, char * Ref)
{
	return Props(Output,Name1[0],Prop1,Name2[0],Prop2,Ref);
}
EXPORT_CODE double CONVENTION Props(char *Output,char Name1, double Prop1, char Name2, double Prop2, char * Ref)
{
	// Go to the std::string, std::string version
	double val = Props(std::string(Output),Name1,Prop1,Name2,Prop2,std::string(Ref));

	/*FILE *fp;
	fp = fopen("c:\\log_Props.txt", "a");
	fprintf(fp,"%s,%c,%g,%c,%g,%s-->%g\n",Output,Name1,Prop1,Name2,Prop2,Ref,val);
	fclose(fp);*/

	return val;
}

EXPORT_CODE double CONVENTION Props1SI(char *Name, char *Output)
{
	long current_unit_system = get_standard_unit_system();
	// Set current unit system to SI (all inputs are already SI)
	set_standard_unit_system(UNIT_SYSTEM_SI);
	double val = Props1(Name, Output);
	set_standard_unit_system(current_unit_system);
	return val;
}

EXPORT_CODE double CONVENTION PropsSI(char *Output, char *Name1, double Prop1, char *Name2, double Prop2, char * Ref)
{
	long current_unit_system = get_standard_unit_system();
	// Set current unit system to SI (all inputs are already SI)
	set_standard_unit_system(UNIT_SYSTEM_SI);
	// Go to the std::string, std::string version
	long i1 = get_param_index(Name1);
	long i2 = get_param_index(Name2);
	long iOutput = get_param_index(Output);
	double val = PropsS(Output, Name1, Prop1, Name2, Prop2, Ref);
	set_standard_unit_system(current_unit_system);
	return val;
}

// All the function interfaces that point to the single-input Props function
EXPORT_CODE double CONVENTION Props1(char *Ref, char *Output)
{
	// Go to the std::string, std::string version
	return Props1(std::string(Ref),std::string(Output));
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

EXPORT_CODE double CONVENTION DerivTerms(char *Term, double T, double rho, char * Ref)
{
	// Go to the std::string, std::string version
	double val = DerivTerms(std::string(Term),T,rho,std::string(Ref));
	return val;
}

EXPORT_CODE double CONVENTION fromSI(char *input, double value, char *new_system)
{
	return convert_from_SI_to_unit_system(input, value, new_system);
}
EXPORT_CODE double CONVENTION toSI(char *input, double value, char *old_system)
{
	return convert_from_unit_system_to_SI(input, value, old_system);
}

EXPORT_CODE int CONVENTION get_debug_level(){return debug_level;}
EXPORT_CODE void CONVENTION set_debug_level(int level){debug_level=level;}
EXPORT_CODE long CONVENTION get_Fluid_index(char * param)
{
	return get_Fluid_index(std::string(param));
}
EXPORT_CODE long CONVENTION get_index_units(long param, char * units)
{
	strcpy(units, (char*)get_index_units(param).c_str());
	return 0;
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
EXPORT_CODE long CONVENTION get_param_index(char * param)
{
	return get_param_index(std::string(param));
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
		catch (std::exception &){
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
		catch (NotImplementedError &)
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

EXPORT_CODE double CONVENTION rhosatL_anc(char* FluidName, double T)
{
	try{
		// Try to load the CoolProp Fluid
		Fluid *pFluid = get_fluid(get_Fluid_index(FluidName));
		return pFluid->rhosatL(T);
	}
	catch(NotImplementedError &){
		return -_HUGE;
	}
	return -_HUGE;
}
EXPORT_CODE double CONVENTION rhosatV_anc(char* FluidName, double T)
{
	try{
		// Try to load the CoolProp Fluid
		Fluid *pFluid = get_fluid(get_Fluid_index(FluidName));
		return pFluid->rhosatV(T);
	}
	catch(NotImplementedError &){
		return -_HUGE;
	}
	return -_HUGE;
}
EXPORT_CODE double CONVENTION psatL_anc(char* FluidName, double T)
{
	try{
		// Try to load the CoolProp Fluid
		Fluid *pFluid = get_fluid(get_Fluid_index(FluidName));
		return pFluid->psatL_anc(T);
	}
	catch(NotImplementedError &){
		return -_HUGE;
	}
	return -_HUGE;
}
EXPORT_CODE double CONVENTION psatV_anc(char* FluidName, double T)
{
	try{
		// Try to load the CoolProp Fluid
		Fluid *pFluid = get_fluid(get_Fluid_index(FluidName));
		return pFluid->psatV_anc(T);
	}
	catch(NotImplementedError &){
		return -_HUGE;
	}
	return -_HUGE;
}
#ifndef SWIG
EXPORT_CODE long CONVENTION get_global_param_string(char *param, char * Output)
{
	strcpy(Output,get_global_param_string(std::string(param)).c_str());
	return 0;
}
EXPORT_CODE long CONVENTION get_fluid_param_string(char *fluid, char *param, char * Output)
{
	strcpy(Output, get_fluid_param_string(std::string(fluid), std::string(param)).c_str());
	return 0;
}
#endif
EXPORT_CODE long CONVENTION Phase(char *Fluid, double T, double p, char *Phase_str)
{
	strcpy(Phase_str,(char*)Phase(std::string(Fluid),T,p).c_str());
	return 0;
}
EXPORT_CODE long CONVENTION Phase_Tp(char *Fluid, double T, double p, char *Phase_str)
{
	strcpy(Phase_str,(char*)Phase(std::string(Fluid),T,p).c_str());
	return 0;
}
EXPORT_CODE long CONVENTION Phase_Trho(char *Fluid, double T, double rho, char *Phase_str)
{
	strcpy(Phase_str,(char*)Phase_Trho(std::string(Fluid),T,rho).c_str());
	return 0;
}


// A function to enforce the state if known
EXPORT_CODE void CONVENTION set_phase(char *Phase_str){
	set_phase(std::string(Phase_str));
}

/// Enable the TTSE for this fluid
EXPORT_CODE bool CONVENTION enable_TTSE_LUT(char *FluidName){
	long iFluid = get_Fluid_index(FluidName); if (iFluid<0){ return false; };
	Fluid *pFluid = get_fluid(iFluid);
	pFluid->enable_TTSE_LUT();
	return true;
};
/// Check if TTSE is enabled
EXPORT_CODE bool CONVENTION isenabled_TTSE_LUT(char *FluidName){
	long iFluid = get_Fluid_index(FluidName); if (iFluid<0){ return false; };
	Fluid *pFluid = get_fluid(iFluid);
	return pFluid->isenabled_TTSE_LUT();
}
/// Disable the TTSE for this fluid
EXPORT_CODE bool CONVENTION disable_TTSE_LUT(char *FluidName){
	long iFluid = get_Fluid_index(FluidName); if (iFluid<0){ return false; };
	Fluid *pFluid = get_fluid(iFluid);
	pFluid->disable_TTSE_LUT();
	return true;
}
/// Enable the writing of TTSE tables to file for this fluid
EXPORT_CODE bool CONVENTION enable_TTSE_LUT_writing(char *FluidName){
	long iFluid = get_Fluid_index(FluidName); if (iFluid<0){ return true;};
	Fluid *pFluid = get_fluid(iFluid);
	pFluid->enable_TTSE_LUT_writing();
	return true;
};
/// Check if the writing of TTSE tables to file is enabled
EXPORT_CODE bool CONVENTION isenabled_TTSE_LUT_writing(char *FluidName){
	long iFluid = get_Fluid_index(FluidName); if (iFluid<0){ return false;};
	Fluid *pFluid = get_fluid(iFluid);
	return pFluid->isenabled_TTSE_LUT_writing();
}
/// Disable the writing of TTSE tables to file for this fluid
EXPORT_CODE bool CONVENTION disable_TTSE_LUT_writing(char *FluidName){
	long iFluid = get_Fluid_index(FluidName); if (iFluid<0){ return false;};
	Fluid *pFluid = get_fluid(iFluid);
	pFluid->disable_TTSE_LUT_writing();
	return true;
}
/// Over-ride the default size of both of the saturation LUT
EXPORT_CODE bool CONVENTION set_TTSESat_LUT_size(char *FluidName, int Nsat){
	long iFluid = get_Fluid_index(FluidName); if (iFluid<0){ return false; };
	Fluid *pFluid = get_fluid(iFluid);
	pFluid->set_TTSESat_LUT_size(Nsat);
	return true;
}
/// Over-ride the default size of the single-phase LUT
EXPORT_CODE bool CONVENTION set_TTSESinglePhase_LUT_size(char *FluidName, int Np, int Nh){
	long iFluid = get_Fluid_index(FluidName); if (iFluid<0){ return false;};
	Fluid *pFluid = get_fluid(iFluid);
	pFluid->set_TTSESinglePhase_LUT_size(Np,Nh);
	return true;
}
/// Over-ride the default range of the single-phase LUT
EXPORT_CODE bool CONVENTION set_TTSESinglePhase_LUT_range(char *FluidName, double hmin, double hmax, double pmin, double pmax){
	long iFluid = get_Fluid_index(FluidName); if (iFluid<0){ return false;};
	Fluid *pFluid = get_fluid(iFluid);
	pFluid->set_TTSESinglePhase_LUT_range(hmin,hmax,pmin,pmax);
	return true;
}
/// Get the current range of the single-phase LUT
EXPORT_CODE bool CONVENTION get_TTSESinglePhase_LUT_range(char *FluidName, double *hmin, double *hmax, double *pmin, double *pmax){
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
EXPORT_CODE int CONVENTION get_standard_unit_system(void){return _get_standard_unit_system();}
/// Sets the flag for the integer flag corresponding to the current set of units
/// @param val The integer value for the current set of units, one of enumerated values UNIT_SYSTEM_SI, UNIT_SYSTEM_KSI (see GlobalConstants.h)
EXPORT_CODE void CONVENTION set_standard_unit_system(int val){return _set_standard_unit_system(val);}
