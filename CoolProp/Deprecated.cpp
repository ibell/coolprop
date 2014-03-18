#if defined(_MSC_VER)
#define _CRTDBG_MAP_ALLOC
#define _CRT_SECURE_NO_WARNINGS
#include <crtdbg.h>
#endif

#include "CoolProp.h"
#include "CPState.h"
#include "FluidClass.h"

#include "Deprecated.h"

EXPORT_CODE double CONVENTION conformal_Trho(const char* FluidName, const char* ReferenceFluidName, double T, double rho, double *Tconform, double *rhoconform)
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

EXPORT_CODE double CONVENTION viscosity_dilute(const char* FluidName, double T)
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
EXPORT_CODE double CONVENTION viscosity_residual(const char* FluidName, double T, double rho)
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

EXPORT_CODE double CONVENTION conductivity_critical(const char* FluidName, double T, double rho)
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
EXPORT_CODE double CONVENTION conductivity_background(const char* FluidName, double T, double rho)
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

EXPORT_CODE double CONVENTION rhosatL_anc(const char* FluidName, double T)
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
EXPORT_CODE double CONVENTION rhosatV_anc(const char* FluidName, double T)
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
EXPORT_CODE double CONVENTION psatL_anc(const char* FluidName, double T)
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
EXPORT_CODE double CONVENTION psatV_anc(const char* FluidName, double T)
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
