
#if defined(_MSC_VER)
#define _CRTDBG_MAP_ALLOC
#define _CRT_SECURE_NO_WARNINGS
#include <stdlib.h>
#include <crtdbg.h>
#else
#include <stdlib.h>
#endif

#include "Units.h"
#include "math.h"
#include "stdio.h"
#include <string.h>
#include "Brine.h"
#include "CoolPropTools.h"
#include "CoolProp.h"
#include "IncompBase.h"
#include "IncompLiquid.h"


/// Base class for simplified brine/solution models
/** Employs the base functions implemented in IncompBase.h and
 *  provides properties as function of temperature, pressure
 *  and composition. */
void IncompressibleLiquid::testInputs(double T_K, double p){
	double result = 0.;

	printf(" %s \n"," ");
	printf("Testing  %s \n",this->get_name().c_str());
	printf("Inputs:  T = %3.3f degC \t p = %2.4f bar \n",T_K-273.15,p/1e5);

    result = this->rho(T_K,p);
    printf("From object:    rho = %4.2f \t kg/m3    \n",result);
    result = this->cp(T_K,p);
    printf("From object:     cp = %1.5f \t kJ/kg-K  \n",result/1e3);
    result = this->h(T_K,p);
	printf("From object:      h = %3.3f \t kJ/kg    \n",result/1e3);
	result = this->s(T_K,p);
	printf("From object:      s = %1.5f \t kJ/kg-K  \n",result/1e3);
	result = this->visc(T_K,p);
	printf("From object:    eta = %1.5f \t 1e-5 Pa-s\n",result*1e5);
	result = this->cond(T_K,p);
	printf("From object: lambda = %1.5f \t W/m-k    \n",result*1e3);
	result = this->u(T_K,p);
	printf("From object:      u = %3.3f \t kJ/kg    \n",result/1e3);
	result = this->psat(T_K);
	printf("From object:   psat = %2.4f \t bar      \n",result/1e5);
}


/*
 * Some more functions to provide a single implementation
 * of important routines.
 * We start with the check functions that can validate input
 * in terms of pressure p and temperature T.
 */

/// Check validity of temperature input.
/** Compares the given temperature T to a stored minimum and
 *  maximum temperature. Enforces the redefinition of Tmin and
 *  Tmax since the default values cause an error. */
bool IncompressibleLiquid::checkT(double T_K){
	if( Tmin < 0. ) {
		throw ValueError("Please specify the minimum temperature.");
	} else if( Tmax < 0.) {
		throw ValueError("Please specify the maximum temperature.");
	} else if ( (Tmin>T_K) || (T_K>Tmax) ) {
		throw ValueError(format("Your temperature %f is not between %f and %f.",T_K,Tmin,Tmax));
	} else {
		return true;
	}
	return false;
}

/// Check validity of pressure input.
/** Compares the given pressure p to the saturation pressure at
 *  temperature T and throws and exception if p is lower than
 *  the saturation conditions.
 *  The default value for psat is -1 yielding true if psat
 *  is not redefined in the subclass.
 *  */
bool IncompressibleLiquid::checkP(double T_K, double p) {
	double ps = psat(T_K);
	if (p<ps) {
		throw ValueError(format("Equations are valid for liquid phase only: %f < %f. ",p,ps));
	} else {
		return true;
	}
}

/// Check validity of temperature and pressure input.
bool IncompressibleLiquid::checkTP(double T, double p) {
	return (checkT(T) && checkP(T,p));
}



/** Handle all the objects in a single list of incompressible
 *  liquids. */
LiquidsContainer::LiquidsContainer() {
	std::vector<IncompressibleLiquid*> tmpVector;

	tmpVector.push_back(new DEBLiquidClass());
	tmpVector.push_back(new HCMLiquidClass());
	tmpVector.push_back(new HFELiquidClass());
	tmpVector.push_back(new PMS1LiquidClass());
	tmpVector.push_back(new PMS2LiquidClass());
	tmpVector.push_back(new SABLiquidClass());
	tmpVector.push_back(new HCBLiquidClass());
	tmpVector.push_back(new TCOLiquidClass());
	// Add new fluids based on data sheets
	tmpVector.push_back(new TherminolD12Class());
	tmpVector.push_back(new TherminolVP1Class());
	tmpVector.push_back(new Therminol72Class());
	tmpVector.push_back(new Therminol66Class());

	tmpVector.push_back(new DowthermJClass());
	tmpVector.push_back(new DowthermQClass());

	tmpVector.push_back(new Texatherm22Class());

	tmpVector.push_back(new NitrateSaltClass());

	tmpVector.push_back(new SylthermXLTClass());

    tmpVector.push_back(new HC10Class());
    tmpVector.push_back(new HC20Class());
    tmpVector.push_back(new HC30Class());
    tmpVector.push_back(new HC40Class());
    tmpVector.push_back(new HC50Class());
    
    // Add new fluids based SecCool software
    tmpVector.push_back(new AS10Class());
    tmpVector.push_back(new AS20Class());
    tmpVector.push_back(new AS30Class());
    tmpVector.push_back(new AS40Class());
    tmpVector.push_back(new AS55Class());
    tmpVector.push_back(new ZS10Class());
    tmpVector.push_back(new ZS25Class());
    tmpVector.push_back(new ZS40Class());
    tmpVector.push_back(new ZS45Class());
    tmpVector.push_back(new ZS55Class());


	// Now we store the vector in the variable
	// and overwrite the map.
	set_liquids(tmpVector);
}

LiquidsContainer::~LiquidsContainer() {
	while (!liquid_list.empty()) {
		delete liquid_list.back();
		liquid_list.pop_back();
	}
}

IncompressibleLiquid * LiquidsContainer::get_liquid(long index){
	return liquid_list[index];
}

IncompressibleLiquid * LiquidsContainer::get_liquid(std::string name){
	std::map<std::string,IncompressibleLiquid*>::iterator it;
	// Try to find using the map if Fluid name is provided
	it = liquid_map.find(name);
	// If it is found the iterator will not be equal to end
	if (it != liquid_map.end() )
	{
		// Return a pointer to the class
		return (*it).second;
	}
	else{
		return NULL;
	}
}

void LiquidsContainer::set_liquids(std::vector<IncompressibleLiquid*> list){
	liquid_list = list;
	// Build the map of fluid names mapping to pointers to the liquid class instances
	for (std::vector<IncompressibleLiquid*>::iterator it = liquid_list.begin(); it != liquid_list.end(); it++)
	{
		// Load up entry in map
		liquid_map[(*it)->get_name()] = *it;
	}
}

LiquidsContainer Liquids = LiquidsContainer();

bool IsIncompressibleLiquid(std::string name)
{
	IncompressibleLiquid *pLiquid = Liquids.get_liquid(name);
	if (pLiquid == NULL){ //Not found since NULL pointer returned
		return false;
	}
	else{
		return true;
	}
}


double pIncompLiquidSI(long iOutput, double T, double p_SI, IncompressibleLiquid *pLiquid)
{
	double out;
	switch (iOutput)
	{
		case iT:
			out = T; break;
		case iP:
			out = p_SI; break;
		case iD:
			out = pLiquid->rho(T,p_SI); break;
		case iC:
			out = pLiquid->cp(T,p_SI); break;
		case iS:
			out = pLiquid->s(T,p_SI); break;
		case iU:
			out = pLiquid->u(T,p_SI); break;
		case iH:
			out = pLiquid->h(T,p_SI); break;
		case iV:
			out = pLiquid->visc(T,p_SI); break;
		case iL:
			out = pLiquid->cond(T,p_SI); break;
		case iTmin:
			out = pLiquid->getTmin(); break;
		case iTmax:
			out = pLiquid->getTmax(); break;
		case iPsat:
			out = pLiquid->psat(T); break;
		default:
			throw ValueError(format("Your index [%d] is invalid for the incompressible liquid %s",iOutput,pLiquid->getName().c_str()));
			out=0; break;
	}
	return out; 
}

double pIncompLiquid(long iOutput, double T, double p, IncompressibleLiquid *pLiquid)
{
    // Convert to SI pressure
    double p_SI = convert_from_unit_system_to_SI(iP,p,get_standard_unit_system());
    
    // Call the function
    double out = pIncompLiquidSI(iOutput,T,p_SI,pLiquid);
    
    // Convert back to current unit system
    return convert_from_SI_to_unit_system(iOutput,out,get_standard_unit_system());
}

double IncompLiquidSI(long iOutput, double T, double p_SI, long iFluid)
{
	return pIncompLiquidSI(iOutput,T,p_SI,Liquids.get_liquid(iFluid));
}

double IncompLiquidSI(long iOutput, double T, double p_SI, std::string name)
{
	return pIncompLiquidSI(iOutput,T,p_SI,Liquids.get_liquid(name));
}


