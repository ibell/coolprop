
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


/** Handle all the objects in a single list of incompressible liquids
 *  and a list of solutions / brines. The singleton pattern assures
 *  there is only one list.
 */
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

double pIncompLiquid(long iOutput, double T, double p, IncompressibleLiquid *pLiquid)
{
	double p_SI = convert_from_unit_system_to_SI(iP,p,get_standard_unit_system());
	double out;
	switch (iOutput)
	{
		case iT:
			out = T; break;
		case iP:
			out = p; break;
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
		default:
			throw ValueError(format("Your index [%d] is invalid for IncompLiquid",iOutput));
			out=0; break;
	}
	return convert_from_SI_to_unit_system(iOutput,out,get_standard_unit_system());
}

double IncompLiquid(long iOutput, double T, double p, std::string name)
{
	return pIncompLiquid(iOutput,T,p,Liquids.get_liquid(name));
}
///
/// @param T Temperature in K
/// @param p Pressure in kPa
/// @param iFluid Index of fluid
double IncompLiquid(long iOutput, double T, double p, long iFluid)
{
	return pIncompLiquid(iOutput,T,p,Liquids.get_liquid(iFluid));
}




///** Handle all the objects in a single list of incompressible liquids
// *  and a list of solutions / brines. The singleton pattern assures
// *  there is only one list.
// */
//IncompressibleLiquid * LiquidsContainer::get_liquid(long index){
//	return liquid_list[index];
//}
//
//IncompressibleLiquid * LiquidsContainer::get_liquid(std::string name){
//	std::map<std::string,IncompressibleLiquid*>::iterator it;
//	// Try to find using the map if Fluid name is provided
//	it = liquid_map.find(name);
//	// If it is found the iterator will not be equal to end
//	if (it != liquid_map.end() )
//	{
//		// Return a pointer to the class
//		return (*it).second;
//	}
//	else{
//		return NULL;
//	}
//}
//
//void LiquidsContainer::set_liquids(std::vector<IncompressibleLiquid*> list){
//	liquid_list = list;
//	// Build the map of fluid names mapping to pointers to the liquid class instances
//	for (std::vector<IncompressibleLiquid*>::iterator it = liquid_list.begin(); it != liquid_list.end(); it++)
//	{
//		// Load up entry in map
//		liquid_map[(*it)->get_name()] = *it;
//	}
//}
//
//void LiquidsContainer::destruct() {
//	delete MInstance;
//	MInstance = 0;
//}
//
//LiquidsContainer* LiquidsContainer::MInstance = 0;
//
//LiquidsContainer::LiquidsContainer() {
//	std::vector<IncompressibleLiquid*> tmpVector;
//
//	tmpVector.push_back(new DEBLiquidClass());
//	tmpVector.push_back(new HCMLiquidClass());
//	tmpVector.push_back(new HFELiquidClass());
//	tmpVector.push_back(new PMS1LiquidClass());
//	tmpVector.push_back(new PMS2LiquidClass());
//	tmpVector.push_back(new SABLiquidClass());
//	tmpVector.push_back(new HCBLiquidClass());
//	tmpVector.push_back(new TCOLiquidClass());
//	// Add new fluids based on data sheets
//	tmpVector.push_back(new TherminolD12Class());
//	tmpVector.push_back(new TherminolVP1Class());
//	tmpVector.push_back(new Therminol72Class());
//	tmpVector.push_back(new Therminol66Class());
//
//	tmpVector.push_back(new DowthermJClass());
//	tmpVector.push_back(new DowthermQClass());
//
//	tmpVector.push_back(new Texatherm22Class());
//
//	tmpVector.push_back(new NitrateSaltClass());
//
//	tmpVector.push_back(new SylthermXLTClass());
//
//	// Now we store the vector in the variable
//	// and overwrite the map.
//	set_liquids(tmpVector);
//	// register the destructor
//	atexit(&destruct);
//}
//
//LiquidsContainer::~LiquidsContainer() {
//	while (!liquid_list.empty()) {
//		delete liquid_list.back();
//		liquid_list.pop_back();
//	}
//}
//
//
////LiquidsContainer* Liquids = LiquidsContainer::Instance();
//
//bool IsIncompressibleLiquid(std::string name)
//{
//	IncompressibleLiquid *pLiquid = LiquidsContainer::Instance()->get_liquid(name);
//	if (pLiquid == NULL){ //Not found since NULL pointer returned
//		return false;
//	}
//	else{
//		return true;
//	}
//}
//
//double pIncompLiquid(long iOutput, double T, double p, IncompressibleLiquid *pLiquid)
//{
//	double p_SI = convert_from_unit_system_to_SI(iP,p,get_standard_unit_system());
//	double out;
//	switch (iOutput)
//	{
//		case iT:
//			out = T; break;
//		case iP:
//			out = p; break;
//		case iD:
//			out = pLiquid->rho(T,p_SI); break;
//		case iC:
//			out = pLiquid->cp(T,p_SI); break;
//		case iS:
//			out = pLiquid->s(T,p_SI); break;
//		case iU:
//			out = pLiquid->u(T,p_SI); break;
//		case iH:
//			out = pLiquid->h(T,p_SI); break;
//		case iV:
//			out = pLiquid->visc(T,p_SI); break;
//		case iL:
//			out = pLiquid->cond(T,p_SI); break;
//		default:
//			throw ValueError(format("Your index [%d] is invalid for IncompLiquid",iOutput));
//			out=0; break;
//	}
//	return convert_from_SI_to_unit_system(iOutput,out,get_standard_unit_system());
//}
//
//double IncompLiquid(long iOutput, double T, double p, std::string name)
//{
//	return pIncompLiquid(iOutput,T,p,LiquidsContainer::Instance()->get_liquid(name));
//}
/////
///// @param T Temperature in K
///// @param p Pressure in kPa
///// @param iFluid Index of fluid
//double IncompLiquid(long iOutput, double T, double p, long iFluid)
//{
//	return pIncompLiquid(iOutput,T,p,LiquidsContainer::Instance()->get_liquid(iFluid));
//}

