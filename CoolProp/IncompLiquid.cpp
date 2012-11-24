
#if defined(_MSC_VER)
#define _CRTDBG_MAP_ALLOC
#define _CRT_SECURE_NO_WARNINGS
#include <stdlib.h>
#include <crtdbg.h>
#else
#include <stdlib.h>
#endif

#include "math.h"
#include "stdio.h"
#include <string.h>
#include "Brine.h"

#include "CoolProp.h"
#include "IncompLiquid.h"

class LiquidsContainer{
private:
	std::vector<IncompressibleLiquid*> LiquidsList;
	std::map<std::string,IncompressibleLiquid*> liquid_map;
public:
	LiquidsContainer(){
		LiquidsList.push_back(new DEBLiquidClass());
		LiquidsList.push_back(new HCMLiquidClass());
		LiquidsList.push_back(new HFELiquidClass());
		LiquidsList.push_back(new PMS1LiquidClass());
		LiquidsList.push_back(new PMS2LiquidClass());
		LiquidsList.push_back(new SABLiquidClass());
		LiquidsList.push_back(new HCBLiquidClass());
		LiquidsList.push_back(new TCOLiquidClass());

		// Build the map of fluid names mapping to pointers to the Fluid class instances
		for (std::vector<IncompressibleLiquid*>::iterator it = LiquidsList.begin(); it != LiquidsList.end(); it++)
		{
			// Load up entry in map
			liquid_map[(*it)->get_name()] = *it;
		}
	};
	~LiquidsContainer(){
		while (!LiquidsList.empty())
		{
			delete LiquidsList.back();
			LiquidsList.pop_back();
		}
	};

	IncompressibleLiquid * get_liquid(long index){
		return LiquidsList[index];
	}

	IncompressibleLiquid * get_liquid(std::string name){
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
};

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
	switch (iOutput)
	{
		case iT:
			return T; break;
		case iP:
			return p; break;
		case iD:
			return pLiquid->rho(T); break;
		case iC:
			return pLiquid->cp(T); break;
		case iS:
			return pLiquid->s(T); break;
		case iU:
			return pLiquid->u(T); break;
		case iH:
			return pLiquid->h(T,p); break;
		case iV:
			return pLiquid->visc(T); break;
		case iL:
			return pLiquid->cond(T); break;
		default:
			throw ValueError(format("Your index [%d] is invalid for IncompLiquid",iOutput)); break;
	}
	return pLiquid->rho(T);
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

