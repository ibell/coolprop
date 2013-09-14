#include "Units.h"
#include "GlobalConstants.h"
#include "CoolProp.h"

double convert_from_unit_system_to_SI(long iInput, double value, int old_system)
{
	if (get_debug_level() > 8)
	{
		std::cout << format("convert_from_unit_system_to_SI(%d,%g,%d)\n",iInput,value,old_system).c_str();
	}
	if (old_system == UNIT_SYSTEM_SI)
	{
		return value;
	}

	switch (iInput)
	{
	case iP:
	case iPtriple:
	case iPcrit:
	case iC:
	case iC0:
	case iS:
	case iG:
	case iO:
	case iH:
	case iU:
	case iL:
		switch (old_system)
		{
		case UNIT_SYSTEM_KSI:
			return value*1000.0;
		case UNIT_SYSTEM_SI:
			return value;
		default:
			throw ValueError(format("Unit System [%d] is undefined",old_system).c_str());
		}
	case iD:
	case iA:
	case iQ:
	case iV:
	case iT:
	case iTcrit:
	case iTtriple: 
	case iTreduce:
	case iTmin:
	case iPHASE_LIQUID:
	case iPHASE_GAS:
	case iPHASE_SUPERCRITICAL:
	case iPHASE_TWOPHASE:
	case iODP:
	case iGWP100:
	case iMM:
	case iRhocrit:
	case iRhoreduce: 
	case iAccentric:
		return value;
	default:
		throw ValueError(format("index [%d] is invalid in convert_from_unit_system_to_SI",iInput).c_str());
	}
}

double convert_from_SI_to_unit_system(long iInput, double value, int new_system)
{
	if (get_debug_level() > 8)
	{
		std::cout << format("convert_from_SI_to_unit_system(%d,%g,%d)\n",iInput,value,new_system).c_str();
	}

	if (new_system == UNIT_SYSTEM_SI)
	{
		return value;
	}

	switch (iInput)
	{
	case iP:
	case iPtriple:
	case iPcrit:
	case iC:
	case iC0:
	case iS:
	case iG:
	case iO:
	case iH:
	case iU:
	case iL:
		switch (new_system)
		{
		case UNIT_SYSTEM_KSI:
			return value/1000.0;
		case UNIT_SYSTEM_SI:
			return value;
		default:
			throw ValueError(format("Unit System [%d] is undefined",new_system).c_str());
		}
	case iD:
	case iQ:
	case iA:
	case iV:
	case iT:
	case iTmin:
	case iTcrit:
	case iTtriple: 
	case iTreduce:
	case iPHASE_LIQUID:
	case iPHASE_GAS:
	case iPHASE_SUPERCRITICAL:
	case iPHASE_TWOPHASE:
	case iODP:
	case iGWP100:
	case iMM:
	case iRhocrit:
	case iRhoreduce: 
	case iAccentric:
		return value;
	default:
		throw ValueError(format("index [%d] is invalid in convert_from_SI_to_unit_system",iInput).c_str());
	}
}