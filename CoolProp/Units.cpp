#include "Units.h"
#include "GlobalConstants.h"

double convert_from_unit_system_to_SI(long iInput, double value, int old_system)
{
	if (old_system == UNIT_SYSTEM_SI)
	{
		return value;
	}

	switch (iInput)
	{
	case iP:
	case iC:
	case iC0:
	case iS:
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
			throw ValueError();
		}
	case iT:
		return value;
	default:
		throw ValueError(format("index [%d] is invalid convert_from_unit_system_to_SI",iInput).c_str());
	}
}

double convert_from_SI_to_unit_system(long iInput, double value, int new_system)
{
	if (new_system == UNIT_SYSTEM_SI)
	{
		return value;
	}

	switch (iInput)
	{
	case iP:
	case iC:
	case iC0:
	case iS:
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
			throw ValueError();
		}
	case iT:
		return value;
	default:
		throw ValueError(format("index [%d] is invalid in convert_from_SI_to_unit_system",iInput).c_str());
	}
}