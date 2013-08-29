#include "Units.h"
#include "GlobalConstants.h"

double convert_from_unit_system_to_SI(long iInput, double value, int old_system)
{
	PressureUnit PU;
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
		return value;
	}
}

double convert_from_SI_to_unit_system(long iInput, double value, int new_system)
{
	PressureUnit PU;
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
		return value;
	}
}