#include "Units.h"
#include "GlobalConstants.h"

double convert_from_unit_system_to_SI(long iInput, double value, int old_system)
{
	PressureUnit PU;
	switch (iInput)
	{
	case iP:
		PU = PressureUnit();
		PU.set_in_unit_system(value, old_system);
		return PU.Pa;
	case iC:
	case iC0:
	case iS:
	case iO:
	case iH:
	case iU:
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

	switch (iInput)
	{
	case iP:
		PU = PressureUnit();
		PU.set_SI(value);
		return PU.get_in_unit_system(new_system);
	case iC:
	case iC0:
	case iS:
	case iO:
	case iH:
	case iU:
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