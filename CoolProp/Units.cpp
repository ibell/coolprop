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
	case iT:
		return value;
	default:
		return value;
	}
}

//double convert_between_unit_systems(long iInput, double value, int old_system, int new_system)
//{
//	if (old_system == new_system)
//	{
//		return value;
//	}
//	else
//	{
//		double SI_val = convert_from_unit_system_to_SI(iInput, value, old_system);
//
//		switch (iInput)
//		{
//		case iP:
//			switch (old_system)
//			{
//			default:
//				return value;
//			}
//		default:
//			return value;
//		}
//	}
//}