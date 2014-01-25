#ifndef UNITS_H
#define UNITS_H

#include "GlobalConstants.h"
#include "CoolPropTools.h"
#include "CPExceptions.h"


class StandardUnit
{
protected:
	double SI_val;

public:
	virtual ~StandardUnit(){};
	virtual void set_SI(double val) = 0;
};

class PressureUnit : public StandardUnit
{
public:
	double kPa, Pa, bar;

	PressureUnit(){};
	PressureUnit(double val, int unit)
	{
		// Convert given units to Pa
		switch (unit)
		{
		case UNIT_PA: 
			break;
		case UNIT_KPA:
			val *= 1000;
			break;
		case UNIT_BAR:
			val *= 1e5;
			break;
		default:
			throw ValueError(format("Your unit [%d] is not a pressure unit",unit).c_str());
			break;
		};
		// Set the internal variable using Pa units
		set_SI(val);
	};
	void set_in_unit_system(double val, int unit_system)
	{
		// Convert given units to Pa
		switch (unit_system)
		{
		case UNIT_SYSTEM_SI:
			break;
		case UNIT_SYSTEM_KSI:
			val *= 1000;
			break;
		default:
			throw ValueError(format("Your unit system [%d] is not valid",unit_system).c_str());
			break;
		};
		set_SI(val);
	};
	double get_in_unit_system(int unit_system)
	{
		// Convert given units to Pa
		switch (unit_system)
		{
		case UNIT_SYSTEM_SI: 
			return Pa;
		case UNIT_SYSTEM_KSI:
			return kPa;
		default:
			throw ValueError(format("Your unit system [%d] is not valid",unit_system).c_str());
			break;
		};
		return _HUGE;
	};

	/// Return the SI system value (in Pa)
	double get_SI(void)	{ return SI_val; };
	/*!
	Set the internal value in the class based on the SI unit system value
	*/
	void set_SI(double SI_val)
	{
		this->SI_val = SI_val;
		Pa = SI_val;
		kPa = SI_val/1000;
		bar = SI_val/1e5;
	};
};

double convert_from_unit_system_to_SI(long iInput, double value, int old_system);
double convert_from_SI_to_unit_system(long iInput, double value, int new_system);
double convert_from_SI_to_unit_system(std::string input, double value, std::string new_system);
double convert_from_unit_system_to_SI(std::string input, double value, std::string old_system);
double convert_from_SI_to_unit_system(char *input, double value, char *new_system);
double convert_from_unit_system_to_SI(char *input, double value, char *old_system);

// Get a conversion factor from SI to the unit system for string
// Example num = "T*T" is units of temperature^2, or "P*P" is pressure^2 or "T/P" is temperature/pressure
double conversion_factor(std::string num);

#endif
