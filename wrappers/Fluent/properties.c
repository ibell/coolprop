/* 
* UDF TO CALCULATE FLUID PROPERTIES BASED ON THE OPEN-SOURCE 
* THERMODYNAMIC LIBRARY COOLPROP
*/

#include "udf.h"

const char* FLUID = "CarbonDioxide";
const real gauge = 101325; //operating pressure in pascal (as defined in fluent)

double Props (char*, char, double, char, double, char*);

/* REAL GAS VISCOSITY
 * Coolprop returns in [Pa s], Fluent admits in [Pa s]
 */
DEFINE_PROPERTY(cell_viscosity, c, t)
{
	real viscosity;
	real temperature = C_T(c, t);
	real pressure = C_P(c, t)+gauge;
	viscosity = Props((char*)"V",'T',temperature,'P',pressure/1000, (char*)FLUID);
	return viscosity;
}

/* REAL GAS THERMAL CONDUCTIVITY
 * Coolprop returns in [kW/(m K)], Fluent admits in [W/(m K)]
*/
DEFINE_PROPERTY(cell_thermalConductivity, c, t)
{
	real thermalConductivity;
	real temperature = C_T(c, t);
	real pressure = C_P(c, t)+gauge;
	thermalConductivity = Props((char*)"L",'T', temperature, 'P', pressure/1000, (char*)FLUID);
	return thermalConductivity*1000; 
}

/* REAL GAS SPECIFIC MASS
 * Coolprop returns in [kg/m³], Fluent admits in [kg/m³]
*/
DEFINE_PROPERTY(cell_density, c, t)
{
	real density;
	real temperature = C_T(c, t);
	real pressure = C_P(c, t)+gauge;
	density = Props((char*)"D",'T', temperature, 'P', pressure/1000, (char*)FLUID);
	return density;
}

/* REAL GAS SPECIFIC HEAT
 * Coolprop returns in [kJ/(kg K)], Fluent admits in [J/kg K]
*/
DEFINE_SPECIFIC_HEAT(cell_specificHeat, temperature, Tref, enthalpy, yi)
{	 
	real pressure;
	pressure = gauge;
	
	/* The following commented code is supposed to get the pressure
	 from the cell to use with Coolprop. Left commented because
	 specific heat depends very little on pressure. Will increase
	 computational time significantly. */
	 
	/*
	* Domain *domain = Get_Domain(1);
	* Thread *t;
	* cell_t c;
	* thread_loop_c(t, domain)
	* {
	* 	begin_c_loop(c, t)
	* 	{
	* 		pressure = C_P(c, t);
	* 	}end_c_loop(c, t)
	* }
	*/
	real specificHeat;
	
	specificHeat = Props((char*)"C",'T', temperature, 'P', pressure/1000, (char*)FLUID)*1000;
	*enthalpy = specificHeat*(temperature-Tref);
	return specificHeat;
}

/* Execute on demand UDF to test if the library was built correctly */
DEFINE_ON_DEMAND(call_coolprop)
{
	real p, t, density, specificHeat, thermalConductivity, enthalpy, viscosity;

	p = 100000.0;
	t = 300.0;

	density = Props((char*)"D",'T',t,'P',p/1000,(char*)FLUID);
	specificHeat = Props((char*)"C",'T',t,'P',p/1000,(char*)FLUID);
	viscosity = Props((char*)"V",'T',t,'P',p/1000,(char*)FLUID);
	thermalConductivity = Props((char*)"L",'T',t,'P',p/1000,(char*)FLUID);
	enthalpy = Props((char*)"H",'T',t,'P',p/1000,(char*)FLUID);

	Message("p = %lf, T = %lf => density = %lf\n", p/1000, t, density);
	Message("p = %lf, T = %lf => specific heat = %lf\n", p/1000, t, specificHeat*1000);
	Message("p = %lf, T = %lf => viscosity = %lf\n", p/1000, t, viscosity);
	Message("p = %lf, T = %lf => thermal conductivity = %lf\n", p/1000, t, thermalConductivity*1000);
	Message("p = %lf, T = %lf => enthalpy = %lf\n", p/1000, t, enthalpy*1000);
}
