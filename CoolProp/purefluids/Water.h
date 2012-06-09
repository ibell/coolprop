#include "FluidClass.h"
#ifndef _WATER_H
#define WATER_H

	class WaterClass : public Fluid{

	public:
		WaterClass();
		~WaterClass(){};
		virtual double conductivity_Trho(double, double);
		virtual double viscosity_Trho(double, double);
		double psat(double);
		double rhosatL(double);
		double rhosatV(double);
		double surface_tension_T(double T);
	};

#endif
