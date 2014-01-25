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
		void ECSParams(double *e_k, double *sigma){
		// Poling
		*e_k = 809.1; *sigma = 0.2641; 
		};
		double surface_tension_T(double T);
	};

#endif
